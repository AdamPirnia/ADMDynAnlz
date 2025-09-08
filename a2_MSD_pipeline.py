#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
a2_MSD Pipeline GUI - Optimized Version

A comprehensive molecular dynamics analysis pipeline for generating Mean Square 
Displacement (MSD) and non-Gaussian parameter analyses from NAMD DCD trajectories.

Features:
- α₂(t) and MSD calculation: Standard non-Gaussian parameter α₂(t) = 3⟨Δr⁴⟩/(5⟨Δr²⟩²) - 1
- α_xz(t) calculation: Directional correlation parameter α_xz(t) = ⟨Δx²·Δz²⟩/(⟨Δx²⟩·⟨Δz²⟩) - 1
- 3-10x performance improvements with parallel processing
- Optimized memory usage and enhanced numerical stability
- Comprehensive error handling and data validation
- Cross-platform support (Linux, macOS, Windows)
- Professional GUI with tooltips and progress monitoring
The pipeline processes DCD trajectories through four main steps:
1. coordinates_extract: Extract raw coordinates using VMD
2. unwrap_coords: Remove periodic boundary artifacts  
3. COM_calc: Calculate center-of-mass trajectories
4. alpha2_MSD or alpha anisotropy: Compute selected statistical parameters

All functions support parallel processing, chunked memory management, and 
provide detailed progress reporting and data quality metrics.
"""
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import json
import multiprocessing as mp
import shutil
import sys

# Configure ttk style for better appearance
def configure_ttk_style():
    """Configure ttk widgets for better visual appearance"""
    style = ttk.Style()
    
    # Configure Entry style
    style.configure('Smooth.TEntry',
                   fieldbackground='white',
                   borderwidth=1,
                   relief='solid',
                   focuscolor='#3498db')
    
    # Configure Checkbutton style  
    style.configure('Smooth.TCheckbutton',
                   background='#f0f0f0',
                   focuscolor='none')
    
    return style

# Path to save & load last inputs
CONFIG_PATH = os.path.expanduser("~/.pipeline_gui_config.json")

class ToolTip:
    """
    Create a tooltip for a given widget
    """
    def __init__(self, widget, text='widget info'):
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.enter)
        self.widget.bind("<Leave>", self.leave)
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def enter(self, event=None):
        self.schedule()

    def leave(self, event=None):
        self.unschedule()
        self.hidetip()

    def schedule(self):
        self.unschedule()
        self.id = self.widget.after(500, self.showtip)  # 500ms delay

    def unschedule(self):
        id = self.id
        self.id = None
        if id:
            self.widget.after_cancel(id)

    def showtip(self, event=None):
        x, y, cx, cy = self.widget.bbox("insert") if hasattr(self.widget, 'bbox') else (0, 0, 0, 0)
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        
        # Creates a toplevel window
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry("+%d+%d" % (x, y))
        
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                      background="#f8f8f8", relief=tk.SOLID, borderwidth=1,
                      font=("Arial", 10, "normal"), wraplength=350, fg='#2c3e50')
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def create_tooltip(widget, text):
    """Helper function to create tooltips easily"""
    return ToolTip(widget, text)

class CalculationsWindow(tk.Toplevel):
    """Separate window for Alpha2/MSD and Dipole calculations"""
    
    def __init__(self, parent):
        super().__init__(parent)
        self.parent = parent
        self.title("Analysis Calculations - Alpha2/MSD & Dipole Moments")
        self.geometry("1100x800")  # Decreased width by 300, reduced height for better fit
        self.configure(bg='#f0f0f0')
        
        # Make window modal
        self.transient(parent)
        self.grab_set()
        
        # Initialize variables
        self.init_variables()
        
        # Create UI
        self.create_ui()
        
        # Load configuration if available
        self.load_calculation_config()
        
        # Center window on parent
        self.center_on_parent()
    
    def init_variables(self):
        """Initialize all calculation variables"""
        # Common parameters (shared with main window)
        self.common_term_var = tk.StringVar(value="")
        self.baseDir_var = tk.StringVar(value="")
        self.num_dcd_var = tk.StringVar(value="")
        self.num_particles_var = tk.StringVar(value="")
        self.max_workers_var = tk.StringVar(value="4")
        self.output_folder_var = tk.StringVar(value="")
        
        # SLURM parameters
        self.sbatch_vars = []
        for i in range(10):  # 10 SLURM parameters
            self.sbatch_vars.append(tk.StringVar(value=""))
        
        # Output files parameters
        self.mainfile_var = tk.StringVar(value="")
        self.submitfile_var = tk.StringVar(value="")
        
        # Alpha2/MSD variables
        self.skip_alpha2 = tk.BooleanVar(value=False)
        self.calc_type_var = tk.StringVar(value="alpha2_msd")
        self.a2_vars = []
        self.a2_chunk_processing = tk.BooleanVar(value=True)
        self.a2_validate = tk.BooleanVar(value=True)
        
        # Dipole calculation variables
        self.skip_dipole = tk.BooleanVar(value=True)  # Skip dipole by default
        self.dipole_calc_type = tk.StringVar(value="individual")  # "individual" or "collective"
        self.dipole_vars = []  # Will be set to individual_dipole_vars or collective_dipole_vars
        self.individual_dipole_vars = []
        self.collective_dipole_vars = []
        self.dipole_parallel = tk.BooleanVar(value=True)
        self.dipole_validate = tk.BooleanVar(value=True)
        self.dipole_chunk_processing = tk.BooleanVar(value=True)
        
        # Velocity extraction variables
        self.skip_velocity = tk.BooleanVar(value=True)  # Skip velocity by default
        self.velocity_vars = []
        self.velocity_parallel = tk.BooleanVar(value=True)
        self.velocity_validate = tk.BooleanVar(value=True)
        
        # DCD selection variables for analysis calculations
        self.alpha2_dcd_selection = tk.StringVar(value="")
        self.individual_dipole_dcd_selection = tk.StringVar(value="")
        self.collective_dipole_dcd_selection = tk.StringVar(value="")
        self.velocity_dcd_selection = tk.StringVar(value="")
    
    def create_ui(self):
        """Create the user interface with scrolling"""
        # Create canvas and scrollbar for scrolling
        self.canvas = tk.Canvas(self, bg='#f0f0f0', highlightthickness=0)
        self.scrollbar = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = tk.Frame(self.canvas, bg='#f0f0f0')
        
        # Configure scrolling
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )
        
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        
        # Pack canvas and scrollbar
        self.canvas.pack(side="left", fill="both", expand=True, padx=(10, 0), pady=10)
        self.scrollbar.pack(side="right", fill="y", pady=10)
        
        # Bind mousewheel to canvas
        self.bind_all("<MouseWheel>", self._on_mousewheel_calc)
        self.bind_all("<Button-4>", self._on_mousewheel_calc)  # Linux scroll up
        self.bind_all("<Button-5>", self._on_mousewheel_calc)  # Linux scroll down
        
        # Header
        header_frame = tk.Frame(self.scrollable_frame, bg='#f0f0f0')
        header_frame.pack(fill="x", pady=(0, 15))
        
        tk.Label(header_frame, text="Analysis Calculations", 
                font=("Arial", 16, "bold"), fg='#2c3e50', bg='#f0f0f0').pack()
        tk.Label(header_frame, text="Configure Alpha2/MSD and Dipole moment calculations", 
                font=("Arial", 10), fg='#7f8c8d', bg='#f0f0f0').pack()
        
        # Common Parameters section
        self.create_common_parameters_section(self.scrollable_frame)
        
        # Create notebook for tabs
        notebook = ttk.Notebook(self.scrollable_frame)
        notebook.pack(fill=tk.BOTH, expand=True, pady=(0, 15))
        
        # Alpha2/MSD tab
        alpha2_frame = tk.Frame(notebook, bg='#f0f0f0')
        notebook.add(alpha2_frame, text="α₂(t) and MSD Calculation")
        self.create_alpha2_section(alpha2_frame)
        
        # Dipole tab
        dipole_frame = tk.Frame(notebook, bg='#f0f0f0')
        notebook.add(dipole_frame, text="Dipole Moments")
        self.create_dipole_section(dipole_frame)
        
        # Velocity Extraction tab
        velocity_frame = tk.Frame(notebook, bg='#f0f0f0')
        notebook.add(velocity_frame, text="Velocity Extraction")
        self.create_velocity_section(velocity_frame)
        
        # SLURM and Output Files section (moved to bottom like main window)
        self.create_slurm_output_section(self.scrollable_frame)
        
        # Buttons frame
        buttons_frame = tk.Frame(self.scrollable_frame, bg='#f0f0f0')
        buttons_frame.pack(fill="x", pady=20)
        
        # Generate button
        generate_btn = tk.Button(buttons_frame, text="Generate Calculation Scripts", 
                                command=self.generate_calculation_files, 
                                bg="#52c77a", fg="black", 
                                font=("Arial", 12, "bold"), padx=20, pady=10,
                                relief='raised', borderwidth=2, cursor='hand2')
        generate_btn.pack(side=tk.LEFT, padx=5)
        
        # Close button
        close_btn = tk.Button(buttons_frame, text="Close", 
                             command=self.close_window, 
                             bg="#e74c3c", fg="white", 
                             font=("Arial", 10), padx=15, pady=10,
                             relief='raised', borderwidth=2, cursor='hand2')
        close_btn.pack(side=tk.RIGHT, padx=5)
    
    def save_calculation_config(self):
        """Save Analysis calculations window configuration"""
        import json
        
        data = {
            # Alpha2/MSD settings
            "skip_alpha2": self.skip_alpha2.get(),
            "calc_type_var": self.calc_type_var.get(),
            "a2_vars": [v.get() for v in self.a2_vars] if hasattr(self, 'a2_vars') else [],
            "a2_chunk_processing": self.a2_chunk_processing.get(),
            "a2_validate": self.a2_validate.get(),
            "alpha2_dcd_selection": self.alpha2_dcd_selection.get(),
            
            # Dipole settings
            "skip_dipole": self.skip_dipole.get(),
            "dipole_calc_type": self.dipole_calc_type.get(),
            "individual_dipole_vars": [v.get() for v in self.individual_dipole_vars] if hasattr(self, 'individual_dipole_vars') else [],
            "collective_dipole_vars": [v.get() for v in self.collective_dipole_vars] if hasattr(self, 'collective_dipole_vars') else [],
            "dipole_parallel": self.dipole_parallel.get(),
            "dipole_validate": self.dipole_validate.get(),
            "dipole_chunk_processing": self.dipole_chunk_processing.get(),
            "individual_dipole_dcd_selection": self.individual_dipole_dcd_selection.get(),
            "collective_dipole_dcd_selection": self.collective_dipole_dcd_selection.get(),
            
            # Velocity extraction settings
            "skip_velocity": self.skip_velocity.get(),
            "velocity_vars": [v.get() for v in self.velocity_vars] if hasattr(self, 'velocity_vars') else [],
            "velocity_parallel": self.velocity_parallel.get(),
            "velocity_validate": self.velocity_validate.get(),
            "velocity_dcd_selection": self.velocity_dcd_selection.get(),
            
            # Common parameters
            "common_term_var": self.common_term_var.get(),
            "baseDir_var": self.baseDir_var.get(),
            "num_dcd_var": self.num_dcd_var.get(),
            "num_particles_var": self.num_particles_var.get(),
            "max_workers_var": self.max_workers_var.get(),
            "output_folder_var": self.output_folder_var.get(),
            "mainfile_var": self.mainfile_var.get(),
            "submitfile_var": self.submitfile_var.get(),
            "sbatch_vars": [v.get() for v in self.sbatch_vars] if hasattr(self, 'sbatch_vars') else []
        }
        
        # Save to a separate config file for calculations
        calc_config_path = os.path.expanduser("~/.pipeline_calculations_config.json")
        try:
            with open(calc_config_path, "w") as f:
                json.dump(data, f)
        except Exception as e:
            print(f"Warning: Could not save calculation config: {e}")

    def close_window(self):
        """Close the calculations window"""
        # Save configuration before closing
        self.save_calculation_config()
        self.destroy()
    
    def _on_mousewheel_calc(self, event):
        """Handle mouse wheel scrolling in calculations window"""
        if event.num == 4 or event.delta > 0:
            self.canvas.yview_scroll(-1, "units")
        elif event.num == 5 or event.delta < 0:
            self.canvas.yview_scroll(1, "units")
    
    def create_common_parameters_section(self, parent):
        """Create the common parameters section"""
        # Common Parameters frame
        common_frame = tk.LabelFrame(parent, text="Common Parameters", 
                                   font=("Arial", 12, "bold"), fg='#2c3e50', bg='#f0f0f0',
                                   relief='groove', borderwidth=2)
        common_frame.pack(fill="x", padx=5, pady=(0, 10))
        
        # Create a grid layout for common parameters
        row = 0
        
        # Base Directory
        tk.Label(common_frame, text="Base Directory:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=row, column=0, sticky="e", padx=5, pady=5)
        baseDir_entry = tk.Entry(common_frame, textvariable=self.baseDir_var, width=80, 
                               font=("Arial", 10), relief='solid', borderwidth=1)
        baseDir_entry.grid(row=row, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
        create_tooltip(baseDir_entry, "Base directory containing your simulation data")
        row += 1
        
        # Number of DCDs and Particles on same row
        tk.Label(common_frame, text="Number of DCDs:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=row, column=0, sticky="e", padx=5, pady=5)
        num_dcd_entry = tk.Entry(common_frame, textvariable=self.num_dcd_var, width=15, 
                               font=("Arial", 10), relief='solid', borderwidth=1)
        num_dcd_entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")
        create_tooltip(num_dcd_entry, "Total number of DCD trajectory files to process")
        
        tk.Label(common_frame, text="Number of Particles:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=row, column=3, sticky="e", padx=5, pady=5)
        num_particles_entry = tk.Entry(common_frame, textvariable=self.num_particles_var, width=15, 
                                     font=("Arial", 10), relief='solid', borderwidth=1)
        num_particles_entry.grid(row=row, column=4, padx=5, pady=5, sticky="w")
        create_tooltip(num_particles_entry, "Number of particles/molecules in your system")
        row += 1
        
        # Max Workers and Output Folder on same row
        tk.Label(common_frame, text="Max Workers:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=row, column=0, sticky="e", padx=5, pady=5)
        max_workers_entry = tk.Entry(common_frame, textvariable=self.max_workers_var, width=15, 
                                   font=("Arial", 10), relief='solid', borderwidth=1)
        max_workers_entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")
        create_tooltip(max_workers_entry, "Number of CPU cores for parallel processing")
        
        tk.Label(common_frame, text="Output Folder:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=row, column=3, sticky="e", padx=5, pady=5)
        output_folder_entry = tk.Entry(common_frame, textvariable=self.output_folder_var, width=20, 
                                     font=("Arial", 10), relief='solid', borderwidth=1)
        output_folder_entry.grid(row=row, column=4, padx=5, pady=5, sticky="w")
        create_tooltip(output_folder_entry, "Name of the output folder for generated scripts")
        row += 1
        
        # Common term expansion
        tk.Label(common_frame, text="Common Term:", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=row, column=0, sticky="e", padx=5, pady=5)
        common_term_entry = tk.Entry(common_frame, textvariable=self.common_term_var, width=80,
                                     font=("Arial", 10), relief='solid', borderwidth=1)
        common_term_entry.grid(row=row, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
        create_tooltip(common_term_entry, "Common path prefix for parameter expansion. Use * as placeholder in other fields (e.g., set this to '/path/to/sim' and use '*/input' in other fields)")
        row += 1
        
        # Sync button
        sync_frame = tk.Frame(common_frame, bg='#f0f0f0')
        sync_frame.grid(row=row, column=0, columnspan=5, pady=10)
        
        tk.Button(sync_frame, text="↓ Load from Main Window", 
                 command=self.load_from_parent, 
                 bg="#3498db", fg="white", font=("Arial", 10),
                 relief='raised', borderwidth=1, cursor='hand2').pack(side=tk.LEFT, padx=5)
        
        tk.Button(sync_frame, text="↑ Update Main Window", 
                 command=self.update_parent, 
                 bg="#e67e22", fg="white", font=("Arial", 10),
                 relief='raised', borderwidth=1, cursor='hand2').pack(side=tk.LEFT, padx=5)
        
        # Load parameters from parent on initialization
        self.load_from_parent()
    
    def load_from_parent(self):
        """Load common parameters from the main window"""
        try:
            self.baseDir_var.set(self.parent.baseDir_var.get())
            self.num_dcd_var.set(self.parent.num_dcd_var.get())
            self.num_particles_var.set(self.parent.num_particles_var.get())
            self.max_workers_var.set(self.parent.max_workers_var.get())
            self.output_folder_var.set(self.parent.output_folder_var.get())
            self.common_term_var.set(self.parent.common_term_var.get() if hasattr(self.parent, 'common_term_var') else "")
            
            # Also load SLURM and output parameters
            self.load_slurm_output_from_parent()
            
            print("✓ Loaded all parameters from main window")
        except Exception as e:
            print(f"Note: Could not load some parameters from main window: {e}")
    
    def update_parent(self):
        """Update main window parameters from calculations window"""
        try:
            self.parent.baseDir_var.set(self.baseDir_var.get())
            self.parent.num_dcd_var.set(self.num_dcd_var.get())
            self.parent.num_particles_var.set(self.num_particles_var.get())
            self.parent.max_workers_var.set(self.max_workers_var.get())
            self.parent.output_folder_var.set(self.output_folder_var.get())
            if hasattr(self.parent, 'common_term_var'):
                self.parent.common_term_var.set(self.common_term_var.get())
            
            # Also update SLURM and output parameters
            self.update_slurm_output_to_parent()
            
            print("✓ Updated all main window parameters")
        except Exception as e:
            print(f"Error updating main window parameters: {e}")
    
    def expand_common_term(self, value):
        """Replace asterisks (*) in value with the common term (same as main window)"""
        common_term = self.common_term_var.get().strip()
        if common_term and '*' in value:
            return value.replace('*', common_term)
        return value
    
    def create_slurm_output_section(self, parent):
        """Create the SLURM and Output Files section"""
        # Container for side-by-side layout
        params_container = tk.Frame(parent, bg='#f0f0f0')
        params_container.pack(fill="x", padx=5, pady=(0, 10))
        
        # SLURM parameters on the left
        slurm_frame = tk.LabelFrame(params_container, text="SLURM Submission Parameters",
                                   font=("Arial", 12, "bold"), fg='#2c3e50', bg='#f0f0f0',
                                   relief='groove', borderwidth=2)
        slurm_frame.grid(row=0, column=0, padx=(0, 5), pady=0, sticky="nw")
        
        slurm_labels = ["Nodes", "Partition", "QOS", "CPUs", "Tasks", "Memory (GB)", "Walltime", "Output prefix", "Email", "Module"]
        slurm_tooltips = [
            "Number of compute nodes to request (usually 1 for single-node jobs)",
            "SLURM partition/queue name (e.g., 'standard', 'gpu', 'high-mem')", 
            "Quality of Service level for job priority",
            "Number of CPU cores per node to request",
            "Number of tasks/processes (usually 1 for serial jobs)",
            "Memory allocation in GB (e.g., 160 for 160GB total memory)",
            "Maximum runtime (format: HH:MM:SS or DD-HH:MM:SS)",
            "Prefix for output log files",
            "Email address for job notifications",
            "Module to load before execution (optional, e.g., 'python/3.9', 'anaconda3')"
        ]
        
        for i, (lbl, tooltip) in enumerate(zip(slurm_labels, slurm_tooltips)):
            # Module field is optional, others are required
            required = "*" if lbl != "Module" else ""
            tk.Label(slurm_frame, text=f"{lbl}:{required}", font=("Arial", 10), 
                    bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=3)
            entry = tk.Entry(slurm_frame, textvariable=self.sbatch_vars[i], width=15, 
                           font=("Arial", 10), relief='solid', borderwidth=1)
            entry.grid(row=i, column=1, sticky="w", padx=5, pady=3)
            create_tooltip(entry, tooltip)
        
        # Output Files on the right
        output_frame = tk.LabelFrame(params_container, text="Output Files", 
                                   font=("Arial", 12, "bold"), fg='#2c3e50', bg='#f0f0f0',
                                   relief='groove', borderwidth=2)
        output_frame.grid(row=0, column=1, padx=(5, 0), pady=0, sticky="nw")
        
        # Output folder name
        tk.Label(output_frame, text="Output folder name:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=0, sticky="e", padx=5, pady=3)
        folder_entry = tk.Entry(output_frame, textvariable=self.output_folder_var, width=25, 
                              font=("Arial", 10), relief='solid', borderwidth=1)
        folder_entry.grid(row=0, column=1, padx=5, pady=3)
        create_tooltip(folder_entry, "Name of the folder to contain all generated files. Will be created in current directory or use Browse to select location.")
        tk.Button(output_frame, text="Browse...", command=self.browse_output_folder_calc,
                 font=("Arial", 9), bg='#ecf0f1', fg='#2c3e50', 
                 relief='raised', borderwidth=1, cursor='hand2').grid(
            row=0, column=2, padx=5, pady=3)
        
        # Main script file
        tk.Label(output_frame, text="Main script file:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=1, column=0, sticky="e", padx=5, pady=3)
        main_entry = tk.Entry(output_frame, textvariable=self.mainfile_var, width=25, 
                            font=("Arial", 10), relief='solid', borderwidth=1)
        main_entry.grid(row=1, column=1, padx=5, pady=3)
        create_tooltip(main_entry, "Python script filename (will be created in output folder)")
        
        # Submit script file
        tk.Label(output_frame, text="Submit script file:*", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=2, column=0, sticky="e", padx=5, pady=3)
        submit_entry = tk.Entry(output_frame, textvariable=self.submitfile_var, width=25, 
                              font=("Arial", 10), relief='solid', borderwidth=1)
        submit_entry.grid(row=2, column=1, padx=5, pady=3)
        create_tooltip(submit_entry, "SLURM batch script filename (.sh extension will be added automatically)")
        
        # Sync buttons for SLURM/Output parameters (spanning both frames)
        sync_slurm_frame = tk.Frame(params_container, bg='#f0f0f0')
        sync_slurm_frame.grid(row=1, column=0, columnspan=2, pady=10)
        
        tk.Button(sync_slurm_frame, text="↓ Load SLURM/Output from Main", 
                 command=self.load_slurm_output_from_parent, 
                 bg="#3498db", fg="white", font=("Arial", 10),
                 relief='raised', borderwidth=1, cursor='hand2').pack(side=tk.LEFT, padx=5)
        
        tk.Button(sync_slurm_frame, text="↑ Update SLURM/Output in Main", 
                 command=self.update_slurm_output_to_parent, 
                 bg="#e67e22", fg="white", font=("Arial", 10),
                 relief='raised', borderwidth=1, cursor='hand2').pack(side=tk.LEFT, padx=5)
        
        # Load SLURM and output parameters from parent on initialization
        self.load_slurm_output_from_parent()
    
    def browse_output_folder_calc(self):
        """Browse for output folder in calculations window"""
        directory = filedialog.askdirectory()
        if directory:
            # Use just the folder name, not the full path
            self.output_folder_var.set(os.path.basename(directory))
    
    def load_slurm_output_from_parent(self):
        """Load SLURM and output parameters from the main window"""
        try:
            # Load SLURM parameters
            if hasattr(self.parent, 'sbatch_vars'):
                for i, var in enumerate(self.sbatch_vars):
                    if i < len(self.parent.sbatch_vars):
                        var.set(self.parent.sbatch_vars[i].get())
            
            # Load output parameters
            self.output_folder_var.set(self.parent.output_folder_var.get())
            self.mainfile_var.set(self.parent.mainfile_var.get())
            self.submitfile_var.set(self.parent.submitfile_var.get())
        except Exception as e:
            print(f"Warning: Could not load SLURM/output parameters: {e}")

    def update_slurm_output_to_parent(self):
        """Update SLURM and output parameters in the main window"""
        try:
            # Update SLURM parameters
            if hasattr(self.parent, 'sbatch_vars'):
                for i, var in enumerate(self.sbatch_vars):
                    if i < len(self.parent.sbatch_vars):
                        self.parent.sbatch_vars[i].set(var.get())
            
            # Update output parameters
            self.parent.output_folder_var.set(self.output_folder_var.get())
            self.parent.mainfile_var.set(self.mainfile_var.get())
            self.parent.submitfile_var.set(self.submitfile_var.get())
        except Exception as e:
            print(f"Warning: Could not update SLURM/output parameters: {e}")
    
    def create_alpha2_section(self, parent):
        """Create the Alpha2/MSD calculation section"""
        # Title frame with skip checkbox
        title_frame = tk.Frame(parent, bg='#f0f0f0')
        title_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Label(title_frame, text="Non-Gaussian Parameter alpha_2(t) and MSD Calculation", 
                font=("Arial", 13, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(side="left")
        skip_check = tk.Checkbutton(title_frame, text="Skip", variable=self.skip_alpha2,
                                   command=lambda: self.toggle_frame(alpha2_content, self.skip_alpha2.get()),
                                   font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50', 
                                   selectcolor='#e8f4fd', activebackground='#f0f0f0')
        skip_check.pack(side="right", padx=10)
        create_tooltip(skip_check, "Check to skip α₂(t)/αₓₓ calculation. Leave unchecked to generate calculation scripts for this step.")
        
        # Content frame
        alpha2_content = tk.LabelFrame(parent, text="", relief='groove', borderwidth=1, bg='#f0f0f0')
        alpha2_content.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Calculation type selection
        calc_type_frame = tk.Frame(alpha2_content, bg='#f0f0f0')
        calc_type_frame.pack(fill="x", padx=10, pady=10)
        
        tk.Label(calc_type_frame, text="Calculation Type:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').pack(side="left")
        
        radio_frame = tk.Frame(calc_type_frame, bg='#f0f0f0')
        radio_frame.pack(side="left", padx=20)
        
        # Unicode fonts for Greek characters
        unicode_fonts = [("Segoe UI", 10), ("Arial Unicode MS", 10), ("DejaVu Sans", 10), 
                        ("Liberation Sans", 10), ("Arial", 10), ("TkDefaultFont", 10)]
        
        unicode_font = ("TkDefaultFont", 10)
        for font in unicode_fonts:
            try:
                test_widget = tk.Label(radio_frame, text="α", font=font)
                unicode_font = font
                test_widget.destroy()
                break
            except:
                continue
        
        alpha2_radio = tk.Radiobutton(radio_frame, text="α₂(t) and MSD", variable=self.calc_type_var, 
                      value="alpha2_msd", font=unicode_font, bg='#f0f0f0', fg='#2c3e50',
                      selectcolor='#e8f4fd', activebackground='#f0f0f0')
        alpha2_radio.pack(side=tk.LEFT, padx=5)
        create_tooltip(alpha2_radio, "Calculate standard non-Gaussian parameter alpha2(t) = 3<Dr^4>/(5<Dr^2>^2) - 1")
        
        alpha_anisotropy_radio = tk.Radiobutton(radio_frame, text="Alpha Anisotropy", variable=self.calc_type_var, 
                      value="alpha_xz", font=unicode_font, bg='#f0f0f0', fg='#2c3e50',
                      selectcolor='#e8f4fd', activebackground='#f0f0f0')
        alpha_anisotropy_radio.pack(side=tk.LEFT, padx=5)
        create_tooltip(alpha_anisotropy_radio, "Calculate directional correlation parameter alpha_xz(t) = <Dx^2*Dz^2>/(<Dx^2>*<Dz^2>) - 1")
        
        # Parameters grid
        params_frame = tk.Frame(alpha2_content, bg='#f0f0f0')
        params_frame.pack(fill="x", padx=10, pady=10)
        
        labels = ["Input File Pattern", "Output File Pattern", "Min frames"]
        tooltips = [
            "Path pattern for input COM files. Use * for common term, {i} for file index. Example: anlz/NVT_*/com_data/com_{i}.dat",
            "Path pattern for output files. Use * for common term, {i} for file index. Example: anlz/NVT_*/analysis/alpha2_{i}.dat",
            "Minimum number of time frames required per trajectory file"
        ]
        
        for i, (lbl, tooltip) in enumerate(zip(labels, tooltips)):
            tk.Label(params_frame, text=f"{lbl}:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
                row=i, column=0, sticky="e", padx=5, pady=5)
            v = tk.StringVar()
            entry = tk.Entry(params_frame, textvariable=v, width=80, font=("Arial", 10), relief='solid', borderwidth=1)
            entry.grid(row=i, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
            create_tooltip(entry, tooltip)
            self.a2_vars.append(v)
        
        # DCD selection for alpha2
        tk.Label(params_frame, text="DCD Selection (optional):", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=len(labels), column=0, sticky="e", padx=5, pady=5)
        alpha2_dcd_entry = tk.Entry(params_frame, textvariable=self.alpha2_dcd_selection, width=80, font=("Arial", 10), relief='solid', borderwidth=1)
        alpha2_dcd_entry.grid(row=len(labels), column=1, columnspan=2, padx=5, pady=5, sticky="ew")
        create_tooltip(alpha2_dcd_entry, "Optional DCD indices to process (e.g., '4-10', '1,3,5', 'range(10,20)'). Leave empty to process all DCDs.")
        
        # Advanced options
        options_frame = tk.Frame(alpha2_content, bg='#f0f0f0')
        options_frame.pack(fill="x", padx=10, pady=10)
        
        tk.Label(options_frame, text="Chunk Processing:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=0, column=0, sticky="e", padx=5, pady=5)
        chunk_check = tk.Checkbutton(options_frame, variable=self.a2_chunk_processing, font=("Arial", 10),
                                    bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        chunk_check.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(chunk_check, "Process data in chunks for better memory efficiency with large datasets")
        
        tk.Label(options_frame, text="Validate Data:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=1, column=0, sticky="e", padx=5, pady=5)
        validate_check = tk.Checkbutton(options_frame, variable=self.a2_validate, font=("Arial", 10),
                                       bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        validate_check.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(validate_check, "Perform data validation checks during processing")
        
        self.alpha2_content_frame = alpha2_content
        self.toggle_frame(alpha2_content, self.skip_alpha2.get())
    
    def create_dipole_section(self, parent):
        """Create the dipole calculation section"""
        # Title frame with skip checkbox
        title_frame = tk.Frame(parent, bg='#f0f0f0')
        title_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Label(title_frame, text="Dipole Moment Calculation", 
                font=("Arial", 13, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(side="left")
        skip_check = tk.Checkbutton(title_frame, text="Skip", variable=self.skip_dipole,
                                   command=lambda: self.toggle_frame(dipole_content, self.skip_dipole.get()),
                                   font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50', 
                                   selectcolor='#e8f4fd', activebackground='#f0f0f0')
        skip_check.pack(side="right", padx=10)
        create_tooltip(skip_check, "Check to skip dipole calculation (default: checked). Uncheck only if you want to generate dipole calculation scripts.")
        
        # Content frame
        dipole_content = tk.LabelFrame(parent, text="", relief='groove', borderwidth=1, bg='#f0f0f0')
        dipole_content.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Dipole calculation type selection
        dipole_type_frame = tk.Frame(dipole_content, bg='#f0f0f0')
        dipole_type_frame.pack(fill="x", padx=10, pady=10)
        
        tk.Label(dipole_type_frame, text="Calculation Method:*", font=("Arial", 10, "bold"), bg='#f0f0f0', fg='#2c3e50').pack(side="left")
        
        # Initialize dipole calculation type variable
        self.dipole_calc_type = tk.StringVar(value="individual")
        
        radio_frame = tk.Frame(dipole_type_frame, bg='#f0f0f0')
        radio_frame.pack(side="left", padx=20)
        
        individual_radio = tk.Radiobutton(radio_frame, text="Individual Dipole Moments", 
                                        variable=self.dipole_calc_type, value="individual",
                                        command=self.update_dipole_fields,
                                        font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50',
                                        selectcolor='#e8f4fd', activebackground='#f0f0f0')
        individual_radio.pack(side=tk.LEFT, padx=5)
        create_tooltip(individual_radio, "Calculate individual dipole moments for each molecule using Python implementation (requires charges and atoms per particle)")
        
        collective_radio = tk.Radiobutton(radio_frame, text="Collective Dipole Moment", 
                                        variable=self.dipole_calc_type, value="collective",
                                        command=self.update_dipole_fields,
                                        font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50',
                                        selectcolor='#e8f4fd', activebackground='#f0f0f0')
        collective_radio.pack(side=tk.LEFT, padx=5)
        create_tooltip(collective_radio, "Calculate collective dipole moment of selected atoms using VMD implementation (requires VMD path and target selection)")
        
        # Create frames for different parameter sets
        self.individual_params_frame = tk.Frame(dipole_content, bg='#f0f0f0')
        self.collective_params_frame = tk.Frame(dipole_content, bg='#f0f0f0')
        
        # Individual dipole parameters
        individual_labels = [
            "Coordinates Pattern", "COM Pattern", "Dipole Vectors Pattern", "Dipole Magnitudes Pattern",
            "Atomic charges", "Atoms per molecule", "Molecules to process (optional)", "Stride (optional)"
        ]
        individual_tooltips = [
            "Path pattern for coordinate files. Use * for common term, {i} for file index. Example: anlz/NVT_*/coords/xyz_{i}.dat",
            "Path pattern for COM files. Use * for common term, {i} for file index. Example: anlz/NVT_*/com_data/com_{i}.dat",
            "Output pattern for dipole vector files. Use * for common term, {i} for file index. Example: anlz/NVT_*/dipole/vectors_{i}.dat",
            "Output pattern for dipole magnitude files. Use * for common term, {i} for file index. Example: anlz/NVT_*/dipole/magnitudes_{i}.dat",
            "Comma-separated atomic charges (e.g., '-0.8,0.4,0.4' for water O H H)",
            "Number of atoms per molecule (e.g., 3 for water)",
            "Number of molecules to process (leave empty to process all)",
            "Processing stride (default: 1, process every frame)"
        ]
        
        self.individual_dipole_vars = []
        for i, (lbl, tooltip) in enumerate(zip(individual_labels, individual_tooltips)):
            required = "*" if "optional" not in lbl.lower() else ""
            tk.Label(self.individual_params_frame, text=f"{lbl}:{required}", font=("Arial", 10), 
                    bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=3)
            v = tk.StringVar()
            entry = tk.Entry(self.individual_params_frame, textvariable=v, width=60, 
                           font=("Arial", 10), relief='solid', borderwidth=1)
            entry.grid(row=i, column=1, columnspan=2, padx=5, pady=3, sticky="ew")
            create_tooltip(entry, tooltip)
            self.individual_dipole_vars.append(v)
        
        # DCD selection for individual dipole
        tk.Label(self.individual_params_frame, text="DCD Selection (optional):", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=len(individual_labels), column=0, sticky="e", padx=5, pady=3)
        individual_dipole_dcd_entry = tk.Entry(self.individual_params_frame, textvariable=self.individual_dipole_dcd_selection, width=60, 
                                             font=("Arial", 10), relief='solid', borderwidth=1)
        individual_dipole_dcd_entry.grid(row=len(individual_labels), column=1, columnspan=2, padx=5, pady=3, sticky="ew")
        create_tooltip(individual_dipole_dcd_entry, "Optional DCD indices to process (e.g., '4-10', '1,3,5', 'range(10,20)'). Leave empty to process all DCDs.")
        
        # Collective dipole parameters  
        collective_labels = [
            "Trajectory Pattern", "PSF Pattern", "DCD Pattern", "Dipole Output Pattern", "Target Selection", "VMD Path"
        ]
        collective_tooltips = [
            "Base path pattern for trajectory files. Use * for common term. Example: anlz/NVT_*/trajectories",
            "Path pattern for PSF files. Use * for common term, {i} for file index. Example: anlz/NVT_*/run_{i}/system.psf",
            "Path pattern for DCD files. Use * for common term, {i} for file index. Example: anlz/NVT_*/run_{i}/traj.dcd",
            "Output pattern for dipole result files. Use * for common term, {i} for file index. Example: anlz/NVT_*/dipole/collective_{i}.dat",
            "VMD atom selection for dipole calculation (e.g., 'water' or 'residue 1 to 100')",
            "Full path to VMD executable"
        ]
        
        self.collective_dipole_vars = []
        for i, (lbl, tooltip) in enumerate(zip(collective_labels, collective_tooltips)):
            tk.Label(self.collective_params_frame, text=f"{lbl}:*", font=("Arial", 10), 
                    bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=3)
            v = tk.StringVar()
            entry = tk.Entry(self.collective_params_frame, textvariable=v, width=60, 
                           font=("Arial", 10), relief='solid', borderwidth=1)
            entry.grid(row=i, column=1, columnspan=2, padx=5, pady=3, sticky="ew")
            create_tooltip(entry, tooltip)
            self.collective_dipole_vars.append(v)
        
        # DCD selection for collective dipole
        tk.Label(self.collective_params_frame, text="DCD Selection (optional):", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=len(collective_labels), column=0, sticky="e", padx=5, pady=3)
        collective_dipole_dcd_entry = tk.Entry(self.collective_params_frame, textvariable=self.collective_dipole_dcd_selection, width=60, 
                                             font=("Arial", 10), relief='solid', borderwidth=1)
        collective_dipole_dcd_entry.grid(row=len(collective_labels), column=1, columnspan=2, padx=5, pady=3, sticky="ew")
        create_tooltip(collective_dipole_dcd_entry, "Optional DCD indices to process (e.g., '4-10', '1,3,5', 'range(10,20)'). Leave empty to process all DCDs.")
        
        # Keep reference to dipole_vars for backward compatibility
        self.dipole_vars = self.individual_dipole_vars
        
        # Advanced options (shared between both calculation types)
        options_frame = tk.Frame(dipole_content, bg='#f0f0f0')
        options_frame.pack(fill="x", padx=10, pady=10)
        
        tk.Label(options_frame, text="Use Parallel Processing:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=0, column=0, sticky="e", padx=5, pady=5)
        parallel_check = tk.Checkbutton(options_frame, variable=self.dipole_parallel, font=("Arial", 10),
                                       bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        parallel_check.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(parallel_check, "Enable parallel processing for faster dipole calculations")
        
        tk.Label(options_frame, text="Chunk Processing:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=1, column=0, sticky="e", padx=5, pady=5)
        chunk_check = tk.Checkbutton(options_frame, variable=self.dipole_chunk_processing, font=("Arial", 10),
                                    bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        chunk_check.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(chunk_check, "Process data in chunks for memory efficiency")
        
        tk.Label(options_frame, text="Validate Data:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=2, column=0, sticky="e", padx=5, pady=5)
        validate_check = tk.Checkbutton(options_frame, variable=self.dipole_validate, font=("Arial", 10),
                                       bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        validate_check.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(validate_check, "Perform data validation checks during processing")
        
        # Store references to content frame and options
        self.dipole_content_frame = dipole_content
        self.dipole_options_frame = options_frame
        
        # Initialize the display based on default selection
        self.update_dipole_fields()
        
        self.toggle_frame(dipole_content, self.skip_dipole.get())
    
    def create_velocity_section(self, parent):
        """Create the velocity extraction section"""
        # Title frame with skip checkbox
        title_frame = tk.Frame(parent, bg='#f0f0f0')
        title_frame.pack(fill="x", padx=10, pady=5)
        
        tk.Label(title_frame, text="Velocity Extraction from VELDCD Files", 
                font=("Arial", 13, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(side="left")
        skip_check = tk.Checkbutton(title_frame, text="Skip", variable=self.skip_velocity,
                                   command=lambda: self.toggle_frame(velocity_content, self.skip_velocity.get()),
                                   font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50', 
                                   selectcolor='#e8f4fd', activebackground='#f0f0f0')
        skip_check.pack(side="right", padx=10)
        create_tooltip(skip_check, "Check to skip velocity extraction (default: checked). Uncheck to generate velocity extraction scripts for NAMD VELDCD files.")
        
        # Content frame
        velocity_content = tk.LabelFrame(parent, text="", relief='groove', borderwidth=1, bg='#f0f0f0')
        velocity_content.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Velocity extraction parameters
        params_frame = tk.Frame(velocity_content, bg='#f0f0f0')
        params_frame.pack(fill="x", padx=10, pady=10)
        
        # Parameter labels and tooltips
        velocity_labels = [
            "PSF Pattern", 
            "VELDCD Pattern", 
            "Output Pattern",
            "Number of Molecules",
            "VMD Executable Path"
        ]
        
        velocity_tooltips = [
            "Pattern for PSF files (e.g., 'path/to/system.psf' or 'run_{i}/system.psf')",
            "Pattern for VELDCD files (e.g., 'path/to/traj.veldcd' or 'run_{i}/traj.veldcd')", 
            "Pattern for output velocity files (e.g., 'velocities/velCOM_{i}.dat')",
            "Number of molecules/residues for COM velocity calculation",
            "Full path to VMD executable (e.g., '/usr/local/bin/vmd')"
        ]
        
        self.velocity_vars = []
        for i, (lbl, tooltip) in enumerate(zip(velocity_labels, velocity_tooltips)):
            required = "*" if i < 5 else ""  # All fields are required
            tk.Label(params_frame, text=f"{lbl}:{required}", font=("Arial", 10), 
                    bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=(5,10), pady=3)
            
            v = tk.StringVar()
            entry = tk.Entry(params_frame, textvariable=v, width=40, 
                           font=("Arial", 10), relief='solid', borderwidth=1)
            entry.grid(row=i, column=1, columnspan=2, padx=5, pady=3, sticky="ew")
            create_tooltip(entry, tooltip)
            self.velocity_vars.append(v)
        
        # DCD selection for velocity
        tk.Label(params_frame, text="DCD Selection (optional):", font=("Arial", 10), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=len(velocity_labels), column=0, sticky="e", padx=(5,10), pady=3)
        velocity_dcd_entry = tk.Entry(params_frame, textvariable=self.velocity_dcd_selection, width=40,
                                     font=("Arial", 10), relief='solid', borderwidth=1)
        velocity_dcd_entry.grid(row=len(velocity_labels), column=1, columnspan=2, padx=5, pady=3, sticky="ew")
        create_tooltip(velocity_dcd_entry, "Optional DCD indices to process (e.g., '4-10', '1,3,5', 'range(10,20)'). Leave empty to process all DCDs.")
        
        # Browse button for VMD executable
        def browse_vmd_velocity():
            f = filedialog.askopenfilename(filetypes=[("Executable","*"), ("All files", "*")])
            if f:
                self.velocity_vars[4].set(f)  # VMD path is at index 4
        
        tk.Button(params_frame, text="Browse VMD", command=browse_vmd_velocity, 
                 bg="#3498db", fg="white", font=("Arial", 9), padx=10, pady=2,
                 relief='raised', borderwidth=1, cursor='hand2').grid(row=4, column=3, padx=5, pady=3)
        
        # Advanced options
        options_frame = tk.Frame(velocity_content, bg='#f0f0f0')
        options_frame.pack(fill="x", padx=10, pady=10)
        
        tk.Label(options_frame, text="Advanced Options:", font=("Arial", 11, "bold"), 
                bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=0, columnspan=3, sticky="w", pady=(5,10))
        
        # Parallel processing checkbox
        tk.Label(options_frame, text="Parallel Processing:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=1, column=0, sticky="e", padx=5, pady=5)
        parallel_check = tk.Checkbutton(options_frame, variable=self.velocity_parallel, font=("Arial", 10),
                                       bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        parallel_check.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(parallel_check, "Enable parallel processing for faster velocity extraction")
        
        # Data validation checkbox
        tk.Label(options_frame, text="Data Validation:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(
            row=2, column=0, sticky="e", padx=5, pady=5)
        validate_check = tk.Checkbutton(options_frame, variable=self.velocity_validate, font=("Arial", 10),
                                       bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        validate_check.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(validate_check, "Perform data validation checks during processing")
        
        # Configure grid weights for proper resizing
        params_frame.columnconfigure(1, weight=1)
        options_frame.columnconfigure(1, weight=1)
        
        # Store references to content frame
        self.velocity_content_frame = velocity_content
        
        self.toggle_frame(velocity_content, self.skip_velocity.get())
    
    def toggle_frame(self, frame, skip):
        """Disable or enable all entries & buttons in a frame"""
        state = "disabled" if skip else "normal"
        for widget in frame.winfo_children():
            if isinstance(widget, (tk.Entry, tk.Button)):
                widget.configure(state=state)
            elif hasattr(widget, 'winfo_children'):
                for child in widget.winfo_children():
                    if isinstance(child, (tk.Entry, tk.Button)):
                        child.configure(state=state)
    
    def center_on_parent(self):
        """Center this window on the parent window"""
        self.update_idletasks()
        parent_x = self.parent.winfo_x()
        parent_y = self.parent.winfo_y()
        parent_width = self.parent.winfo_width()
        parent_height = self.parent.winfo_height()
        
        x = parent_x + (parent_width - self.winfo_width()) // 2
        y = parent_y + (parent_height - self.winfo_height()) // 2
        
        self.geometry(f"+{x}+{y}")
    
    def load_calculation_config(self):
        """Load calculation configuration from saved file"""
        import json
        
        calc_config_path = os.path.expanduser("~/.pipeline_calculations_config.json")
        if not os.path.exists(calc_config_path):
            return
        
        try:
            with open(calc_config_path) as f:
                data = json.load(f)
            
            # Load Alpha2/MSD settings
            if "skip_alpha2" in data:
                self.skip_alpha2.set(data["skip_alpha2"])
            if "calc_type_var" in data:
                self.calc_type_var.set(data["calc_type_var"])
            if "a2_chunk_processing" in data:
                self.a2_chunk_processing.set(data["a2_chunk_processing"])
            if "a2_validate" in data:
                self.a2_validate.set(data["a2_validate"])
            if "alpha2_dcd_selection" in data:
                self.alpha2_dcd_selection.set(data["alpha2_dcd_selection"])
            
            # Load a2_vars if they exist
            if "a2_vars" in data and hasattr(self, 'a2_vars'):
                for i, value in enumerate(data["a2_vars"]):
                    if i < len(self.a2_vars):
                        self.a2_vars[i].set(value)
            
            # Load Dipole settings
            if "skip_dipole" in data:
                self.skip_dipole.set(data["skip_dipole"])
            if "dipole_calc_type" in data:
                self.dipole_calc_type.set(data["dipole_calc_type"])
            if "dipole_parallel" in data:
                self.dipole_parallel.set(data["dipole_parallel"])
            if "dipole_validate" in data:
                self.dipole_validate.set(data["dipole_validate"])
            if "dipole_chunk_processing" in data:
                self.dipole_chunk_processing.set(data["dipole_chunk_processing"])
            if "individual_dipole_dcd_selection" in data:
                self.individual_dipole_dcd_selection.set(data["individual_dipole_dcd_selection"])
            if "collective_dipole_dcd_selection" in data:
                self.collective_dipole_dcd_selection.set(data["collective_dipole_dcd_selection"])
            
            # Load dipole vars if they exist
            if "individual_dipole_vars" in data and hasattr(self, 'individual_dipole_vars'):
                for i, value in enumerate(data["individual_dipole_vars"]):
                    if i < len(self.individual_dipole_vars):
                        self.individual_dipole_vars[i].set(value)
            
            if "collective_dipole_vars" in data and hasattr(self, 'collective_dipole_vars'):
                for i, value in enumerate(data["collective_dipole_vars"]):
                    if i < len(self.collective_dipole_vars):
                        self.collective_dipole_vars[i].set(value)
            
            # Load Velocity extraction settings
            if "skip_velocity" in data:
                self.skip_velocity.set(data["skip_velocity"])
            if "velocity_parallel" in data:
                self.velocity_parallel.set(data["velocity_parallel"])
            if "velocity_validate" in data:
                self.velocity_validate.set(data["velocity_validate"])
            if "velocity_dcd_selection" in data:
                self.velocity_dcd_selection.set(data["velocity_dcd_selection"])
            
            # Load velocity vars if they exist
            if "velocity_vars" in data and hasattr(self, 'velocity_vars'):
                for i, value in enumerate(data["velocity_vars"]):
                    if i < len(self.velocity_vars):
                        self.velocity_vars[i].set(value)
            
            # Load Common parameters
            if "common_term_var" in data:
                self.common_term_var.set(data["common_term_var"])
            if "baseDir_var" in data:
                self.baseDir_var.set(data["baseDir_var"])
            if "num_dcd_var" in data:
                self.num_dcd_var.set(data["num_dcd_var"])
            if "num_particles_var" in data:
                self.num_particles_var.set(data["num_particles_var"])
            if "max_workers_var" in data:
                self.max_workers_var.set(data["max_workers_var"])
            if "output_folder_var" in data:
                self.output_folder_var.set(data["output_folder_var"])
            if "mainfile_var" in data:
                self.mainfile_var.set(data["mainfile_var"])
            if "submitfile_var" in data:
                self.submitfile_var.set(data["submitfile_var"])
            
            # Load sbatch_vars if they exist
            if "sbatch_vars" in data and hasattr(self, 'sbatch_vars'):
                for i, value in enumerate(data["sbatch_vars"]):
                    if i < len(self.sbatch_vars):
                        self.sbatch_vars[i].set(value)
                        
        except Exception as e:
            print(f"Warning: Could not load calculation config: {e}")
    
    def generate_calculation_files(self):
        """Generate separate script files for the selected calculations"""
        try:
            # Get common parameters from local variables (can be different from main window)
            baseDir = self.baseDir_var.get().strip()
            if not baseDir:
                raise ValueError("Base Directory is required")
            
            num_dcd = int(self.num_dcd_var.get())
            num_particles = int(self.num_particles_var.get())
            max_workers = int(self.max_workers_var.get())
            
            # Get output folder from local variables
            output_folder = self.output_folder_var.get().strip()
            if not output_folder:
                raise ValueError("Output folder name is required")
            
            # Use current directory for output path
            output_path = os.path.join(os.getcwd(), output_folder)
            
            # Create the output directory if it doesn't exist
            os.makedirs(output_path, exist_ok=True)
            
            # Copy main_functions directory from current working directory to output folder
            current_dir = os.getcwd()
            main_functions_src = os.path.join(current_dir, "main_functions")
            main_functions_dst = os.path.join(output_path, "main_functions")
            
            if os.path.exists(main_functions_src) and not os.path.exists(main_functions_dst):
                try:
                    # Copy the entire main_functions folder
                    self._copy_compiled_functions(main_functions_src, main_functions_dst)
                    print(f"✓ Copied main_functions folder from {main_functions_src} to {main_functions_dst}")
                except Exception as e:
                    print(f"Warning: Could not copy main_functions: {e}")
            elif not os.path.exists(main_functions_src):
                print(f"Warning: main_functions directory not found in current directory: {current_dir}")
            else:
                print(f"✓ main_functions already exists in {main_functions_dst}")
            
            files_created = []
            
            # Generate Alpha2/MSD calculation script
            if not self.skip_alpha2.get():
                alpha2_file = self.generate_alpha2_script(baseDir, num_dcd, num_particles, max_workers, output_path)
                files_created.append(alpha2_file)
            
            # Generate Dipole calculation script
            if not self.skip_dipole.get():
                dipole_file = self.generate_dipole_script(baseDir, num_dcd, num_particles, max_workers, output_path)
                files_created.append(dipole_file)
            
            # Generate Velocity extraction script
            if not self.skip_velocity.get():
                velocity_file = self.generate_velocity_script(baseDir, num_dcd, num_particles, max_workers, output_path)
                files_created.append(velocity_file)
            
            # Generate SLURM submission script if we have calculations to run
            if files_created:
                slurm_file = self.generate_slurm_script(files_created, output_path)
                files_created.append(slurm_file)
            
            if not files_created:
                alpha2_skipped = self.skip_alpha2.get()
                dipole_skipped = self.skip_dipole.get()
                
                if alpha2_skipped and dipole_skipped:
                    messagebox.showwarning("No calculations selected", 
                                         "Please enable at least one calculation type by unchecking the 'Skip' checkbox:\n\n"
                                         "• Uncheck 'Skip' for α₂(t)/αₓₓ calculation\n"
                                         "• Uncheck 'Skip' for Dipole calculation\n\n"
                                         "At least one calculation must be enabled.")
                else:
                    messagebox.showwarning("Generation failed", 
                                         "No calculation scripts were generated. Please check that all required parameters are provided for the enabled calculations.")
                return
            
            message = f"Successfully generated calculation scripts:\n\n"
            for file in files_created:
                message += f"• {os.path.basename(file)}\n"
            message += f"\nLocation: {output_path}\n\n"
            message += "These scripts can be run independently or as part of your pipeline."
            
            # Save configuration
            self.save_calculation_config()
            messagebox.showinfo("Success", message)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate calculation files:\n{e}")
    
    def generate_slurm_script(self, calculation_files, output_path):
        """Generate SLURM submission script for the calculations"""
        # Get SLURM parameters
        sb = [v.get().strip() for v in self.sbatch_vars]
        nodes, partition, qos, cpus, tasks, memory, walltime, output_prefix, email, module = sb
        
        # Get output file names and expand common term
        main_file = self.mainfile_var.get().strip()
        submit_file = self.submitfile_var.get().strip()
        
        # Apply common term expansion to file names
        main_file = self.expand_common_term(main_file) if main_file else ""
        submit_file = self.expand_common_term(submit_file) if submit_file else "submit_calculations"
        
        # Ensure .sh extension for submit file
        if submit_file and not submit_file.endswith('.sh'):
            submit_file += '.sh'
        
        # Auto-determine SLURM output prefix based on common term vs output folder
        common_term = self.common_term_var.get().strip()
        output_folder = self.output_folder_var.get().strip()
        
        if common_term:
            slurm_output_name = f"sbatch_{common_term}"
        elif output_folder:
            slurm_output_name = f"sbatch_{output_folder}"
        else:
            # Fallback to submit file name (without .sh) 
            slurm_output_name = submit_file[:-3] if submit_file.endswith('.sh') else submit_file
        
        script_path = os.path.join(output_path, submit_file)
        
        # Create SLURM submission script content
        script_content = f'''#!/bin/bash
# SLURM Submission Script for Analysis Calculations
# Generated by MD Analysis Pipeline GUI v2.0

'''
        
        # Add SLURM directives
        if nodes:
            script_content += f"#SBATCH -N {nodes}\n"
        if partition:
            script_content += f"#SBATCH -p {partition}\n"
        if qos:
            script_content += f"#SBATCH -q {qos}\n"
        if cpus:
            script_content += f"#SBATCH -c {cpus}\n"
        if tasks:
            script_content += f"#SBATCH -n {tasks}\n"
        if memory:
            script_content += f"#SBATCH --mem={memory}G\n"
        if walltime:
            script_content += f"#SBATCH -t {walltime}\n"
        # Use submit script filename for SLURM output (not the output_prefix parameter)
        script_content += f"#SBATCH -o {slurm_output_name}.log\n"
        if email:
            script_content += f"#SBATCH --mail-type=ALL\n"
            script_content += f"#SBATCH --mail-user={email}\n"
        
        script_content += f"#SBATCH --export=NONE\n\n"
        
        # Add module loading section
        if module:
            script_content += f'''# Load required module
module load {module}

########################################################
# Analysis Calculations Execution
########################################################

echo 'Starting analysis calculations...'
echo 'Job started at:' $(date)

'''
        else:
            script_content += '''# No module specified
# Add module load commands here if needed, e.g.:
# module load python/3.8
# module load conda

########################################################
# Analysis Calculations Execution
########################################################

echo 'Starting analysis calculations...'
echo 'Job started at:' $(date)

'''
        
        # Add execution commands for each calculation script
        for calc_file in calculation_files:
            if calc_file.endswith('.py'):  # Only run Python scripts
                script_name = os.path.basename(calc_file)
                script_content += f"echo 'Running {script_name}...'\n"
                script_content += f"python {script_name}\n"
                script_content += f"echo 'Completed {script_name}'\n\n"
        
        script_content += '''echo 'All calculations completed at:' $(date)
'''
        
        # Write the script
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Make it executable
        import stat
        os.chmod(script_path, os.stat(script_path).st_mode | stat.S_IEXEC)
        
        return script_path
    
    def generate_alpha2_script(self, baseDir, num_dcd, num_particles, max_workers, output_path):
        """Generate Alpha2/MSD calculation script"""
        in4, out4, minf = [v.get().strip() for v in self.a2_vars]
        # Get common term for script generation
        common_term = self.common_term_var.get().strip()
        # Apply common term expansion to alpha2/MSD variables
        in4, out4, minf = [self.expand_common_term(var) for var in [in4, out4, minf]]
        
        calc_type = self.calc_type_var.get()
        chunk_processing = self.a2_chunk_processing.get()
        validate_data = self.a2_validate.get()
        
        # Get DCD selection
        alpha2_dcd_selection = self.alpha2_dcd_selection.get().strip()
        
        if not all([in4, out4, minf]):
            missing_params = []
            if not in4: missing_params.append("Input File Pattern")
            if not out4: missing_params.append("Output File Pattern") 
            if not minf: missing_params.append("Min frames")
            
            calc_name = "α₂(t) and MSD" if calc_type == "alpha2_msd" else "αₓₓ(t)"
            raise ValueError(f"Missing required {calc_name} parameters: {', '.join(missing_params)}")
        
        # Use user-specified main file name with common term expansion (exact name, no suffixes)
        main_file = self.mainfile_var.get().strip()
        if main_file:
            script_name = self.expand_common_term(main_file)
            # Ensure .py extension
            if not script_name.endswith('.py'):
                script_name += '.py'
        else:
            script_name = "alpha2_calculation.py" if calc_type == "alpha2_msd" else "alpha_xz_calculation.py"
        
        script_path = os.path.join(output_path, script_name)
        
        function_name = "a2_MSD" if calc_type == "alpha2_msd" else "alpha_xz"
        module_name = "alpha2_MSD" if calc_type == "alpha2_msd" else "axz"
        calculation_type = "α₂(t) and MSD" if calc_type == "alpha2_msd" else "α_xz(t)"
        
        script_content = f'''#!/usr/bin/env python3
"""
{calculation_type} Calculation Script
Generated by MD Analysis Pipeline GUI v2.0
"""

import sys
import time
import os

# Add main_functions to Python path for imports
main_functions_path = os.path.join(os.getcwd(), 'main_functions')
if main_functions_path not in sys.path:
    sys.path.insert(0, main_functions_path)

def main():
    print('='*60)
    print('{calculation_type.upper()} CALCULATION')
    print('='*60)
    start_time = time.time()
    
    try:
        # Import Python function directly
        import main_functions.{module_name}
        {function_name} = main_functions.{module_name}.{function_name}

        # Parse DCD selection if provided
        dcd_selection = "{alpha2_dcd_selection}"
        if dcd_selection:
            # Parse DCD selection (e.g., "4-10", "1,3,5", "range(10,20)")
            try:
                if 'range(' in dcd_selection:
                    # Handle range() expressions
                    dcd_indices = list(eval(dcd_selection))
                elif '-' in dcd_selection and ',' in dcd_selection:
                    # Handle mixed format like "1,3,5-8"
                    dcd_indices = []
                    for part in dcd_selection.split(','):
                        if '-' in part:
                            start, end = map(int, part.split('-'))
                            dcd_indices.extend(range(start, end + 1))
                        else:
                            dcd_indices.append(int(part))
                elif '-' in dcd_selection:
                    # Handle range format like "4-10"
                    start, end = map(int, dcd_selection.split('-'))
                    dcd_indices = list(range(start, end + 1))
                else:
                    # Handle comma-separated format like "1,3,5"
                    dcd_indices = [int(x.strip()) for x in dcd_selection.split(',')]
                
                print(f'Using DCD selection: {{dcd_indices}}')
                actual_num_dcd = len(dcd_indices)
            except Exception as e:
                print(f'Warning: Invalid DCD selection "{{dcd_selection}}", using all DCDs: {{e}}')
                dcd_indices = None
                actual_num_dcd = {num_dcd}
        else:
            dcd_indices = None
            actual_num_dcd = {num_dcd}

        results = {function_name}(
            baseDir={repr(baseDir)},
            input_pattern={repr(in4)},
            output_pattern={repr(out4)},
            num_dcd=actual_num_dcd,
            partcl_num={num_particles},
            numFrames=int({minf}),
            chunk_processing={chunk_processing},
            validate_data={validate_data},
            common_term={repr(common_term)},
            dcd_indices=dcd_indices
        )
        
        total_time = time.time() - start_time
        
        print(f'\\n✓ {calculation_type} calculation completed!')
        print(f'  Successful files: {{results["success"]}}/{num_dcd}')
        print(f'  Total time: {{total_time:.2f}}s')
        if 'data_quality' in results:
            print(f'  Data quality: {{results["data_quality"]}}')
        
    except Exception as e:
        import sys  # Ensure sys is available in exception handler
        print(f'✗ {calculation_type} calculation failed: {{e}}')
        sys.exit(1)

if __name__ == '__main__':
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def generate_dipole_script(self, baseDir, num_dcd, num_particles, max_workers, output_path):
        """Generate dipole calculation script based on selected calculation type"""
        calc_type = self.dipole_calc_type.get()
        
        if calc_type == "individual":
            return self.generate_individual_dipole_script(baseDir, num_dcd, num_particles, max_workers, output_path)
        else:
            return self.generate_collective_dipole_script(baseDir, num_dcd, num_particles, max_workers, output_path)
    
    def generate_individual_dipole_script(self, baseDir, num_dcd, num_particles, max_workers, output_path):
        """Generate individual dipole moment calculation script"""
        coords_in, com_in, vectors_out, magnitudes_out, charges_str, atoms_per_particle_str, molecules_to_process_str, stride_str = [v.get().strip() for v in self.individual_dipole_vars]
        # Apply common term expansion to individual dipole variables  
        common_term = self.common_term_var.get().strip()
        coords_in, com_in, vectors_out, magnitudes_out = [self.expand_common_term(var) for var in [coords_in, com_in, vectors_out, magnitudes_out]]
        
        # Get DCD selection
        individual_dipole_dcd_selection = self.individual_dipole_dcd_selection.get().strip()
        
        if not all([coords_in, com_in, vectors_out, magnitudes_out, charges_str, atoms_per_particle_str]):
            raise ValueError("All required individual dipole parameters must be provided")
        
        # Parse charges
        try:
            charges = [float(x.strip()) for x in charges_str.split(',') if x.strip()]
        except ValueError:
            raise ValueError("Charges must be comma-separated numbers (e.g., -0.8476,0.4238,0.4238)")
        
        try:
            atoms_per_particle = int(atoms_per_particle_str)
        except ValueError:
            raise ValueError("Atoms per particle must be an integer")
        
        # Parse molecules to process (optional)
        molecules_to_process = None
        if molecules_to_process_str:
            try:
                molecules_to_process = int(molecules_to_process_str)
            except ValueError:
                raise ValueError("Molecules to process must be an integer")
        
        stride = 1
        if stride_str:
            try:
                stride = int(stride_str)
            except ValueError:
                raise ValueError("Stride must be an integer")
        
        if len(charges) != atoms_per_particle:
            raise ValueError(f"Number of charges ({len(charges)}) must match atoms per particle ({atoms_per_particle})")
        
        use_parallel = self.dipole_parallel.get()
        chunk_processing = self.dipole_chunk_processing.get()
        validate_data = self.dipole_validate.get()
        
        # Use user-specified main file name with common term expansion (exact name, no suffixes)
        main_file = self.mainfile_var.get().strip()
        if main_file:
            script_name = self.expand_common_term(main_file)
            # Ensure .py extension
            if not script_name.endswith('.py'):
                script_name += '.py'
        else:
            script_name = "individual_dipole_calculation.py"
        
        script_path = os.path.join(output_path, script_name)
        
        script_content = f'''#!/usr/bin/env python3
"""
Individual Dipole Moment Calculation Script
Generated by MD Analysis Pipeline GUI v2.0
"""

import sys
import time
import os

# Add main_functions to Python path for imports
main_functions_path = os.path.join(os.getcwd(), 'main_functions')
if main_functions_path not in sys.path:
    sys.path.insert(0, main_functions_path)

def main():
    print('='*60)
    print('INDIVIDUAL DIPOLE MOMENT CALCULATION')
    print('='*60)
    start_time = time.time()
    
    try:
        # Import Python dipole function directly
        import main_functions.dipole_function
        dipole_functions = main_functions.dipole_function.dipole_functions
        
        # Parse DCD selection if provided
        dcd_selection = "{individual_dipole_dcd_selection}"
        if dcd_selection:
            # Parse DCD selection (e.g., "4-10", "1,3,5", "range(10,20)")
            try:
                if 'range(' in dcd_selection:
                    # Handle range() expressions
                    dcd_indices = list(eval(dcd_selection))
                elif '-' in dcd_selection and ',' in dcd_selection:
                    # Handle mixed format like "1,3,5-8"
                    dcd_indices = []
                    for part in dcd_selection.split(','):
                        if '-' in part:
                            start, end = map(int, part.split('-'))
                            dcd_indices.extend(range(start, end + 1))
                        else:
                            dcd_indices.append(int(part))
                elif '-' in dcd_selection:
                    # Handle range format like "4-10"
                    start, end = map(int, dcd_selection.split('-'))
                    dcd_indices = list(range(start, end + 1))
                else:
                    # Handle comma-separated format like "1,3,5"
                    dcd_indices = [int(x.strip()) for x in dcd_selection.split(',')]
                
                print(f'Using DCD selection: {{dcd_indices}}')
                actual_num_dcd = len(dcd_indices)
            except Exception as e:
                print(f'Warning: Invalid DCD selection "{{dcd_selection}}", using all DCDs: {{e}}')
                dcd_indices = None
                actual_num_dcd = {num_dcd}
        else:
            dcd_indices = None
            actual_num_dcd = {num_dcd}
        
        results = dipole_functions(
            baseDir={repr(baseDir)},
            coords_pattern={repr(coords_in)},
            com_pattern={repr(com_in)},
            output_pattern={repr(vectors_out)},
            magnitudes_pattern={repr(magnitudes_out)},
            Charges={charges},
            num_dcds=actual_num_dcd,
            num_particles={num_particles},
            atoms_per_particle={atoms_per_particle},
            stride={stride},
            max_workers={max_workers if use_parallel else 1},
            chunk_processing={chunk_processing},
            validate_data={validate_data},
            molecules_to_process={molecules_to_process},
            common_term={repr(common_term)},
            dcd_indices=dcd_indices
        )
        
        total_time = time.time() - start_time
        
        print(f'\\n✓ Individual dipole moment calculation completed!')
        print(f'  Successful files: {{results["success"]}}/{num_dcd}')
        print(f'  Total time: {{total_time:.2f}}s')
        if 'overall_mean_magnitude' in results:
            print(f'  Average dipole magnitude: {{results["overall_mean_magnitude"]:.3f}} ± {{results["overall_std_magnitude"]:.3f}} D')
        
    except Exception as e:
        import sys  # Ensure sys is available in exception handler
        print(f'✗ Individual dipole moment calculation failed: {{e}}')
        sys.exit(1)

if __name__ == '__main__':
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def generate_collective_dipole_script(self, baseDir, num_dcd, num_particles, max_workers, output_path):
        """Generate collective dipole moment calculation script using VMD"""
        traj_in, psf_base, dcd_base, dipole_out, target_sel, vmd_path = [v.get().strip() for v in self.collective_dipole_vars]
        # Apply common term expansion to collective dipole variables
        common_term = self.common_term_var.get().strip()
        traj_in, psf_base, dcd_base, dipole_out, vmd_path = [self.expand_common_term(var) for var in [traj_in, psf_base, dcd_base, dipole_out, vmd_path]]
        
        # Get DCD selection
        collective_dipole_dcd_selection = self.collective_dipole_dcd_selection.get().strip()
        
        if not all([traj_in, psf_base, dcd_base, dipole_out, target_sel, vmd_path]):
            raise ValueError("All collective dipole parameters must be provided")
        
        use_parallel = self.dipole_parallel.get()
        
        # Use user-specified main file name with common term expansion (exact name, no suffixes)
        main_file = self.mainfile_var.get().strip()
        if main_file:
            script_name = self.expand_common_term(main_file)
            # Ensure .py extension
            if not script_name.endswith('.py'):
                script_name += '.py'
        else:
            script_name = "collective_dipole_calculation.py"
        
        script_path = os.path.join(output_path, script_name)
        
        script_content = f'''#!/usr/bin/env python3
"""
Collective Dipole Moment Calculation Script (VMD)
Generated by MD Analysis Pipeline GUI v2.0
"""

import sys
import time
import os

# Add main_functions to Python path for imports
main_functions_path = os.path.join(os.getcwd(), 'main_functions')
if main_functions_path not in sys.path:
    sys.path.insert(0, main_functions_path)

def main():
    print('='*60)
    print('COLLECTIVE DIPOLE MOMENT CALCULATION (VMD)')
    print('='*60)
    start_time = time.time()
    
    try:
        # Import Python VMD dipole function directly
        import main_functions.vmd_dipole
        vmd_dipole_collective = main_functions.vmd_dipole.vmd_dipole_collective
        
        # Parse DCD selection if provided
        dcd_selection = "{collective_dipole_dcd_selection}"
        if dcd_selection:
            # Parse DCD selection (e.g., "4-10", "1,3,5", "range(10,20)")
            try:
                if 'range(' in dcd_selection:
                    # Handle range() expressions
                    dcd_indices = list(eval(dcd_selection))
                elif '-' in dcd_selection and ',' in dcd_selection:
                    # Handle mixed format like "1,3,5-8"
                    dcd_indices = []
                    for part in dcd_selection.split(','):
                        if '-' in part:
                            start, end = map(int, part.split('-'))
                            dcd_indices.extend(range(start, end + 1))
                        else:
                            dcd_indices.append(int(part))
                elif '-' in dcd_selection:
                    # Handle range format like "4-10"
                    start, end = map(int, dcd_selection.split('-'))
                    dcd_indices = list(range(start, end + 1))
                else:
                    # Handle comma-separated format like "1,3,5"
                    dcd_indices = [int(x.strip()) for x in dcd_selection.split(',')]
                
                print(f'Using DCD selection: {{dcd_indices}}')
                actual_num_dcd = len(dcd_indices)
            except Exception as e:
                print(f'Warning: Invalid DCD selection "{{dcd_selection}}", using all DCDs: {{e}}')
                dcd_indices = None
                actual_num_dcd = {num_dcd}
        else:
            dcd_indices = None
            actual_num_dcd = {num_dcd}
        
        # Note: vmd_dipole_collective function needs updating to pattern format
        # For now, use legacy format for collective dipole
        results = vmd_dipole_collective(
            baseDir={repr(baseDir)},
            INdir={repr(traj_in)},
            OUTdir={repr(dipole_out)},
            num_dcd=actual_num_dcd,
            num_particles={num_particles},
            psf_pattern={repr(psf_base)},
            dcd_pattern={repr(dcd_base)},
            target_selection={repr(target_sel)},
            vmd_path={repr(vmd_path)},
            common_term={repr(common_term)},
            max_workers={max_workers if use_parallel else 1},
            dcd_indices=dcd_indices
        )
        
        total_time = time.time() - start_time
        
        print(f'\\n✓ Collective dipole moment calculation completed!')
        print(f'  Successful files: {{results["success"]}}/{num_dcd}')
        print(f'  Total time: {{total_time:.2f}}s')
        print(f'  Target selection: {target_sel}')
        
    except Exception as e:
        import sys  # Ensure sys is available in exception handler
        print(f'✗ Collective dipole moment calculation failed: {{e}}')
        sys.exit(1)

if __name__ == '__main__':
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def generate_velocity_script(self, baseDir, num_dcd, num_particles, max_workers, output_path):
        """Generate velocity extraction script"""
        psf_pattern, veldcd_pattern, output_pattern, num_molecules_str, vmd_path = [v.get().strip() for v in self.velocity_vars]
        
        # Apply common term expansion to velocity variables
        common_term = self.common_term_var.get().strip()
        psf_pattern, veldcd_pattern, output_pattern, vmd_path = [self.expand_common_term(var) for var in [psf_pattern, veldcd_pattern, output_pattern, vmd_path]]
        
        # Validate required fields
        if not all([psf_pattern, veldcd_pattern, output_pattern, num_molecules_str, vmd_path]):
            raise ValueError("All velocity extraction parameters are required")
        
        try:
            num_molecules = int(num_molecules_str)
        except ValueError:
            raise ValueError("Number of molecules must be a valid integer")
        
        # Get DCD selection
        velocity_dcd_selection = self.velocity_dcd_selection.get().strip()
        
        # Get processing options
        parallel_processing = self.velocity_parallel.get()
        validate_data = self.velocity_validate.get()
        
        # Use user-defined main file name or fallback to default
        main_file = self.mainfile_var.get().strip()
        if main_file:
            script_name = self.expand_common_term(main_file)
            # Ensure .py extension
            if not script_name.endswith('.py'):
                script_name += '.py'
        else:
            script_name = "velocity_extraction.py"
        
        script_path = os.path.join(output_path, script_name)
        
        script_content = f'''#!/usr/bin/env python3
"""
Velocity Extraction Script
Generated by MD Analysis Pipeline GUI v2.0
"""

import sys
import time
import os

# Add main_functions to Python path for imports
main_functions_path = os.path.join(os.getcwd(), 'main_functions')
if main_functions_path not in sys.path:
    sys.path.insert(0, main_functions_path)

def main():
    print('='*60)
    print('VELOCITY EXTRACTION')
    print('='*60)
    start_time = time.time()
    
    try:
        # Import Python function directly
        import main_functions.velocity_extract
        extract_velocities = main_functions.velocity_extract.extract_velocities

        # Parse DCD selection if provided
        dcd_selection = "{velocity_dcd_selection}"
        if dcd_selection:
            # Parse DCD selection (e.g., "4-10", "1,3,5", "range(10,20)")
            try:
                if 'range(' in dcd_selection:
                    # Handle range() expressions
                    dcd_indices = list(eval(dcd_selection))
                elif '-' in dcd_selection and ',' in dcd_selection:
                    # Handle mixed format like "1,3,5-8"
                    dcd_indices = []
                    for part in dcd_selection.split(','):
                        if '-' in part:
                            start, end = map(int, part.split('-'))
                            dcd_indices.extend(range(start, end + 1))
                        else:
                            dcd_indices.append(int(part))
                elif '-' in dcd_selection:
                    # Handle range format like "4-10"
                    start, end = map(int, dcd_selection.split('-'))
                    dcd_indices = list(range(start, end + 1))
                else:
                    # Handle comma-separated format like "1,3,5"
                    dcd_indices = [int(x.strip()) for x in dcd_selection.split(',')]
                
                print(f'Using DCD selection: {{dcd_indices}}')
                actual_num_dcd = len(dcd_indices)
            except Exception as e:
                print(f'Warning: Invalid DCD selection "{{dcd_selection}}", using all DCDs: {{e}}')
                dcd_indices = None
                actual_num_dcd = {num_dcd}
        else:
            dcd_indices = None
            actual_num_dcd = {num_dcd}

        results = extract_velocities(
            baseDir={repr(baseDir)},
            psf_pattern={repr(psf_pattern)},
            veldcd_pattern={repr(veldcd_pattern)},
            output_pattern={repr(output_pattern)},
            num_dcd=actual_num_dcd,
            num_molecules={num_molecules},
            vmd={repr(vmd_path)},
            max_workers={max_workers},
            dcd_indices=dcd_indices,
            common_term={repr(common_term)}
        )
        
        total_time = time.time() - start_time
        
        print(f'\\n✓ Velocity extraction completed!')
        print(f'  Successful files: {{results["successful"]}}/{num_dcd}')
        print(f'  Failed files: {{results["failed"]}}')
        print(f'  Total time: {{total_time:.2f}}s')
        
    except Exception as e:
        import sys  # Ensure sys is available in exception handler
        print(f'✗ Velocity extraction failed: {{e}}')
        sys.exit(1)

if __name__ == '__main__':
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def update_dipole_fields(self):
        """Update the displayed dipole parameter fields based on calculation type"""
        calc_type = self.dipole_calc_type.get()
        
        # Hide both frames first
        self.individual_params_frame.pack_forget()
        self.collective_params_frame.pack_forget()
        
        # Show the appropriate frame
        if calc_type == "individual":
            self.individual_params_frame.pack(fill="x", padx=10, pady=10)
            # Update the dipole_vars reference for backward compatibility
            self.dipole_vars = self.individual_dipole_vars
        else:  # collective
            self.collective_params_frame.pack(fill="x", padx=10, pady=10)
            # For collective mode, we'll handle parameters differently in script generation
            
        # Force the window to update its layout
        self.update_idletasks()
    
    def browse_vmd_dipole(self, var):
        """Browse for VMD executable for dipole calculations"""
        from tkinter import filedialog
        f = filedialog.askopenfilename(
            title="Select VMD Executable",
            filetypes=[("Executable","*"), ("All files", "*")]
        )
        if f:
            var.set(f)
    
    def close_window(self):
        """Close the calculations window"""
        self.grab_release()
        self.destroy()

class PipelineGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MD Dynamics Analysis Pipeline - Data Preprocessing")
        self.geometry("1600x1000")  # Wider and shorter for better layout
        self.configure(bg='#f0f0f0')  # Light gray background for modern look
        
        # Configure ttk styling for better appearance
        self.style = configure_ttk_style()
        
        # Initialize output path variable
        self.output_full_path = None

        # Create main frame and scrollbar
        self.main_frame = tk.Frame(self, bg='#f0f0f0')
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # Create canvas and scrollbars
        self.canvas = tk.Canvas(self.main_frame, bg='#f0f0f0', highlightthickness=0)
        self.v_scrollbar = ttk.Scrollbar(self.main_frame, orient="vertical", command=self.canvas.yview)
        self.h_scrollbar = ttk.Scrollbar(self.main_frame, orient="horizontal", command=self.canvas.xview)
        self.scrollable_frame = tk.Frame(self.canvas, bg='#f0f0f0')

        # Configure scrollable frame
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )

        # Create window in canvas
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.v_scrollbar.set, xscrollcommand=self.h_scrollbar.set)

        # Pack canvas and scrollbars
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.v_scrollbar.grid(row=0, column=1, sticky="ns")
        self.h_scrollbar.grid(row=1, column=0, sticky="ew")
        
        # Configure grid weights for proper resizing
        self.main_frame.grid_rowconfigure(0, weight=1)
        self.main_frame.grid_columnconfigure(0, weight=1)
        
        # Configure scrollable frame grid for the new layout
        self.scrollable_frame.grid_columnconfigure(0, weight=1)
        
        # Initialize all essential variables first
        self.baseDir_var = tk.StringVar(value="")
        self.num_dcd_var = tk.IntVar(value=1)
        self.num_particles_var = tk.IntVar(value=1000)
        self.common_term_var = tk.StringVar(value="")
        self.output_folder_var = tk.StringVar(value="")
        self.mainfile_var = tk.StringVar(value="")
        self.submitfile_var = tk.StringVar(value="")
        self.max_workers_var = tk.IntVar(value=min(4, mp.cpu_count()))
        
        # Initialize lists for step variables
        self.coords_vars = []
        self.unwrap_vars = []
        self.com_vars = []
        self.sbatch_vars = []
        
        # Initialize additional variables
        self.unwrap_chunk_var = tk.StringVar(value="5000")
        self.coords_dcd_selection = tk.StringVar(value="")
        self.unwrap_dcd_selection = tk.StringVar(value="")
        self.com_dcd_selection = tk.StringVar(value="")
        self.unwrap_opt = [tk.StringVar(value=""), tk.StringVar(value="")]
        
        # Initialize trajectory characteristics for intelligent sizing
        self.traj_file_size_mb = tk.StringVar(value="")
        self.traj_num_frames = tk.StringVar(value="")
        self.traj_num_atoms = tk.StringVar(value="")
        self.available_memory_gb = tk.StringVar(value="")
        self.coord_precision = tk.StringVar(value="single")
        
        # Add callback to clear browse path when user manually types folder name
        self.output_folder_var.trace('w', self.on_folder_name_changed)

        # Bind mouse wheel to canvas (Windows)
        self.canvas.bind("<MouseWheel>", self._on_mousewheel)
        self.canvas.bind("<Shift-MouseWheel>", self._on_horizontal_mousewheel)
        self.bind_all("<MouseWheel>", self._on_mousewheel)
        self.bind_all("<Shift-MouseWheel>", self._on_horizontal_mousewheel)
        
        # Bind mouse wheel for Linux
        self.canvas.bind("<Button-4>", self._on_mousewheel)
        self.canvas.bind("<Button-5>", self._on_mousewheel)
        self.canvas.bind("<Shift-Button-4>", self._on_horizontal_mousewheel)
        self.canvas.bind("<Shift-Button-5>", self._on_horizontal_mousewheel)
        self.bind_all("<Button-4>", self._on_mousewheel)
        self.bind_all("<Button-5>", self._on_mousewheel)
        self.bind_all("<Shift-Button-4>", self._on_horizontal_mousewheel)
        self.bind_all("<Shift-Button-5>", self._on_horizontal_mousewheel)
        
        # Bind canvas click to focus
        self.canvas.bind("<Button-1>", lambda e: self.canvas.focus_set())

        # * fields are required
        header_frame = tk.Frame(self.scrollable_frame, bg='#f0f0f0')
        header_frame.grid(row=0, column=0, sticky="ew", padx=15, pady=(15,0))
        
        tk.Label(header_frame, text="MD Data Preprocessing Pipeline", 
                font=("Arial", 14, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(anchor='w')
        tk.Label(header_frame, text="Steps 1-3: Extract → Unwrap → COM | * Required fields | v2.0.2", 
                font=("Arial", 8), fg='#7f8c8d', bg='#f0f0f0').pack(anchor='w')

        # --- Common Parameters ---
        common = tk.LabelFrame(self.scrollable_frame, text="Common Parameters", 
                              font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0',
                              relief='groove', borderwidth=1)
        common.grid(row=1, column=0, padx=15, pady=10, sticky="ew")

        tk.Label(common, text="Base Directory:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=0, sticky="e", padx=5, pady=5)
        base_entry = tk.Entry(common, textvariable=self.baseDir_var, width=70, font=("Arial", 10), 
                             relief='solid', borderwidth=1)
        base_entry.grid(row=0, column=1, padx=5, pady=5)
        create_tooltip(base_entry, "Root directory containing all simulation data and analysis folders")
        tk.Label(common, text="Number of DCDs:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=3, sticky="e", padx=5, pady=5)
        dcd_entry = tk.Entry(common, textvariable=self.num_dcd_var, width=15, font=("Arial", 10),
                            relief='solid', borderwidth=1)
        dcd_entry.grid(row=0, column=4, sticky="w", padx=5, pady=5)
        create_tooltip(dcd_entry, "Total number of trajectory DCD files to process (e.g., 100 for com_0.dat through com_99.dat)")

        tk.Label(common, text="Number of Particles:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=2, column=0, sticky="e", padx=5, pady=5)
        particles_entry = tk.Entry(common, textvariable=self.num_particles_var, width=15, font=("Arial", 10),
                                  relief='solid', borderwidth=1)
        particles_entry.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(particles_entry, "Total number of molecules/particles in each trajectory file")

        # Common term placeholder
        tk.Label(common, text="Common Term {*}:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=2, column=3, sticky="e", padx=5, pady=5)
        common_term_entry = tk.Entry(common, textvariable=self.common_term_var, width=30, font=("Arial", 10),
                                     relief='solid', borderwidth=1)
        common_term_entry.grid(row=2, column=4, sticky="w", padx=5, pady=5)
        create_tooltip(common_term_entry, "Common term to replace '*' in all path/name fields. Example: '/scratch/user/project' - then use '*' in paths like '*/input' → '/scratch/user/project/input'")

        # --- Trajectory Characteristics for Intelligent Optimization ---
        traj_info = tk.LabelFrame(self.scrollable_frame, text="Trajectory Characteristics (for Smart Optimization)", 
                                 font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0',
                                 relief='groove', borderwidth=1)
        traj_info.grid(row=2, column=0, padx=15, pady=10, sticky="ew")

        # First row - file size and frames
        tk.Label(traj_info, text="Single DCD file size (MB):", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=0, sticky="e", padx=5, pady=5)
        file_size_entry = tk.Entry(traj_info, textvariable=self.traj_file_size_mb, width=15, font=("Arial", 10),
                                  relief='solid', borderwidth=1)
        file_size_entry.grid(row=0, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(file_size_entry, "Size of a single DCD file in MB (e.g., 500). Used to estimate memory requirements and optimize chunk sizes.")

        tk.Label(traj_info, text="Frames per DCD:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=2, sticky="e", padx=5, pady=5)
        frames_entry = tk.Entry(traj_info, textvariable=self.traj_num_frames, width=15, font=("Arial", 10),
                               relief='solid', borderwidth=1)
        frames_entry.grid(row=0, column=3, sticky="w", padx=5, pady=5)
        create_tooltip(frames_entry, "Number of frames in each DCD file (e.g., 25000). Used to calculate memory per frame.")

        # Second row - atoms and precision
        tk.Label(traj_info, text="Atoms being extracted:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=1, column=0, sticky="e", padx=5, pady=5)
        atoms_extracted_entry = tk.Entry(traj_info, textvariable=self.traj_num_atoms, width=15, font=("Arial", 10),
                                        relief='solid', borderwidth=1)
        atoms_extracted_entry.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(atoms_extracted_entry, "Number of atoms being extracted in Step 1 (e.g., 3000 for 1000 molecules × 3 atoms each). Used to calculate actual output file sizes.")

        tk.Label(traj_info, text="Coordinate precision:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=1, column=2, sticky="e", padx=5, pady=5)
        self.coord_precision = tk.StringVar(value="single")
        precision_combo = ttk.Combobox(traj_info, textvariable=self.coord_precision, width=12, font=("Arial", 10),
                                      values=["single", "double"], state="readonly")
        precision_combo.grid(row=1, column=3, sticky="w", padx=5, pady=5)
        create_tooltip(precision_combo, "Coordinate precision: 'single' (4 bytes, ~50% smaller files) or 'double' (8 bytes, full precision)")

        # Third row - available memory and max workers
        tk.Label(traj_info, text="Available memory (GB):", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=4, sticky="e", padx=5, pady=5)
        memory_entry = tk.Entry(traj_info, textvariable=self.available_memory_gb, width=15, font=("Arial", 10),
                               relief='solid', borderwidth=1)
        memory_entry.grid(row=0, column=5, sticky="w", padx=5, pady=5)
        create_tooltip(memory_entry, "Available system memory in GB (e.g., 240). Used to optimize chunk sizes and worker counts.")

        tk.Label(traj_info, text="Max Workers:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=1, column=4, sticky="e", padx=5, pady=5)
        workers_entry = tk.Entry(traj_info, textvariable=self.max_workers_var, width=15, font=("Arial", 10),
                                relief='solid', borderwidth=1)
        workers_entry.grid(row=1, column=5, sticky="w", padx=5, pady=5)
        create_tooltip(workers_entry, "Number of CPU cores for parallel processing. Used in chunk size optimization: more workers = faster but smaller chunk sizes.")
        tk.Label(traj_info, text=f"(auto: {mp.cpu_count()})", 
                font=("Arial", 9), bg='#f0f0f0', fg='#7f8c8d').grid(row=1, column=6, sticky="w", padx=5)

        # Calculate button
        calc_btn = tk.Button(traj_info, text="🧠 Calculate Optimal Settings", 
                            command=self.calculate_optimal_settings, 
                            bg="#3498db", fg="white", 
                            font=("Arial", 10, "bold"), padx=15, pady=5,
                            relief='raised', borderwidth=2, cursor='hand2')
        calc_btn.grid(row=3, column=0, columnspan=2, padx=5, pady=10, sticky="w")
        create_tooltip(calc_btn, "Automatically calculate optimal chunk sizes, worker counts, and memory settings based on your trajectory characteristics")

        # Results label
        self.optimization_results_label = tk.Label(traj_info, text="Enter trajectory info and click 'Calculate Optimal Settings' for recommendations", 
                                                  font=("Arial", 9), bg='#f0f0f0', fg='#7f8c8d', wraplength=600)
        self.optimization_results_label.grid(row=3, column=2, columnspan=2, sticky="w", padx=5, pady=10)

        # --- Steps container ---
        steps = tk.Frame(self.scrollable_frame, bg='#f0f0f0')
        steps.grid(row=3, column=0, padx=10, pady=5, sticky="ew")

        # ----- STEP 1: Coordinate Extraction -----
        self.skip1 = tk.BooleanVar(value=False)
        coords_frame = tk.Frame(steps, bg='#f0f0f0')
        coords_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nw")
        
        # Create title frame with skip checkbox
        title_frame = tk.Frame(coords_frame, bg='#f0f0f0')
        title_frame.pack(fill="x", padx=2, pady=2)
        
        tk.Label(title_frame, text="Step 1: coordinates_extract", 
                font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(side="left")
        skip1_check = tk.Checkbutton(title_frame, text="Skip", variable=self.skip1,
                                    command=lambda: self.toggle_frame(coords, self.skip1.get()),
                                    font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50', 
                                    selectcolor='#e8f4fd', activebackground='#f0f0f0')
        skip1_check.pack(side="right", padx=10)
        create_tooltip(skip1_check, "Skip coordinate extraction step if already completed")
        
        coords = tk.LabelFrame(coords_frame, text="", relief='groove', borderwidth=1, bg='#f0f0f0')
        coords.pack(fill="both", expand=True, padx=2, pady=2)

        labels1 = ["PSF Pattern","DCD Pattern","Output Pattern","Target Selections","VMD path"]
        tooltips1 = [
            "Path pattern for PSF files. Use * for common term, {i} for file index. Example: trajectories_*/run_{i}/system.psf",
            "Path pattern for DCD files. Use * for common term, {i} for file index. Example: trajectories_*/run_{i}/traj.dcd",
            "Path pattern for output files. Use * for common term, {i} for file index. Example: path/to/xyz_{i}.dat",
            "VMD atom selection (e.g., 'resname WAT and residue 0 to 999', 'water', 'protein')",
            "Full path to VMD executable"
        ]
        self.coords_vars = []
        for i, (lbl, tooltip) in enumerate(zip(labels1, tooltips1), start=0):
            tk.Label(coords, text=f"{lbl}:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=5)
            v = tk.StringVar()
            ent = tk.Entry(coords, textvariable=v, width=50, font=("Arial", 10), relief='solid', borderwidth=1)
            ent.grid(row=i, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
            create_tooltip(ent, tooltip)
            if lbl=="Target Selections":
                tk.Label(coords, text='e.g., "resname WAT"', font=("Arial", 7), bg='#f0f0f0', fg='#7f8c8d').grid(
                    row=i, column=3, sticky="w", padx=5
                )
            self.coords_vars.append(v)
        
        # DCD selection for coordinates_extract
        tk.Label(coords, text="DCD Selection (optional):", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=len(labels1), column=0, sticky="e", padx=5, pady=5)
        self.coords_dcd_selection = tk.StringVar()
        coords_dcd_entry = tk.Entry(coords, textvariable=self.coords_dcd_selection, width=25, font=("Arial", 10), relief='solid', borderwidth=1)
        coords_dcd_entry.grid(row=len(labels1), column=1, padx=5, pady=5)
        create_tooltip(coords_dcd_entry, "Optional: Specify which DCDs to process (e.g., '4-10' for DCDs 4 to 10, '4-6,8-10' for DCDs 4-6 and 8-10, or 'range(10,20)' for Python range). Leave empty to process all DCDs.")
        tk.Label(coords, text='e.g., "4-10", "4-6,8-10", or "range(10,20)"', font=("Arial", 7), bg='#f0f0f0', fg='#7f8c8d').grid(
            row=len(labels1), column=2, sticky="w", padx=5
        )

        # Advanced options for coordinates_extract
        tk.Label(coords, text="Use Parallel VMD:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=len(labels1)+1, column=0, sticky="e", padx=5, pady=5)
        self.coords_parallel = tk.BooleanVar(value=True)
        coords_parallel_check = tk.Checkbutton(coords, variable=self.coords_parallel, font=("Arial", 10), 
                                              bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        coords_parallel_check.grid(row=len(labels1)+1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(coords_parallel_check, "Enable parallel VMD processing for faster coordinate extraction")
        
        self.toggle_frame(coords, self.skip1.get())

        # ----- STEP 2: Unwrap Coordinates -----
        self.skip2 = tk.BooleanVar(value=False)
        unwrap_frame = tk.Frame(steps, bg='#f0f0f0')
        unwrap_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nw")
        
        # Create title frame with skip checkbox
        title_frame2 = tk.Frame(unwrap_frame, bg='#f0f0f0')
        title_frame2.pack(fill="x", padx=2, pady=2)
        
        tk.Label(title_frame2, text="Step 2: unwrap_coords", 
                font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(side="left")
        skip2_check = tk.Checkbutton(title_frame2, text="Skip", variable=self.skip2,
                                    command=lambda: self.toggle_frame(unwrap, self.skip2.get()),
                                    font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50', 
                                    selectcolor='#e8f4fd', activebackground='#f0f0f0')
        skip2_check.pack(side="right", padx=10)
        create_tooltip(skip2_check, "Skip coordinate unwrapping step if already completed")
        
        unwrap = tk.LabelFrame(unwrap_frame, text="", relief='groove', borderwidth=1, bg='#f0f0f0')
        unwrap.pack(fill="both", expand=True, padx=2, pady=2)

        labels2 = ["Input Pattern", "Output Pattern", "XSC file", "Num atoms"]
        tooltips2 = [
            "Path pattern for input coordinate files. Use * for common term, {i} for file index. Example: anlz/NVT_*/wrapped/xyz_{i}.dat",
            "Path pattern for output unwrapped files. Use * for common term, {i} for file index. Example: anlz/NVT_*/unwrapped/unwrapped_xyz_{i}.dat",
            "Path to XSC file containing box dimensions (e.g., restart_equil.xsc)",
            "Total number of atoms for which coordinates will be processed"
        ]
        self.unwrap_vars = []
        for i, (lbl, tooltip) in enumerate(zip(labels2, tooltips2), start=0):
            tk.Label(unwrap, text=f"{lbl}:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=5)
            v = tk.StringVar()
            ent = tk.Entry(unwrap, textvariable=v, width=50, font=("Arial", 10), relief='solid', borderwidth=1)
            ent.grid(row=i, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
            create_tooltip(ent, tooltip)
            self.unwrap_vars.append(v)
        
        # Optional interval & stride
        self.unwrap_opt = []
        for j,lbl in enumerate(["Interval (optional)","Stride (optional)"], start=len(labels2)):
            tk.Label(unwrap, text=f"{lbl}:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=j, column=0, sticky="e", padx=5, pady=5)
            v = tk.StringVar()
            ent = tk.Entry(unwrap, textvariable=v, width=25, font=("Arial", 10), relief='solid', borderwidth=1)
            ent.grid(row=j, column=1, padx=5, pady=5)
            self.unwrap_opt.append(v)
        
        # DCD selection for unwrap_coords
        row_offset = len(labels2) + len(self.unwrap_opt)
        tk.Label(unwrap, text="DCD Selection (optional):", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=row_offset, column=0, sticky="e", padx=5, pady=5)
        self.unwrap_dcd_selection = tk.StringVar()
        unwrap_dcd_entry = tk.Entry(unwrap, textvariable=self.unwrap_dcd_selection, width=25, font=("Arial", 10), relief='solid', borderwidth=1)
        unwrap_dcd_entry.grid(row=row_offset, column=1, padx=5, pady=5)
        create_tooltip(unwrap_dcd_entry, "Optional: Specify which DCDs to process (e.g., '4-10' for DCDs 4 to 10, '4-6,8-10' for DCDs 4-6 and 8-10, or 'range(10,20)' for Python range). Leave empty to process all DCDs.")
        tk.Label(unwrap, text='e.g., "4-10", "4-6,8-10", or "range(10,20)"', font=("Arial", 7), bg='#f0f0f0', fg='#7f8c8d').grid(
            row=row_offset, column=2, sticky="w", padx=5
        )

        # Advanced options for unwrap_coords
        tk.Label(unwrap, text="Chunk Size:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=row_offset+1, column=0, sticky="e", padx=5, pady=5)
        self.unwrap_chunk_var = tk.StringVar(value="5000")
        chunk_entry = tk.Entry(unwrap, textvariable=self.unwrap_chunk_var, width=15, font=("Arial", 10), relief='solid', borderwidth=1)
        chunk_entry.grid(row=row_offset+1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(chunk_entry, "Memory chunk size for processing large files. Smaller values (e.g., 5000) use less memory but take longer. Use 'auto' for automatic sizing, or specify number of frames.")
        
        tk.Label(unwrap, text="Use Parallel:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=row_offset+2, column=0, sticky="e", padx=5, pady=5)
        self.unwrap_parallel = tk.BooleanVar(value=True)
        parallel_check = tk.Checkbutton(unwrap, variable=self.unwrap_parallel, font=("Arial", 10),
                                       bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        parallel_check.grid(row=row_offset+2, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(parallel_check, "Enable parallel processing for faster unwrapping of coordinates")
        
        self.toggle_frame(unwrap, self.skip2.get())

        # ----- STEP 3: COM Calculation -----
        self.skip3 = tk.BooleanVar(value=False)
        com_frame = tk.Frame(steps, bg='#f0f0f0')
        com_frame.grid(row=0, column=2, padx=5, pady=5, sticky="nw")
        
        # Create title frame with skip checkbox
        title_frame3 = tk.Frame(com_frame, bg='#f0f0f0')
        title_frame3.pack(fill="x", padx=2, pady=2)
        
        tk.Label(title_frame3, text="Step 3: COM_calc", 
                font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0').pack(side="left")
        skip3_check = tk.Checkbutton(title_frame3, text="Skip", variable=self.skip3,
                                    command=lambda: self.toggle_frame(com, self.skip3.get()),
                                    font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50', 
                                    selectcolor='#e8f4fd', activebackground='#f0f0f0')
        skip3_check.pack(side="right", padx=10)
        create_tooltip(skip3_check, "Skip COM calculation step if already completed")
        
        com = tk.LabelFrame(com_frame, text="", relief='groove', borderwidth=1, bg='#f0f0f0')
        com.pack(fill="both", expand=True, padx=2, pady=2)

        labels3 = ["Input Pattern", "Output Pattern", "Atoms per particle", "Mass list"]
        tooltips3 = [
            "Path pattern for input unwrapped coordinate files. Use * for common term, {i} for file index. Example: anlz/NVT_*/unwrapped/unwrapped_xyz_{i}.dat",
            "Path pattern for output COM files. Use * for common term, {i} for file index. Example: anlz/NVT_*/com_data/com_{i}.dat",
            "Number of atoms per molecule (e.g., 3 for water)",
            "Comma-separated atomic masses (e.g., '16.0,1.008,1.008' for water O H H)"
        ]
        self.com_vars = []
        for i, (lbl, tooltip) in enumerate(zip(labels3, tooltips3), start=0):
            tk.Label(com, text=f"{lbl}:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=5)
            v = tk.StringVar()
            ent = tk.Entry(com, textvariable=v, width=50, font=("Arial", 10), relief='solid', borderwidth=1)
            ent.grid(row=i, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
            create_tooltip(ent, tooltip)
            if lbl=="Mass list":
                tk.Label(com, text='e.g., "16.0,1.008,1.008"', font=("Arial", 7), bg='#f0f0f0', fg='#7f8c8d').grid(
                    row=i, column=3, sticky="w", padx=5
                )
            self.com_vars.append(v)
        
        # DCD selection for COM_calc
        row_offset = len(labels3)
        tk.Label(com, text="DCD Selection (optional):", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=row_offset, column=0, sticky="e", padx=5, pady=5)
        self.com_dcd_selection = tk.StringVar()
        com_dcd_entry = tk.Entry(com, textvariable=self.com_dcd_selection, width=25, font=("Arial", 10), relief='solid', borderwidth=1)
        com_dcd_entry.grid(row=row_offset, column=1, padx=5, pady=5)
        create_tooltip(com_dcd_entry, "Optional: Specify which DCDs to process (e.g., '4-10' for DCDs 4 to 10, '4-6,8-10' for DCDs 4-6 and 8-10, or 'range(10,20)' for Python range). Leave empty to process all DCDs.")
        tk.Label(com, text='e.g., "4-10", "4-6,8-10", or "range(10,20)"', font=("Arial", 7), bg='#f0f0f0', fg='#7f8c8d').grid(
            row=row_offset, column=2, sticky="w", padx=5
        )

        # Advanced options for COM_calc
        tk.Label(com, text="Use Parallel:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=row_offset+1, column=0, sticky="e", padx=5, pady=5)
        self.com_parallel = tk.BooleanVar(value=True)
        com_parallel_check = tk.Checkbutton(com, variable=self.com_parallel, font=("Arial", 10),
                                           bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        com_parallel_check.grid(row=row_offset+1, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(com_parallel_check, "Enable parallel processing for faster COM calculations")
        
        tk.Label(com, text="Use Memory Map:", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=row_offset+2, column=0, sticky="e", padx=5, pady=5)
        self.com_memmap = tk.BooleanVar(value=True)
        com_memmap_check = tk.Checkbutton(com, variable=self.com_memmap, font=("Arial", 10),
                                         bg='#f0f0f0', fg='#2c3e50', selectcolor='#e8f4fd', activebackground='#f0f0f0')
        com_memmap_check.grid(row=row_offset+2, column=1, sticky="w", padx=5, pady=5)
        create_tooltip(com_memmap_check, "Use memory mapping for very large files to reduce memory usage")
        tk.Label(com, text="(for very large files)", font=("Arial", 7), bg='#f0f0f0', fg='#7f8c8d').grid(row=row_offset+2, column=2, sticky="w", padx=5)
        
        self.toggle_frame(com, self.skip3.get())



        # --- SLURM parameters and Output Files side by side ---
        params_container = tk.Frame(self.scrollable_frame, bg='#f0f0f0')
        params_container.grid(row=4, column=0, padx=15, pady=10, sticky="ew")
        
        # SLURM parameters on the left
        sb = tk.LabelFrame(params_container, text="SLURM Submission Parameters",
                          font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0',
                          relief='groove', borderwidth=1)
        sb.grid(row=0, column=0, padx=(0, 10), pady=0, sticky="nw")
        
        sb_labels = ["Nodes","Partition","QOS","CPUs","Tasks","Memory (GB)","Walltime","Output prefix","Email","Module"]
        sb_tooltips = [
            "Number of compute nodes to request (usually 1 for single-node jobs)",
            "SLURM partition/queue name (e.g., 'standard', 'gpu', 'high-mem')", 
            "Quality of Service level for job priority",
            "Number of CPU cores per node to request",
            "Number of tasks/processes (usually 1 for serial jobs)",
            "Memory allocation in GB (e.g., 160 for 160GB total memory)",
            "Maximum runtime (format: HH:MM:SS or DD-HH:MM:SS)",
            "Prefix for output log files",
            "Email address for job notifications",
            "Module to load before execution (optional, e.g., 'python/3.9', 'anaconda3')"
        ]
        
        self.sbatch_vars = []
        for i, (lbl, tooltip) in enumerate(zip(sb_labels, sb_tooltips)):
            # Module field is optional, others are required
            required = "*" if lbl != "Module" else ""
            tk.Label(sb, text=f"{lbl}:{required}", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=i, column=0, sticky="e", padx=5, pady=3)
            v = tk.StringVar()
            if lbl == "CPUs":
                v.set(str(self.max_workers_var.get()))  # Default to max_workers
            elif lbl == "Tasks":
                v.set("1")  # Usually 1 for this type of job
            entry = tk.Entry(sb, textvariable=v, width=18, font=("Arial", 10))
            entry.grid(row=i, column=1, sticky="w", padx=5, pady=3)
            create_tooltip(entry, tooltip)
            self.sbatch_vars.append(v)

        # --- Output Files on the right ---
        files = tk.LabelFrame(params_container, text="Output Files", 
                             font=("Arial", 11, "bold"), fg='#2c3e50', bg='#f0f0f0',
                             relief='groove', borderwidth=1)
        files.grid(row=0, column=1, padx=(10, 0), pady=0, sticky="nw")
        
        tk.Label(files, text="Output folder name:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=0, column=0, sticky="e", padx=5, pady=3)
        folder_entry = tk.Entry(files, textvariable=self.output_folder_var, width=60, font=("Arial", 10), relief='solid', borderwidth=1)
        folder_entry.grid(row=0, column=1, columnspan=2, padx=5, pady=3, sticky="ew")
        create_tooltip(folder_entry, "Name of the folder to contain all generated files. Will be created in current directory.")
        
        # Show output path
        self.output_path_label = tk.Label(files, text="", font=("Arial", 8), fg='#7f8c8d', bg='#f0f0f0', wraplength=400)
        self.output_path_label.grid(row=1, column=0, columnspan=3, sticky="w", padx=5)
        self.update_output_path_label()
        
        # Main script file
        tk.Label(files, text="Main script file:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=2, column=0, sticky="e", padx=5, pady=3)
        main_entry = tk.Entry(files, textvariable=self.mainfile_var, width=60, font=("Arial", 10), relief='solid', borderwidth=1)
        main_entry.grid(row=2, column=1, columnspan=2, padx=5, pady=3, sticky="ew")
        create_tooltip(main_entry, "Python script filename (will be created in output folder)")
        
        # Submit script file
        tk.Label(files, text="Submit script file:*", font=("Arial", 10), bg='#f0f0f0', fg='#2c3e50').grid(row=3, column=0, sticky="e", padx=5, pady=3)
        submit_entry = tk.Entry(files, textvariable=self.submitfile_var, width=60, font=("Arial", 10), relief='solid', borderwidth=1)
        submit_entry.grid(row=3, column=1, columnspan=2, padx=5, pady=3, sticky="ew")
        create_tooltip(submit_entry, "SLURM batch script filename (.sh extension will be added automatically)")

        # Generate buttons
        generate_frame = tk.Frame(files, bg='#f0f0f0')
        generate_frame.grid(row=4, column=0, columnspan=3, pady=10)
        
        generate_btn = tk.Button(generate_frame, text="Generate Files",
                                command=self.generate_files, bg="#27ae60", fg="white",
                                font=("Arial", 11, "bold"), relief='raised', borderwidth=3,
                                cursor='hand2', padx=20, pady=5)
        generate_btn.pack(side=tk.LEFT, padx=10)
        create_tooltip(generate_btn, "Generate all pipeline scripts with optimized settings for your specified configuration")
        
        calc_window_btn = tk.Button(generate_frame, text="Open Analysis Calculations", 
                                   command=self.open_calculations_window, 
                                   bg="#3498db", fg="white", 
                                   font=("Arial", 11, "bold"), padx=20, pady=5,
                                   relief='raised', borderwidth=3, cursor='hand2')
        calc_window_btn.pack(side=tk.LEFT, padx=10)
        create_tooltip(calc_window_btn, "Open separate window for Alpha2/MSD and Dipole moment calculations")
        
        mult_btn = tk.Button(generate_frame, text="Multi-Run Analysis",
                           command=self.generate_mult_run, bg="#e67e22", fg="white",
                           font=("Arial", 11, "bold"), relief='raised', borderwidth=3,
                           cursor='hand2', padx=20, pady=5)
        mult_btn.pack(side=tk.LEFT, padx=10)
        create_tooltip(mult_btn, "Generate scripts for multi-temperature or multi-condition analysis workflows")

        # Load last inputs if any
        self.load_config()
        
        # Window initialization complete
        pass

    def _on_mousewheel(self, event):
        """Handle mouse wheel scrolling."""
        # Cross-platform mouse wheel handling
        if event.delta:
            # Windows
            self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        else:
            # Linux
            if event.num == 4:
                self.canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                self.canvas.yview_scroll(1, "units")

    def _on_horizontal_mousewheel(self, event):
        """Handle horizontal mouse wheel scrolling (Shift + wheel)."""
        # Cross-platform horizontal mouse wheel handling
        if event.delta:
            # Windows
            self.canvas.xview_scroll(int(-1*(event.delta/120)), "units")
        else:
            # Linux
            if event.num == 4:
                self.canvas.xview_scroll(-1, "units")
            elif event.num == 5:
                self.canvas.xview_scroll(1, "units")


    def toggle_frame(self, frame, skip):
        """Disable or enable all entries & buttons in a step frame,
           but leave the skip‐checkbox itself active."""
        state = "disabled" if skip else "normal"
        for w in frame.winfo_children():
            # only disable entries and normal buttons
            if isinstance(w, (tk.Entry, tk.Button)):
                w.configure(state=state)
    

    
    def parse_dcd_selection(self, selection_str, num_dcd):
        """
        Parse DCD selection string and return list of DCD indices.
        
        Parameters:
        -----------
        selection_str : str
            Selection string like "4-10", "4-6,8-10", "range(10,20)", or empty string
        num_dcd : int
            Total number of DCDs available
            
        Returns:
        --------
        list : List of DCD indices to process
        
        Examples:
        ---------
        >>> parse_dcd_selection("4-10", 20)
        [4, 5, 6, 7, 8, 9, 10]
        >>> parse_dcd_selection("4-6,8-10", 20)
        [4, 5, 6, 8, 9, 10]
        >>> parse_dcd_selection("range(10,20)", 30)
        [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        >>> parse_dcd_selection("", 5)
        [0, 1, 2, 3, 4]
        """
        if not selection_str.strip():
            # Empty string means process all DCDs
            return list(range(num_dcd))
        
        dcd_indices = []
        selection_str = selection_str.strip()
        
        try:
            # Check if it's a Python range expression
            if selection_str.startswith('range(') and selection_str.endswith(')'):
                # Safely evaluate the range expression
                try:
                    # Extract the range parameters
                    range_content = selection_str[6:-1]  # Remove 'range(' and ')'
                    # Split by comma and evaluate
                    parts = [p.strip() for p in range_content.split(',')]
                    
                    if len(parts) == 1:
                        # range(stop)
                        stop = int(parts[0])
                        dcd_indices = list(range(stop))
                    elif len(parts) == 2:
                        # range(start, stop)
                        start = int(parts[0])
                        stop = int(parts[1])
                        dcd_indices = list(range(start, stop))
                    elif len(parts) == 3:
                        # range(start, stop, step)
                        start = int(parts[0])
                        stop = int(parts[1])
                        step = int(parts[2])
                        dcd_indices = list(range(start, stop, step))
                    else:
                        raise ValueError(f"Invalid range format: {selection_str}")
                    
                    # Validate indices are within bounds
                    for idx in dcd_indices:
                        if idx < 0 or idx >= num_dcd:
                            raise ValueError(f"Range {selection_str} contains index {idx} outside valid DCD range (0-{num_dcd-1})")
                    
                except (ValueError, IndexError) as e:
                    raise ValueError(f"Invalid range format: {selection_str}. Use range(start, stop) or range(start, stop, step)")
                
            else:
                # Original parsing logic for ranges like "4-10" or "4-6,8-10"
                ranges = [r.strip() for r in selection_str.split(',')]
                
                for range_str in ranges:
                    if '-' in range_str:
                        # Range like "4-10"
                        start, end = range_str.split('-', 1)
                        start = int(start.strip())
                        end = int(end.strip())
                        
                        if start > end:
                            raise ValueError(f"Invalid range: {range_str} (start > end)")
                        if start < 0 or end >= num_dcd:
                            raise ValueError(f"Range {range_str} is outside valid DCD range (0-{num_dcd-1})")
                        
                        dcd_indices.extend(range(start, end + 1))
                    else:
                        # Single number like "5"
                        index = int(range_str.strip())
                        if index < 0 or index >= num_dcd:
                            raise ValueError(f"DCD index {index} is outside valid range (0-{num_dcd-1})")
                        dcd_indices.append(index)
            
            # Remove duplicates and sort
            dcd_indices = sorted(list(set(dcd_indices)))
            
            if not dcd_indices:
                raise ValueError("No valid DCD indices found")
                
            return dcd_indices
            
        except ValueError as e:
            if "invalid literal for int()" in str(e):
                raise ValueError(f"Invalid DCD selection format: '{selection_str}'. Use formats like '4-10', '4-6,8-10', or 'range(10,20)'")
            else:
                raise e


    # —— Browsers & Save dialogs ——
    def browse_baseDir(self):
        d = filedialog.askdirectory()
        if d:
            self.baseDir_var.set(d)
            self.update_output_path_label()

    def browse_vmd(self):
        f = filedialog.askopenfilename(filetypes=[("Executable","*"), ("All files", "*")])
        if f:
            # last var of step1
            self.coords_vars[-1].set(f)

    def browse_xsc(self):
        f = filedialog.askopenfilename(filetypes=[("XSC files","*.xsc"),("All","*")])
        if f:
            self.unwrap_vars[2].set(f)

    def browse_output_folder(self):
        d = filedialog.askdirectory()
        if d:
            # Store the full path in a new variable
            self.output_full_path = d
            # Still show just the folder name in the UI
            folder_name = os.path.basename(d)
            self.output_folder_var.set(folder_name)
            self.update_output_path_label()
    
    def update_output_path_label(self):
        """Update the label showing the current output path."""
        if hasattr(self, 'output_full_path') and self.output_full_path:
            path_text = f"Output location: {self.output_full_path}"
        else:
            current_dir = os.getcwd()
            path_text = f"Output location: {current_dir}/<folder_name>"
        
        if hasattr(self, 'output_path_label'):
            self.output_path_label.config(text=path_text)
    
    def on_folder_name_changed(self, *args):
        """Clear the browse path when user manually types a folder name."""
        # If the current folder name doesn't match the browse path basename, clear the browse path
        if hasattr(self, 'output_full_path') and self.output_full_path:
            current_folder = self.output_folder_var.get().strip()
            browse_folder = os.path.basename(self.output_full_path)
            if current_folder != browse_folder:
                self.output_full_path = None
                self.update_output_path_label()

    def reset_output_path(self):
        """Reset the output path state (useful after testing or to clear stale paths)."""
        self.output_full_path = None
        self.update_output_path_label()

    def open_calculations_window(self):
        """Open the calculations window for Alpha2/MSD and Dipole calculations"""
        try:
            calc_window = CalculationsWindow(self)
            calc_window.focus()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open calculations window:\n{e}")

    def save_mainfile(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".py", filetypes=[("Python","*.py"),("All","*")]
        )
        if path:
            self.mainfile_var.set(os.path.basename(path))

    def save_submitfile(self):
        path = filedialog.asksaveasfilename(
            defaultextension=".sh", filetypes=[("Shell Script","*.sh"),("All","*")]
        )
        if path:
            self.submitfile_var.set(os.path.basename(path))


    # —— Persistence ——
    def save_config(self):
        data = {
            "baseDir": self.baseDir_var.get(),
            "num_dcd": self.num_dcd_var.get(),
            "num_particles": self.num_particles_var.get(),
            "common_term": self.common_term_var.get(),
            "max_workers": self.max_workers_var.get(),
            "coords": [v.get() for v in self.coords_vars],
            "coords_parallel": self.coords_parallel.get(),
            "coords_dcd_selection": self.coords_dcd_selection.get(),
            "unwrap": [v.get() for v in self.unwrap_vars],
            "unwrap_opt": [v.get() for v in self.unwrap_opt],
            "unwrap_chunk": self.unwrap_chunk_var.get(),
            "unwrap_parallel": self.unwrap_parallel.get(),
            "unwrap_dcd_selection": self.unwrap_dcd_selection.get(),
            "com": [v.get() for v in self.com_vars],
            "com_parallel": self.com_parallel.get(),
            "com_memmap": self.com_memmap.get(),
            "com_dcd_selection": self.com_dcd_selection.get(),
            "sbatch": [v.get() for v in self.sbatch_vars],
            "output_folder": self.output_folder_var.get(),
            "output_full_path": getattr(self, 'output_full_path', None),
            "mainfile": self.mainfile_var.get(),
            "submitfile": self.submitfile_var.get(),
            "traj_file_size_mb": self.traj_file_size_mb.get(),
            "traj_num_frames": self.traj_num_frames.get(),
            "traj_num_atoms": self.traj_num_atoms.get(),
            "available_memory_gb": self.available_memory_gb.get(),
            "coord_precision": self.coord_precision.get()
        }
        with open(CONFIG_PATH, "w") as f:
            json.dump(data, f)
    
    def save_calculation_config(self):
        """Save Analysis calculations window configuration"""
        data = {
            # Alpha2/MSD settings
            "skip_alpha2": self.skip_alpha2.get(),
            "calc_type_var": self.calc_type_var.get(),
            "a2_vars": [v.get() for v in self.a2_vars],
            "a2_chunk_processing": self.a2_chunk_processing.get(),
            "a2_validate": self.a2_validate.get(),
            "alpha2_dcd_selection": self.alpha2_dcd_selection.get(),
            
            # Dipole settings
            "skip_dipole": self.skip_dipole.get(),
            "dipole_calc_type": self.dipole_calc_type.get(),
            "individual_dipole_vars": [v.get() for v in self.individual_dipole_vars],
            "collective_dipole_vars": [v.get() for v in self.collective_dipole_vars],
            "dipole_parallel": self.dipole_parallel.get(),
            "dipole_validate": self.dipole_validate.get(),
            "dipole_chunk_processing": self.dipole_chunk_processing.get(),
            "individual_dipole_dcd_selection": self.individual_dipole_dcd_selection.get(),
            "collective_dipole_dcd_selection": self.collective_dipole_dcd_selection.get(),
            
            # Velocity extraction settings
            "skip_velocity": self.skip_velocity.get(),
            "velocity_vars": [v.get() for v in self.velocity_vars] if hasattr(self, 'velocity_vars') else [],
            "velocity_parallel": self.velocity_parallel.get(),
            "velocity_validate": self.velocity_validate.get(),
            "velocity_dcd_selection": self.velocity_dcd_selection.get(),
            
            # Common parameters
            "common_term_var": self.common_term_var.get(),
            "baseDir_var": self.baseDir_var.get(),
            "num_dcd_var": self.num_dcd_var.get(),
            "num_particles_var": self.num_particles_var.get(),
            "max_workers_var": self.max_workers_var.get(),
            "output_folder_var": self.output_folder_var.get(),
            "mainfile_var": self.mainfile_var.get(),
            "submitfile_var": self.submitfile_var.get(),
            "sbatch_vars": [v.get() for v in self.sbatch_vars] if hasattr(self, 'sbatch_vars') else []
        }
        
        # Save to a separate config file for calculations
        calc_config_path = os.path.expanduser("~/.pipeline_calculations_config.json")
        with open(calc_config_path, "w") as f:
            json.dump(data, f)

    def load_config(self):
        if not os.path.exists(CONFIG_PATH):
            return
        try:
            with open(CONFIG_PATH) as f:
                data = json.load(f)
        except Exception:
            return

        self.baseDir_var.set(data.get("baseDir",""))
        self.num_dcd_var.set(data.get("num_dcd",1))
        self.num_particles_var.set(data.get("num_particles",1000))
        self.common_term_var.set(data.get("common_term", ""))
        self.max_workers_var.set(data.get("max_workers", min(4, mp.cpu_count())))

        for v,val in zip(self.coords_vars, data.get("coords",[])):
            v.set(val)
        self.coords_parallel.set(data.get("coords_parallel", True))
        self.coords_dcd_selection.set(data.get("coords_dcd_selection", ""))

        for v,val in zip(self.unwrap_vars, data.get("unwrap",[])):
            v.set(val)
        for v,val in zip(self.unwrap_opt, data.get("unwrap_opt",[])):
            v.set(val)
        self.unwrap_chunk_var.set(data.get("unwrap_chunk", "auto"))
        self.unwrap_parallel.set(data.get("unwrap_parallel", True))
        self.unwrap_dcd_selection.set(data.get("unwrap_dcd_selection", ""))

        for v,val in zip(self.com_vars, data.get("com",[])):
            v.set(val)
        self.com_parallel.set(data.get("com_parallel", True))
        self.com_memmap.set(data.get("com_memmap", False))
        self.com_dcd_selection.set(data.get("com_dcd_selection", ""))

        for v,val in zip(self.sbatch_vars, data.get("sbatch",[])):
            v.set(val)

        self.output_folder_var.set(data.get("output_folder", "pipeline_run"))
        self.output_full_path = data.get("output_full_path", None)
        self.mainfile_var.set(data.get("mainfile",""))
        self.submitfile_var.set(data.get("submitfile",""))
        
        # Load trajectory characteristics
        self.traj_file_size_mb.set(data.get("traj_file_size_mb", ""))
        self.traj_num_frames.set(data.get("traj_num_frames", ""))
        self.traj_num_atoms.set(data.get("traj_num_atoms", ""))
        self.available_memory_gb.set(data.get("available_memory_gb", ""))
        self.coord_precision.set(data.get("coord_precision", "single"))

        # Note: Skip checkboxes are NOT loaded from config as requested
    
    def load_calculation_config(self):
        """Load Analysis calculations window configuration"""
        calc_config_path = os.path.expanduser("~/.pipeline_calculations_config.json")
        if not os.path.exists(calc_config_path):
            return
        try:
            with open(calc_config_path) as f:
                data = json.load(f)
        except Exception:
            return

        # Alpha2/MSD settings
        self.skip_alpha2.set(data.get("skip_alpha2", False))
        self.calc_type_var.set(data.get("calc_type_var", "alpha2_msd"))
        for v, val in zip(self.a2_vars, data.get("a2_vars", [])):
            v.set(val)
        self.a2_chunk_processing.set(data.get("a2_chunk_processing", True))
        self.a2_validate.set(data.get("a2_validate", True))
        self.alpha2_dcd_selection.set(data.get("alpha2_dcd_selection", ""))
        
        # Dipole settings
        self.skip_dipole.set(data.get("skip_dipole", True))
        self.dipole_calc_type.set(data.get("dipole_calc_type", "individual"))
        for v, val in zip(self.individual_dipole_vars, data.get("individual_dipole_vars", [])):
            v.set(val)
        for v, val in zip(self.collective_dipole_vars, data.get("collective_dipole_vars", [])):
            v.set(val)
        self.dipole_parallel.set(data.get("dipole_parallel", True))
        self.dipole_validate.set(data.get("dipole_validate", True))
        self.dipole_chunk_processing.set(data.get("dipole_chunk_processing", True))
        self.individual_dipole_dcd_selection.set(data.get("individual_dipole_dcd_selection", ""))
        self.collective_dipole_dcd_selection.set(data.get("collective_dipole_dcd_selection", ""))
        
        # Velocity extraction settings
        self.skip_velocity.set(data.get("skip_velocity", True))
        for v, val in zip(self.velocity_vars, data.get("velocity_vars", [])):
            v.set(val)
        self.velocity_parallel.set(data.get("velocity_parallel", True))
        self.velocity_validate.set(data.get("velocity_validate", True))
        self.velocity_dcd_selection.set(data.get("velocity_dcd_selection", ""))
        
        # Common parameters
        self.common_term_var.set(data.get("common_term_var", ""))
        self.baseDir_var.set(data.get("baseDir_var", ""))
        self.num_dcd_var.set(data.get("num_dcd_var", ""))
        self.num_particles_var.set(data.get("num_particles_var", ""))
        self.max_workers_var.set(data.get("max_workers_var", "4"))
        self.output_folder_var.set(data.get("output_folder_var", ""))
        self.mainfile_var.set(data.get("mainfile_var", ""))
        self.submitfile_var.set(data.get("submitfile_var", ""))
        
        # SLURM parameters
        if hasattr(self, 'sbatch_vars'):
            for v, val in zip(self.sbatch_vars, data.get("sbatch_vars", [])):
                v.set(val)

    def generate_mult_run(self):
        """Generate a multi_run.sh script that submits all .sh files in the output directory using sbatch."""
        try:
            # Get output folder path
            output_folder = self.output_folder_var.get().strip()
            if not output_folder:
                raise ValueError("Output folder name is required")
            
            # Use the full path if set by browse AND the folder name matches, otherwise use current directory
            if (hasattr(self, 'output_full_path') and self.output_full_path and 
                os.path.basename(self.output_full_path) == output_folder):
                output_path = self.output_full_path
            else:
                # Check if folder exists in current directory
                output_path = os.path.join(os.getcwd(), output_folder)
            
            if not os.path.exists(output_path):
                raise ValueError(f"Output directory does not exist: {output_path}")
            
            # Find all .sh files in the output directory
            sh_files = []
            for file in os.listdir(output_path):
                if file.endswith('.sh'):
                    sh_files.append(file)
            
            if not sh_files:
                raise ValueError(f"No .sh files found in directory: {output_path}")
            
            # Sort the files for consistent order
            sh_files.sort()
            
            # Generate multi_run.sh content
            multi_run_content = [
                "#!/bin/bash",
                "# Multi-run script generated by MD Analysis Pipeline GUI v2.0",
                "# This script submits all .sh files in the directory using sbatch",
                "",
                f"# Found {len(sh_files)} .sh files to submit:",
            ]
            
            # Add comment listing all files
            for sh_file in sh_files:
                multi_run_content.append(f"# - {sh_file}")
            
            multi_run_content.extend([
                "",
                "echo 'Starting multi-run submission...'",
                "echo 'Submitting jobs:'",
                ""
            ])
            
            # Add sbatch commands for each .sh file
            for sh_file in sh_files:
                multi_run_content.extend([
                    f"echo 'Submitting {sh_file}...'",
                    f"sbatch {sh_file}",
                    ""
                ])
            
            multi_run_content.extend([
                "echo 'All jobs submitted!'",
                "echo 'Check job status with: squeue -u $USER'"
            ])
            
            # Write multi_run.sh file
            multi_run_file = os.path.join(output_path, "multi_run.sh")
            with open(multi_run_file, 'w') as f:
                f.write('\n'.join(multi_run_content))
            
            # Make the script executable
            os.chmod(multi_run_file, 0o755)
            
            # Show success message
            message = f"""Successfully generated multi_run.sh script!

Location: {multi_run_file}

Found {len(sh_files)} .sh files:
{chr(10).join(['• ' + f for f in sh_files])}

Usage:
1. Upload to your supercomputer
2. Run: ./multi_run.sh

The script will submit all .sh files using sbatch commands."""

            messagebox.showinfo("Multi-run Script Generated", message)
                              
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate multi-run script:\n{e}")



    def _copy_compiled_functions(self, src_dir, dest_dir):
        """Copy Python function files"""
        import os
        import shutil
        import glob
        
        # Remove destination directory if it exists
        if os.path.exists(dest_dir):
            shutil.rmtree(dest_dir)
        
        # Create destination directory
        os.makedirs(dest_dir, exist_ok=True)
        
        # Copy .py files
        py_files = glob.glob(os.path.join(src_dir, "*.py"))
        # Exclude __init__.py and other special files
        py_files = [f for f in py_files if not os.path.basename(f).startswith('__')]
        
        if not py_files:
            raise ValueError(f"No Python .py files found in {src_dir}.")
        
        copied_count = 0
        for py_file in py_files:
            filename = os.path.basename(py_file)
            dest_file = os.path.join(dest_dir, filename)
            shutil.copy2(py_file, dest_file)
            copied_count += 1
            print(f"  ✓ Copied Python function: {filename}")
        
        print(f"Successfully copied {copied_count} Python function files")

    # —— Generate driver & submission scripts ——
    def generate_files(self):
        try:
            # --- Gather common inputs ---
            bd = self.baseDir_var.get().strip()
            if not bd:
                raise ValueError("Base Directory is required")
            nd = int(self.num_dcd_var.get())
            max_workers = int(self.max_workers_var.get())
            
            # Get common term for use in script generation
            common_term = self.common_term_var.get().strip()
            # Apply common term expansion to base directory
            bd = self.expand_common_term(bd)
            
            # Create output folder
            output_folder = self.output_folder_var.get().strip()
            if not output_folder:
                raise ValueError("Output folder name is required")
            
            # Use the full path if set by browse AND the folder name matches, otherwise use current directory
            if (hasattr(self, 'output_full_path') and self.output_full_path and 
                os.path.basename(self.output_full_path) == output_folder):
                output_path = self.output_full_path
            else:
                # Create folder in current directory
                output_path = os.path.join(os.getcwd(), output_folder)
            
            # Create the output directory
            os.makedirs(output_path, exist_ok=True)
            
            # Copy main_functions folder to output directory if it doesn't exist
            main_functions_dest = os.path.join(output_path, "main_functions")
            if not os.path.exists(main_functions_dest):
                # Look for main_functions folder with multiple fallback paths
                script_dir = os.path.dirname(os.path.abspath(__file__))
                current_dir = os.getcwd()
                
                # Try multiple potential locations for main_functions
                potential_paths = [
                    os.path.join(script_dir, "main_functions"),          # Same dir as script
                    os.path.join(current_dir, "main_functions"),         # Current working dir
                    os.path.join(os.path.dirname(script_dir), "main_functions"),  # Parent dir
                    "main_functions"  # Relative path
                ]
                
                main_functions_src = None
                for path in potential_paths:
                    if os.path.exists(path):
                        main_functions_src = path
                        break
                
                if main_functions_src:
                    try:
                        # Copy the entire main_functions folder
                        self._copy_compiled_functions(main_functions_src, main_functions_dest)
                        print(f"Copied main_functions folder from {main_functions_src} to {main_functions_dest}")
                    except Exception as e:
                        raise ValueError(f"Failed to copy main_functions folder: {e}")
                else:
                    # Provide detailed error information
                    error_msg = f"main_functions folder not found. Searched in:\n"
                    error_msg += f"  Script directory: {script_dir}\n"
                    error_msg += f"  Current directory: {current_dir}\n"
                    error_msg += f"  Checked paths:\n"
                    for path in potential_paths:
                        exists = "✓" if os.path.exists(path) else "✗"
                        error_msg += f"    {exists} {path}\n"
                    error_msg += f"Please ensure the main_functions folder exists in one of these locations."
                    raise ValueError(error_msg)

            # --- Build the lines of the Python driver ---
            lines = [
                "#!/usr/bin/env python3",
                '"""',
                "MD Analysis Pipeline - Generated by GUI v2.0",
                "This script uses the optimized pipeline functions with file patterns",
                "for enhanced performance, parallel processing, and robust error handling.",
                '"""',
                "",
                "import sys",
                "import time",
                "import os",
                "",
                "def main():",
                "    print('='*60)",
                "    print('MD ANALYSIS PIPELINE')",
                "    print('='*60)",
                "    start_time = time.time()",
                "    all_results = {}",
                ""
            ]

            # Step 1: Coordinate Extraction
            if not self.skip1.get():
                psf_pattern, dcd_pattern, output_pattern, target_selection, vmd = [v.get().strip() for v in self.coords_vars]
                # Apply common term expansion
                psf_pattern, dcd_pattern, output_pattern, target_selection, vmd = [self.expand_common_term(var) for var in [psf_pattern, dcd_pattern, output_pattern, target_selection, vmd]]
                
                use_parallel = self.coords_parallel.get()
                coords_dcd_selection = self.coords_dcd_selection.get().strip()
                
                # Parse DCD selection
                try:
                    coords_dcd_list = self.parse_dcd_selection(coords_dcd_selection, nd)
                except ValueError as e:
                    raise ValueError(f"Invalid DCD selection for coordinates extraction: {e}")
                
                lines.extend([
                    "    # Step 1: Coordinate Extraction",
                    "    print('\\nStep 1: Extracting coordinates from DCD files...')",
                    f"    selected_dcds = {coords_dcd_list}",
                    f"    print(f'Processing DCDs: {{selected_dcds}}')",
                    "    try:",
                    "        from main_functions.coordinates_extract import raw_coords",
                    f"        results_coords = raw_coords(",
                    f"            baseDir={repr(bd)},",
                    f"            psf_pattern={repr(psf_pattern)},",
                    f"            dcd_pattern={repr(dcd_pattern)},",
                    f"            output_pattern={repr(output_pattern)},",
                    f"            num_dcd=len(selected_dcds),",
                    f"            target_selection={repr(target_selection)},",
                    f"            vmd={repr(vmd)},",
                    f"            max_workers={max_workers if use_parallel else 1},",
                    f"            dcd_indices=selected_dcds,",
                    f"            common_term={repr(common_term)}",
                    "        )",
                    "        all_results['coordinates_extract'] = results_coords",
                    f"        expected_files = len(selected_dcds)",
                    f"        if results_coords['successful'] and len(results_coords['successful']) >= expected_files:",
                    "            print('✓ All coordinate extraction completed successfully')",
                    "        else:",
                    f"            print(f'Warning: Only {{len(results_coords.get(\"successful\", []))}} out of {{expected_files}} coordinate files extracted successfully')",
                    f"            if results_coords.get('failed'):",
                    f"                print(f'Failed DCDs: {{results_coords[\"failed\"]}}')",
                    "    except Exception as e:",
                    "        print(f'✗ Coordinate extraction failed: {e}')",
                    "        sys.exit(1)",
                    ""
                ])

            # Step 2: Unwrap Coordinates  
            if not self.skip2.get():
                input_pattern, output_pattern, xsc_pattern, na = [v.get().strip() for v in self.unwrap_vars]
                # Apply common term expansion  
                input_pattern, output_pattern, xsc_pattern = [self.expand_common_term(var) for var in [input_pattern, output_pattern, xsc_pattern]]
                
                iv = self.unwrap_opt[0].get().strip() or "slice(None)"
                st = self.unwrap_opt[1].get().strip() or "1"
                chunk_size = self.unwrap_chunk_var.get().strip()
                use_parallel = self.unwrap_parallel.get()
                unwrap_dcd_selection = self.unwrap_dcd_selection.get().strip()
                
                # Parse DCD selection
                try:
                    unwrap_dcd_list = self.parse_dcd_selection(unwrap_dcd_selection, nd)
                except ValueError as e:
                    raise ValueError(f"Invalid DCD selection for unwrap coordinates: {e}")
                
                chunk_param = "None" if chunk_size == "auto" else chunk_size
                
                lines.extend([
                    "    # Step 2: Unwrap Coordinates",
                    "    print('\\nStep 2: Unwrapping periodic boundary conditions...')",
                    f"    selected_dcds = {unwrap_dcd_list}",
                    f"    print(f'Processing DCDs: {{selected_dcds}}')",
                    "    try:",
                    "        from main_functions.unwrap_coords import unwrapper",
                    f"        results_unwrap = unwrapper(",
                    f"            baseDir={repr(bd)},",
                    f"            input_pattern={repr(input_pattern)},",
                    f"            output_pattern={repr(output_pattern)},",
                    f"            xsc_pattern={repr(xsc_pattern)},",
                    f"            num_dcd=len(selected_dcds),",
                    f"            num_atoms=int({na}),",
                    f"            interval={iv},",
                    f"            stride={st},",
                    f"            max_workers={max_workers if use_parallel else 1},",
                    f"            chunk_size={chunk_param},",
                    f"            dcd_indices=selected_dcds,",
                    f"            common_term={repr(common_term)}",
                    "        )",
                    "        all_results['unwrap_coords'] = results_unwrap",
                    f"        expected_files = len(selected_dcds)",
                    f"        if results_unwrap['success'] >= expected_files:",
                    "            print('✓ All coordinate unwrapping completed successfully')",
                    "        else:",
                    f"            print(f'Warning: Only {{results_unwrap[\"success\"]}} out of {{expected_files}} unwrap files processed successfully')",
                    f"            if results_unwrap.get('failed'):",
                    f"                print(f'Failed DCDs: {{results_unwrap[\"failed\"]}}')",
                    "    except Exception as e:",
                    "        print(f'✗ Coordinate unwrapping failed: {e}')",
                    "        sys.exit(1)",
                    ""
                ])

            # Step 3: COM Calculation
            if not self.skip3.get():
                input_pattern, output_pattern, ap, ml = [v.get().strip() for v in self.com_vars]
                # Apply common term expansion
                input_pattern, output_pattern = [self.expand_common_term(var) for var in [input_pattern, output_pattern]]
                
                use_memmap = self.com_memmap.get()
                num_particles = int(self.num_particles_var.get())
                com_dcd_selection = self.com_dcd_selection.get().strip()
                
                # Parse DCD selection
                try:
                    com_dcd_list = self.parse_dcd_selection(com_dcd_selection, nd)
                except ValueError as e:
                    raise ValueError(f"Invalid DCD selection for COM calculation: {e}")
                
                lines.extend([
                    "    # Step 3: Center-of-Mass Calculation", 
                    "    print('\\nStep 3: Computing center-of-mass coordinates...')",
                    f"    selected_dcds = {com_dcd_list}",
                    f"    print(f'Processing DCDs: {{selected_dcds}}')",
                    "    try:",
                    "        from main_functions.COM_calc import coms",
                    f"        results_com = coms(",
                    f"            baseDir={repr(bd)},",
                                         f"            input_pattern={repr(input_pattern)},",
                     f"            output_pattern={repr(output_pattern)},",
                    f"            num_dcd=len(selected_dcds),",
                    f"            prtcl_num={num_particles},",
                    f"            prtcl_atoms=int({ap}),",
                    f"            particl_mass={ml.split(',')},",
                    f"            max_workers={max_workers if self.com_parallel.get() else 1},",
                    f"            use_memmap={use_memmap},",
                    f"            dcd_indices=selected_dcds,",
                    f"            common_term={repr(common_term)}",
                    "        )",
                    "        all_results['COM_calc'] = results_com",
                    f"        expected_files = len(selected_dcds)",
                    f"        if results_com['success'] >= expected_files:",
                    "            print('✓ All COM calculations completed successfully')",
                    "        else:",
                    f"            print(f'Warning: Only {{results_com[\"success\"]}} out of {{expected_files}} COM files processed successfully')",
                    f"            if results_com.get('failed'):",
                    f"                print(f'Failed DCDs: {{results_com[\"failed\"]}}')",
                    "    except Exception as e:",
                    "        print(f'✗ COM calculation failed: {e}')",
                    "        sys.exit(1)",
                    ""
                ])

            lines.extend([
                "    # Summary",
                "    total_time = time.time() - start_time",
                "    print('\\n' + '='*60)",
                "    print('PIPELINE EXECUTION SUMMARY')",
                "    print('='*60)",
                "    for step, results in all_results.items():",
                "        if 'total_time' in results:",
                "            print(f'{step:20s}: {results[\"total_time\"]:.2f}s')",
                "    print(f'{\"Total time\":20s}: {total_time:.2f}s')",
                "    print('='*60)",
                "    print('Pipeline completed successfully!')",
                "",
                "if __name__ == '__main__':",
                "    main()"
            ])

            # --- Write the main .py file ---
            main_base = self.mainfile_var.get().strip()
            main_base = self.expand_common_term(main_base)
            if not main_base:
                main_base = "main_analysis_pipeline"
            
            main_fn = main_base if main_base.endswith(".py") else main_base + ".py"
            main_full_path = os.path.join(output_path, main_fn)

            with open(main_full_path, "w") as f:
                f.write("\n".join(lines))

            # --- Write the submission .sh file ---
            sub_base = self.submitfile_var.get().strip()
            sub_base = self.expand_common_term(sub_base)
            if not sub_base:
                sub_base = "submit_analysis_pipeline"
            
            sub_fn = sub_base if sub_base.endswith(".sh") else sub_base + ".sh"
            sub_full_path = os.path.join(output_path, sub_fn)

            sb = [v.get().strip() for v in self.sbatch_vars]
            # sb order: Nodes, Partition, QOS, CPUs, Tasks, Memory (GB), Walltime, Output prefix, Email, Module
            
            # Auto-determine SLURM output prefix based on common term vs output folder
            common_term = self.common_term_var.get().strip()
            output_folder = self.output_folder_var.get().strip()
            
            if common_term:
                slurm_output_name = f"sbatch_{common_term}"
            elif output_folder:
                slurm_output_name = f"sbatch_{output_folder}"
            else:
                # Fallback to submit file name (without .sh)
                slurm_output_name = sub_base[:-3] if sub_base.endswith('.sh') else sub_base

            with open(sub_full_path, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("# Generated by Optimized MD Analysis Pipeline GUI v2.0\n")
                f.write("# This script uses the optimized pipeline with file patterns\n\n")
                f.write(f"#SBATCH -N {sb[0]}\n")
                f.write(f"#SBATCH -p {sb[1]}\n")
                f.write(f"#SBATCH -q {sb[2]}\n")
                f.write(f"#SBATCH -c {sb[3]}\n")
                f.write(f"#SBATCH -n {sb[4]}\n")
                if sb[5].strip():  # Only add memory if specified
                    f.write(f"#SBATCH --mem={sb[5]}G\n")
                f.write(f"#SBATCH -t {sb[6]}\n")
                f.write(f"#SBATCH -o {slurm_output_name}.log\n")
                f.write("#SBATCH --mail-type=ALL\n")
                email = sb[8]
                if "@" not in email:
                    email = email + "@asu.edu"
                f.write(f"#SBATCH --mail-user={email}\n")
                f.write("#SBATCH --export=NONE\n\n")
                
                # Add module loading section
                if len(sb) > 9 and sb[9]:  # Check if module parameter exists and is not empty
                    f.write("# Load required module\n")
                    f.write(f"module load {sb[9]}\n\n")
                else:
                    f.write("# Load any required modules (modify as needed)\n")
                    f.write("# module load python/3.8\n")
                    f.write("# module load vmd\n\n")
                f.write("########################################################\n")
                f.write("# Optimized MD Analysis Pipeline Execution\n")
                f.write("########################################################\n\n")
                f.write(f"echo 'Starting optimized pipeline with {max_workers} workers...'\n")
                f.write("echo 'Job started at:' $(date)\n")
                f.write(f"python {main_fn}\n")
                f.write("echo 'Job completed at:' $(date)\n")

            # --- Save & notify ---
            self.save_config()
            
            message = f"""Successfully generated optimized pipeline files in folder:

Output Folder: {output_path}
Main Script: {main_fn}
SLURM Script: {sub_fn}

Key Features:
• All files organized in a single folder: {output_folder}
• main_functions folder included for portability
• Uses new file pattern interface (no more INdir/OUTdir)
• Parallel processing with {max_workers} workers
• Enhanced error handling and progress reporting
• Memory-efficient processing for large datasets
• Simplified and cleaner code structure

The generated scripts use the updated pipeline functions
with file patterns for improved flexibility and performance."""

            messagebox.showinfo("Success", message)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate files:\n{e}")

    def load_application_icon(self):
        """Load the application icon with Ubuntu/built-app compatibility.
        
        Comprehensive approach for Ubuntu taskbar icon visibility:
        1. System icon installation for built apps
        2. Multiple icon setting methods simultaneously  
        3. Desktop file creation for proper integration
        4. Fallback methods for different Ubuntu configurations
        """
        # Get the directory where the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Try different icon files in order of preference (PNG first for Ubuntu compatibility)
        icon_files = [
            "icon.png",   # PNG format (works best on Ubuntu with PIL)
            "icon.ico",   # Windows ICO format - fallback
            "icon.gif",   # GIF format (native Tkinter support)
            "pipeline_icon.png", 
            "pipeline_icon.ico",
            "pipeline_icon.gif",
            "app_icon.png",
            "app_icon.ico",
            "app_icon.gif"
        ]
        
        icon_loaded = False
        
        # First, try to install icon to system locations for built apps
        self.install_system_icon(script_dir)
        
        for icon_file in icon_files:
            icon_path = os.path.join(script_dir, icon_file)
            
            if os.path.exists(icon_path):
                try:
                    if icon_file.endswith('.ico'):
                        # Try to load ICO file - but modern ICO files often fail on Linux
                        ico_success = False
                        try:
                            self.iconbitmap(icon_path)
                            print(f"✓ Loaded ICO icon with iconbitmap: {icon_file}")
                            ico_success = True
                        except Exception as ico_err:
                            print(f"iconbitmap failed for {icon_file}: {ico_err}")
                        
                        # If iconbitmap failed, try to convert ICO to usable format
                        if not ico_success:
                            try:
                                from PIL import Image, ImageTk
                                # Load the ICO and convert to standard sizes
                                image = Image.open(icon_path)
                                
                                # Create a simple 32x32 ICO that Linux Tkinter can handle
                                simple_ico_path = icon_path + ".simple.ico"
                                simple_image = image.resize((32, 32), Image.Resampling.LANCZOS)
                                # Convert to RGB to avoid transparency issues
                                if simple_image.mode == 'RGBA':
                                    # Create white background and composite
                                    background = Image.new('RGB', simple_image.size, (255, 255, 255))
                                    background.paste(simple_image, mask=simple_image.split()[-1])
                                    simple_image = background
                                simple_image.save(simple_ico_path, format='ICO', sizes=[(32, 32)])
                                
                                # Try the simplified ICO
                                try:
                                    self.iconbitmap(simple_ico_path)
                                    os.unlink(simple_ico_path)  # Clean up
                                    print(f"✓ Loaded ICO icon with simplified conversion: {icon_file}")
                                    ico_success = True
                                except Exception as simple_err:
                                    try:
                                        os.unlink(simple_ico_path)  # Clean up failed file
                                    except:
                                        pass
                                    print(f"Simplified ICO also failed: {simple_err}")
                                
                                # If ICO methods fail, try iconphoto with the original image
                                if not ico_success:
                                    try:
                                        resized = image.resize((48, 48), Image.Resampling.LANCZOS)
                                        photo = ImageTk.PhotoImage(resized)
                                        self.iconphoto(True, photo)  # type: ignore
                                        self._icon_photo = photo
                                        print(f"✓ Loaded ICO as PNG icon: {icon_file}")
                                        ico_success = True
                                    except Exception as png_err:
                                        print(f"PNG conversion from ICO failed: {png_err}")
                                        
                            except Exception as pil_err:
                                print(f"PIL processing of ICO failed: {pil_err}")
                        
                        if ico_success:
                            icon_loaded = True
                            break
                    elif icon_file.endswith(('.png', '.gif')):
                                                    # Ubuntu-optimized PNG/GIF loading - exact diagnostic approach
                        try:
                            # Try PIL method exactly as it worked in diagnostic
                            from PIL import Image, ImageTk
                            image = Image.open(icon_path)
                            
                            # Create a simple window test first to ensure compatibility
                            test_window = tk.Toplevel(self)
                            test_window.withdraw()  # Hide it immediately
                            
                            try:
                                # Test the method that worked in diagnostic
                                resized = image.resize((32, 32), Image.Resampling.LANCZOS)
                                photo = ImageTk.PhotoImage(resized)
                                
                                # Test on the hidden window first
                                test_window.iconphoto(True, photo)  # type: ignore
                                test_window.destroy()  # Clean up test window
                                
                                # Apply multiple icon methods simultaneously for maximum compatibility
                                success_methods = []
                                
                                # Method 1: Standard iconphoto
                                try:
                                    self.iconphoto(True, photo)  # type: ignore
                                    success_methods.append("iconphoto")
                                except:
                                    pass
                                
                                # Method 2: Multiple sizes for different Ubuntu configurations
                                try:
                                    sizes = [16, 22, 24, 32, 48, 64]
                                    for size in sizes:
                                        sized_img = image.resize((size, size), Image.Resampling.LANCZOS)
                                        sized_photo = ImageTk.PhotoImage(sized_img)
                                        self.iconphoto(True, sized_photo)  # type: ignore
                                    success_methods.append("multi-size")
                                except:
                                    pass
                                
                                # Method 3: Try iconbitmap with converted icon for built apps
                                try:
                                    import tempfile
                                    temp_ico = tempfile.NamedTemporaryFile(suffix='.ico', delete=False)
                                    temp_ico.close()
                                    
                                    # Convert to simple ICO for system compatibility
                                    ico_image = image.resize((32, 32), Image.Resampling.LANCZOS)
                                    if ico_image.mode != 'RGB':
                                        ico_image = ico_image.convert('RGB')
                                    ico_image.save(temp_ico.name, 'ICO')
                                    
                                    self.iconbitmap(temp_ico.name)
                                    os.unlink(temp_ico.name)
                                    success_methods.append("iconbitmap")
                                except:
                                    pass
                                
                                # Method 4: System integration method
                                try:
                                    # This relies on the system icon we installed
                                    self.wm_iconname("md-analysis-pipeline")  # type: ignore
                                    success_methods.append("system-icon")
                                except:
                                    pass
                                
                                self._icon_photo = photo
                                methods_str = "+".join(success_methods) if success_methods else "fallback"
                                print(f"✓ Icon loaded ({methods_str}): {icon_file}")
                                icon_loaded = True
                                break
                                
                            except Exception as test_err:
                                test_window.destroy()  # Clean up test window
                                print(f"PIL iconphoto test failed: {test_err}")
                                # Try direct approach without test
                                try:
                                    resized = image.resize((48, 48), Image.Resampling.LANCZOS)
                                    photo = ImageTk.PhotoImage(resized)
                                    self._icon_photo = photo  # Store reference first
                                    print(f"✓ Icon loaded (PIL fallback mode): {icon_file}")
                                    icon_loaded = True
                                    break
                                except Exception as direct_err:
                                    print(f"PIL direct approach failed: {direct_err}")
                                    # Continue to native Tkinter
                                
                        except ImportError:
                            # Fallback to native Tkinter (diagnostic showed this works for PNG!)
                            if icon_file.endswith('.png'):
                                # Native PNG support works on Ubuntu as shown by diagnostic
                                try:
                                    photo = tk.PhotoImage(file=icon_path)
                                    self.iconphoto(True, photo)  # type: ignore
                                    self._icon_photo = photo
                                    print(f"✓ Loaded PNG icon with native Tkinter: {icon_file}")
                                    icon_loaded = True
                                    break
                                except Exception as png_native_err:
                                    print(f"Native PNG loading failed: {png_native_err}")
                            elif icon_file.endswith('.gif'):
                                try:
                                    photo = tk.PhotoImage(file=icon_path)
                                    self.iconphoto(True, photo)  # type: ignore
                                    # Keep a reference to prevent garbage collection
                                    self._icon_photo = photo
                                    print(f"✓ Loaded GIF icon: {icon_file}")
                                    icon_loaded = True
                                    break
                                except Exception as gif_err:
                                    print(f"GIF loading failed: {gif_err}")
                        except Exception as e:
                            print(f"Failed to load {icon_file}: {e}")
                            continue
                except Exception as e:
                    print(f"Failed to load {icon_file}: {e}")
                    continue
        
        if not icon_loaded:
            # Create a simple default icon programmatically
            try:
                self.create_default_icon()
                print("✓ Created default application icon")
            except Exception as e:
                print(f"Could not create default icon: {e}")
                print("Application will run without a custom icon")
        
        # Additional Linux compatibility measures
        self.setup_linux_window_properties()
        
        # Ubuntu/Linux desktop feedback
        if not icon_loaded:
            print("Note: Icon loading failed. For Ubuntu desktop, try:")
            print("1. sudo apt-get install python3-pil python3-tk")
            print("2. Right-click taskbar icon → 'Keep in Launcher'")
            print("3. Check ~/.local/share/applications/ for desktop file")
        else:
            print("📌 Ubuntu tip: Right-click taskbar icon → 'Keep in Launcher' for permanent access")
            print("📁 Desktop file installed for system integration")
    
    def setup_linux_window_properties(self):
        """Set up Unity/Ubuntu-specific window properties for taskbar integration."""
        try:
            # Unity-specific window class - very important for taskbar recognition
            try:
                self.wm_class("md-analysis-pipeline", "MD-Analysis-Pipeline")  # type: ignore
            except AttributeError:
                pass
            
            # Set window name and icon name for Unity
            try:
                self.wm_iconname("MD Analysis Pipeline")  # type: ignore
                self.title("MD Dynamics Analysis Pipeline")  # Ensure title is set
            except AttributeError:
                pass
            
            # Unity-specific window manager hints
            try:
                # Set window type hint for Unity
                self.wm_attributes('-type', 'normal')  # type: ignore
            except (AttributeError, tk.TclError):
                pass
                
            try:
                # Set Unity-specific properties
                self.wm_attributes('-topmost', False)  # type: ignore
                self.wm_attributes('-zoomed', False)   # type: ignore
            except (AttributeError, tk.TclError):
                pass
            
            # Create a temporary .desktop file for proper Unity integration
            self.create_unity_desktop_entry()
            
            # Force window manager recognition
            self.update_idletasks()
            self.lift()
            self.focus_force()
            
            # Additional Unity compatibility
            try:
                self.wm_protocol("WM_DELETE_WINDOW", self.on_closing)
            except:
                pass
            
        except Exception as e:
            # If any setup fails, continue silently
            pass
    
    def create_unity_desktop_entry(self):
        """Create a temporary .desktop file for Unity taskbar integration."""
        try:
            import tempfile
            import os
            
            # Create a temporary .desktop file that Unity can recognize
            desktop_content = f"""[Desktop Entry]
Version=1.0
Type=Application
Name=MD Analysis Pipeline
Comment=Molecular Dynamics Analysis Pipeline
Exec=python3 {os.path.abspath(__file__)}
Icon={os.path.abspath("icon.png") if os.path.exists("icon.png") else "applications-science"}
Terminal=false
Categories=Science;Education;
StartupWMClass=MD-Analysis-Pipeline
"""
            
            # Write to a temporary location that Unity might check
            temp_desktop = os.path.expanduser("~/.local/share/applications/md-analysis-pipeline-temp.desktop")
            os.makedirs(os.path.dirname(temp_desktop), exist_ok=True)
            
            with open(temp_desktop, 'w') as f:
                f.write(desktop_content)
            
            # Store the path for cleanup
            self.temp_desktop_file = temp_desktop
            
        except Exception:
            # If desktop file creation fails, continue without it
            pass
    
    def calculate_optimal_settings(self):
        """Calculate optimal processing settings based on trajectory characteristics"""
        try:
            # Get trajectory characteristics
            file_size_mb = float(self.traj_file_size_mb.get().strip()) if self.traj_file_size_mb.get().strip() else 0
            num_frames = int(self.traj_num_frames.get().strip()) if self.traj_num_frames.get().strip() else 0
            num_atoms_extracted = int(self.traj_num_atoms.get().strip()) if self.traj_num_atoms.get().strip() else 0
            available_memory_gb = float(self.available_memory_gb.get().strip()) if self.available_memory_gb.get().strip() else 0
            precision = self.coord_precision.get()
            max_workers = int(self.max_workers_var.get())  # Get user-specified max workers
            
            if not all([file_size_mb, num_frames, num_atoms_extracted, available_memory_gb]):
                messagebox.showwarning("Missing Information", 
                                     "Please fill in all trajectory characteristics fields for optimal calculations.")
                return
            
            # Calculate memory requirements and optimal settings
            recommendations = self._calculate_memory_optimization(file_size_mb, num_frames, num_atoms_extracted, available_memory_gb, precision, max_workers)
            
            # Update GUI with recommendations (chunk size only, since max_workers is user-controlled)
            self._apply_recommendations(recommendations)
            
            # Display results
            results_text = f"✓ Chunk size optimized for {recommendations['max_workers']} workers: {recommendations['chunk_size']:,} frames, SLURM Mem: {recommendations['slurm_memory_gb']}GB"
            self.optimization_results_label.config(text=results_text, fg='#27ae60')
            
            # Show detailed popup
            self._show_optimization_details(recommendations)
            
        except ValueError as e:
            messagebox.showerror("Invalid Input", f"Please enter valid numbers for all fields.\nError: {e}")
        except Exception as e:
            messagebox.showerror("Calculation Error", f"Error calculating optimal settings:\n{e}")
    
    def _calculate_memory_optimization(self, file_size_mb, num_frames, num_atoms_extracted, available_memory_gb, precision, max_workers):
        """Calculate optimal memory settings based on trajectory characteristics and real-world data"""
        
        # Real-world calibration constants (based on user data: 24GB DCD → 42GB output, 200k frames, 3k atoms)
        REFERENCE_DCD_SIZE_GB = 24.0
        REFERENCE_OUTPUT_SIZE_GB = 42.0
        REFERENCE_FRAMES = 200000
        REFERENCE_ATOMS = 3000
        DCD_TO_OUTPUT_RATIO = REFERENCE_OUTPUT_SIZE_GB / REFERENCE_DCD_SIZE_GB  # 1.75
        
        # Precision factors
        BYTES_PER_COORD_DOUBLE = 8  # Double precision
        BYTES_PER_COORD_SINGLE = 4  # Single precision
        bytes_per_coord = BYTES_PER_COORD_SINGLE if precision == "single" else BYTES_PER_COORD_DOUBLE
        precision_factor = 0.5 if precision == "single" else 1.0
        
        # System overhead constants
        VMD_OVERHEAD_FACTOR = 2.5  # VMD memory overhead during processing
        PYTHON_OVERHEAD_FACTOR = 1.8  # Python/NumPy object overhead
        SAFETY_FACTOR = 0.75  # Use 75% of available memory for safety
        BASE_MEMORY_GB = 3  # Base system overhead
        
        # Calculate expected output file size after Step 1 (coordinate extraction)
        # Scale from reference data and adjust for precision
        scaling_factor = (num_frames / REFERENCE_FRAMES) * (num_atoms_extracted / REFERENCE_ATOMS)
        expected_output_size_gb = (file_size_mb / 1024) * DCD_TO_OUTPUT_RATIO * scaling_factor * precision_factor
        expected_output_size_mb = expected_output_size_gb * 1024
        
        # Calculate memory per frame (more accurate based on actual output scaling)
        memory_per_frame_mb = (expected_output_size_mb / num_frames) * PYTHON_OVERHEAD_FACTOR
        
        # Calculate VMD memory usage during Step 1
        vmd_memory_per_file_mb = file_size_mb * VMD_OVERHEAD_FACTOR
        
        # Available memory for processing (in MB)
        available_memory_mb = (available_memory_gb - BASE_MEMORY_GB) * 1024 * SAFETY_FACTOR
        
        # Calculate optimal chunk size for Step 2 (unwrapping) based on user-specified max_workers
        # This is the most memory-intensive step
        # Each worker needs memory for: chunk processing + file I/O overhead + working memory
        available_memory_per_worker_mb = available_memory_mb / max_workers
        memory_overhead_factor = 1.8  # 1.8x for processing overhead
        
        # Calculate optimal chunk size that fits within memory constraints
        max_chunk_size_by_memory = int(available_memory_per_worker_mb / (memory_per_frame_mb * memory_overhead_factor))
        optimal_chunk_size = max(1000, min(max_chunk_size_by_memory, num_frames))  # Clamp between 1000 and num_frames
        
        # Use the user-specified number of workers
        optimal_workers = max_workers
        
        # Calculate actual memory per worker with optimized chunk size
        memory_per_worker_mb = optimal_chunk_size * memory_per_frame_mb * memory_overhead_factor
        
        # Calculate recommended SLURM memory (with buffer for all steps)
        # Account for VMD memory + parallel processing + output files
        peak_memory_step1_gb = (vmd_memory_per_file_mb + expected_output_size_mb) / 1024
        peak_memory_step2_gb = (optimal_workers * memory_per_worker_mb + expected_output_size_mb) / 1024
        peak_memory_step3_gb = expected_output_size_gb * 0.8  # COM calculation is lighter
        
        max_memory_needed_gb = max(peak_memory_step1_gb, peak_memory_step2_gb, peak_memory_step3_gb)
        total_estimated_memory_gb = max_memory_needed_gb + BASE_MEMORY_GB
        recommended_slurm_memory = max(16, int(total_estimated_memory_gb * 1.4))  # 40% buffer, minimum 16GB
        
        # Calculate processing time estimates (improved)
        estimated_vmd_time_hours = (file_size_mb / 1024) * 0.5  # ~30 minutes per GB of DCD
        estimated_unwrap_time_hours = (num_frames * num_atoms_extracted) / (600000 * optimal_workers)  # Empirical rate
        estimated_com_time_hours = num_frames / (100000 * optimal_workers)  # COM is faster
        total_time_per_file_hours = estimated_vmd_time_hours + estimated_unwrap_time_hours + estimated_com_time_hours
        
        # Determine recommended DCD batch size based on memory and time
        if total_estimated_memory_gb > available_memory_gb * 0.8:  # Close to memory limit
            recommended_batch_size = 1
        elif total_time_per_file_hours > 8:  # Very long processing time
            recommended_batch_size = 2
        elif total_time_per_file_hours > 4:  # Long processing time
            recommended_batch_size = 3
        elif available_memory_gb < 100:  # Low memory systems
            recommended_batch_size = 3
        else:
            recommended_batch_size = 4
        
        return {
            'chunk_size': optimal_chunk_size,
            'max_workers': optimal_workers,
            'slurm_memory_gb': recommended_slurm_memory,
            'recommended_batch_size': recommended_batch_size,
            'memory_per_frame_mb': memory_per_frame_mb,
            'memory_per_worker_mb': memory_per_worker_mb,
            'estimated_time_per_file_hours': total_time_per_file_hours,
            'vmd_memory_mb': vmd_memory_per_file_mb,
            'total_estimated_memory_gb': total_estimated_memory_gb,
            'expected_output_size_gb': expected_output_size_gb,
            'precision': precision,
            'precision_factor': precision_factor,
            'peak_memory_step1_gb': peak_memory_step1_gb,
            'peak_memory_step2_gb': peak_memory_step2_gb,
            'peak_memory_step3_gb': peak_memory_step3_gb
        }
    
    def _apply_recommendations(self, recommendations):
        """Apply calculated recommendations to GUI fields"""
        # Update chunk size (max_workers is now user-controlled in GUI)
        self.unwrap_chunk_var.set(str(recommendations['chunk_size']))
        
        # Update SLURM memory if SLURM fields exist
        if hasattr(self, 'sbatch_vars') and len(self.sbatch_vars) > 5:
            self.sbatch_vars[5].set(str(recommendations['slurm_memory_gb']))  # Memory field
        
        # Update CPU count in SLURM to match user-specified workers
        if hasattr(self, 'sbatch_vars') and len(self.sbatch_vars) > 3:
            self.sbatch_vars[3].set(str(recommendations['max_workers']))  # CPUs field
    
    def _show_optimization_details(self, recommendations):
        """Show detailed optimization results in a popup window"""
        details_window = tk.Toplevel(self)
        details_window.title("🧠 Optimization Results")
        details_window.geometry("600x500")
        details_window.configure(bg='#f0f0f0')
        details_window.transient(self)
        details_window.grab_set()
        
        # Header
        header_frame = tk.Frame(details_window, bg='#3498db')
        header_frame.pack(fill="x")
        tk.Label(header_frame, text="🧠 Intelligent Optimization Results", 
                font=("Arial", 14, "bold"), fg='white', bg='#3498db', pady=10).pack()
        
        # Content frame with scrolling
        content_frame = tk.Frame(details_window, bg='#f0f0f0')
        content_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        # Results text
        precision_text = f"{recommendations['precision']} precision" + (f" (~{(1-recommendations['precision_factor'])*100:.0f}% smaller files)" if recommendations['precision_factor'] < 1.0 else "")
        
        results_text = f"""📊 MEMORY ANALYSIS:
Expected output file size: {recommendations['expected_output_size_gb']:.1f} GB ({precision_text})
Memory per frame: {recommendations['memory_per_frame_mb']:.2f} MB
Memory per worker: {recommendations['memory_per_worker_mb']:.1f} MB
VMD memory estimate: {recommendations['vmd_memory_mb']:.1f} MB

📈 PEAK MEMORY BY STEP:
Step 1 (VMD extraction): {recommendations['peak_memory_step1_gb']:.1f} GB
Step 2 (Unwrapping): {recommendations['peak_memory_step2_gb']:.1f} GB
Step 3 (COM calculation): {recommendations['peak_memory_step3_gb']:.1f} GB
Total estimated: {recommendations['total_estimated_memory_gb']:.1f} GB

⚙️ RECOMMENDED SETTINGS:
Optimal chunk size: {recommendations['chunk_size']:,} frames
Max workers: {recommendations['max_workers']}
SLURM memory request: {recommendations['slurm_memory_gb']} GB
Coordinate precision: {recommendations['precision']}

⏱️ TIME ESTIMATES:
Estimated time per DCD: {recommendations['estimated_time_per_file_hours']:.1f} hours
Recommended batch size: {recommendations['recommended_batch_size']} DCDs at once

💡 OPTIMIZATION TIPS:
• Process {recommendations['recommended_batch_size']} DCDs at a time using DCD Selection
• Use {recommendations['precision']} precision for {"~50% smaller" if recommendations['precision'] == "single" else "full precision"} files
• Step 2 (unwrapping) is most memory-intensive - chunk size optimized for this
• Monitor memory usage during first run to fine-tune settings
• Use memory mapping for Step 3 (COM calculation)

🎯 APPLIED TO GUI:
✓ Chunk Size optimized to {recommendations['chunk_size']:,} frames (based on {recommendations['max_workers']} workers)
✓ SLURM Memory updated to {recommendations['slurm_memory_gb']} GB
✓ SLURM CPUs updated to {recommendations['max_workers']}
✓ Coordinate Precision: {recommendations['precision']}

💡 NOTE: Max Workers ({recommendations['max_workers']}) is user-controlled in Trajectory Characteristics.
Chunk size is optimized based on your worker setting for best memory efficiency."""

        text_widget = tk.Text(content_frame, wrap=tk.WORD, font=("Consolas", 10), 
                             bg='white', fg='#2c3e50', relief='sunken', borderwidth=2)
        text_widget.pack(fill="both", expand=True)
        text_widget.insert("1.0", results_text)
        text_widget.config(state='disabled')
        
        # Close button
        tk.Button(details_window, text="Close", command=details_window.destroy,
                 bg="#e74c3c", fg="white", font=("Arial", 12, "bold"), 
                 padx=20, pady=5).pack(pady=10)
        
        # Center on parent
        details_window.update_idletasks()
        x = self.winfo_x() + (self.winfo_width() - details_window.winfo_width()) // 2
        y = self.winfo_y() + (self.winfo_height() - details_window.winfo_height()) // 2
        details_window.geometry(f"+{x}+{y}")

    def on_closing(self):
        """Clean up when window is closed."""
        try:
            # Clean up temporary desktop file
            if hasattr(self, 'temp_desktop_file') and os.path.exists(self.temp_desktop_file):
                os.unlink(self.temp_desktop_file)
        except:
            pass
        
        # Close the application
        self.destroy()
    
    def install_system_icon(self, script_dir):
        """Install icon to system locations for Ubuntu taskbar recognition."""
        try:
            icon_installed = False
            
            # Find the best icon file
            for icon_name in ["icon.png", "icon.ico"]:
                icon_path = os.path.join(script_dir, icon_name)
                if os.path.exists(icon_path):
                    
                    # Install to user's local icon directory
                    try:
                        local_icons_dir = os.path.expanduser("~/.local/share/icons/hicolor")
                        sizes = ["16x16", "22x22", "24x24", "32x32", "48x48", "64x64", "128x128"]
                        
                        from PIL import Image
                        original_image = Image.open(icon_path)
                        
                        for size_str in sizes:
                            size = int(size_str.split('x')[0])
                            size_dir = os.path.join(local_icons_dir, size_str, "apps")
                            os.makedirs(size_dir, exist_ok=True)
                            
                            # Create resized icon
                            resized = original_image.resize((size, size), Image.Resampling.LANCZOS)
                            icon_dest = os.path.join(size_dir, "md-analysis-pipeline.png")
                            resized.save(icon_dest, "PNG")
                        
                        # Update icon cache
                        try:
                            import subprocess
                            subprocess.run(["gtk-update-icon-cache", local_icons_dir], 
                                         capture_output=True, timeout=5)
                        except:
                            pass
                        
                        print(f"✓ System icon installed from {icon_name}")
                        icon_installed = True
                        break
                        
                    except Exception as install_err:
                        print(f"Icon installation failed: {install_err}")
                        continue
            
            # Create persistent desktop file
            if icon_installed:
                self.create_persistent_desktop_file(script_dir)
                
        except Exception as e:
            # If system installation fails, continue with normal icon loading
            pass
    
    def create_persistent_desktop_file(self, script_dir):
        """Create a persistent .desktop file for proper Ubuntu integration."""
        try:
            desktop_dir = os.path.expanduser("~/.local/share/applications")
            os.makedirs(desktop_dir, exist_ok=True)
            
            desktop_file = os.path.join(desktop_dir, "md-analysis-pipeline.desktop")
            
            # Determine the executable path (script or built app)
            if getattr(sys, 'frozen', False):
                # Built application
                exec_path = sys.executable
            else:
                # Python script
                exec_path = f"python3 {os.path.abspath(__file__)}"
            
            desktop_content = f"""[Desktop Entry]
Version=1.0
Type=Application
Name=MD Analysis Pipeline
Comment=Molecular Dynamics Analysis Pipeline
Exec={exec_path}
Icon=md-analysis-pipeline
Terminal=false
Categories=Science;Education;
StartupWMClass=MD-Analysis-Pipeline
MimeType=application/x-md-analysis;
Keywords=molecular;dynamics;analysis;MSD;science;
"""
            
            with open(desktop_file, 'w') as f:
                f.write(desktop_content)
            
            # Make executable
            os.chmod(desktop_file, 0o755)
            
            # Update desktop database
            try:
                import subprocess
                subprocess.run(["update-desktop-database", desktop_dir], 
                             capture_output=True, timeout=5)
            except:
                pass
                
            print(f"✓ Desktop file created: {desktop_file}")
            
        except Exception as e:
            # If desktop file creation fails, continue
            pass
    
    def create_default_icon(self):
        """Create a simple default icon using Tkinter with multiple sizes for better compatibility."""
        try:
            # Create multiple icon sizes for better cross-platform support
            sizes = [16, 32, 48, 64]
            photos = []
            
            for size in sizes:
                # Create icon with molecular theme
                photo = tk.PhotoImage(width=size, height=size)
                
                # Fill background with gradient-like effect
                center = size // 2
                for x in range(size):
                    for y in range(size):
                        # Create a subtle gradient effect
                        distance = ((x - center) ** 2 + (y - center) ** 2) ** 0.5
                        if distance < center * 0.8:
                            photo.put("#3498db", (x, y))  # Light blue
                        else:
                            photo.put("#2980b9", (x, y))  # Darker blue
                
                # Add molecular structure scaled to size
                atom_size = max(1, size // 16)
                bond_thickness = max(1, size // 32)
                
                # Central molecule (larger)
                central_size = max(2, size // 8)
                for dx in range(central_size):
                    for dy in range(central_size):
                        photo.put("#ffffff", (center - central_size//2 + dx, center - central_size//2 + dy))
                
                # Surrounding atoms (scaled positions)
                atom_positions = [
                    (center - size//4, center - size//4),  # Top-left
                    (center + size//4, center - size//4),  # Top-right
                    (center - size//4, center + size//4),  # Bottom-left
                    (center + size//4, center + size//4),  # Bottom-right
                    (center, center - size//3),           # Top
                    (center - size//3, center),           # Left
                    (center + size//3, center),           # Right
                    (center, center + size//3),           # Bottom
                ]
                
                for ax, ay in atom_positions:
                    if 0 <= ax < size and 0 <= ay < size:
                        for dx in range(atom_size):
                            for dy in range(atom_size):
                                if 0 <= ax + dx < size and 0 <= ay + dy < size:
                                    photo.put("#2c3e50", (ax + dx, ay + dy))
                
                photos.append(photo)
                
                # Try to set icon for each size (helps with different DE requirements)
                try:
                    self.iconphoto(True, photo)  # type: ignore
                    icon_set_success = True
                except Exception as e:
                    # iconphoto might not work on this system
                    if "not a photo image" in str(e).lower():
                        # Common Linux issue - continue anyway
                        pass
                    else:
                        print(f"Icon setting failed: {e}")
            
            # Keep references to prevent garbage collection
            self._icon_photos = photos
            
            # Try alternative icon setting methods if iconphoto failed
            try:
                # Create a simple bitmap fallback
                import tempfile
                temp_ico = tempfile.NamedTemporaryFile(suffix='.ico', delete=False)
                temp_ico.close()
                
                # Create a simple ICO file programmatically
                self.create_simple_ico_file(temp_ico.name)
                self.iconbitmap(temp_ico.name)
                os.unlink(temp_ico.name)
                print(f"✓ Created default molecular icon via ICO fallback (sizes: {', '.join(map(str, sizes))})")
            except Exception:
                print(f"✓ Created default molecular icon (sizes: {', '.join(map(str, sizes))}) - limited compatibility")
            
        except Exception as e:
            # If default icon creation fails, just continue without icon
            print(f"Could not create default icon: {e}")
            pass
    
    def create_simple_ico_file(self, filename):
        """Create a simple ICO file for systems where PhotoImage iconphoto doesn't work."""
        try:
            # Create a simple 32x32 ICO file manually
            # This is a minimal ICO file structure
            import struct
            
            # ICO header: 6 bytes
            ico_header = struct.pack('<HHH', 0, 1, 1)  # Reserved, Type, Count
            
            # ICO directory entry: 16 bytes  
            ico_dir = struct.pack('<BBBBHHII', 
                32, 32,     # Width, Height
                0,          # Color count (0 = no palette)
                0,          # Reserved
                1,          # Color planes
                32,         # Bits per pixel
                32*32*4,    # Size of image data
                22          # Offset to image data
            )
            
            # Create simple RGBA image data (32x32 pixels)
            image_data = bytearray()
            center = 16
            
            for y in range(32):
                for x in range(32):
                    # Create a simple molecular pattern
                    dist = ((x - center) ** 2 + (y - center) ** 2) ** 0.5
                    
                    if dist < 3:  # Center atom
                        # White center
                        image_data.extend([255, 255, 255, 255])  # BGRA
                    elif 8 <= dist <= 10:  # Ring of atoms
                        # Dark atoms
                        image_data.extend([64, 64, 64, 255])   # BGRA
                    elif dist < 14:  # Background
                        # Light blue background
                        image_data.extend([219, 119, 52, 255])  # BGRA (#3477db)
                    else:  # Outer edge
                        # Darker blue
                        image_data.extend([185, 128, 41, 255])  # BGRA (#2980b9)
            
            # Write ICO file
            with open(filename, 'wb') as f:
                f.write(ico_header)
                f.write(ico_dir)
                f.write(image_data)
                
        except Exception as e:
            # If ICO creation fails, create an even simpler fallback
            print(f"ICO creation failed: {e}")
            # Just create an empty file so iconbitmap doesn't crash
            with open(filename, 'wb') as f:
                f.write(b'')

    def expand_common_term(self, value):
        """Replace asterisks (*) in value with the common term"""
        common_term = self.common_term_var.get().strip()
        if common_term and '*' in value:
            return value.replace('*', common_term)
        return value

    def update_dipole_fields(self):
        """Update the displayed dipole parameter fields based on calculation type"""
        calc_type = self.dipole_calc_type.get()


if __name__=="__main__":
    PipelineGUI().mainloop()
