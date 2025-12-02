# ADMDynAnlz

A comprehensive molecular dynamics analysis pipeline for NAMD simulations. 

---

## Contents

- **Data Processing Pipeline**: 
  - Coordinate extraction from DCD files
  - Coordinate unwrapping (periodic boundary condition removal)
  - Center-of-mass calculation
  
- **Analysis Calculations**:
  - **α₂(t) and MSD**: Standard non-Gaussian parameter α₂(t) = 3⟨Δr⁴⟩/(5⟨Δr⟩²)² - 1 with mean square displacement
  - **α_xz(t)**: Directional correlation parameter α_ij(t) = ⟨Δi²·Δj²⟩/(⟨Δi²⟩·⟨Δj²⟩) - 1
  - **Individual Dipole Moments**: Calculate and save individual dipole moments (vector and magnitude) of molecules
  - **Collective Dipole Moments**: Extract collective dipole moments of system components (e.g., protein, specific residues)
  - **Velocity extraction**: Extract velocities of center of masses for individual molecules from VELDCD files produced by NAMD.

- **GUI Interface**: User-friendly graphical interface for pipeline configuration and execution
- **Executable Version**: Standalone executable for Linux systems

---

## What Can Be Done With It

This pipeline enables researchers to:

1. **Analyze molecular motion** through mean square displacement calculations
2. **Quantify non-Gaussian behavior** using the α₂(t) parameter to detect deviations from normal diffusion
3. **Study directional correlations** with α_xz(t) to understand anisotropic motion
4. **Calculate dipole moments** for both individual molecules and collective system properties
5. **Process large trajectory datasets** from molecular dynamics simulations
6. **Generate publication-ready data** for diffusion, dynamics, and electrostatic analyses

---

## How It Works

The pipeline processes NAMD DCD trajectory files through a systematic workflow:

### **Data Processing Steps:**
1. **Coordinate Extraction**: Uses VMD to extract atomic coordinates from DCD files based on user-defined selections
2. **Coordinate Unwrapping**: Removes periodic boundary condition artifacts to obtain continuous molecular trajectories
3. **Center-of-Mass Calculation**: Computes molecular centers of mass for multi-atom molecules

### **Analysis Calculations:**
- **MSD/α₂(t)**: Analyzes displacement statistics to calculate mean square displacements and detect non-Gaussian diffusion behavior
- **α_xz(t)**: Computes directional correlation functions to study anisotropic motion
- **Dipole Analysis**: Calculates molecular dipole moments using atomic charges and coordinates, with options for individual molecules or collective system properties

The pipeline handles the mathematical calculations, file I/O, and data organization automatically, allowing researchers to focus on scientific interpretation of results.

---

## How to Use

### **Option 1: Graphical Interface (Recommended)**

```bash
# Clone the repository
git clone https://github.com/AdamPirnia/ADMDynAnlz.git
cd ADMDynAnlz

# Install dependencies (if needed)
sudo apt update && sudo apt install python3-tk -y  # Ubuntu/Debian

# Run the GUI
python3 ADMDynAnlz.py
```

**Using the GUI:**
1. **Set Common Parameters**: Specify base directory, number of DCD files, and number of particles
2. **Configure Data Processing**: Set up coordinate extraction, unwrapping, and COM calculation steps
3. **Choose Analysis Type**: Select α₂(t)/MSD, α_xz(t), or dipole moment calculations
4. **Generate Scripts**: Create customized analysis scripts for your system
5. **Execute**: Run the generated scripts locally or on computing clusters

### **Option 2: Standalone Executable (Linux Only)**

```bash
# Download and run the executable
chmod +x ADMDynAnlz_Linux
./ADMDynAnlz_Linux
```

### **Requirements**
- **Python 3.6+** with NumPy and tkinter
- **VMD** (Visual Molecular Dynamics) for coordinate extraction
- **Input files**: NAMD DCD trajectories, PSF structure files, atomic charges (for dipole calculations)

---
![Alt text](example.png)
---

## Platform Availability

- **Current Version**: Linux executable available (`ADMDynAnlz_Linux`)
- **Source Code**: Compatible with Linux, macOS, and Windows
- **Future Releases**: Executables for macOS and Windows will be provided soon

---

## FAQ

**Q: What file formats does the pipeline support?**  
A: The pipeline works with NAMD DCD trajectory files and PSF structure files. Coordinates are output as plain text files.

**Q: How long does a typical analysis take?**  
A: Processing time depends on trajectory size and system complexity. Small systems (< 1000 particles) typically process in minutes, while large systems may take hours.

**Q: Can I analyze only part of my trajectory?**  
A: Yes, the GUI allows you to select specific DCD files or frame ranges for processing.

**Q: What units are used for the outputs?**  
A: MSD values are in Ų (square Angstroms), dipole moments are in Debye units, and α₂(t)/α_xz(t) are dimensionless parameters.

**Q: Do I need programming experience to use this?**  
A: No, the graphical interface guides you through the setup process. However, basic familiarity with molecular dynamics simulations is helpful.

**Q: Can I run this on a computing cluster?**  
A: Yes, the pipeline can generate SLURM scripts for high-performance computing environments.

**Q: How do I cite this software?**  
A: Please refer to the license file for citation requirements and contact the author for publication guidelines.

**Q: What if I encounter errors or need help?**  
A: Check that all file paths are correct, VMD is properly installed, and input files are valid. Contact the author for technical support.

---

## License

© 2025 Adam Pirnia — All rights reserved.

This software is provided under a proprietary license. Please read the `LICENSE.txt` file for complete terms and conditions, including usage rights, restrictions, and attribution requirements.
