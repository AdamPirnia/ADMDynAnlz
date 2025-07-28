import os
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

def vmd_dipole_collective(baseDir, INdir, OUTdir, psf, dcd, num_dcd, target, vmd, max_workers=None, dcd_indices=None, output_pattern=None):
    """
    Calculates collective dipole moments from a series of DCD trajectories using VMD in parallel.
    
    This function computes the collective dipole moment of selected atoms using VMD's 
    measure dipole command, which is useful for analyzing overall molecular orientations
    and dipole fluctuations in MD simulations.

    For each trajectory chunk (0→1 ns, 1→2 ns, …), this function:
      1. Creates output directories based on output pattern or OUTdir.
      2. Writes a VMD Tcl script that:
         - Loads the PSF and DCD for that chunk.
         - Iterates over every frame.
         - Selects atoms according to the target selection.
         - Calculates collective dipole moment using VMD's measure dipole.
         - Writes results using output pattern or to `OUTdir/dipole_{i}.dat`.
      3. Invokes VMD in text mode (`-dispdev text`) in parallel processes.
      4. Captures stdout/stderr to `logs/log_{i}.lgo`.

    Parameters
    ----------
    baseDir : str
        Root directory containing the `INdir` subfolders named
        "0to1ns", "1to2ns", … up to `num_dcd-1`to`num_dcd`ns.
    INdir : str
        Name of the subfolder under `baseDir` where all trajectory chunks live.
    OUTdir : str
        Name of the output subdirectory under `baseDir` in which results will be 
        created as the `.dat` files. OUTdir will be created if it does not exist.
        (Used when output_pattern is None for backward compatibility)
    psf : str
        Base filename (+path between OUTdir and the file without extension) of the 
        topology file in each chunk (e.g. `"system"` if your files are `system.psf`).
    dcd : str
        Base filename (+path between OUTdir and the file without extension) of the 
        trajectory file in each chunk (e.g. `"traj"` if your files are `traj.dcd`).
    output_pattern : str, optional
        Path pattern for output files. Can contain * (common term) and {i} (file index).
        Example: "path/to/dipole_{i}.dat". If None, uses OUTdir for backward compatibility.
    num_dcd : int
        Number of trajectory chunks to process.  Generates scripts for
        indices `0` through `num_dcd-1`.
    target : str
        VMD atom selection string (e.g., "water", "resname WAT", "protein", 
        "segname PROT and backbone", "index 0 to 999").
    vmd : str
        The path to the VMD executable.
    max_workers : int, optional
        Maximum number of parallel workers. Defaults to min(num_dcd, CPU count).
    dcd_indices : list, optional
        List of DCD indices to process (e.g., [0, 1, 4, 5] to process only DCDs 0, 1, 4, 5).
        If None, processes all DCDs from 0 to num_dcd-1. Default is None.
    
    Side Effects
    ------------
    - Creates directories:
        - `{baseDir}/{OUTdir}`
        - `writenCodes`
        - `logs`
    - Writes Tcl scripts to `writenCodes/dipole_calculator_{i}.tcl`.
    - Runs VMD on each script in parallel, logging to `logs/log_{i}.lgo`.
    - Produces dipole files `dipole_{i}.dat` in `{baseDir}/{OUTdir}`.

    Returns
    -------
    dict
        Summary of processing results with success/failure counts and timing.
        
    Example
    -------
    >>> results = vmd_dipole_collective(
    ...     baseDir="/home/user/sim",
    ...     INdir="trajectories", 
    ...     OUTdir="dipoles",
    ...     psf="system",
    ...     dcd="traj",
    ...     num_dcd=6,
    ...     target="resname WAT",
    ...     vmd="/usr/local/bin/vmd",
    ...     max_workers=4
    ... )
    # Creates:
    #   Output files according to output_pattern or /home/user/sim/dipoles/dipole_0.dat, …, dipole_5.dat
    #   writenCodes/dipole_calculator_0.tcl, …, dipole_calculator_5.tcl
    #   logs/log_0.lgo, …, log_5.lgo
    """
    
    start_time = time.time()
    
    # Create output directories - handle both pattern and legacy formats
    if output_pattern:
        # New pattern-based approach
        from path_utils import expand_path_pattern
        # Create directory from first output file pattern
        first_output = expand_path_pattern(output_pattern, "", 0)
        output_dir = os.path.dirname(os.path.join(baseDir, first_output))
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
    else:
        # Legacy approach for backward compatibility
        os.makedirs(f'{baseDir}/{OUTdir}', exist_ok=True)
    
    os.makedirs('writenCodes', exist_ok=True)
    os.makedirs('logs', exist_ok=True)
    
    # Determine which DCDs to process
    if dcd_indices is None:
        dcd_list = list(range(num_dcd))
    else:
        dcd_list = dcd_indices
    
    # Set up parallel processing
    if max_workers is None:
        max_workers = min(len(dcd_list), mp.cpu_count())
    
    print(f"Processing DCDs: {dcd_list}")
    print(f"Using {max_workers} parallel workers for {len(dcd_list)} trajectory chunks...")
    
    # Generate all TCL scripts first
    for i in dcd_list:
        _write_dipole_tcl_script(i, baseDir, INdir, OUTdir, psf, dcd, target, output_pattern)
    
    # Process VMD scripts in parallel
    results = {'success': 0, 'failed': [], 'total_time': 0}
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_index = {
            executor.submit(_run_vmd_dipole_script, i, vmd): i 
            for i in dcd_list
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                success, stdout, stderr = future.result()
                if success:
                    results['success'] += 1
                    print(f"✓ Completed dipole calculation for chunk {index}")
                else:
                    results['failed'].append(index)
                    print(f"✗ Failed dipole calculation for chunk {index}")
                    if stderr:
                        print(f"  Error: {stderr}")
                        
            except Exception as exc:
                results['failed'].append(index)
                print(f"✗ Exception in dipole calculation for chunk {index}: {exc}")
    
    results['total_time'] = time.time() - start_time
    
    print(f"\n{'='*50}")
    print(f"COLLECTIVE DIPOLE CALCULATION SUMMARY")
    print(f"{'='*50}")
    print(f"Total chunks: {num_dcd}")
    print(f"Successful: {results['success']}")
    print(f"Failed: {len(results['failed'])}")
    if results['failed']:
        print(f"Failed indices: {results['failed']}")
    print(f"Total time: {results['total_time']:.2f} seconds")
    print(f"Parallel efficiency: {max_workers} workers")
    print(f"{'='*50}\n")
    
    return results


def _write_dipole_tcl_script(i, baseDir, INdir, OUTdir, psf, dcd, target, output_pattern=None):
    """Write optimized TCL script for collective dipole moment calculation."""
    
    # Determine output file path
    if output_pattern:
        from path_utils import expand_path_pattern
        output_path = expand_path_pattern(output_pattern, "", i)
        output_full_path = os.path.join(baseDir, output_path)
    else:
        # Legacy format for backward compatibility
        output_full_path = f"{baseDir}/{OUTdir}/dipole_{i}.dat"
    
    with open(f"writenCodes/dipole_calculator_{i}.tcl", "w") as f:
        f.write(f"""# Optimized collective dipole moment calculation script for chunk {i}
# Generated by vmd_dipole.py

puts "Starting collective dipole moment calculation for chunk {i}"
puts "Timestamp: [clock format [clock seconds]]"

# Set paths
set baseDir "{baseDir}"
set dir "{INdir}"

# Open output file
set outfile [open "{output_full_path}" w]

# Load trajectory
set PSF "${{baseDir}}/${{dir}}/{i}to{i+1}ns/{psf}.psf"
set DCD "${{baseDir}}/${{dir}}/{i}to{i+1}ns/{dcd}.dcd"

puts "Loading PSF: $PSF"
puts "Loading DCD: $DCD"

set molid [mol load psf $PSF dcd $DCD]
puts "Molecule loaded with ID: $molid"

# Get trajectory info
set nframes [molinfo $molid get numframes]
puts "Number of frames: $nframes"

# Create selection (done once for efficiency)
set sel [atomselect $molid "{target}"]
set natoms [$sel num]
puts "Selected $natoms atoms with selection: {target}"

# Verify selection is valid
if {{$natoms == 0}} {{
    puts "ERROR: No atoms selected with selection '{target}'"
    exit 1
}}

# Write header to output file
puts $outfile "# Collective dipole moment calculation"
puts $outfile "# Selection: {target} ($natoms atoms)"
puts $outfile "# Format: frame dipole_x dipole_y dipole_z magnitude_debye"

# Process each frame efficiently
for {{set frame 0}} {{$frame < $nframes}} {{incr frame}} {{
    # Update selection to current frame
    $sel frame $frame
    
    # Calculate dipole moment vector and magnitude
    set mu_vector [measure dipole $sel]
    set mu_magnitude [veclength [measure dipole $sel -debye]]
    
    # Extract individual components
    set mu_x [lindex $mu_vector 0]
    set mu_y [lindex $mu_vector 1] 
    set mu_z [lindex $mu_vector 2]
    
    # Write frame data
    puts $outfile "$frame $mu_x $mu_y $mu_z $mu_magnitude"
    
    # Progress indicator for large trajectories
    if {{$frame % 100 == 0}} {{
        puts "Processed frame $frame/$nframes"
    }}
}}

# Cleanup
$sel delete
mol delete $molid
close $outfile

puts "Collective dipole moment calculation completed for chunk {i}"
puts "Output saved to: {output_full_path}"

# Force garbage collection and exit
gc
exit
""")


def _run_vmd_dipole_script(index, vmd_path):
    """Run a single VMD dipole script and return results."""
    
    script_path = f"writenCodes/dipole_calculator_{index}.tcl"
    log_path = f"logs/log_{index}.lgo"
    
    # Build VMD command
    command = f'"{vmd_path}" -dispdev text -nt 1 -e "{script_path}"'
    
    try:
        # Run VMD process
        process = subprocess.run(
            command, 
            shell=True, 
            capture_output=True, 
            text=True,
            timeout=3600  # 1 hour timeout per chunk
        )
        
        # Save log file
        with open(log_path, 'w') as log_file:
            log_file.write(f"Command: {command}\n")
            log_file.write(f"Return code: {process.returncode}\n")
            log_file.write(f"STDOUT:\n{process.stdout}\n")
            log_file.write(f"STDERR:\n{process.stderr}\n")
        
        success = process.returncode == 0
        return success, process.stdout, process.stderr
        
    except subprocess.TimeoutExpired:
        error_msg = f"VMD process for chunk {index} timed out after 1 hour"
        with open(log_path, 'w') as log_file:
            log_file.write(f"ERROR: {error_msg}\n")
        return False, "", error_msg
        
    except Exception as e:
        error_msg = f"Exception running VMD for chunk {index}: {str(e)}"
        with open(log_path, 'w') as log_file:
            log_file.write(f"ERROR: {error_msg}\n")
        return False, "", error_msg



    
