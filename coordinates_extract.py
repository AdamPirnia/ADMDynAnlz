import os
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

def raw_coords(baseDir, INdir, OUTdir, psf, dcd, num_dcd, particles, resnam, vmd, max_workers=None):
    """
    Extracts raw XYZ coordinates from a series of DCD trajectories using VMD in parallel,
    by generating per-segment Tcl scripts, running them in batch, and saving
    the results as text files.

    For each trajectory chunk (0→1 ns, 1→2 ns, …), this function:
      1. Creates output directories (`OUTdir/data`, `writenCodes`, `logs`).
      2. Writes a VMD Tcl script that:
         - Loads the PSF and DCD for that chunk.
         - Iterates over every frame.
         - Selects residues 0 through `particles` (e.g. water molecules).
         - Extracts their x,y,z coordinates efficiently.
         - Writes one line per frame to `OUTdir/data/xyz_{i}.dat`.
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
        Name of the output subdirectory under `baseDir` in which `data/` will be 
        created to hold the `.dat` files. OUTdir should exist.
    psf : str
        Base filename (without extension) of the topology file in each chunk
        (e.g. `"system"` if your files are `system.psf`).
    dcd : str
        Base filename (without extension) of the trajectory file in each chunk
        (e.g. `"traj"` if your files are `traj.dcd`).
    num_dcd : int
        Number of trajectory chunks to process.  Generates scripts for
        indices `0` through `num_dcd-1`.
    particles : str
        Particle selection string for VMD (e.g., "0 to 999" for residues 0-999).
    resnam : str
        Residue name for selection (e.g., "WAT" for water molecules).
    vmd : str
        The path to the VMD executable.
    max_workers : int, optional
        Maximum number of parallel workers. Defaults to min(num_dcd, CPU count).
    
    Side Effects
    ------------
    - Creates directories:
        - `{baseDir}/{OUTdir}/data`
        - `writenCodes`
        - `logs`
    - Writes Tcl scripts to `writenCodes/coords_{i}.tcl`.
    - Runs VMD on each script in parallel, logging to `logs/log_{i}.lgo`.
    - Produces coordinate files `xyz_{i}.dat` in `{baseDir}/{OUTdir}/data`.

    Returns
    -------
    dict
        Summary of processing results with success/failure counts and timing.

    Example
    -------
    >>> results = raw_coords(
    ...     baseDir="/home/user/sim",
    ...     INdir="trajectories", 
    ...     OUTdir="extracted",
    ...     psf="system",
    ...     dcd="traj",
    ...     num_dcd=6,
    ...     particles="0 to 999",
    ...     resnam="WAT",
    ...     vmd="/usr/local/bin/vmd",
    ...     max_workers=4
    ... )
    # Creates:
    #   /home/user/sim/extracted/data/xyz_0.dat, …, xyz_5.dat
    #   writenCodes/coords_0.tcl, …, coords_5.tcl
    #   logs/log_0.lgo, …, log_5.lgo
    """
    
    start_time = time.time()
    
    # Create output directories
    os.makedirs(f'{baseDir}/{OUTdir}/data', exist_ok=True)
    os.makedirs('writenCodes', exist_ok=True)
    os.makedirs('logs', exist_ok=True)
    
    # Set up parallel processing
    if max_workers is None:
        max_workers = min(num_dcd, mp.cpu_count())
    
    print(f"Processing {num_dcd} trajectory chunks using {max_workers} parallel workers...")
    
    # Generate all TCL scripts first
    for i in range(num_dcd):
        _write_tcl_script(i, baseDir, INdir, OUTdir, psf, dcd, particles, resnam)
    
    # Process VMD scripts in parallel
    results = {'success': 0, 'failed': [], 'total_time': 0}
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_index = {
            executor.submit(_run_vmd_script, i, vmd): i 
            for i in range(num_dcd)
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                success, stdout, stderr = future.result()
                if success:
                    results['success'] += 1
                    print(f"✓ Completed trajectory chunk {index}")
                else:
                    results['failed'].append(index)
                    print(f"✗ Failed trajectory chunk {index}")
                    if stderr:
                        print(f"  Error: {stderr}")
                        
            except Exception as exc:
                results['failed'].append(index)
                print(f"✗ Exception in trajectory chunk {index}: {exc}")
    
    results['total_time'] = time.time() - start_time
    
    print(f"\n{'='*50}")
    print(f"COORDINATE EXTRACTION SUMMARY")
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


def _write_tcl_script(i, baseDir, INdir, OUTdir, psf, dcd, particles, resnam):
    """Write optimized TCL script for a single trajectory chunk."""
    
    with open(f"writenCodes/coords_{i}.tcl", "w") as f:
        f.write(f"""# Optimized coordinate extraction script for chunk {i}
# Generated by coordinates_extract.py

puts "Starting coordinate extraction for chunk {i}"
puts "Timestamp: [clock format [clock seconds]]"

# Set paths
set baseDir "{baseDir}"
set resDir "${{baseDir}}/{OUTdir}"
set dir "{INdir}"

# Open output file
set outfile [open "${{resDir}}/data/xyz_{i}.dat" w]

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
set sel [atomselect $molid "resname {resnam} and residue {particles}"]
set natoms [$sel num]
puts "Selected $natoms atoms"

# Process each frame efficiently
for {{set frame 0}} {{$frame < $nframes}} {{incr frame}} {{
    # Update selection to current frame
    $sel frame $frame
    
    # Get coordinates as flat list
    set coords [$sel get {{x y z}}]
    
    # Write coordinates efficiently (one puts call per frame)
    set line ""
    foreach coord $coords {{
        append line [format " %.6f %.6f %.6f" [lindex $coord 0] [lindex $coord 1] [lindex $coord 2]]
    }}
    puts $outfile $line
    
    # Progress indicator for large trajectories
    if {{$frame % 100 == 0}} {{
        puts "Processed frame $frame/$nframes"
    }}
}}

# Cleanup
$sel delete
mol delete $molid
close $outfile

puts "Coordinate extraction completed for chunk {i}"
puts "Output saved to: ${{resDir}}/data/xyz_{i}.dat"

# Force garbage collection and exit
gc
exit
""")


def _run_vmd_script(index, vmd_path):
    """Run a single VMD script and return results."""
    
    script_path = f"writenCodes/coords_{index}.tcl"
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


# Legacy function for backward compatibility
def raw_coords_legacy(baseDir, INdir, OUTdir, psf, dcd, num_dcd, particles, resnam, vmd):
    """
    Legacy version of raw_coords - DEPRECATED due to critical bugs.
    Use raw_coords() instead for fixed implementation with parallel processing.
    """
    import warnings
    warnings.warn(
        "raw_coords_legacy contains critical bugs and is deprecated. "
        "Use raw_coords() instead.", 
        DeprecationWarning, 
        stacklevel=2
    )
    
    # Original implementation with bugs - kept for reference
    os.makedirs(f'{baseDir}/{OUTdir}', exist_ok=True)
    os.makedirs(f'writenCodes', exist_ok=True)
    os.makedirs(f'logs', exist_ok=True)

    for i in range(num_dcd):
        with open(f"writenCodes/coords_{i}.tcl", "w") as f:
            f.write(f"""
# Print the start message with date and time
puts "Starting Analysis - date/time"                                            
date                                                                            


# Set the base directory path to the trajectories
set baseDir {baseDir}
set resDir {baseDir}/{OUTdir}

set dir {INdir}

set outfile [open  ${{resDir}}/xyz_{i}.dat w]
set PSF "${{baseDir}}/${{dir}}/{i}to{i+1}ns/{psf}.psf"
set DCD "${{baseDir}}/${{dir}}/{i}to{i+1}ns/{dcd}.dcd"

# Load the velocity trajectory file (PSF and DCD files)
set traj [mol load psf $PSF dcd $DCD]
puts "Trajectory is loaded"

# Get the number of frames in the trajectory
set nf [molinfo $traj get numframes]  
puts $nf 

# Loop over each frame in the trajectory
for {{set frame 0}} {{$frame < $nf}} {{incr frame}} {{
    # Go to the current frame
    animate goto $frame

    # Select water molecules
    set sel [atomselect top "resname {resnam} and residue {particles}"]

    set coor [$sel get {{x y z}}]

    foreach {{x y z}} $coords {{
        puts -nonewline $outfile [format " %.6f %.6f %.6f" $x $y $z]
    }}
    # terminate the line
    puts $outfile ""

    # clean up the selection
    $sel delete

    set frameData ""
    # Loop over each set of 3 coordinates in the coor list
    for {{set j 0}} {{$j < [llength $coor]}} {{incr j}} {{
        
        set x [expr {{double([lindex [lindex $coor $j] 0])}}]
        set y [expr {{double([lindex [lindex $coor $j] 1])}}]
        set z [expr {{double([lindex [lindex $coor $j] 2])}}]


        # Print the coordinates with 3 decimal places into the output file
        append frameData [format " %.6f %.6f %.6f " $x $y $z]
    }}
    puts $outfile $frameData
}}

mol delete all
close $outfile

#garbage collection
gc
exit
""")
        
        print("after writing file.")
        # Run the bash command using subprocess
        command1 = f"{vmd} -dispdev text -nt 64 -e writenCodes/coords_{i}.tcl > logs/log_{i}.lgo"
        print("after executing the command.")
        process = subprocess.run(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        stdout = process.stdout.decode('utf-8')
        stderr = process.stderr.decode('utf-8')

        print(f"STDOUT:\n{stdout}")
        print(f"STDERR:\n{stderr}")
    
