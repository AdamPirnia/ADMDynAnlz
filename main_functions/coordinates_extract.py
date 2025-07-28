import os
import subprocess
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
from path_utils import expand_path_pattern, validate_path_pattern

def raw_coords(baseDir, psf_pattern, dcd_pattern, output_pattern=None, num_dcd=None, target_selection=None, vmd=None, max_workers=None, dcd_indices=None, common_term=""):
    """
    Extracts raw XYZ coordinates from a series of DCD trajectories using VMD in parallel,
    by generating per-segment Tcl scripts, running them in batch, and saving
    the results as text files.

    For each trajectory chunk, this function:
      1. Creates output directories based on output pattern.
      2. Writes a VMD Tcl script that:
         - Loads the PSF and DCD for that chunk using the expanded patterns.
         - Iterates over every frame.
         - Selects atoms using the target_selection criteria.
         - Extracts their x,y,z coordinates efficiently.
         - Writes results using the output pattern.
      3. Invokes VMD in text mode (`-dispdev text`) in parallel processes.
      4. Captures stdout/stderr to `logs/log_{i}.lgo`.

    Parameters
    ----------
    baseDir : str
        Root directory for the simulation.
    psf_pattern : str
        Path pattern for PSF files. Can contain * (common term) and {i} (file index).
        Example: "trajectories_*/run_{i}/system.psf"
    dcd_pattern : str
        Path pattern for DCD files. Can contain * (common term) and {i} (file index).
        Example: "trajectories_*/run_{i}/traj.dcd"
    output_pattern : str
        Path pattern for output files. Can contain * (common term) and {i} (file index).
        Example: "path/to/xyz_{i}.dat"
    num_dcd : int
        Number of trajectory chunks to process.  Generates scripts for
        indices `0` through `num_dcd-1`.
    target_selection : str
        VMD atom selection string (e.g., "resname WAT and residue 0 to 999", "water", "protein").
    vmd : str
        The path to the VMD executable.
    max_workers : int, optional
        Maximum number of parallel workers. Defaults to min(num_dcd, CPU count).
    dcd_indices : list, optional
        List of DCD indices to process (e.g., [0, 1, 4, 5] to process only DCDs 0, 1, 4, 5).
        If None, processes all DCDs from 0 to num_dcd-1. Default is None.
    common_term : str, optional
        Value to replace * placeholders in patterns. Default is "".
    
    Side Effects
    ------------
    - Creates directories:
        - `{baseDir}/{output_pattern}`
        - `writenCodes`
        - `logs`
    - Writes Tcl scripts to `writenCodes/coords_{i}.tcl`.
    - Runs VMD on each script in parallel, logging to `logs/log_{i}.lgo`.
    - Produces coordinate files `xyz_{i}.dat` in `{baseDir}/{output_pattern}`.

    Returns
    -------
    dict
        Summary of processing results with success/failure counts and timing.

    Example
    -------
    >>> results = raw_coords(
    ...     baseDir="/home/user/sim",
    ...     psf_pattern="trajectories_*/run_{i}/system.psf",
    ...     dcd_pattern="trajectories_*/run_{i}/traj.dcd",
    ...     output_pattern="path/to/xyz_{i}.dat",
    ...     num_dcd=6,
    ...     target_selection="resname WAT and residue 0 to 999",
    ...     vmd="/usr/local/bin/vmd",
    ...     common_term="240"
    ... )
    # Creates:
    #   /home/user/sim/path/to/xyz_0.dat, …, path/to/xyz_5.dat
    #   writenCodes/coords_0.tcl, …, coords_5.tcl
    #   logs/log_0.lgo, …, log_5.lgo
    """
    
    start_time = time.time()
    
    # Validate path patterns
    is_valid, error_msg = validate_path_pattern(psf_pattern)
    if not is_valid:
        raise ValueError(f"Invalid PSF pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(dcd_pattern)
    if not is_valid:
        raise ValueError(f"Invalid DCD pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(output_pattern)
    if not is_valid:
        raise ValueError(f"Invalid output pattern: {error_msg}")
    
    print(f"{'='*50}")
    print(f"COORDINATE EXTRACTION")
    print(f"{'='*50}")
    print(f"Base directory: {baseDir}")
    print(f"PSF pattern: {psf_pattern}")
    print(f"DCD pattern: {dcd_pattern}")
    print(f"Output pattern: {output_pattern}")
    print(f"Common term: {common_term}")
    print(f"Number of DCDs: {num_dcd}")
    print(f"Target selection: {target_selection}")
    print(f"VMD executable: {vmd}")
    
    # Determine indices to process
    if dcd_indices is None:
        dcd_list = list(range(num_dcd))
    else:
        dcd_list = dcd_indices
        print(f"Processing selected DCDs: {dcd_list}")
    
    if max_workers is None:
        max_workers = min(len(dcd_list), mp.cpu_count())
    
    print(f"Using {max_workers} parallel workers")
    
    # Create necessary directories
    # Create output directory from first output file pattern
    if dcd_list:
        first_output = expand_path_pattern(output_pattern, common_term, dcd_list[0])
        output_full_path = os.path.join(baseDir, first_output)
        output_dir = os.path.dirname(output_full_path)
        
        print(f"Output pattern: {output_pattern}")
        print(f"Expanded pattern: {first_output}")
        print(f"Full output path: {output_full_path}")
        print(f"Output directory: {output_dir}")
        
        os.makedirs(output_dir, exist_ok=True)
    os.makedirs("writenCodes", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    
    # Note: With flexible target_selection, we can't easily estimate output size
    # The actual number of selected atoms will be determined by VMD at runtime
    
    # Validate that input files exist for first DCD (sanity check)
    if dcd_list:
        test_psf = os.path.join(baseDir, expand_path_pattern(psf_pattern, common_term, dcd_list[0]))
        test_dcd = os.path.join(baseDir, expand_path_pattern(dcd_pattern, common_term, dcd_list[0]))
        
        if not os.path.exists(test_psf):
            raise FileNotFoundError(f"PSF file not found: {test_psf}")
        if not os.path.exists(test_dcd):
            raise FileNotFoundError(f"DCD file not found: {test_dcd}")
        
        print(f"✓ Input files validated for index {dcd_list[0]}")
    
    # Results tracking
    results = {
        'successful': [],
        'failed': [],
        'total_time': 0,
        'num_processed': len(dcd_list)
    }
    
    # Generate all TCL scripts first
    for i in dcd_list:
        _write_tcl_script(i, baseDir, psf_pattern, dcd_pattern, output_pattern, target_selection, common_term)
    
    # Process VMD scripts in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_index = {
            executor.submit(_run_vmd_script, i, vmd, common_term): i 
            for i in dcd_list
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_index):
            index = future_to_index[future]
            try:
                success, stdout, stderr = future.result()
                if success:
                    results['successful'].append(index)
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
    print(f"Successful: {len(results['successful'])}")
    print(f"Failed: {len(results['failed'])}")
    if results['failed']:
        print(f"Failed indices: {results['failed']}")
    print(f"Total time: {results['total_time']:.2f} seconds")
    print(f"Parallel efficiency: {max_workers} workers")
    print(f"{'='*50}\n")
    
    # Note: Output shape validation removed since expected count depends on VMD selection
    
    return results



def _write_tcl_script(i, baseDir, psf_pattern, dcd_pattern, output_pattern, target_selection, common_term=""):
    """Write optimized TCL script for a single trajectory chunk."""
    
    # Expand patterns to get actual file paths
    psf_path = expand_path_pattern(psf_pattern, common_term, i)
    dcd_path = expand_path_pattern(dcd_pattern, common_term, i)
    output_path = expand_path_pattern(output_pattern, common_term, i)
    
    # Create absolute paths
    psf_full_path = os.path.join(baseDir, psf_path)
    dcd_full_path = os.path.join(baseDir, dcd_path)
    output_full_path = os.path.join(baseDir, output_path)
    
    print(f"DEBUG TCL Script {i}:")
    print(f"  Output pattern: {output_pattern}")
    print(f"  Expanded output: {output_path}")
    print(f"  Full output path: {output_full_path}")
    print(f"  Output directory: {os.path.dirname(output_full_path)}")
    
    # Ensure output directory exists before creating TCL script
    output_dir = os.path.dirname(output_full_path)
    if output_dir and output_dir != baseDir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Create more specific TCL filename with common term
    tcl_filename = f"coords_{common_term}_{i}.tcl" if common_term else f"coords_{i}.tcl"
    
    with open(f"writenCodes/{tcl_filename}", "w") as f:
        f.write(f"""# Optimized coordinate extraction script for chunk {i}
# Generated by coordinates_extract.py

puts "Starting coordinate extraction for chunk {i}"
puts "Timestamp: [clock format [clock seconds]]"

# Set paths
set baseDir "{baseDir}"

# Open output file
set outfile [open "{output_full_path}" w]

# Load trajectory
set PSF "{psf_full_path}"
set DCD "{dcd_full_path}"

puts "Loading PSF: $PSF"
puts "Loading DCD: $DCD"

# Debug output - check file existence
puts "About to load: PSF=$PSF, DCD=$DCD"
if {{![file exists $PSF]}} {{
    puts "ERROR: PSF file does not exist: $PSF"
    exit 1
}}
if {{![file exists $DCD]}} {{
    puts "ERROR: DCD file does not exist: $DCD"
    exit 1
}}

set molid [mol load psf "$PSF" dcd "$DCD"]
puts "Molecule loaded with ID: $molid"

# Get number of frames
set num_frames [molinfo $molid get numframes]
puts "Total frames: $num_frames"

# Create selection
set sel [atomselect $molid "{target_selection}"]
set num_atoms [$sel num]
puts "Selected $num_atoms atoms"

if {{$num_atoms == 0}} {{
    puts "ERROR: No atoms selected with criteria '{target_selection}'"
    exit 1
}}

# Extract coordinates for all frames
puts "Extracting coordinates..."
set frame_count 0
for {{set frame 0}} {{$frame < $num_frames}} {{incr frame}} {{
    $sel frame $frame
    set coords [$sel get {{x y z}}]
    
    # Write coordinates in a single line
    set line ""
    foreach coord $coords {{
        append line "[lindex $coord 0] [lindex $coord 1] [lindex $coord 2] "
    }}
    puts $outfile [string trimright $line]
    
    incr frame_count
    if {{$frame_count % 100 == 0}} {{
        puts "Processed $frame_count frames..."
    }}
}}

$sel delete
close $outfile
mol delete $molid

puts "Extraction complete: $frame_count frames processed"
puts "Output written to: {output_full_path}"
puts "Timestamp: [clock format [clock seconds]]"

exit 0
""")


def _run_vmd_script(index, vmd_path, common_term=""):
    """Run a single VMD script and return results."""
    
    # Use the same naming convention as in _write_tcl_script
    tcl_filename = f"coords_{common_term}_{index}.tcl" if common_term else f"coords_{index}.tcl"
    script_path = f"writenCodes/{tcl_filename}"
    log_path = f"logs/log_{index}.lgo"
    
    # Build VMD command - ensure proper quoting for paths with spaces
    command = f'"{vmd_path}" -dispdev text -nt 1 -e "{script_path}"'
    
    print(f"Command: {command}")
    print(f"Script path: {script_path}")
    print(f"Log path: {log_path}")
    
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



    
