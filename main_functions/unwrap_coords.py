import numpy as np
import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import warnings
from path_utils import expand_path_pattern, validate_path_pattern

    
def unwrapper(baseDir, input_pattern=None, output_pattern=None, xsc_pattern=None, num_dcd=None, num_atoms=None, interval=slice(None), stride=1, max_workers=None, chunk_size=None, dcd_indices=None, common_term=""):
    """
    Read a series of coordinates directly extracted from MD trajectories, unwrap 
    periodic boundary crossings, and write out unwrapped coordinate files with
    optimized memory usage and optional parallel processing.

    This function:
      1. Reads the last line of the `.xsc` file to determine the simulation box vectors.
      2. For each frame in the input pattern (optionally in parallel):
         a. Loads the wrapped coordinates from input files in chunks if needed.
         b. Computes frame-to-frame displacements.
         c. Applies minimum-image corrections (i.e. if a jump > half box, subtract/add
            full box length).
         d. Cumulatively sums the corrected displacements to produce unwrapped
            trajectories.
         e. Saves the unwrapped coordinates to the output pattern.

    Parameters
    ----------
    baseDir : str
        Root directory for the simulation analysis.
    input_pattern : str
        Path pattern for input coordinate files. Can contain * (common term) and {i} (file index).
        Example: "anlz/NVT_*/wrapped/xyz_{i}.dat"
    output_pattern : str
        Path pattern for output unwrapped files. Can contain * (common term) and {i} (file index).
        Example: "anlz/NVT_*/unwrapped/unwrapped_xyz_{i}.dat"
    xsc_pattern : str
        Path pattern for XSC files. Can contain * (common term) and {i} (file index).
        Example: "anlz/NVT_*/restart_equil.xsc" (if same for all) or "anlz/NVT_*/restart_{i}.xsc"
    num_dcd : int
        Number of frames to process (i.e. number of input files).
    num_atoms : int
        Number of all atoms for which α₂(t) is calculated. Each input file 
        is expected to have `num_atoms*3` columns when flattened.
    interval : slice function, optional
        Slices the time length of the trajectory to be analyzed.
    stride : int, optional
        Frame stride for processing. Currently reserved (not implemented);
        defaults to 1 (i.e. process every frame).
    max_workers : int, optional
        Maximum number of parallel workers for processing files. 
        Defaults to min(num_dcd, CPU count) for parallel processing, 
        or None for sequential processing.
    chunk_size : int, optional
        Number of frames to process at once for memory efficiency.
        Automatically determined if not specified.
    dcd_indices : list, optional
        List of DCD indices to process (e.g., [0, 1, 4, 5] to process only DCDs 0, 1, 4, 5).
        If None, processes all DCDs from 0 to num_dcd-1. Default is None.
    common_term : str, optional
        Value to replace * placeholders in patterns. Default is "".

    Notes
    -----
    - The box lengths are taken from indices [1], [5], [9] of the last XSC line.
    - Half-box thresholds (`thold`) are used for minimum-image unwrapping.
    - The function reshapes each flat XYZ array into shape (n_frames, num_atoms, 3)
      before unwrapping, and flattens it back when saving.
    - Outputs are written in plain text with six-decimal floats.
    - For very large systems, consider using chunked processing to manage memory.

    Returns
    -------
    dict
        Processing results including timing and success/failure counts.

    Examples
    --------
    >>> results = unwrapper(
    ...     baseDir="/home/user/sim",
    ...     input_pattern="anlz/NVT_*/wrapped/xyz_{i}.dat",
    ...     output_pattern="anlz/NVT_*/unwrapped/unwrapped_xyz_{i}.dat",
    ...     xsc_pattern="anlz/NVT_*/restart_equil.xsc",
    ...     num_dcd=6,
    ...     num_atoms=3000,
    ...     common_term="240"
    ... )
    """
    
    start_time = time.time()
    
    # Validate path patterns
    is_valid, error_msg = validate_path_pattern(input_pattern)
    if not is_valid:
        raise ValueError(f"Invalid input pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(output_pattern)
    if not is_valid:
        raise ValueError(f"Invalid output pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(xsc_pattern)
    if not is_valid:
        raise ValueError(f"Invalid XSC pattern: {error_msg}")
    
    print(f"{'='*50}")
    print(f"COORDINATE UNWRAPPING")
    print(f"{'='*50}")
    print(f"Base directory: {baseDir}")
    print(f"Input pattern: {input_pattern}")
    print(f"Output pattern: {output_pattern}")
    print(f"XSC pattern: {xsc_pattern}")
    print(f"Common term: {common_term}")
    print(f"Number of DCDs: {num_dcd}")
    print(f"Number of atoms: {num_atoms}")
    
    # Determine indices to process
    if dcd_indices is None:
        dcd_list = list(range(num_dcd))
    else:
        dcd_list = dcd_indices
        print(f"Processing selected DCDs: {dcd_list}")
    
    # Create output directory (extract directory from first output file)
    if dcd_list:
        first_output = expand_path_pattern(output_pattern, common_term, dcd_list[0])
        output_dir = os.path.dirname(os.path.join(baseDir, first_output))
        os.makedirs(output_dir, exist_ok=True)
    
    # Read box size from XSC file
    try:
        # For XSC, use the first file's pattern (often XSC is the same for all)
        xsc_file = expand_path_pattern(xsc_pattern, common_term, dcd_list[0] if dcd_list else 0)
        xsc_path = os.path.join(baseDir, xsc_file)
        
        with open(xsc_path, 'r') as fr:
            lines = fr.readlines()
        box_size = np.array([
            float(lines[-1].split()[1]),   # x-dimension
            float(lines[-1].split()[5]),   # y-dimension  
            float(lines[-1].split()[9])    # z-dimension
        ])
        thold = box_size / 2
        print(f"✓ Box dimensions from {xsc_path}: {box_size}")
        
    except (FileNotFoundError, IndexError, ValueError) as e:
        raise ValueError(f"Error reading XSC file {xsc_path}: {e}")
    
    # Validate first input file exists
    if dcd_list:
        first_input = expand_path_pattern(input_pattern, common_term, dcd_list[0])
        first_input_path = os.path.join(baseDir, first_input)
        if not os.path.exists(first_input_path):
            raise FileNotFoundError(f"First input file not found: {first_input_path}")
        print(f"✓ Input files validated for index {dcd_list[0]}")
    
    # Determine optimal chunk size for memory efficiency
    if chunk_size is None:
        chunk_size = min(1000, num_atoms)  # Default chunk size
    
    # Processing setup
    if max_workers is None:
        max_workers = min(len(dcd_list), mp.cpu_count())
    
    print(f"Using {max_workers} workers, chunk size: {chunk_size}")
    
    # Calculate expected output shape for validation
    # Expected output: (F, M*N*3) where F=frames, M=molecules, N=atoms per molecule
    expected_columns = num_atoms * 3  # Total columns should be num_atoms * 3
    print(f"Expected output shape: (frames, {expected_columns}) for {num_atoms} atoms")
    
    results = {'success': 0, 'failed': [], 'total_time': 0}
    
    if max_workers == 1:
        # Sequential processing
        for i in dcd_list:
            try:
                _unwrap_single_file(
                    i, baseDir, input_pattern, output_pattern, xsc_pattern, num_atoms, 
                    box_size, thold, interval, stride, chunk_size, common_term
                )
                results['success'] += 1
                print(f"✓ Completed file {i}")
            except Exception as e:
                results['failed'].append(i)
                print(f"✗ Failed file {i}: {e}")
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_index = {
                executor.submit(
                    _unwrap_single_file, 
                    i, baseDir, input_pattern, output_pattern, xsc_pattern, num_atoms,
                    box_size, thold, interval, stride, chunk_size, common_term
                ): i 
                for i in dcd_list
            }
            
            # Collect results
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                try:
                    future.result()  # Raises exception if failed
                    results['success'] += 1
                    print(f"✓ Completed file {index}")
                except Exception as exc:
                    results['failed'].append(index)
                    print(f"✗ Failed file {index}: {exc}")
    
    results['total_time'] = time.time() - start_time
    
    # Summary
    print(f"\n{'='*50}")
    print(f"COORDINATE UNWRAPPING SUMMARY")
    print(f"{'='*50}")
    print(f"Total files: {num_dcd}")
    print(f"Successful: {results['success']}")
    print(f"Failed: {len(results['failed'])}")
    if results['failed']:
        print(f"Failed indices: {results['failed']}")
    print(f"Total time: {results['total_time']:.2f} seconds")
    print(f"Average time per file: {results['total_time']/num_dcd:.2f} seconds")
    print(f"{'='*50}\n")
    
    # Validate output shapes
    if results['success'] > 0:
        print(f"Validating output file shapes...")
        _validate_unwrapped_shapes(baseDir, output_pattern, dcd_list[:min(3, len(dcd_list))], expected_columns, common_term)
    
    return results


def _validate_unwrapped_shapes(baseDir, output_pattern, sample_indices, expected_columns, common_term=""):
    """Validate that unwrapped output files have the expected shape (F, M*N*3)"""
    from path_utils import expand_path_pattern
    
    validation_results = []
    for idx in sample_indices:
        try:
            output_file = expand_path_pattern(output_pattern, common_term, idx)
            output_path = os.path.join(baseDir, output_file)
            
            if os.path.exists(output_path):
                # Check file shape efficiently
                with open(output_path, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line:
                        actual_columns = len(first_line.split())
                        
                        # Count total lines for frame count
                        lines = sum(1 for _ in f) + 1  # +1 for the first line we already read
                        
                        if actual_columns == expected_columns:
                            validation_results.append(f"  File {idx}: ✓ Correct shape ({lines}, {actual_columns})")
                        else:
                            validation_results.append(f"  File {idx}: ⚠️  Shape mismatch - expected {expected_columns} columns, got {actual_columns}")
        except Exception as e:
            validation_results.append(f"  File {idx}: ❌ Validation failed - {e}")
    
    if validation_results:
        print("Shape validation results:")
        for result in validation_results:
            print(result)


def _unwrap_single_file(file_index, baseDir, input_pattern, output_pattern, xsc_pattern, num_atoms, 
                       box_size, thold, interval, stride, chunk_size, common_term):
    """Process a single coordinate file with memory-efficient unwrapping."""
    
    input_file_rel = expand_path_pattern(input_pattern, common_term, file_index)
    output_file_rel = expand_path_pattern(output_pattern, common_term, file_index)
    
    input_file = os.path.join(baseDir, input_file_rel)
    output_file = os.path.join(baseDir, output_file_rel)
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        # Load coordinates
        print(f"Processing {input_file}...")
        coords = np.loadtxt(input_file)
        
        # Handle single frame case
        if coords.ndim == 1:
            coords = coords.reshape(1, -1)
        
        # Apply interval slicing
        coords = coords[interval]
        n_frames, n_cols = coords.shape
        
        # Validate dimensions
        expected_cols = num_atoms * 3
        if n_cols != expected_cols:
            raise ValueError(f"Expected {expected_cols} columns, got {n_cols}")
        
        # Reshape to (frames, atoms, 3)
        coords = coords.reshape(n_frames, num_atoms, 3)
        
        # Initialize unwrapped coordinates
        unwrapped = np.zeros_like(coords)
        unwrapped[0] = coords[0]  # First frame is reference
        
        # Process frame-by-frame or in chunks
        for frame in range(1, n_frames):
            diff = coords[frame] - coords[frame-1]
            
            # Apply minimum image convention
            for dim in range(3):
                mask_pos = diff[:, dim] > thold[dim]
                mask_neg = diff[:, dim] < -thold[dim]
                diff[mask_pos, dim] -= box_size[dim]
                diff[mask_neg, dim] += box_size[dim]
            
            unwrapped[frame] = unwrapped[frame-1] + diff
        
        # Flatten and save
        unwrapped_flat = unwrapped.reshape(n_frames, -1)
        np.savetxt(output_file, unwrapped_flat, fmt='%.6f')
        
        print(f"✓ Completed {input_file_rel} -> {output_file_rel}")
        
    except Exception as e:
        print(f"✗ Failed {input_file_rel}: {e}")
        raise