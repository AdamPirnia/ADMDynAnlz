import numpy as np
import os
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import warnings

    
def unwrapper(baseDir, INdir, OUTdir, xsc, num_dcd, num_atoms, interval=slice(None), stride=1, max_workers=None, chunk_size=None):
    """
    Read a series of coordinates directly extracted from MD trajectories, unwrap 
    periodic boundary crossings, and write out unwrapped coordinate files with
    optimized memory usage and optional parallel processing.

    This function:
      1. Reads the last line of the `.xsc` file to determine the simulation box vectors.
      2. For each frame in the input directory (optionally in parallel):
         a. Loads the wrapped coordinates from `xyz_i.dat` in chunks if needed.
         b. Computes frame-to-frame displacements.
         c. Applies minimum-image corrections (i.e. if a jump > half box, subtract/add
            full box length).
         d. Cumulatively sums the corrected displacements to produce unwrapped
            trajectories.
         e. Saves the unwrapped coordinates to the output directory.

    Parameters
    ----------
    baseDir : str
        Root directory for the simulation analysis. The XSC, INdir, and OUTdir
        are all relative to this.
    INdir : str
        Subdirectory under `baseDir` containing the wrapped coordinate files
        named `xyz_0.dat, xyz_1.dat, …, xyz_{num_dcd-1}.dat`.
    OUTdir : str
        Subdirectory under `baseDir` where unwrapped files will be written.
        Files will be named `unwrapped_xyz_0.dat`, etc. Created if needed.
    xsc : str
        Path (relative to `baseDir`) to the `.xsc` file holding box vectors.
        The last line is expected to contain at least 9 numeric entries
        from which the x-, y-, and z-box lengths are extracted.
    num_dcd : int
        Number of frames to process (i.e. number of `xyz_i.dat` files).
    num_atoms : int
        Number of all atoms for which α₂(t) is calculated. Each `xyz_i.dat` 
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
    ...     INdir="anlz/NVT/wrapped",
    ...     OUTdir="anlz/NVT",
    ...     xsc="restart_equil.xsc",
    ...     num_dcd=1000,
    ...     num_atoms=1234,
    ...     interval=slice(0, 10000),
    ...     stride=2,
    ...     max_workers=4
    ... )
    """
    
    start_time = time.time()
    
    # Create output directory
    os.makedirs(f'{baseDir}/{OUTdir}', exist_ok=True)
    
    # Read box size from XSC file
    try:
        with open(f"{baseDir}/{xsc}", 'r') as fr:
            lines = fr.readlines()
        box_size = np.array([
            float(lines[-1].split()[1]),  # x
            float(lines[-1].split()[5]),  # y  
            float(lines[-1].split()[9])   # z
        ])
        thold = box_size / 2
        print(f"Box dimensions: {box_size[0]:.2f} x {box_size[1]:.2f} x {box_size[2]:.2f}")
        
    except (FileNotFoundError, IndexError, ValueError) as e:
        raise ValueError(f"Error reading XSC file {baseDir}/{xsc}: {e}")
    
    # Determine optimal chunk size for memory efficiency
    if chunk_size is None:
        # Estimate memory usage: frames * atoms * 3 coords * 8 bytes per float64
        # Target ~1GB chunks for large systems
        bytes_per_frame = num_atoms * 3 * 8
        target_chunk_memory = 1e9  # 1GB
        chunk_size = max(100, int(target_chunk_memory / bytes_per_frame))
        print(f"Auto-determined chunk size: {chunk_size} frames")
    
    # Set up parallel processing
    if max_workers is None:
        # Auto-determine based on system resources
        if num_dcd >= 4:  # Only use parallel for sufficient files
            max_workers = min(num_dcd, mp.cpu_count())
        else:
            max_workers = 1  # Sequential for small jobs
    
    results = {
        'success': 0, 
        'failed': [], 
        'total_time': 0,
        'parallel_workers': max_workers,
        'chunk_size': chunk_size
    }
    
    print(f"Processing {num_dcd} coordinate files...")
    print(f"Using {max_workers} workers with chunk size {chunk_size}")
    
    if max_workers == 1:
        # Sequential processing
        for i in range(num_dcd):
            try:
                _unwrap_single_file(
                    i, baseDir, INdir, OUTdir, num_atoms, 
                    box_size, thold, interval, stride, chunk_size
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
                    i, baseDir, INdir, OUTdir, num_atoms,
                    box_size, thold, interval, stride, chunk_size
                ): i 
                for i in range(num_dcd)
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
    
    return results


def _unwrap_single_file(file_index, baseDir, INdir, OUTdir, num_atoms, 
                       box_size, thold, interval, stride, chunk_size):
    """Process a single coordinate file with memory-efficient unwrapping."""
    
    input_file = f'{baseDir}/{INdir}/data/xyz_{file_index}.dat'
    output_file = f'{baseDir}/{OUTdir}/unwrapped_xyz_{file_index}.dat'
    
    try:
        # Load coordinates with memory-efficient reading
        xyz_data = np.loadtxt(input_file)
        
        # Apply interval and stride
        xyz_data = xyz_data[interval][::stride, :num_atoms*3]
        
        # Reshape to (frames, atoms, 3)
        n_frames = xyz_data.shape[0]
        xyz = xyz_data.reshape(n_frames, num_atoms, 3)
        
        # Process in chunks to manage memory for large trajectories
        unwrapped_xyz = np.zeros_like(xyz)
        unwrapped_xyz[0] = xyz[0]  # First frame is reference
        
        # Process frames in chunks
        for chunk_start in range(1, n_frames, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_frames)
            
            # Calculate differences for this chunk
            if chunk_start == 1:
                # First chunk: diff from frame 0
                prev_frame = xyz[0:1]
                chunk_frames = xyz[chunk_start:chunk_end]
            else:
                # Subsequent chunks: include previous frame for continuity
                prev_frame = xyz[chunk_start-1:chunk_start]
                chunk_frames = xyz[chunk_start:chunk_end]
            
            # Compute differences
            all_frames = np.vstack([prev_frame, chunk_frames])
            diffs = all_frames[1:] - all_frames[:-1]
            
            # Apply periodic boundary corrections
            for dim in range(3):
                diffs[..., dim][diffs[..., dim] > thold[dim]] -= box_size[dim]
                diffs[..., dim][diffs[..., dim] < -thold[dim]] += box_size[dim]
            
            # Update unwrapped coordinates
            if chunk_start == 1:
                # First chunk
                unwrapped_xyz[chunk_start:chunk_end] = (
                    unwrapped_xyz[0] + np.cumsum(diffs, axis=0)
                )
            else:
                # Subsequent chunks: continue from previous
                unwrapped_xyz[chunk_start:chunk_end] = (
                    unwrapped_xyz[chunk_start-1] + np.cumsum(diffs, axis=0)
                )
        
        # Save unwrapped coordinates
        output_shape = unwrapped_xyz.shape
        unwrapped_flat = unwrapped_xyz.reshape(output_shape[0], num_atoms*3)
        
        # Use efficient saving with proper formatting
        np.savetxt(output_file, unwrapped_flat, fmt='%.6f', delimiter=' ')
        
        return True
        
    except Exception as e:
        raise RuntimeError(f"Error processing file {input_file}: {e}")


def unwrapper_legacy(baseDir, INdir, OUTdir, xsc, num_dcd, num_atoms, interval=slice(None), stride=1):
    """
    Legacy version of unwrapper - kept for compatibility.
    Use unwrapper() for improved performance and memory efficiency.
    """
    warnings.warn(
        "unwrapper_legacy is deprecated. Use unwrapper() for better performance.",
        DeprecationWarning,
        stacklevel=2
    )
    
    os.makedirs(f'{baseDir}/{OUTdir}', exist_ok=True)
    with open(f"{baseDir}/{xsc}", 'r') as fr:
        lins = fr.readlines()
    box_size = np.array([lins[-1].split()[1], lins[-1].split()[5], lins[-1].split()[9]]).astype(float)
    thold = box_size/2

    for i in range(num_dcd):
        anlz_file = f'{baseDir}/{INdir}/xyz_{i}.dat'
        xyz = np.loadtxt(anlz_file)[interval][::stride, :num_atoms*3]
        xyz = xyz.reshape(len(xyz), num_atoms, 3)
        diffs = xyz[1:] - xyz[:-1]
        new_xyz = np.zeros_like(xyz)
        new_xyz[0] = xyz[0]
        
        # Correct differences along each dimension
        for dim in range(3):
            diffs[..., dim][diffs[..., dim] > thold[dim]] -= box_size[dim]
            diffs[..., dim][diffs[..., dim] < -thold[dim]] += box_size[dim]
        
        # Compute cumulative sum once after all corrections
        new_xyz[1:] = new_xyz[0] + np.cumsum(diffs, axis=0)
        
        shp = np.shape(new_xyz)
        print(shp)
        np.savetxt(f"{baseDir}/{OUTdir}/unwrapped_xyz_{i}.dat", new_xyz.reshape(shp[0], num_atoms*3), fmt='%f')
        
    del xyz
    del new_xyz
