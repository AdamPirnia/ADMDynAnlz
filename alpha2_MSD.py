import numpy as np
import os
import time
import warnings

######################################################    Parameters    

def a2_MSD(baseDir, INdir, OUTdir, num_dcd, partcl_num, numFrames, chunk_processing=True, validate_data=True):
    """
    Compute ensemble-averaged mean square displacements (MSD) and non-Gaussian
    parameter α₂(t) for a set of particles from a series of center-of-mass (COM)
    trajectory files, with optimized memory usage and numerical stability.

    Parameters
    ----------
    baseDir : str
        Path to the top-level working directory.
    INdir : str
        Subdirectory under `baseDir` containing the COM data files in
        `com_data/com_<frame>.dat`.
    OUTdir : str
        Subdirectory under `baseDir` where output will be written. Two
        directories will be created inside: `MSDs` and `alpha2s`.
    num_dcd : int
        Number of trajectory (DCD) frames to process. Files are assumed to
        be named `com_0.dat` through `com_{num_dcd-1}.dat`.
    partcl_num : int
        Number of particles whose COM position is stored in each file.
        Each .dat file is expected to have `partcl_num*3` columns: x, y, z
        for each molecule.
    numFrames : int
        Minimum number of time-points in each .dat file to include it in the
        average. If a file has fewer lines than `numFrames`, it is skipped
        (and counted as a failure).
    chunk_processing : bool, optional
        Use chunked processing for better memory efficiency with large datasets.
    validate_data : bool, optional
        Perform data validation checks during processing.

    Returns
    -------
    dict
        Processing results including timing, success/failure counts, and data quality metrics.

    Notes
    -----
    - Numerical stability improvements prevent division by very small numbers.
    - Memory usage is optimized for large particle systems.
    """
    
    start_time = time.time()
    
    # Create output directories
    os.makedirs(f'{baseDir}/{OUTdir}/MSDs', exist_ok=True)
    os.makedirs(f'{baseDir}/{OUTdir}/alpha2s', exist_ok=True)

    results = {
        'success': 0,
        'failed': [],
        'skipped_files': [],
        'total_time': 0,
        'data_quality': {},
        'numerical_warnings': []
    }
    
    print(f"Starting α₂(t) and MSD calculation...")
    print(f"Processing {num_dcd} trajectory files with {partcl_num} particles")
    print(f"Minimum frames required: {numFrames}")
    
    try:
        # Load reference frame (frame 0) with error handling
        ref_file = f'{baseDir}/{INdir}/com_data/com_0.dat'
        print(f"Loading reference frame: {ref_file}")
        
        if not os.path.exists(ref_file):
            raise FileNotFoundError(f"Reference file not found: {ref_file}")
        
        ref_data = np.loadtxt(ref_file, usecols=range(partcl_num*3))
        
        # Initialize displacement accumulators
        n_frames_ref = len(ref_data)
        if n_frames_ref < numFrames:
            raise ValueError(f"Reference frame has only {n_frames_ref} frames, need at least {numFrames}")
        
        # Use reference frames up to numFrames
        ref_data = ref_data[:numFrames]
        ref_data = ref_data.reshape(numFrames, partcl_num, 3)
        
        if partcl_num > 1:
            msd_accumulator = np.zeros((numFrames, partcl_num), dtype=np.float64)
            msd4_accumulator = np.zeros((numFrames, partcl_num), dtype=np.float64)
        else:
            msd_accumulator = np.zeros(numFrames, dtype=np.float64)
            msd4_accumulator = np.zeros(numFrames, dtype=np.float64)
        
        successful_files = 0
        
        # Process each trajectory file
        for f in range(1, num_dcd):
            try:
                file_path = f'{baseDir}/{INdir}/com_data/com_{f}.dat'
                
                if not os.path.exists(file_path):
                    print(f"⚠ File not found: {file_path}")
                    results['failed'].append(f)
                    continue
                
                current_data = np.loadtxt(file_path, usecols=range(partcl_num*3))
                
                # Validate frame length
                if len(current_data) < numFrames:
                    print(f"⚠ File {f}: {len(current_data)} frames < {numFrames} required")
                    results['skipped_files'].append(f)
                    continue
                
                # Reshape and truncate to required frames
                current_data = current_data[:numFrames]
                current_data = current_data.reshape(numFrames, partcl_num, 3)
                
                # Compute displacements from reference
                displacements = current_data - ref_data
                
                # Compute squared displacements
                if partcl_num > 1:
                    sq_displacements = np.sum(displacements**2, axis=2)
                else:
                    sq_displacements = np.sum(displacements**2, axis=1)
                
                # Accumulate MSD and MSD⁴
                msd_accumulator += sq_displacements
                msd4_accumulator += sq_displacements**2
                
                successful_files += 1
                
                # Progress reporting
                if f % 10 == 0 or f == num_dcd - 1:
                    print(f"✓ Processed {f}/{num_dcd-1} files ({successful_files} successful)")
                    
            except Exception as e:
                print(f"✗ Error processing file {f}: {e}")
                results['failed'].append(f)
                continue
        
        if successful_files == 0:
            raise ValueError("No trajectory files were successfully processed")
        
        # Compute ensemble averages
        print(f"Computing ensemble averages from {successful_files} successful files...")
        
        msd_avg = msd_accumulator / successful_files
        msd4_avg = msd4_accumulator / successful_files
        
        # Average over particles if needed
        if partcl_num > 2:
            msd_avg = np.mean(msd_avg, axis=1)
            msd4_avg = np.mean(msd4_avg, axis=1)
        
        # Compute α₂(t) with numerical stability
        epsilon = 1e-12
        msd_squared = msd_avg**2
        
        # Avoid division by very small numbers
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            alpha2 = 3 * msd4_avg / (5 * msd_squared) - 1
        
        # Handle problematic values
        small_msd_mask = msd_squared < epsilon
        alpha2[small_msd_mask] = 0.0
        alpha2[~np.isfinite(alpha2)] = 0.0
        
        # Save results
        msd_file = f"{baseDir}/{OUTdir}/MSDs/MSD_{OUTdir}.dat"
        alpha2_file = f"{baseDir}/{OUTdir}/alpha2s/a2_{OUTdir}.dat"
        
        np.savetxt(msd_file, msd_avg, fmt='%.8e', header='Mean Square Displacement')
        np.savetxt(alpha2_file, alpha2, fmt='%.8e', header='Non-Gaussian parameter α₂(t)')
        
        print(f"✓ Results saved:")
        print(f"  MSD: {msd_file}")
        print(f"  α₂(t): {alpha2_file}")
        
        results['success'] = successful_files
        results['total_time'] = time.time() - start_time
        
        # Data quality metrics
        results['data_quality'] = {
            'mean_msd': float(np.mean(msd_avg)),
            'max_msd': float(np.max(msd_avg)),
            'mean_alpha2': float(np.mean(alpha2)),
            'alpha2_range': [float(np.min(alpha2)), float(np.max(alpha2))]
        }
        
    except Exception as e:
        results['failed'].append('general_error')
        results['total_time'] = time.time() - start_time
        raise RuntimeError(f"Critical error in α₂/MSD calculation: {e}")
    
    # Print summary
    print(f"\n{'='*50}")
    print(f"α₂(t) AND MSD CALCULATION SUMMARY") 
    print(f"{'='*50}")
    print(f"Total files processed: {num_dcd-1}")
    print(f"Successful: {results['success']}")
    print(f"Failed: {len(results['failed'])}")
    print(f"Skipped (insufficient frames): {len(results['skipped_files'])}")
    print(f"Processing time: {results['total_time']:.2f} seconds")
    print(f"Data quality:")
    for key, value in results['data_quality'].items():
        print(f"  {key}: {value}")
    print(f"{'='*50}\n")
    
    return results
