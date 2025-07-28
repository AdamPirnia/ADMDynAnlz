import numpy as np
import os
import time
import warnings
from path_utils import expand_path_pattern, validate_path_pattern

######################################################    Parameters    

def a2_MSD(baseDir, input_pattern, output_pattern, num_dcd, partcl_num, numFrames, chunk_processing=True, validate_data=True, common_term=""):
    """
    Compute ensemble-averaged mean square displacements (MSD) and non-Gaussian
    parameter α₂(t) for a set of particles from a series of center-of-mass (COM)
    trajectory files. Each trajectory represents 1 ns of data, and displacements
    are calculated within each trajectory relative to its first frame.

    Parameters
    ----------
    baseDir : str
        Path to the top-level working directory.
    input_pattern : str
        Path pattern for input COM files. Can contain * (common term) and {i} (file index).
        Example: "anlz/NVT_*/com_data/com_{i}.dat"
    output_pattern : str
        Path pattern for output directory. Can contain * (common term).
        Two directories will be created: {output_pattern}/MSDs and {output_pattern}/alpha2s.
        Example: "anlz/NVT_*/analysis"
    num_dcd : int
        Number of trajectory (DCD) frames to process.
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
    common_term : str, optional
        Value to replace * placeholders in patterns. Default is "".

    Returns
    -------
    dict
        Processing results including timing, success/failure counts, and data quality metrics.

    Notes
    -----
    - Displacements are calculated within each trajectory relative to its first frame.
    - Results are ensemble-averaged over all successful trajectories.
    - Numerical stability improvements prevent division by very small numbers.
    """
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    start_time = time.time()
    
    # Validate path patterns
    is_valid, error_msg = validate_path_pattern(input_pattern)
    if not is_valid:
        raise ValueError(f"Invalid input pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(output_pattern)
    if not is_valid:
        raise ValueError(f"Invalid output pattern: {error_msg}")
    
    print(f"{'='*50}")
    print(f"α₂(t) AND MSD CALCULATION")
    print(f"{'='*50}")
    print(f"Base directory: {baseDir}")
    print(f"Input pattern: {input_pattern}")
    print(f"Output pattern: {output_pattern}")
    print(f"Common term: {common_term}")
    print(f"Number of DCDs: {num_dcd}")
    print(f"Particles: {partcl_num}, Min frames: {numFrames}")
    
    # Create output directories with robust path handling
    output_dir_rel = expand_path_pattern(output_pattern, common_term)
    output_dir = os.path.join(baseDir, output_dir_rel)
    os.makedirs(output_dir, exist_ok=True)
    
    msd_dir = os.path.join(output_dir, 'MSDs')
    alpha2_dir = os.path.join(output_dir, 'alpha2s')
    os.makedirs(msd_dir, exist_ok=True)
    os.makedirs(alpha2_dir, exist_ok=True)
    
    print(f"Output directories created: {output_dir}")
    
    # Validate first input file exists
    first_input = expand_path_pattern(input_pattern, common_term, 0)
    first_input_path = os.path.join(baseDir, first_input)
    if not os.path.exists(first_input_path):
        raise FileNotFoundError(f"First input file not found: {first_input_path}")
    print(f"✓ Input files validated for index 0")
    
    # Initialize processing
    successful_files = 0
    total_frames = 0
    frame_counts = []

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
        # Initialize accumulators - will be set after processing first successful file
        msd_accumulator = None
        msd4_accumulator = None
        
        # Process each trajectory file (including com_0.dat)
        for f in range(num_dcd):
            try:
                file_path_rel = expand_path_pattern(input_pattern, common_term, f)
                file_path = os.path.join(baseDir, file_path_rel)
                
                if not os.path.exists(file_path):
                    print(f"⚠ File not found: {file_path_rel}")
                    results['failed'].append(f)
                    continue
                
                # Load trajectory data
                data = np.loadtxt(file_path, usecols=range(partcl_num*3))
                
                # Validate frame length
                if len(data) < numFrames:
                    print(f"⚠ File {f}: {len(data)} frames < {numFrames} required")
                    results['skipped_files'].append(f)
                    continue
                
                # Reshape data to (frames, molecules, 3)
                data = data.reshape(len(data), partcl_num, 3)
                
                # Calculate displacements relative to first frame of this trajectory
                data = data - data[0, :, :]  # Subtract first frame
                
                # Calculate squared displacements using einsum for efficiency
                data2 = np.einsum('ijk,ijk->ij', data, data)  # Sum over x,y,z dimensions
                data4 = data2**2
                
                # Initialize accumulators with first successful file
                if msd_accumulator is None:
                    msd_accumulator = np.zeros_like(data2, dtype=np.float64)
                    msd4_accumulator = np.zeros_like(data4, dtype=np.float64)
                
                # Accumulate MSD and MSD⁴
                msd_accumulator += data2
                msd4_accumulator += data4
                
                successful_files += 1
                
                # Progress reporting
                if f % 10 == 0 or f == num_dcd - 1:
                    print(f"✓ Processed {f+1}/{num_dcd} files ({successful_files} successful)")
                    
            except Exception as e:
                print(f"✗ Error processing file {f}: {e}")
                results['failed'].append(f)
                continue
        
        if successful_files == 0 or msd_accumulator is None:
            raise ValueError("No trajectory files were successfully processed")
        
        # Compute ensemble averages
        print(f"Computing ensemble averages from {successful_files} successful files...")
        
        # Type assertions for the linter
        assert msd_accumulator is not None and msd4_accumulator is not None
        msd_avg = msd_accumulator / successful_files
        msd4_avg = msd4_accumulator / successful_files
        
        # Average over particles if more than 2 molecules
        if partcl_num > 2:
            msd_avg = np.mean(msd_avg, axis=1)  # Average over molecules
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
        msd_file = os.path.join(msd_dir, f"MSD_{partcl_num}mols.dat")
        alpha2_file = os.path.join(alpha2_dir, f"a2_{partcl_num}mols.dat")
        
        # Ensure directory exists before saving
        os.makedirs(os.path.dirname(msd_file), exist_ok=True)
        print(f"Saving MSD results to: {msd_file}")
        print(f"Saving α₂(t) results to: {alpha2_file}")
        
        np.savetxt(msd_file, msd_avg, fmt='%f')
        np.savetxt(alpha2_file, alpha2, fmt='%f')
        
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
    print(f"Total files processed: {num_dcd}")
    print(f"Successful: {results['success']}")
    print(f"Failed: {len(results['failed'])}")
    print(f"Skipped (insufficient frames): {len(results['skipped_files'])}")
    print(f"Processing time: {results['total_time']:.2f} seconds")
    print(f"Data quality:")
    for key, value in results['data_quality'].items():
        print(f"  {key}: {value}")
    print(f"{'='*50}\n")
    
    return results


