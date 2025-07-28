#!/usr/bin/env python3
"""
Optimized Dipole Moment Calculation Module

Calculates molecular dipole moments from coordinate and COM data with enhanced
performance, parallel processing, and robust error handling.
"""
import numpy as np
import os
import time
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import warnings

# Debye conversion constant  
DEBYE_CONVERSION = 0.2081943  # e*Å to Debye

def dipoleM(coords, charges, com, atoms_per_mol, num_molecules_to_calc):
    """
    Calculate dipole moments for a set of molecular coordinates.
    
    Parameters:
    -----------
    coords : numpy.ndarray
        Coordinate array of shape (n_frames, n_particles * atoms_per_mol, 3)
    charges : numpy.ndarray
        Charge array of shape (atoms_per_mol,)
    com : numpy.ndarray
        Center of mass array of shape (n_frames, n_particles, 3)
    atoms_per_mol : int
        Number of atoms per molecule
    num_molecules_to_calc : int
        Number of molecules to calculate (first N molecules)
        
    Returns:
    --------
    tuple : (dipole_vectors, dipole_magnitudes)
    """
    n_frames = coords.shape[0]
    
    # Limit coordinates to only the first num_molecules_to_calc molecules
    coords_limited = coords[:, :num_molecules_to_calc * atoms_per_mol * 3]
    com_limited = com[:, :num_molecules_to_calc, :]
    
    # Reshape coordinates to (n_frames, num_molecules_to_calc, atoms_per_mol, 3)
    coords_reshaped = coords_limited.reshape(n_frames, num_molecules_to_calc, atoms_per_mol, 3)
    
    # Reshape COM to (n_frames, num_molecules_to_calc, 1, 3) for broadcasting
    com_reshaped = com_limited.reshape(n_frames, num_molecules_to_calc, 1, 3)
    
    # Calculate relative coordinates (center of charge)
    relative_coords = coords_reshaped - com_reshaped
    
    # Calculate dipole moments: sum over atoms in each molecule
    # dipole = sum(charge_i * r_i) for each molecule
    dipole_vectors = np.sum(
        relative_coords * charges[np.newaxis, np.newaxis, :, np.newaxis], 
        axis=2
    ) / DEBYE_CONVERSION
    
    # Calculate dipole magnitudes
    dipole_magnitudes = np.linalg.norm(dipole_vectors, axis=2)
    
    return dipole_vectors, dipole_magnitudes

def _process_single_dipole_file(file_idx, baseDir, coords_pattern, com_pattern, output_pattern, 
                                magnitudes_pattern, charges_array, atoms_per_particle, 
                                effective_molecules, stride, common_term):
    """Process a single trajectory file for dipole calculation - at module level for multiprocessing."""
    try:
        from path_utils import expand_path_pattern
        
        # File paths
        coord_file_rel = expand_path_pattern(coords_pattern, common_term, file_idx)
        com_file_rel = expand_path_pattern(com_pattern, common_term, file_idx)
        
        coord_file = os.path.join(baseDir, coord_file_rel)
        com_file = os.path.join(baseDir, com_file_rel)
        
        # Debug: Print file paths being checked
        print(f"Processing file {file_idx}:")
        print(f"  Coord file: {coord_file}")
        print(f"  COM file: {com_file}")
        
        # Check if files exist
        if not os.path.exists(coord_file):
            error_msg = f'Coordinate file not found: {coord_file}'
            print(f"  ERROR: {error_msg}")
            return {'success': False, 'error': error_msg, 'file_idx': file_idx}
        if not os.path.exists(com_file):
            error_msg = f'COM file not found: {com_file}'
            print(f"  ERROR: {error_msg}")
            return {'success': False, 'error': error_msg, 'file_idx': file_idx}
            
        # Load data
        print(f"  Loading data from files...")
        coord_data = np.loadtxt(coord_file)
        com_data = np.loadtxt(com_file)
        print(f"  Loaded coords: {coord_data.shape}, COM: {com_data.shape}")
        
        # Data validation and reshaping (similar to original logic)
        if coord_data.ndim == 1:
            coord_data = coord_data.reshape(1, -1)
        if com_data.ndim == 1:
            com_data = com_data.reshape(1, -1)
        
        # Apply stride if specified
        if stride > 1:
            coord_data = coord_data[::stride]
            com_data = com_data[::stride]
        
        # Reshape COM data to (n_frames, n_particles, 3) for dipoleM function
        n_frames = com_data.shape[0]
        num_particles = com_data.shape[1] // 3  # Total particles in system
        com_reshaped = com_data.reshape(n_frames, num_particles, 3)
        
        print(f"  Processing {n_frames} frames, {num_particles} particles, calculating {effective_molecules} molecules")
        
        # Use the actual dipoleM function to calculate dipole moments
        print(f"  Calling dipoleM with shapes: coords {coord_data.shape}, charges {charges_array.shape}, com {com_reshaped.shape}")
        dipole_vectors, dipole_magnitudes = dipoleM(coord_data, charges_array, com_reshaped, atoms_per_particle, effective_molecules)
        print(f"  dipoleM returned: vectors {dipole_vectors.shape}, magnitudes {dipole_magnitudes.shape}")
        
        # Save results
        if magnitudes_pattern:
            dipole_file = expand_path_pattern(output_pattern, common_term, file_idx)
            magnitude_file = expand_path_pattern(magnitudes_pattern, common_term, file_idx)
            
            dipole_path = os.path.join(baseDir, dipole_file)
            magnitude_path = os.path.join(baseDir, magnitude_file)
            
            os.makedirs(os.path.dirname(dipole_path), exist_ok=True)
            if magnitude_path != dipole_path:
                os.makedirs(os.path.dirname(magnitude_path), exist_ok=True)
        else:
            output_path = os.path.join(baseDir, expand_path_pattern(output_pattern, common_term, file_idx))
            dipole_path = os.path.join(output_path, f'dipoles_{file_idx}.dat')
            magnitude_path = os.path.join(output_path, f'dipole_magnitudes_{file_idx}.dat')
            os.makedirs(output_path, exist_ok=True)
        
        # Reshape dipole vectors to match expected output format
        dipole_vectors_flat = dipole_vectors.reshape(n_frames, -1)
        
        print(f"  Saving to: {dipole_path}")
        print(f"  Saving to: {magnitude_path}")
        np.savetxt(dipole_path, dipole_vectors_flat, fmt='%.8f')
        np.savetxt(magnitude_path, dipole_magnitudes, fmt='%.8f')
        
        print(f"  SUCCESS: File {file_idx} processed, {n_frames} frames")
        return {
            'success': True, 
            'file_idx': file_idx, 
            'frames_processed': n_frames,
            'mean_magnitude': np.mean(dipole_magnitudes),
            'std_magnitude': np.std(dipole_magnitudes)
        }
        
    except Exception as e:
        import traceback
        error_details = f"File {file_idx} failed: {str(e)}\nTraceback: {traceback.format_exc()}"
        print(f"  ERROR: {error_details}")
        return {'success': False, 'error': error_details, 'file_idx': file_idx}

def dipole_functions(baseDir, coords_pattern, com_pattern, output_pattern, Charges, num_dcds, num_particles, 
                    atoms_per_particle=3, stride=1, max_workers=1, chunk_processing=True, 
                    validate_data=True, progress_callback=None, molecules_to_process=None, common_term="",
                    magnitudes_pattern=None, dcd_indices=None):
    """
    Calculate molecular dipole moments from trajectory data.
    
    Parameters:
    -----------
    baseDir : str
        Base directory path
    coords_pattern : str  
        Path pattern for coordinate files. Can contain * (common term) and {i} (file index).
        Example: "anlz/NVT_*/unwrapped/continued_xyz_{i}.dat"
    com_pattern : str
        Path pattern for COM files. Can contain * (common term) and {i} (file index).
        Example: "anlz/NVT_*/com_data/com_{i}.dat"
    output_pattern : str
        Path pattern for output vector files OR output directory (backward compatibility).
        If magnitudes_pattern is None, treated as directory. Otherwise, treated as vector file pattern.
        Example: "anlz/NVT_*/dipole/vectors_{i}.dat"
    magnitudes_pattern : str, optional
        Path pattern for magnitude output files. If provided, output_pattern is treated as vector pattern.
        Example: "anlz/NVT_*/dipole/magnitudes_{i}.dat"
    dcd_indices : list, optional
        List of specific DCD indices to process. If None, processes all DCDs.
    Charges : list
        List of atomic charges for each atom in a molecule
    num_dcds : int
        Number of trajectory files to process
    num_particles : int
        Total number of molecules/particles in each file
    atoms_per_particle : int, optional
        Number of atoms per molecule (default: 3)
    stride : int, optional
        Frame stride for processing (default: 1)
    max_workers : int, optional
        Number of parallel workers (default: 1)
    chunk_processing : bool, optional
        Enable chunked processing for memory efficiency (default: True)
    validate_data : bool, optional
        Perform data validation checks (default: True)
    progress_callback : callable, optional
        Callback function for progress updates
    molecules_to_process : int, optional
        Number of molecules to process (default: all)
    common_term : str, optional
        Common term for path expansion (default: "")
    """
    
    start_time = time.time()
    
    # Validate path patterns
    from path_utils import expand_path_pattern, validate_path_pattern
    
    is_valid, error_msg = validate_path_pattern(coords_pattern)
    if not is_valid:
        raise ValueError(f"Invalid coords pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(com_pattern)
    if not is_valid:
        raise ValueError(f"Invalid COM pattern: {error_msg}")
        
    is_valid, error_msg = validate_path_pattern(output_pattern)
    if not is_valid:
        raise ValueError(f"Invalid output pattern: {error_msg}")
    
    # Determine effective number of molecules to process
    if molecules_to_process is None:
        effective_molecules = num_particles
    else:
        effective_molecules = min(molecules_to_process, num_particles)
        if molecules_to_process > num_particles:
            print(f"Warning: molecules_to_process ({molecules_to_process}) > num_particles ({num_particles}). Processing all {num_particles} molecules.")
    
    # Input validation
    if not isinstance(Charges, (list, tuple, np.ndarray)):
        raise ValueError("Charges must be a list, tuple, or numpy array")
    
    charges_array = np.array(Charges, dtype=float)
    if len(charges_array) != atoms_per_particle:
        raise ValueError(f"Number of charges ({len(charges_array)}) must match atoms_per_particle ({atoms_per_particle})")
    
    # Convert to debye units (charge * distance conversion factor)
    DEBYE_CONVERSION = 0.2081943  # e*Å to Debye
    
    # Create full paths and determine output mode
    vectors_as_files = magnitudes_pattern is not None
    
    if vectors_as_files:
        # New mode: separate file patterns for vectors and magnitudes
        print(f"Output mode: Separate file patterns")
        print(f"Vectors pattern: {output_pattern}")
        print(f"Magnitudes pattern: {magnitudes_pattern}")
    else:
        # Legacy mode: directory pattern with hardcoded filenames
        output_dir_rel = expand_path_pattern(output_pattern, common_term)
        output_path = os.path.join(baseDir, output_dir_rel)
        os.makedirs(output_path, exist_ok=True)
        print(f"Output mode: Directory pattern (legacy)")
        print(f"Output directory: {output_path}")
    
    print(f"{'='*50}")
    print(f"DIPOLE MOMENT CALCULATION")
    print(f"{'='*50}")
    print(f"Base directory: {baseDir}")
    print(f"Coords pattern: {coords_pattern}")
    print(f"COM pattern: {com_pattern}")
    print(f"Common term: {common_term}")
    print(f"Processing {effective_molecules} out of {num_particles} total molecules")
    print(f"Using {max_workers} workers, atoms per particle: {atoms_per_particle}")
    print(f"Charges: {charges_array.tolist()}")
    
    # Handle DCD selection
    if dcd_indices is not None:
        actual_dcd_list = dcd_indices
        actual_num_dcds = len(dcd_indices)
        print(f"Using DCD selection: {dcd_indices}")
    else:
        actual_dcd_list = list(range(num_dcds))
        actual_num_dcds = num_dcds
        print(f"Processing all DCDs: 0 to {num_dcds-1}")
    




    # Process files
    results = []
    successful_files = 0
    total_frames = 0
    quality_metrics = {'mean_magnitudes': [], 'std_magnitudes': []}
    
    print(f"Processing {actual_num_dcds} trajectory files for dipole moment calculation...")
    print(f"Using {max_workers} workers, atoms per particle: {atoms_per_particle}")
    print(f"Processing {effective_molecules} out of {num_particles} total molecules")
    print(f"Charges: {charges_array.tolist()}")
    
    if max_workers > 1:
        # Try parallel processing with fallback to single-threaded
        try:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Submit all tasks
                future_to_idx = {executor.submit(_process_single_dipole_file, i, baseDir, coords_pattern, 
                                                 com_pattern, output_pattern, magnitudes_pattern, 
                                                 charges_array, atoms_per_particle, effective_molecules, 
                                                 stride, common_term): i for i in actual_dcd_list}
                
                # Process completed tasks
                for future in as_completed(future_to_idx):
                    result = future.result()
                    results.append(result)
        except (TypeError, AttributeError) as e:
            if 'pickle' in str(e).lower() or 'local object' in str(e):
                print(f"⚠️  Multiprocessing failed due to function pickling: {e}")
                print("   Falling back to single-threaded processing...")
                max_workers = 1  # Force single-threaded fallback
                # Process files sequentially
                for i in actual_dcd_list:
                    result = _process_single_dipole_file(i, baseDir, coords_pattern, com_pattern, 
                                                       output_pattern, magnitudes_pattern, charges_array, 
                                                       atoms_per_particle, effective_molecules, stride, common_term)
                    results.append(result)
            else:
                raise e  # Re-raise other errors
                
                if result['success']:
                    successful_files += 1
                    total_frames += result['frames_processed']
                    quality_metrics['mean_magnitudes'].append(result['mean_magnitude'])
                    quality_metrics['std_magnitudes'].append(result['std_magnitude'])
                    print(f"✓ File {result['file_idx']}: {result['frames_processed']} frames, avg magnitude: {result['mean_magnitude']:.3f} D")
                else:
                    print(f"✗ File {result['file_idx']}: {result['error']}")
                
                if progress_callback:
                    progress_callback(len(results), actual_num_dcds)
    else:
        # Sequential processing
        for i in actual_dcd_list:
            result = _process_single_dipole_file(i, baseDir, coords_pattern, com_pattern, 
                                               output_pattern, magnitudes_pattern, charges_array, 
                                               atoms_per_particle, effective_molecules, stride, common_term)
            results.append(result)
            
            if result['success']:
                successful_files += 1
                total_frames += result['frames_processed']
                quality_metrics['mean_magnitudes'].append(result['mean_magnitude'])
                quality_metrics['std_magnitudes'].append(result['std_magnitude'])
                print(f"✓ File {result['file_idx']}: {result['frames_processed']} frames, avg magnitude: {result['mean_magnitude']:.3f} D")
            else:
                print(f"✗ File {result['file_idx']}: {result['error']}")
            
            if progress_callback:
                progress_callback(i + 1, actual_num_dcds)
    
    # Calculate summary statistics
    total_time = time.time() - start_time
    
    summary = {
        'success': successful_files,
        'total': actual_num_dcds,
        'total_time': total_time,
        'total_frames': total_frames,
        'avg_time_per_file': total_time / actual_num_dcds if actual_num_dcds > 0 else 0,
        'data_quality': 'Good' if successful_files == actual_num_dcds else f'Partial ({successful_files}/{actual_num_dcds})',
        'output_directory': output_path if not vectors_as_files else None # Only show if not vectors_as_files
    }
    
    if quality_metrics['mean_magnitudes']:
        summary.update({
            'overall_mean_magnitude': np.mean(quality_metrics['mean_magnitudes']),
            'overall_std_magnitude': np.mean(quality_metrics['std_magnitudes']),
            'magnitude_range': [np.min(quality_metrics['mean_magnitudes']), np.max(quality_metrics['mean_magnitudes'])]
        })
    
    print(f"\nDipole calculation completed:")
    print(f"  Successful files: {successful_files}/{actual_num_dcds}")
    print(f"  Total frames processed: {total_frames}")
    print(f"  Total time: {total_time:.2f}s")
    if quality_metrics['mean_magnitudes']:
        print(f"  Average dipole magnitude: {summary['overall_mean_magnitude']:.3f} ± {summary['overall_std_magnitude']:.3f} D")
    
    return summary

# Alias for compatibility
dipole_calculation = dipole_functions


