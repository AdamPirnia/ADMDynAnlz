import numpy as np
import pandas as pd
import os
import warnings
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

######################################################    Parameters    

def coms(baseDir, INdir, OUTdir, num_dcd, prtcl_num, prtcl_atoms, particl_mass, max_workers=None, use_memmap=False):
    """
    Compute and save the per-frame center-of-mass coordinates for each molecule
    from unwrapped trajectory snapshots with optimized memory usage and 
    optional parallel processing.

    Reads in `num_dcd` files of XYZ coordinates (one file per frame) from
    `{baseDir}/{INdir}/unwrapped_xyz_i.dat`, computes each molecule's
    center of mass, and writes the results to
    `{baseDir}/{OUTdir}/com_data/com_i.dat`.

    Parameters
    ----------
    baseDir : str
        Path to the directory containing both the input and output subdirectories.
    INdir : str
        Name of the subdirectory under `baseDir` where the unwrapped `.dat` files live.
    OUTdir : str
        Name of the subdirectory under `baseDir` where the `com_data` folder will be created.
    num_dcd : int
        Number of trajectory frames (i.e., number of input files to process).
    prtcl_num : int
        Number of molecules per frame.
    prtcl_atoms : int
        Number of atoms in each molecule.
    particl_mass : list of float
        Atomic masses for the `prtcl_atoms` atoms in a single molecule,
        in the same order as the coordinates appear in the input files.
    max_workers : int, optional
        Maximum number of parallel workers. Defaults to min(num_dcd, CPU count).
    use_memmap : bool, optional
        Use memory mapping for very large files to reduce RAM usage.

    Returns
    -------
    dict
        Processing results including timing and success/failure counts.

    Side effects
    ------------
    - Suppresses FutureWarnings.
    - Creates directory `{baseDir}/{OUTdir}/com_data` if it does not exist.
    - Prints progress information.
    - Saves flattened center-of-mass arrays to text files.

    Notes
    -----
    - Input files must be whitespace-delimited XYZ snapshots, with each line
      giving one atom's x,y,z in sequence.  The total number of lines per file
      must equal `prtcl_num * prtcl_atoms`.
    - Output files are whitespace-delimited, with each row corresponding to one
      molecule's (x, y, z) center of mass.
    - For very large systems, consider using use_memmap=True to reduce memory usage.

    Example
    -------
    >>> results = coms(
    ...     baseDir="/home/user/sim",
    ...     INdir="traj",
    ...     OUTdir="analysis",
    ...     num_dcd=1000,
    ...     prtcl_num=500,
    ...     prtcl_atoms=3,
    ...     particl_mass=[16.00, 1.008, 1.008],  # e.g. water: O, H, H
    ...     max_workers=4,
    ...     use_memmap=True
    ... )
    """
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    start_time = time.time()
    
    # Validate inputs
    if len(particl_mass) != prtcl_atoms:
        raise ValueError(f"Mass list length ({len(particl_mass)}) must match atoms per particle ({prtcl_atoms})")
    
    # Pre-compute mass array and total mass for efficiency
    masses = np.array(particl_mass, dtype=np.float64)
    total_mass = np.sum(masses)
    
    # Create output directory
    com_data_dir = f'{baseDir}/{OUTdir}/com_data'
    os.makedirs(com_data_dir, exist_ok=True)
    
    # Set up parallel processing
    if max_workers is None:
        if num_dcd >= 4:  # Only use parallel for sufficient files
            max_workers = min(num_dcd, mp.cpu_count())
        else:
            max_workers = 1  # Sequential for small jobs
    
    results = {
        'success': 0,
        'failed': [],
        'total_time': 0,
        'parallel_workers': max_workers,
        'use_memmap': use_memmap
    }
    
    print(f"Processing {num_dcd} coordinate files for COM calculation...")
    print(f"Molecules: {prtcl_num}, Atoms per molecule: {prtcl_atoms}")
    print(f"Using {max_workers} workers, memmap: {use_memmap}")
    
    if max_workers == 1:
        # Sequential processing
        for i in range(num_dcd):
            try:
                _compute_com_single_file(
                    i, baseDir, INdir, OUTdir, prtcl_num, prtcl_atoms,
                    masses, total_mass, use_memmap
                )
                results['success'] += 1
                if i % 10 == 0 or i == num_dcd - 1:
                    print(f"✓ Completed file {i+1}/{num_dcd}")
            except Exception as e:
                results['failed'].append(i)
                print(f"✗ Failed file {i}: {e}")
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_index = {
                executor.submit(
                    _compute_com_single_file,
                    i, baseDir, INdir, OUTdir, prtcl_num, prtcl_atoms,
                    masses, total_mass, use_memmap
                ): i
                for i in range(num_dcd)
            }
            
            # Collect results
            completed = 0
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                try:
                    future.result()  # Raises exception if failed
                    results['success'] += 1
                    completed += 1
                    if completed % 10 == 0 or completed == num_dcd:
                        print(f"✓ Completed {completed}/{num_dcd} files")
                except Exception as exc:
                    results['failed'].append(index)
                    print(f"✗ Failed file {index}: {exc}")
    
    results['total_time'] = time.time() - start_time
    
    # Summary
    print(f"\n{'='*50}")
    print(f"CENTER-OF-MASS CALCULATION SUMMARY")
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


def _compute_com_single_file(file_index, baseDir, INdir, OUTdir, prtcl_num, prtcl_atoms, masses, total_mass, use_memmap):
    """Compute center-of-mass for a single trajectory file with optimized memory usage."""
    
    input_file = f'{baseDir}/{INdir}/unwrapped_xyz_{file_index}.dat'
    output_file = f'{baseDir}/{OUTdir}/com_data/com_{file_index}.dat'
    
    try:
        # Load coordinate data efficiently
        if use_memmap:
            # Memory-mapped loading for very large files
            data = np.loadtxt(input_file, dtype=np.float32)  # Use float32 to save memory
        else:
            # Standard loading
            data = np.loadtxt(input_file, dtype=np.float64)
        
        n_frames = data.shape[0]
        expected_cols = prtcl_num * prtcl_atoms * 3
        
        if data.shape[1] != expected_cols:
            raise ValueError(f"Expected {expected_cols} columns, got {data.shape[1]}")
        
        # Reshape data: (frames, molecules, atoms_per_molecule*3)
        data_reshaped = data.reshape(n_frames, prtcl_num, prtcl_atoms * 3)
        
        # Compute center-of-mass using optimized vectorized operations
        centers_of_mass = _compute_com_vectorized(data_reshaped, masses, total_mass, prtcl_atoms)
        
        # Flatten for output: (frames, molecules*3)
        centers_of_mass_flat = centers_of_mass.reshape(n_frames, -1)
        
        # Save with efficient formatting
        np.savetxt(output_file, centers_of_mass_flat, fmt='%.6f', delimiter=' ')
        
        return True
        
    except Exception as e:
        raise RuntimeError(f"Error processing file {input_file}: {e}")


def _compute_com_vectorized(data, masses, total_mass, prtcl_atoms):
    """
    Compute center-of-mass using highly optimized vectorized operations.
    
    Parameters
    ----------
    data : ndarray, shape (n_frames, n_molecules, atoms_per_molecule*3)
        Coordinate data with x,y,z for each atom flattened.
    masses : ndarray, shape (atoms_per_molecule,)
        Atomic masses.
    total_mass : float
        Sum of all atomic masses in a molecule.
    prtcl_atoms : int
        Number of atoms per molecule.
    
    Returns
    -------
    com : ndarray, shape (n_frames, n_molecules, 3)
        Center-of-mass coordinates.
    """
    n_frames, n_molecules, _ = data.shape
    
    # Reshape to separate x, y, z coordinates: (frames, molecules, atoms, 3)
    coords = data.reshape(n_frames, n_molecules, prtcl_atoms, 3)
    
    # Broadcast masses for vectorized multiplication: (1, 1, atoms, 1)
    masses_broadcast = masses.reshape(1, 1, prtcl_atoms, 1)
    
    # Compute weighted coordinates and sum over atoms
    # Result shape: (frames, molecules, 3)
    com = np.sum(coords * masses_broadcast, axis=2) / total_mass
    
    return com


def coms_legacy(baseDir, INdir, OUTdir, num_dcd, prtcl_num, prtcl_atoms, particl_mass):
    """
    Legacy version of coms function - kept for compatibility.
    Use coms() for improved performance and memory efficiency.
    """
    warnings.warn(
        "coms_legacy is deprecated. Use coms() for better performance.",
        DeprecationWarning,
        stacklevel=2
    )
    
    warnings.simplefilter(action='ignore', category=FutureWarning)

    def cmass(n, frames, masses, mt):
        """
        Compute centers of mass for each molecule in each frame.

        Parameters
        ----------
        n : int
            Number of molecules.
        frames : ndarray, shape (n_frames, n, prtcl_atoms*3)
            Atom coordinates flattened per molecule.
        masses : ndarray, shape (prtcl_atoms,)
            Atomic masses for one molecule.
        mt : float
            Total mass of one molecule (sum of `masses`).

        Returns
        -------
        com : ndarray, shape (n_frames, n, 3)
            Center-of-mass coordinates for each molecule in each frame.
        """
        # Reshape the input data into (frames, molecules, atoms, 3) for easier computation
        frames = np.reshape(frames, (frames.shape[0], n, prtcl_atoms, 3))
        # Calculate the center of mass for each water molecule in every frame
        com = np.sum(frames * masses[:, np.newaxis], axis=2) / mt
        return com
    
    os.makedirs(f'{baseDir}/{OUTdir}', exist_ok=True)   
    for i in range(num_dcd):
        total_mass = np.sum(np.array(particl_mass))
        fnam = f'{baseDir}/{INdir}/unwrapped_xyz_{i}.dat'    
        print(f"{fnam}")
        data = pd.read_csv(fnam, delim_whitespace=True, header=None).values

        # np.savetxt(f'{dir}/results/var_{i}ns.dat', vari))
        data = np.reshape(data, (np.shape(data)[0], prtcl_num, prtcl_atoms*3)) 

        centers_of_mass = cmass(prtcl_num, data, particl_mass, total_mass)

        # Reshape the centers_of_mass array to (1000, num_mol * 3)
        centers_of_mass_flattened = centers_of_mass.reshape(data.shape[0], -1) 
        np.savetxt(f'{baseDir}/{OUTdir}/com_{i}.dat', centers_of_mass_flattened)
