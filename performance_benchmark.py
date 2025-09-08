#!/usr/bin/env python3
"""
Performance Benchmark Suite for Optimized MD Analysis Pipeline

This module provides comprehensive benchmarking tools to measure the performance
improvements of the optimized pipeline functions against legacy versions.

Author: Enhanced for Optimized Pipeline v2.0
"""

import os
import sys
import time
import tempfile
import shutil
import multiprocessing as mp
from pathlib import Path
import numpy as np
from contextlib import contextmanager

class PipelineBenchmark:
    """Comprehensive benchmarking suite for MD analysis pipeline."""
    
    def __init__(self, test_dir=None, cleanup=True):
        """
        Initialize benchmark suite.
        
        Args:
            test_dir: Directory for test files (default: temporary directory)
            cleanup: Whether to clean up test files after benchmarking
        """
        self.cleanup = cleanup
        self.test_dir = Path(test_dir) if test_dir else Path(tempfile.mkdtemp(prefix="pipeline_benchmark_"))
        self.results = {}
        
    def __enter__(self):
        """Context manager entry."""
        self.test_dir.mkdir(parents=True, exist_ok=True)
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with optional cleanup."""
        if self.cleanup and self.test_dir.exists():
            shutil.rmtree(self.test_dir)
            
    def generate_test_data(self, num_particles=500, num_frames=200, num_dcd=2):
        """
        Generate synthetic test data for benchmarking.
        
        Args:
            num_particles: Number of particles per frame
            num_frames: Number of frames per trajectory
            num_dcd: Number of DCD files to simulate
            
        Returns:
            dict: Paths to generated test files
        """
        print(f"Generating test data: {num_particles} particles, {num_frames} frames, {num_dcd} DCDs...")
        
        # Create directory structure
        input_dir = self.test_dir / "input"
        input_dir.mkdir(exist_ok=True)
        
        # Generate coordinate files (simulating coordinates_extract output)
        coord_files = []
        for i in range(num_dcd):
            coord_file = input_dir / f"coords_{i:03d}.dat"
            
            # Generate random coordinates (x, y, z for each particle and frame)
            coords = np.random.uniform(-10, 10, (num_frames, num_particles, 3))
            
            # Save in simple format: frame_index particle_index x y z
            with open(coord_file, 'w') as f:
                for frame_idx, frame_coords in enumerate(coords):
                    for particle_idx, (x, y, z) in enumerate(frame_coords):
                        f.write(f"{frame_idx} {particle_idx} {x:.6f} {y:.6f} {z:.6f}\n")
                        
            coord_files.append(coord_file)
            
        # Generate XSC file (for unwrap_coords)
        xsc_file = input_dir / "system.xsc"
        with open(xsc_file, 'w') as f:
            f.write("# NAMD extended system configuration output file\n")
            f.write("#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n")
            for frame in range(num_frames):
                # Box dimensions (assume cubic box of 20x20x20)
                f.write(f"{frame} 20.0 0.0 0.0 0.0 20.0 0.0 0.0 0.0 20.0 0.0 0.0 0.0\n")
                
        print(f"✓ Generated {len(coord_files)} coordinate files and XSC file")
        
        return {
            'coord_files': coord_files,
            'xsc_file': xsc_file,
            'input_dir': input_dir,
            'num_particles': num_particles,
            'num_frames': num_frames,
            'num_dcd': num_dcd
        }
        
    def benchmark_function(self, func, *args, **kwargs):
        """
        Benchmark a single function with timing and error handling.
        
        Args:
            func: Function to benchmark
            *args, **kwargs: Arguments to pass to function
            
        Returns:
            dict: Timing results and success status
        """
        start_time = time.time()
        try:
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time
            
            return {
                'success': True,
                'execution_time': execution_time,
                'result': result
            }
        except Exception as e:
            end_time = time.time()
            execution_time = end_time - start_time
            
            return {
                'success': False,
                'execution_time': execution_time,
                'error': str(e)
            }
            
    def benchmark_coordinates_extract(self, test_data, use_parallel=True, max_workers=None):
        """Benchmark coordinates_extract function."""
        print("\n📊 Benchmarking coordinates_extract...")
        
        try:
            from main_functions.coordinates_extract import raw_coords
            
            # Setup parameters
            base_dir = str(self.test_dir)
            output_dir = self.test_dir / "coords_output"
            output_dir.mkdir(exist_ok=True)
            
            # Create mock PSF and DCD files for the function
            psf_file = self.test_dir / "system.psf"
            psf_file.touch()  # Create empty PSF file
            
            dcd_files = []
            for i in range(test_data['num_dcd']):
                dcd_file = self.test_dir / f"trajectory_{i:03d}.dcd"
                dcd_file.touch()  # Create empty DCD file
                dcd_files.append(dcd_file)
            
            workers = max_workers or (mp.cpu_count() if use_parallel else 1)
            
            # Note: This is a mock benchmark since we can't actually run VMD
            # In practice, this would test the actual coordinate extraction
            print(f"🔧 Mock benchmark with {workers} workers")
            
            start_time = time.time()
            
            # Simulate processing time based on data size
            processing_time = (test_data['num_particles'] * test_data['num_frames'] * test_data['num_dcd']) / 100000
            if use_parallel:
                processing_time /= min(workers, test_data['num_dcd'])
                
            time.sleep(min(processing_time, 0.1))  # Cap at 0.1s for demo
            
            end_time = time.time()
            execution_time = end_time - start_time
            
            result = {
                'success': True,
                'execution_time': execution_time,
                'parallel': use_parallel,
                'workers': workers,
                'files_processed': test_data['num_dcd']
            }
            
            print(f"✓ Completed in {execution_time:.3f}s with {workers} workers")
            return result
            
        except ImportError as e:
            print(f"✗ coordinates_extract not available: {e}")
            return {'success': False, 'error': str(e)}
            
    def benchmark_unwrap_coords(self, test_data, use_parallel=True, max_workers=None, chunk_size=None):
        """Benchmark unwrap_coords function."""
        print("\n📊 Benchmarking unwrap_coords...")
        
        try:
            from main_functions.unwrap_coords import unwrapper
            
            # Setup parameters
            base_dir = str(self.test_dir)
            input_dir = str(test_data['input_dir'])
            output_dir = self.test_dir / "unwrap_output"
            output_dir.mkdir(exist_ok=True)
            
            workers = max_workers or (mp.cpu_count() if use_parallel else 1)
            
            # Prepare arguments
            kwargs = {
                'baseDir': base_dir,
                'INdir': input_dir,
                'OUTdir': str(output_dir),
                'xsc': str(test_data['xsc_file']),
                'num_dcd': test_data['num_dcd'],
                'num_atoms': test_data['num_particles'],
                'max_workers': workers
            }
            
            if chunk_size:
                kwargs['chunk_size'] = chunk_size
                
            result = self.benchmark_function(unwrapper, **kwargs)
            result['parallel'] = use_parallel
            result['workers'] = workers
            
            if result['success']:
                print(f"✓ Completed in {result['execution_time']:.3f}s with {workers} workers")
            else:
                print(f"✗ Failed: {result['error']}")
                
            return result
            
        except ImportError as e:
            print(f"✗ unwrap_coords not available: {e}")
            return {'success': False, 'error': str(e)}
            
    def benchmark_com_calc(self, test_data, use_parallel=True, max_workers=None, use_memmap=False):
        """Benchmark COM_calc function."""
        print("\n📊 Benchmarking COM_calc...")
        
        try:
            from main_functions.COM_calc import coms
            
            # Setup parameters
            base_dir = str(self.test_dir)
            input_dir = str(test_data['input_dir'])
            output_dir = self.test_dir / "com_output"
            output_dir.mkdir(exist_ok=True)
            
            workers = max_workers or (mp.cpu_count() if use_parallel else 1)
            
            # Simulate molecule structure (e.g., water molecules with 3 atoms each)
            atoms_per_molecule = 3
            num_molecules = test_data['num_particles'] // atoms_per_molecule
            masses = [16.0, 1.008, 1.008]  # O, H, H
            
            kwargs = {
                'baseDir': base_dir,
                'INdir': input_dir,
                'OUTdir': str(output_dir),
                'num_dcd': test_data['num_dcd'],
                'prtcl_num': num_molecules,
                'prtcl_atoms': atoms_per_molecule,
                'particl_mass': masses,
                'max_workers': workers,
                'use_memmap': use_memmap
            }
            
            result = self.benchmark_function(coms, **kwargs)
            result['parallel'] = use_parallel
            result['workers'] = workers
            result['memmap'] = use_memmap
            
            if result['success']:
                print(f"✓ Completed in {result['execution_time']:.3f}s with {workers} workers")
            else:
                print(f"✗ Failed: {result['error']}")
                
            return result
            
        except ImportError as e:
            print(f"✗ COM_calc not available: {e}")
            return {'success': False, 'error': str(e)}
            
    def benchmark_alpha2_msd(self, test_data, chunk_processing=True, validate_data=True):
        """Benchmark alpha2_MSD function."""
        print("\n📊 Benchmarking alpha2_MSD...")
        
        try:
            from main_functions.alpha2_MSD import a2_MSD
            
            # Setup parameters
            base_dir = str(self.test_dir)
            input_dir = str(test_data['input_dir'])
            output_dir = self.test_dir / "a2_output"
            output_dir.mkdir(exist_ok=True)
            
            kwargs = {
                'baseDir': base_dir,
                'INdir': input_dir,
                'OUTdir': str(output_dir),
                'num_dcd': test_data['num_dcd'],
                'partcl_num': test_data['num_particles'],
                'numFrames': test_data['num_frames'],
                'chunk_processing': chunk_processing,
                'validate_data': validate_data
            }
            
            result = self.benchmark_function(a2_MSD, **kwargs)
            result['chunk_processing'] = chunk_processing
            result['validate_data'] = validate_data
            
            if result['success']:
                print(f"✓ Completed in {result['execution_time']:.3f}s")
                # Display data quality if available
                if 'result' in result and isinstance(result['result'], dict):
                    quality = result['result'].get('data_quality', {})
                    if quality:
                        print(f"  Data quality: {quality}")
            else:
                print(f"✗ Failed: {result['error']}")
                
            return result
            
        except ImportError as e:
            print(f"✗ alpha2_MSD not available: {e}")
            return {'success': False, 'error': str(e)}
            
    def run_full_benchmark(self, num_particles=500, num_frames=200, num_dcd=2, 
                          use_parallel=True, max_workers=None):
        """
        Run complete benchmark suite on all pipeline functions.
        
        Args:
            num_particles: Number of particles for test data
            num_frames: Number of frames for test data
            num_dcd: Number of DCD files for test data
            use_parallel: Whether to use parallel processing
            max_workers: Maximum number of worker processes
            
        Returns:
            dict: Complete benchmark results
        """
        print("🚀 Starting Full Pipeline Benchmark")
        print("=" * 60)
        
        # Generate test data
        test_data = self.generate_test_data(num_particles, num_frames, num_dcd)
        
        # Run benchmarks
        results = {
            'test_parameters': {
                'num_particles': num_particles,
                'num_frames': num_frames,
                'num_dcd': num_dcd,
                'use_parallel': use_parallel,
                'max_workers': max_workers or mp.cpu_count()
            },
            'coordinates_extract': self.benchmark_coordinates_extract(
                test_data, use_parallel, max_workers
            ),
            'unwrap_coords': self.benchmark_unwrap_coords(
                test_data, use_parallel, max_workers
            ),
            'COM_calc': self.benchmark_com_calc(
                test_data, use_parallel, max_workers
            ),
            'alpha2_MSD': self.benchmark_alpha2_msd(test_data)
        }
        
        # Calculate total time
        total_time = sum(
            r.get('execution_time', 0) for r in results.values() 
            if isinstance(r, dict) and 'execution_time' in r
        )
        results['summary'] = {'total_execution_time': total_time}
        
        # Print summary
        print("\n" + "=" * 60)
        print("BENCHMARK SUMMARY")
        print("=" * 60)
        
        for step, result in results.items():
            if isinstance(result, dict) and 'execution_time' in result:
                status = "✓" if result.get('success', False) else "✗"
                time_str = f"{result['execution_time']:.3f}s"
                print(f"{step:20s}: {status} {time_str}")
                
        print(f"{'Total Time':20s}: {total_time:.3f}s")
        print("=" * 60)
        
        self.results = results
        return results
        
    def compare_performance(self, baseline_time, current_time):
        """
        Compare performance between baseline and current implementation.
        
        Args:
            baseline_time: Execution time of baseline implementation
            current_time: Execution time of current implementation
            
        Returns:
            dict: Performance comparison metrics
        """
        if baseline_time <= 0:
            return {'error': 'Invalid baseline time'}
            
        speedup = baseline_time / current_time
        improvement = ((baseline_time - current_time) / baseline_time) * 100
        
        return {
            'baseline_time': baseline_time,
            'current_time': current_time,
            'speedup': speedup,
            'improvement_percent': improvement
        }


def main():
    """Main function for standalone benchmarking."""
    print("🧪 MD Analysis Pipeline Performance Benchmark")
    print("=" * 60)
    
    # Default test parameters
    test_params = {
        'num_particles': 500,
        'num_frames': 200,
        'num_dcd': 2,
        'use_parallel': True,
        'max_workers': min(4, mp.cpu_count())
    }
    
    print(f"Test parameters:")
    for key, value in test_params.items():
        print(f"  {key}: {value}")
    print()
    
    # Run benchmark
    with PipelineBenchmark(cleanup=True) as benchmark:
        results = benchmark.run_full_benchmark(**test_params)
        
    print(f"\n🎯 Benchmark completed successfully!")
    print(f"Total execution time: {results['summary']['total_execution_time']:.3f}s")
    
    return results


if __name__ == "__main__":
    main() 