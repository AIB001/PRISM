#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Management System Demonstration

This script demonstrates the comprehensive data management capabilities
of the PRISM system including storage, processing, caching, and backup.
"""

import sys
import time
import numpy as np
from pathlib import Path

# Add PRISM to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from prism.data import (
    DataStorageManager, 
    StorageTier, 
    CompressionType,
    create_stream_processor,
    create_batch_processor,
    TrajectoryProcessor,
    EnergyProcessor,
    ProcessingMode,
    ProcessingJob,
    create_cache_manager,
    CachePolicy,
    TrajectoryHandler,
    EnergyHandler,
    PMFDataHandler,
    FrameData,
    EnergyData,
    PMFResult,
    create_backup_manager,
    BackupType
)
from prism.utils.logging_system import setup_logging


def demonstrate_storage_system():
    """Demonstrate the data storage system"""
    print("\n" + "="*50)
    print("DEMONSTRATING STORAGE SYSTEM")
    print("="*50)
    
    # Create storage manager
    storage_path = Path("demo_storage")
    storage = DataStorageManager(storage_path, enable_hierarchical=True, enable_compression=True)
    
    # Store different types of data
    print("\n1. Storing various data types:")
    
    # Store JSON data in hot tier
    config_data = {
        "simulation_type": "PMF",
        "temperature": 300.0,
        "steps": 1000000,
        "settings": {"dt": 0.002, "output_freq": 1000}
    }
    
    config_path = storage.store_data(
        "simulation_config", 
        config_data, 
        storage_tier=StorageTier.HOT,
        custom_metadata={"type": "configuration", "version": "1.0"}
    )
    print(f"  ✓ Stored configuration data: {config_path}")
    
    # Store large binary data in warm tier
    trajectory_data = np.random.random((1000, 100, 3))  # Mock trajectory
    traj_path = storage.store_data(
        "trajectory_001", 
        trajectory_data,
        content_type="trajectory",
        storage_tier=StorageTier.WARM,
        compression=CompressionType.LZMA
    )
    print(f"  ✓ Stored trajectory data: {traj_path}")
    
    # Store analysis results
    analysis_results = {
        "binding_energy": -8.5,
        "rmsd_values": np.random.random(1000).tolist(),
        "convergence": True
    }
    
    analysis_path = storage.store_data(
        "pmf_analysis", 
        analysis_results,
        storage_tier=StorageTier.HOT
    )
    print(f"  ✓ Stored analysis results: {analysis_path}")
    
    # Retrieve data
    print("\n2. Retrieving stored data:")
    
    retrieved_config, config_meta = storage.retrieve_data("simulation_config")
    print(f"  ✓ Retrieved config: {retrieved_config['simulation_type']}")
    print(f"    Access count: {config_meta.access_count}")
    
    retrieved_traj, traj_meta = storage.retrieve_data("trajectory_001")
    print(f"  ✓ Retrieved trajectory: shape {retrieved_traj.shape}")
    print(f"    Compression ratio: {traj_meta.size_bytes / len(retrieved_traj.tobytes()):.2f}x")
    
    # Show storage statistics
    stats = storage.get_storage_stats()
    print(f"\n3. Storage Statistics:")
    print(f"  Total objects: {stats['total_objects']}")
    print(f"  Total size: {stats['total_size_mb']:.1f} MB")
    
    for tier, tier_stats in stats['tier_statistics'].items():
        print(f"  {tier}: {tier_stats['objects']} objects, {tier_stats['size_bytes']/(1024*1024):.1f} MB")
    
    return storage


def demonstrate_processing_system():
    """Demonstrate the data processing system"""
    print("\n" + "="*50)
    print("DEMONSTRATING PROCESSING SYSTEM")
    print("="*50)
    
    # Create sample data
    n_frames = 100
    n_atoms = 50
    trajectory_data = np.random.random((n_frames, n_atoms, 3)) * 5.0
    
    print(f"\n1. Processing trajectory data ({n_frames} frames, {n_atoms} atoms)")
    
    # Stream processing
    print("\n  Stream Processing:")
    stream_processor = create_stream_processor(buffer_size=10, max_workers=2)
    traj_processor = TrajectoryProcessor()
    
    stream_processor.start_processing(traj_processor, ProcessingMode.THREADED)
    
    # Submit jobs
    jobs_submitted = 0
    for i in range(5):
        # Create frame subset
        frame_subset = trajectory_data[i*10:(i+1)*10]
        
        job = ProcessingJob(
            job_id=f"traj_job_{i}",
            input_data=frame_subset,
            processor_kwargs={'dt': 0.002, 'masses': np.ones(n_atoms)}
        )
        
        if stream_processor.submit_job(job):
            jobs_submitted += 1
    
    print(f"    Submitted {jobs_submitted} jobs")
    
    # Collect results
    results_collected = 0
    start_time = time.time()
    
    while results_collected < jobs_submitted and time.time() - start_time < 10:
        result = stream_processor.get_result(timeout=1.0)
        if result:
            results_collected += 1
            print(f"    ✓ Job {result.job_id}: {result.status} "
                 f"({result.processing_time:.3f}s)")
            if result.output_data:
                print(f"      RMSD range: {min(result.output_data['rmsd']):.3f} - "
                     f"{max(result.output_data['rmsd']):.3f}")
    
    stream_processor.shutdown()
    
    # Batch processing
    print("\n  Batch Processing:")
    batch_processor = create_batch_processor(max_workers=2, chunk_size=10)
    
    # Create energy data for processing
    energy_data_list = []
    for i in range(10):
        energy_dict = {
            'Potential': np.random.normal(-50000, 1000, 100),
            'Kinetic': np.random.normal(25000, 500, 100),
            'Temperature': np.random.normal(300, 10, 100)
        }
        energy_data_list.append(energy_dict)
    
    energy_processor = EnergyProcessor()
    
    def progress_callback(completed, total):
        print(f"    Progress: {completed}/{total} ({100*completed/total:.1f}%)")
    
    batch_results = batch_processor.process_batch(
        energy_data_list, 
        energy_processor, 
        ProcessingMode.THREADED,
        progress_callback
    )
    
    successful_results = [r for r in batch_results if r.status == "completed"]
    print(f"    ✓ Processed {len(successful_results)}/{len(batch_results)} energy datasets")
    
    if successful_results:
        sample_result = successful_results[0].output_data
        print(f"    Sample result - Potential mean: {sample_result['Potential']['mean']:.1f}")
    
    return stream_processor, batch_processor


def demonstrate_caching_system():
    """Demonstrate the caching system"""
    print("\n" + "="*50)
    print("DEMONSTRATING CACHING SYSTEM")
    print("="*50)
    
    # Create multi-level cache
    cache_manager = create_cache_manager(l1_size_mb=10, l2_size_mb=50)
    
    print("\n1. Caching different data types:")
    
    # Cache simulation parameters
    sim_params = {
        "force_field": "amber99sb-ildn",
        "water_model": "tip3p",
        "temperature": 300.0,
        "pressure": 1.0,
        "box_size": [5.0, 5.0, 5.0]
    }
    
    cache_manager.put("simulation_parameters", sim_params, ttl=3600)
    print("  ✓ Cached simulation parameters")
    
    # Cache computed results
    pmf_results = {
        "binding_energy": -8.5,
        "pmf_profile": np.random.random(50).tolist(),
        "error_estimates": np.random.random(50).tolist()
    }
    
    cache_manager.put("pmf_calculation", pmf_results, ttl=1800)
    print("  ✓ Cached PMF results")
    
    # Cache large trajectory subset
    large_data = np.random.random((500, 100, 3))
    cache_manager.put("trajectory_subset", large_data, ttl=7200)
    print("  ✓ Cached trajectory subset")
    
    print("\n2. Cache retrieval and performance:")
    
    # Retrieve cached data
    start_time = time.time()
    retrieved_params = cache_manager.get("simulation_parameters")
    retrieval_time = time.time() - start_time
    
    if retrieved_params:
        print(f"  ✓ Retrieved parameters in {retrieval_time*1000:.2f}ms")
        print(f"    Force field: {retrieved_params['force_field']}")
    
    # Test cache miss
    start_time = time.time()
    missing_data = cache_manager.get("nonexistent_key")
    miss_time = time.time() - start_time
    
    print(f"  ✓ Cache miss handled in {miss_time*1000:.2f}ms")
    
    # Show cache statistics
    stats = cache_manager.get_stats()
    print(f"\n3. Cache Statistics:")
    print(f"  Total hit ratio: {stats['total_hit_ratio']:.2%}")
    
    for level_name, level_stats in stats['levels'].items():
        print(f"  {level_name}: {level_stats['entry_count']} entries, "
             f"{level_stats['size_mb']:.1f}MB, "
             f"{level_stats['hit_ratio']:.2%} hit ratio")
    
    return cache_manager


def demonstrate_format_handlers():
    """Demonstrate data format handlers"""
    print("\n" + "="*50)
    print("DEMONSTRATING FORMAT HANDLERS")
    print("="*50)
    
    # Create demo directory
    demo_dir = Path("demo_formats")
    demo_dir.mkdir(exist_ok=True)
    
    print("\n1. Creating sample data files:")
    
    # Create sample trajectory data
    trajectory_handler = TrajectoryHandler()
    
    # Generate sample frames
    frames = []
    for i in range(5):
        coords = np.random.random((10, 3)) * 5.0
        frame = FrameData(
            frame_number=i,
            time_ps=i * 2.0,
            coordinates=coords,
            box_vectors=np.eye(3) * 5.0
        )
        frames.append(frame)
    
    # Write GRO file
    gro_file = demo_dir / "sample.gro"
    success = trajectory_handler.write(
        frames, 
        gro_file,
        title="Demo trajectory",
        atom_names=["CA"] * 10,
        residue_names=["ALA"] * 10
    )
    
    if success:
        print(f"  ✓ Created GRO file: {gro_file}")
    
    # Create sample energy data
    energy_handler = EnergyHandler()
    
    time_points = np.arange(0, 10, 0.1)
    energy_data = EnergyData(
        time_ps=time_points.tolist(),
        potential_energy=(np.random.normal(-50000, 1000, len(time_points))).tolist(),
        kinetic_energy=(np.random.normal(25000, 500, len(time_points))).tolist(),
        total_energy=(np.random.normal(-25000, 800, len(time_points))).tolist(),
        temperature=(np.random.normal(300, 10, len(time_points))).tolist(),
        pressure=(np.random.normal(1.0, 0.1, len(time_points))).tolist(),
        volume=(np.random.normal(125, 5, len(time_points))).tolist(),
        density=(np.random.normal(1000, 20, len(time_points))).tolist()
    )
    
    # Write XVG file
    xvg_file = demo_dir / "energy.xvg"
    success = energy_handler.write(energy_data, xvg_file)
    
    if success:
        print(f"  ✓ Created XVG file: {xvg_file}")
    
    # Create sample PMF data
    pmf_handler = PMFDataHandler()
    
    reaction_coord = np.linspace(0.5, 4.0, 30)
    pmf_values = -8.5 + 2.0 * (reaction_coord - 1.5)**2 + np.random.normal(0, 0.5, len(reaction_coord))
    
    pmf_result = PMFResult(
        reaction_coordinate=reaction_coord.tolist(),
        pmf_values=pmf_values.tolist(),
        error_estimates=np.random.uniform(0.5, 2.0, len(reaction_coord)).tolist(),
        binding_energy=float(np.min(pmf_values)),
        convergence_status=True,
        method="WHAM",
        temperature=300.0,
        metadata={"simulation_time_ns": 100, "windows": 30}
    )
    
    # Write JSON file
    pmf_json_file = demo_dir / "pmf_results.json"
    success = pmf_handler.write(pmf_result, pmf_json_file)
    
    if success:
        print(f"  ✓ Created PMF JSON file: {pmf_json_file}")
    
    print("\n2. Reading back data files:")
    
    # Read trajectory
    if gro_file.exists():
        frames_read = list(trajectory_handler.read(gro_file))
        print(f"  ✓ Read {len(frames_read)} frames from GRO file")
        if frames_read:
            print(f"    First frame: {frames_read[0].coordinates.shape[0]} atoms")
    
    # Read energy data
    if xvg_file.exists():
        energy_read = energy_handler.read(xvg_file)
        print(f"  ✓ Read energy data: {energy_read.n_frames} time points")
        print(f"    Temperature range: {min(energy_read.temperature):.1f} - {max(energy_read.temperature):.1f} K")
    
    # Read PMF data
    if pmf_json_file.exists():
        pmf_read = pmf_handler.read(pmf_json_file)
        print(f"  ✓ Read PMF data: {pmf_read.n_points} points")
        print(f"    Binding energy: {pmf_read.binding_energy:.2f} kJ/mol")
    
    return demo_dir


def demonstrate_backup_system():
    """Demonstrate the backup system"""
    print("\n" + "="*50)
    print("DEMONSTRATING BACKUP SYSTEM")
    print("="*50)
    
    # Create demo data directory
    data_dir = Path("demo_data")
    data_dir.mkdir(exist_ok=True)
    
    # Create sample files
    print("\n1. Creating sample data files:")
    
    files_created = []
    for i in range(5):
        file_path = data_dir / f"simulation_{i:03d}.dat"
        with open(file_path, 'w') as f:
            f.write(f"Simulation {i} data\n")
            f.write("Time: " + " ".join([str(j) for j in range(100)]) + "\n")
            f.write("Energy: " + " ".join([str(np.random.random()) for j in range(100)]) + "\n")
        files_created.append(file_path)
        print(f"  ✓ Created {file_path}")
    
    # Create backup manager
    backup_dir = Path("demo_backups")
    backup_manager = create_backup_manager(backup_dir)
    
    print(f"\n2. Creating backups:")
    
    def backup_progress(phase, count, item):
        if count % 10 == 0 or count <= 5:
            print(f"  {phase}: {count} items processed ({item})")
    
    # Create full backup
    full_backup = backup_manager.create_backup(
        data_dir,
        BackupType.FULL,
        CompressionType.GZIP,
        exclude_patterns=["*.tmp", "*.log"],
        progress_callback=backup_progress
    )
    
    print(f"  ✓ Full backup created: {full_backup.backup_id}")
    print(f"    Files: {full_backup.file_count}")
    print(f"    Size: {full_backup.total_size_bytes} -> {full_backup.compressed_size_bytes} bytes")
    print(f"    Compression: {full_backup.compression_ratio:.2f}x")
    
    # Modify some files (simulate changes)
    time.sleep(1)  # Ensure different timestamp
    for i in range(2):
        file_path = data_dir / f"simulation_{i:03d}.dat"
        with open(file_path, 'a') as f:
            f.write(f"Additional data at {time.time()}\n")
    
    # Create incremental backup
    incremental_backup = backup_manager.create_backup(
        data_dir,
        BackupType.INCREMENTAL,
        CompressionType.GZIP,
        progress_callback=backup_progress
    )
    
    print(f"  ✓ Incremental backup created: {incremental_backup.backup_id}")
    print(f"    Files: {incremental_backup.file_count}")
    
    print(f"\n3. Backup management:")
    
    # List backups
    all_backups = backup_manager.list_backups()
    print(f"  Total backups: {len(all_backups)}")
    for backup in all_backups:
        print(f"    {backup.backup_id}: {backup.backup_type.value}, "
             f"{backup.file_count} files, {backup.compressed_size_bytes} bytes")
    
    # Show backup statistics
    stats = backup_manager.get_backup_stats()
    print(f"\n  Statistics:")
    print(f"    Total backups: {stats['total_backups']}")
    print(f"    Total size: {stats['total_size_mb']:.1f} MB")
    print(f"    Average compression: {stats['average_compression_ratio']:.2f}x")
    
    print(f"\n4. Backup restoration:")
    
    # Create restore directory
    restore_dir = Path("demo_restore")
    restore_dir.mkdir(exist_ok=True)
    
    def restore_progress(phase, count, item):
        if count % 10 == 0 or count <= 5:
            print(f"  {phase}: {count} items processed")
    
    # Restore full backup
    success = backup_manager.restore_backup(
        full_backup.backup_id,
        restore_dir,
        overwrite_existing=True,
        progress_callback=restore_progress
    )
    
    if success:
        print(f"  ✓ Backup restored to {restore_dir}")
        
        # Verify restored files
        restored_files = list(restore_dir.glob("*.dat"))
        print(f"    Restored {len(restored_files)} files")
        
        # Compare content
        original_file = data_dir / "simulation_000.dat"
        restored_file = restore_dir / "simulation_000.dat"
        
        if original_file.exists() and restored_file.exists():
            original_content = original_file.read_text()
            restored_content = restored_file.read_text()
            
            # Compare first part (before modifications)
            original_lines = original_content.split('\n')[:3]
            restored_lines = restored_content.split('\n')[:3]
            
            if original_lines == restored_lines:
                print(f"    ✓ File content verification passed")
            else:
                print(f"    ⚠ File content differs (expected for incremental)")
    
    return backup_manager, data_dir, backup_dir


def main():
    """Run comprehensive data management demonstration"""
    print("PRISM Data Management System Demonstration")
    print("="*50)
    
    # Setup logging
    setup_logging(level="INFO", enable_monitoring=False)
    
    try:
        # Demonstrate each component
        storage = demonstrate_storage_system()
        
        processors = demonstrate_processing_system()
        
        cache_manager = demonstrate_caching_system()
        
        format_demo_dir = demonstrate_format_handlers()
        
        backup_manager, data_dir, backup_dir = demonstrate_backup_system()
        
        print("\n" + "="*50)
        print("DEMONSTRATION COMPLETED SUCCESSFULLY")
        print("="*50)
        
        print("\nGenerated demo files:")
        print(f"  Storage data: demo_storage/")
        print(f"  Format examples: {format_demo_dir}/")
        print(f"  Data files: {data_dir}/")
        print(f"  Backups: {backup_dir}/")
        
        print("\nFeatures demonstrated:")
        print("  ✓ Hierarchical storage with compression")
        print("  ✓ Stream and batch data processing")
        print("  ✓ Multi-level caching system")
        print("  ✓ Format handlers (GRO, XVG, JSON)")
        print("  ✓ Backup and recovery system")
        
        print("\nCleanup:")
        cleanup = input("Delete demo files? (y/N): ").lower().strip()
        if cleanup == 'y':
            import shutil
            
            for path in [Path("demo_storage"), Path("demo_formats"), 
                        Path("demo_data"), Path("demo_backups"), Path("demo_restore")]:
                if path.exists():
                    shutil.rmtree(path)
                    print(f"  ✓ Removed {path}")
        
    except Exception as e:
        print(f"\nDemonstration failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()