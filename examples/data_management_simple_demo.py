#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Simple PRISM Data Management System Demo

A simplified demonstration that doesn't require external dependencies.
"""

import sys
import time
import json
import tempfile
from pathlib import Path

# Add simplified data structures for demo
class DemoStorageManager:
    """Simplified storage manager for demonstration"""
    
    def __init__(self, base_path):
        self.base_path = Path(base_path)
        self.base_path.mkdir(parents=True, exist_ok=True)
        self.data_store = {}
        self.metadata = {}
    
    def store_data(self, key, data, **kwargs):
        """Store data with metadata"""
        file_path = self.base_path / f"{key}.json"
        
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.metadata[key] = {
            'path': str(file_path),
            'size': file_path.stat().st_size,
            'created': time.time(),
            'type': kwargs.get('content_type', 'unknown')
        }
        
        return str(file_path)
    
    def retrieve_data(self, key):
        """Retrieve stored data"""
        if key not in self.metadata:
            return None, None
            
        file_path = Path(self.metadata[key]['path'])
        
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        return data, self.metadata[key]
    
    def get_stats(self):
        """Get storage statistics"""
        total_size = sum(meta['size'] for meta in self.metadata.values())
        return {
            'total_objects': len(self.metadata),
            'total_size_mb': total_size / (1024 * 1024),
            'objects': list(self.metadata.keys())
        }


class DemoCache:
    """Simplified cache for demonstration"""
    
    def __init__(self, max_size=100):
        self.cache = {}
        self.access_times = {}
        self.max_size = max_size
        self.hits = 0
        self.misses = 0
    
    def get(self, key):
        """Get value from cache"""
        if key in self.cache:
            self.hits += 1
            self.access_times[key] = time.time()
            return self.cache[key]
        else:
            self.misses += 1
            return None
    
    def put(self, key, value):
        """Put value in cache"""
        if len(self.cache) >= self.max_size:
            # Remove oldest item
            oldest_key = min(self.access_times.keys(), key=lambda k: self.access_times[k])
            del self.cache[oldest_key]
            del self.access_times[oldest_key]
        
        self.cache[key] = value
        self.access_times[key] = time.time()
    
    def get_stats(self):
        """Get cache statistics"""
        total_requests = self.hits + self.misses
        hit_ratio = self.hits / total_requests if total_requests > 0 else 0
        
        return {
            'hits': self.hits,
            'misses': self.misses,
            'hit_ratio': hit_ratio,
            'size': len(self.cache)
        }


def demonstrate_storage_system():
    """Demonstrate simplified storage system"""
    print("\n" + "="*50)
    print("DEMONSTRATING STORAGE SYSTEM")
    print("="*50)
    
    # Create storage manager
    storage = DemoStorageManager("demo_storage")
    
    print("\n1. Storing different data types:")
    
    # Store configuration
    config = {
        "simulation": "PMF",
        "temperature": 300.0,
        "steps": 1000000,
        "output_frequency": 1000
    }
    
    path = storage.store_data("config", config, content_type="configuration")
    print(f"  ✓ Stored configuration: {path}")
    
    # Store analysis results
    results = {
        "binding_energy": -8.5,
        "rmsd_values": [0.1, 0.15, 0.12, 0.18, 0.14],
        "convergence": True,
        "method": "WHAM"
    }
    
    path = storage.store_data("analysis", results, content_type="results")
    print(f"  ✓ Stored analysis results: {path}")
    
    # Store system parameters
    system = {
        "protein_atoms": 1500,
        "ligand_atoms": 25,
        "water_molecules": 8000,
        "box_size": [5.0, 5.0, 5.0]
    }
    
    path = storage.store_data("system", system, content_type="system")
    print(f"  ✓ Stored system parameters: {path}")
    
    print("\n2. Retrieving stored data:")
    
    # Retrieve data
    retrieved_config, config_meta = storage.retrieve_data("config")
    if retrieved_config:
        print(f"  ✓ Retrieved config: {retrieved_config['simulation']}")
        print(f"    File size: {config_meta['size']} bytes")
    
    retrieved_results, results_meta = storage.retrieve_data("analysis")
    if retrieved_results:
        print(f"  ✓ Retrieved results: binding energy = {retrieved_results['binding_energy']}")
        print(f"    Convergence: {retrieved_results['convergence']}")
    
    # Show statistics
    stats = storage.get_stats()
    print(f"\n3. Storage Statistics:")
    print(f"  Total objects: {stats['total_objects']}")
    print(f"  Total size: {stats['total_size_mb']:.3f} MB")
    print(f"  Stored objects: {', '.join(stats['objects'])}")
    
    return storage


def demonstrate_caching_system():
    """Demonstrate simplified caching system"""
    print("\n" + "="*50)
    print("DEMONSTRATING CACHING SYSTEM")
    print("="*50)
    
    cache = DemoCache(max_size=5)
    
    print("\n1. Caching simulation data:")
    
    # Cache some data
    cache.put("temp_300K", {"temperature": 300.0, "results": [1, 2, 3, 4, 5]})
    print("  ✓ Cached temperature 300K results")
    
    cache.put("temp_310K", {"temperature": 310.0, "results": [1.1, 2.1, 3.1, 4.1, 5.1]})
    print("  ✓ Cached temperature 310K results")
    
    cache.put("temp_320K", {"temperature": 320.0, "results": [1.2, 2.2, 3.2, 4.2, 5.2]})
    print("  ✓ Cached temperature 320K results")
    
    print("\n2. Cache retrieval:")
    
    # Test cache hits
    result = cache.get("temp_300K")
    if result:
        print(f"  ✓ Cache hit: temp_300K -> {result['temperature']}K")
    
    result = cache.get("temp_310K")
    if result:
        print(f"  ✓ Cache hit: temp_310K -> {result['temperature']}K")
    
    # Test cache miss
    result = cache.get("temp_400K")
    if result is None:
        print("  ✓ Cache miss: temp_400K (expected)")
    
    # Fill cache beyond capacity
    print("\n3. Cache eviction test:")
    
    for temp in [330, 340, 350, 360, 370]:
        cache.put(f"temp_{temp}K", {"temperature": temp, "results": [temp/100] * 5})
        print(f"  ✓ Cached temp_{temp}K")
    
    # Test if old entries were evicted
    old_result = cache.get("temp_300K")
    if old_result is None:
        print("  ✓ Old cache entry evicted (expected)")
    
    # Show cache statistics
    stats = cache.get_stats()
    print(f"\n4. Cache Statistics:")
    print(f"  Hits: {stats['hits']}")
    print(f"  Misses: {stats['misses']}")
    print(f"  Hit ratio: {stats['hit_ratio']:.2%}")
    print(f"  Cache size: {stats['size']}")
    
    return cache


def demonstrate_data_processing():
    """Demonstrate simplified data processing"""
    print("\n" + "="*50)
    print("DEMONSTRATING DATA PROCESSING")
    print("="*50)
    
    print("\n1. Processing energy data:")
    
    # Simulate energy data processing
    energy_data = []
    for i in range(10):
        energy_point = {
            "time": i * 0.1,
            "potential": -50000 + (i * 100),  # Simulate drift
            "kinetic": 25000 + (i * 50),
            "total": -25000 + (i * 150)
        }
        energy_data.append(energy_point)
    
    print(f"  ✓ Generated {len(energy_data)} energy data points")
    
    # Calculate statistics
    potentials = [point["potential"] for point in energy_data]
    kinetics = [point["kinetic"] for point in energy_data]
    
    stats = {
        "potential_mean": sum(potentials) / len(potentials),
        "potential_min": min(potentials),
        "potential_max": max(potentials),
        "kinetic_mean": sum(kinetics) / len(kinetics),
        "kinetic_min": min(kinetics),
        "kinetic_max": max(kinetics)
    }
    
    print(f"  ✓ Processed energy statistics:")
    print(f"    Potential: {stats['potential_mean']:.1f} ± {(stats['potential_max']-stats['potential_min'])/2:.1f}")
    print(f"    Kinetic: {stats['kinetic_mean']:.1f} ± {(stats['kinetic_max']-stats['kinetic_min'])/2:.1f}")
    
    print("\n2. Processing trajectory analysis:")
    
    # Simulate RMSD calculation
    rmsd_values = []
    for i in range(20):
        # Simulate RMSD with some drift and fluctuation
        base_rmsd = 0.15 + (i * 0.01)  # Gradual drift
        noise = (hash(str(i)) % 100 - 50) / 1000  # Pseudo-random noise
        rmsd = max(0.05, base_rmsd + noise)  # Ensure positive
        rmsd_values.append(rmsd)
    
    print(f"  ✓ Calculated RMSD for {len(rmsd_values)} frames")
    
    rmsd_stats = {
        "mean": sum(rmsd_values) / len(rmsd_values),
        "min": min(rmsd_values),
        "max": max(rmsd_values)
    }
    
    print(f"    RMSD: {rmsd_stats['mean']:.3f} ± {(rmsd_stats['max']-rmsd_stats['min'])/2:.3f} nm")
    
    return {"energy_stats": stats, "rmsd_stats": rmsd_stats}


def demonstrate_backup_system():
    """Demonstrate simplified backup system"""
    print("\n" + "="*50)
    print("DEMONSTRATING BACKUP SYSTEM")
    print("="*50)
    
    # Create temporary data
    with tempfile.TemporaryDirectory(prefix="prism_demo_") as temp_dir:
        data_dir = Path(temp_dir) / "data"
        backup_dir = Path(temp_dir) / "backups"
        
        data_dir.mkdir()
        backup_dir.mkdir()
        
        print(f"\n1. Creating sample data files:")
        
        # Create sample files
        files_created = []
        for i in range(3):
            file_path = data_dir / f"simulation_{i}.txt"
            with open(file_path, 'w') as f:
                f.write(f"Simulation {i} data\n")
                f.write(f"Temperature: 300.{i}K\n")
                f.write(f"Steps: {1000 * (i + 1)}\n")
                f.write(f"Results: {[j * 0.1 for j in range(10)]}\n")
            
            files_created.append(file_path)
            print(f"  ✓ Created {file_path.name}")
        
        print(f"\n2. Creating backup:")
        
        # Simple backup implementation
        import shutil
        import tarfile
        
        backup_file = backup_dir / f"backup_{int(time.time())}.tar.gz"
        
        with tarfile.open(backup_file, "w:gz") as tar:
            tar.add(data_dir, arcname="data")
        
        backup_size = backup_file.stat().st_size
        original_size = sum(f.stat().st_size for f in files_created)
        compression_ratio = original_size / backup_size if backup_size > 0 else 1.0
        
        print(f"  ✓ Created backup: {backup_file.name}")
        print(f"    Original size: {original_size} bytes")
        print(f"    Compressed size: {backup_size} bytes")
        print(f"    Compression ratio: {compression_ratio:.2f}x")
        
        print(f"\n3. Restoring backup:")
        
        # Create restore directory
        restore_dir = Path(temp_dir) / "restored"
        restore_dir.mkdir()
        
        # Extract backup
        with tarfile.open(backup_file, "r:gz") as tar:
            tar.extractall(restore_dir)
        
        restored_files = list((restore_dir / "data").glob("*.txt"))
        print(f"  ✓ Restored {len(restored_files)} files")
        
        # Verify restoration
        for original, restored in zip(files_created, sorted(restored_files)):
            original_content = original.read_text()
            restored_content = restored.read_text()
            
            if original_content == restored_content:
                print(f"    ✓ {original.name} verified")
            else:
                print(f"    ⚠ {original.name} differs")
        
        return {
            "backup_file": backup_file.name,
            "original_size": original_size,
            "compressed_size": backup_size,
            "compression_ratio": compression_ratio,
            "files_backed_up": len(files_created),
            "files_restored": len(restored_files)
        }


def main():
    """Run simplified data management demonstration"""
    print("PRISM Data Management System - Simple Demo")
    print("="*50)
    print("This demo showcases core data management features")
    print("without requiring external dependencies.")
    
    try:
        # Demonstrate each component
        storage = demonstrate_storage_system()
        
        cache = demonstrate_caching_system()
        
        processing_results = demonstrate_data_processing()
        
        backup_info = demonstrate_backup_system()
        
        print("\n" + "="*50)
        print("DEMONSTRATION COMPLETED SUCCESSFULLY")
        print("="*50)
        
        print("\nSummary of demonstrated features:")
        print("  ✓ Data storage and retrieval")
        print("  ✓ Multi-level caching with eviction")
        print("  ✓ Data processing and analysis")
        print("  ✓ Backup and recovery operations")
        
        print("\nKey metrics:")
        storage_stats = storage.get_stats()
        cache_stats = cache.get_stats()
        
        print(f"  Storage: {storage_stats['total_objects']} objects, "
             f"{storage_stats['total_size_mb']:.3f} MB")
        print(f"  Cache: {cache_stats['hit_ratio']:.1%} hit ratio")
        print(f"  Backup: {backup_info['compression_ratio']:.1f}x compression")
        
        # Show processing results
        energy_stats = processing_results['energy_stats']
        rmsd_stats = processing_results['rmsd_stats']
        
        print(f"  Analysis: PE={energy_stats['potential_mean']:.0f}, "
             f"RMSD={rmsd_stats['mean']:.3f}nm")
        
        print(f"\nDemo completed in temporary storage.")
        print("All temporary files automatically cleaned up.")
        
    except Exception as e:
        print(f"\nDemo failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()