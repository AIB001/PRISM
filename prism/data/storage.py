#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Storage Management System

Advanced data storage management with hierarchical storage, compression,
and efficient data organization for molecular dynamics simulations.
"""

import os
import json
import pickle
import gzip
import bz2
import lzma
import hashlib
import shutil
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Tuple, Iterator
from dataclasses import dataclass, asdict
from enum import Enum
import threading
import time
from contextlib import contextmanager

from ..utils.logging_system import PrismLogger


class StorageTier(Enum):
    """Storage tier levels for hierarchical storage"""
    HOT = "hot"           # Frequently accessed, fast storage (SSD)
    WARM = "warm"         # Occasionally accessed, standard storage
    COLD = "cold"         # Rarely accessed, slow/archival storage
    FROZEN = "frozen"     # Long-term archival, compressed storage


class CompressionType(Enum):
    """Supported compression algorithms"""
    NONE = "none"
    GZIP = "gzip"
    BZIP2 = "bzip2"
    LZMA = "lzma"
    ZSTD = "zstd"


@dataclass
class StorageMetadata:
    """Metadata for stored data objects"""
    object_id: str
    file_path: str
    size_bytes: int
    compression: CompressionType
    storage_tier: StorageTier
    created_timestamp: float
    accessed_timestamp: float
    access_count: int
    checksum: str
    content_type: str
    custom_metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.custom_metadata is None:
            self.custom_metadata = {}
    
    def to_dict(self) -> Dict[str, Any]:
        data = asdict(self)
        data['compression'] = self.compression.value
        data['storage_tier'] = self.storage_tier.value
        return data
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'StorageMetadata':
        data['compression'] = CompressionType(data['compression'])
        data['storage_tier'] = StorageTier(data['storage_tier'])
        return cls(**data)


class CompressionManager:
    """Manages data compression and decompression"""
    
    def __init__(self):
        self.logger = PrismLogger("compression_manager")
        
        # Compression handlers
        self.compressors = {
            CompressionType.GZIP: self._gzip_compress,
            CompressionType.BZIP2: self._bzip2_compress,
            CompressionType.LZMA: self._lzma_compress
        }
        
        self.decompressors = {
            CompressionType.GZIP: self._gzip_decompress,
            CompressionType.BZIP2: self._bzip2_decompress,
            CompressionType.LZMA: self._lzma_decompress
        }
    
    def compress_data(self, data: bytes, compression_type: CompressionType) -> bytes:
        """Compress data using specified algorithm"""
        if compression_type == CompressionType.NONE:
            return data
        
        if compression_type not in self.compressors:
            raise ValueError(f"Unsupported compression type: {compression_type}")
        
        start_time = time.time()
        compressed_data = self.compressors[compression_type](data)
        compression_time = time.time() - start_time
        
        compression_ratio = len(data) / len(compressed_data) if len(compressed_data) > 0 else 1
        
        self.logger.info(f"Compressed {len(data)} bytes to {len(compressed_data)} bytes "
                        f"({compression_ratio:.2f}x) in {compression_time:.3f}s using {compression_type.value}")
        
        return compressed_data
    
    def decompress_data(self, data: bytes, compression_type: CompressionType) -> bytes:
        """Decompress data using specified algorithm"""
        if compression_type == CompressionType.NONE:
            return data
        
        if compression_type not in self.decompressors:
            raise ValueError(f"Unsupported compression type: {compression_type}")
        
        return self.decompressors[compression_type](data)
    
    def _gzip_compress(self, data: bytes) -> bytes:
        return gzip.compress(data, compresslevel=6)
    
    def _gzip_decompress(self, data: bytes) -> bytes:
        return gzip.decompress(data)
    
    def _bzip2_compress(self, data: bytes) -> bytes:
        return bz2.compress(data, compresslevel=6)
    
    def _bzip2_decompress(self, data: bytes) -> bytes:
        return bz2.decompress(data)
    
    def _lzma_compress(self, data: bytes) -> bytes:
        return lzma.compress(data, preset=6)
    
    def _lzma_decompress(self, data: bytes) -> bytes:
        return lzma.decompress(data)
    
    def choose_best_compression(self, data: bytes, 
                              target_compression_ratio: float = 2.0) -> CompressionType:
        """Choose best compression algorithm based on data characteristics"""
        data_size = len(data)
        
        # For small data, compression overhead may not be worth it
        if data_size < 1024:  # 1KB
            return CompressionType.NONE
        
        # Test different compression algorithms on a sample
        sample_size = min(10000, data_size)  # 10KB sample or full data
        sample_data = data[:sample_size]
        
        results = {}
        
        for comp_type in [CompressionType.GZIP, CompressionType.BZIP2, CompressionType.LZMA]:
            try:
                start_time = time.time()
                compressed = self.compress_data(sample_data, comp_type)
                compression_time = time.time() - start_time
                
                compression_ratio = len(sample_data) / len(compressed)
                speed = len(sample_data) / compression_time / (1024 * 1024)  # MB/s
                
                # Score based on compression ratio and speed
                score = compression_ratio * 0.7 + speed * 0.3 / 100  # Normalize speed
                
                results[comp_type] = {
                    'ratio': compression_ratio,
                    'speed': speed,
                    'score': score
                }
                
            except Exception as e:
                self.logger.warning(f"Failed to test {comp_type}: {e}")
                results[comp_type] = {'score': 0}
        
        # Choose best algorithm
        best_algorithm = max(results.keys(), key=lambda k: results[k]['score'])
        
        # If compression ratio is too low, don't compress
        if results[best_algorithm]['ratio'] < target_compression_ratio:
            return CompressionType.NONE
        
        return best_algorithm


class HierarchicalStorage:
    """Hierarchical storage management with automatic tier migration"""
    
    def __init__(self, storage_config: Dict[str, Path]):
        self.storage_config = storage_config
        self.logger = PrismLogger("hierarchical_storage")
        
        # Ensure storage directories exist
        for tier, path in storage_config.items():
            path.mkdir(parents=True, exist_ok=True)
        
        # Access tracking
        self.access_tracker = {}
        self.lock = threading.RLock()
        
        # Migration policies
        self.migration_policies = {
            StorageTier.HOT: {
                'max_age_days': 7,
                'min_access_count': 10,
                'max_size_gb': 100
            },
            StorageTier.WARM: {
                'max_age_days': 30,
                'min_access_count': 3,
                'max_size_gb': 500
            },
            StorageTier.COLD: {
                'max_age_days': 180,
                'min_access_count': 1,
                'max_size_gb': 1000
            }
        }
    
    def get_storage_path(self, tier: StorageTier) -> Path:
        """Get storage path for specified tier"""
        tier_key = tier.value
        if tier_key not in self.storage_config:
            raise ValueError(f"Storage tier {tier} not configured")
        return self.storage_config[tier_key]
    
    def should_migrate(self, metadata: StorageMetadata) -> Optional[StorageTier]:
        """Determine if data should be migrated to different tier"""
        current_tier = metadata.storage_tier
        current_age_days = (time.time() - metadata.created_timestamp) / 86400
        
        # Check migration from hot to warm
        if current_tier == StorageTier.HOT:
            policy = self.migration_policies[StorageTier.HOT]
            if (current_age_days > policy['max_age_days'] or 
                metadata.access_count < policy['min_access_count']):
                return StorageTier.WARM
        
        # Check migration from warm to cold
        elif current_tier == StorageTier.WARM:
            policy = self.migration_policies[StorageTier.WARM]
            if (current_age_days > policy['max_age_days'] or 
                metadata.access_count < policy['min_access_count']):
                return StorageTier.COLD
        
        # Check migration from cold to frozen
        elif current_tier == StorageTier.COLD:
            policy = self.migration_policies[StorageTier.COLD]
            if (current_age_days > policy['max_age_days'] or 
                metadata.access_count < policy['min_access_count']):
                return StorageTier.FROZEN
        
        return None
    
    def migrate_data(self, object_id: str, from_metadata: StorageMetadata, 
                    to_tier: StorageTier) -> StorageMetadata:
        """Migrate data to different storage tier"""
        from_path = Path(from_metadata.file_path)
        to_path = self.get_storage_path(to_tier) / from_path.name
        
        # Move file
        shutil.move(str(from_path), str(to_path))
        
        # Update metadata
        from_metadata.file_path = str(to_path)
        from_metadata.storage_tier = to_tier
        
        self.logger.info(f"Migrated {object_id} from {from_metadata.storage_tier.value} to {to_tier.value}")
        
        return from_metadata
    
    def track_access(self, object_id: str):
        """Track data access for migration decisions"""
        with self.lock:
            current_time = time.time()
            if object_id not in self.access_tracker:
                self.access_tracker[object_id] = {
                    'count': 0,
                    'last_access': current_time
                }
            
            self.access_tracker[object_id]['count'] += 1
            self.access_tracker[object_id]['last_access'] = current_time


class DataStorageManager:
    """Main data storage management system"""
    
    def __init__(self, base_storage_path: Path, 
                 enable_hierarchical: bool = True,
                 enable_compression: bool = True):
        self.base_path = base_storage_path
        self.base_path.mkdir(parents=True, exist_ok=True)
        
        self.logger = PrismLogger("data_storage_manager")
        self.compression_manager = CompressionManager() if enable_compression else None
        
        # Metadata storage
        self.metadata_file = self.base_path / ".metadata.json"
        self.metadata: Dict[str, StorageMetadata] = self._load_metadata()
        
        # Hierarchical storage setup
        if enable_hierarchical:
            storage_config = {
                StorageTier.HOT.value: self.base_path / "hot",
                StorageTier.WARM.value: self.base_path / "warm", 
                StorageTier.COLD.value: self.base_path / "cold",
                StorageTier.FROZEN.value: self.base_path / "frozen"
            }
            self.hierarchical_storage = HierarchicalStorage(storage_config)
        else:
            self.hierarchical_storage = None
        
        # Thread safety
        self.lock = threading.RLock()
        
        # Auto-cleanup thread
        self._start_maintenance_thread()
    
    def store_data(self, object_id: str, data: Any, 
                  content_type: str = "binary",
                  storage_tier: StorageTier = StorageTier.HOT,
                  compression: Optional[CompressionType] = None,
                  custom_metadata: Optional[Dict[str, Any]] = None) -> str:
        """Store data object with metadata"""
        
        # Serialize data
        if isinstance(data, (dict, list)):
            serialized_data = json.dumps(data).encode('utf-8')
            content_type = "json"
        elif isinstance(data, str):
            serialized_data = data.encode('utf-8')
            content_type = "text"
        elif isinstance(data, bytes):
            serialized_data = data
        else:
            # Use pickle for other objects
            serialized_data = pickle.dumps(data)
            content_type = "pickle"
        
        # Choose compression if not specified
        if compression is None and self.compression_manager:
            compression = self.compression_manager.choose_best_compression(serialized_data)
        elif compression is None:
            compression = CompressionType.NONE
        
        # Compress data
        if self.compression_manager and compression != CompressionType.NONE:
            compressed_data = self.compression_manager.compress_data(serialized_data, compression)
        else:
            compressed_data = serialized_data
            compression = CompressionType.NONE
        
        # Generate file path
        if self.hierarchical_storage:
            storage_path = self.hierarchical_storage.get_storage_path(storage_tier)
        else:
            storage_path = self.base_path
        
        file_name = f"{object_id}.dat"
        file_path = storage_path / file_name
        
        # Calculate checksum
        checksum = hashlib.sha256(compressed_data).hexdigest()
        
        # Write data to file
        with open(file_path, 'wb') as f:
            f.write(compressed_data)
        
        # Create metadata
        current_time = time.time()
        metadata = StorageMetadata(
            object_id=object_id,
            file_path=str(file_path),
            size_bytes=len(compressed_data),
            compression=compression,
            storage_tier=storage_tier,
            created_timestamp=current_time,
            accessed_timestamp=current_time,
            access_count=1,
            checksum=checksum,
            content_type=content_type,
            custom_metadata=custom_metadata or {}
        )
        
        # Store metadata
        with self.lock:
            self.metadata[object_id] = metadata
            self._save_metadata()
        
        self.logger.info(f"Stored {object_id}: {len(compressed_data)} bytes "
                        f"({compression.value} compression, {storage_tier.value} tier)")
        
        return str(file_path)
    
    def retrieve_data(self, object_id: str) -> Tuple[Any, StorageMetadata]:
        """Retrieve data object with metadata"""
        
        with self.lock:
            if object_id not in self.metadata:
                raise KeyError(f"Object {object_id} not found in storage")
            
            metadata = self.metadata[object_id]
        
        # Track access
        if self.hierarchical_storage:
            self.hierarchical_storage.track_access(object_id)
        
        # Update access metadata
        current_time = time.time()
        metadata.accessed_timestamp = current_time
        metadata.access_count += 1
        
        # Read data from file
        file_path = Path(metadata.file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"Data file not found: {file_path}")
        
        with open(file_path, 'rb') as f:
            compressed_data = f.read()
        
        # Verify checksum
        actual_checksum = hashlib.sha256(compressed_data).hexdigest()
        if actual_checksum != metadata.checksum:
            raise ValueError(f"Data corruption detected for {object_id}")
        
        # Decompress data
        if self.compression_manager and metadata.compression != CompressionType.NONE:
            serialized_data = self.compression_manager.decompress_data(compressed_data, metadata.compression)
        else:
            serialized_data = compressed_data
        
        # Deserialize data
        if metadata.content_type == "json":
            data = json.loads(serialized_data.decode('utf-8'))
        elif metadata.content_type == "text":
            data = serialized_data.decode('utf-8')
        elif metadata.content_type == "pickle":
            data = pickle.loads(serialized_data)
        else:
            data = serialized_data
        
        # Update metadata
        with self.lock:
            self.metadata[object_id] = metadata
            self._save_metadata()
        
        return data, metadata
    
    def delete_data(self, object_id: str) -> bool:
        """Delete data object and metadata"""
        
        with self.lock:
            if object_id not in self.metadata:
                return False
            
            metadata = self.metadata[object_id]
            file_path = Path(metadata.file_path)
            
            # Delete file
            if file_path.exists():
                file_path.unlink()
            
            # Remove metadata
            del self.metadata[object_id]
            self._save_metadata()
        
        self.logger.info(f"Deleted {object_id}")
        return True
    
    def list_objects(self, storage_tier: Optional[StorageTier] = None,
                    content_type: Optional[str] = None) -> List[str]:
        """List stored objects with optional filtering"""
        
        with self.lock:
            objects = []
            for object_id, metadata in self.metadata.items():
                if storage_tier and metadata.storage_tier != storage_tier:
                    continue
                if content_type and metadata.content_type != content_type:
                    continue
                objects.append(object_id)
        
        return objects
    
    def get_storage_stats(self) -> Dict[str, Any]:
        """Get storage usage statistics"""
        
        with self.lock:
            total_objects = len(self.metadata)
            total_size = sum(meta.size_bytes for meta in self.metadata.values())
            
            # Statistics by tier
            tier_stats = {}
            for tier in StorageTier:
                tier_objects = [m for m in self.metadata.values() if m.storage_tier == tier]
                tier_stats[tier.value] = {
                    'objects': len(tier_objects),
                    'size_bytes': sum(m.size_bytes for m in tier_objects)
                }
            
            # Statistics by compression
            compression_stats = {}
            for comp_type in CompressionType:
                comp_objects = [m for m in self.metadata.values() if m.compression == comp_type]
                compression_stats[comp_type.value] = {
                    'objects': len(comp_objects),
                    'size_bytes': sum(m.size_bytes for m in comp_objects)
                }
        
        return {
            'total_objects': total_objects,
            'total_size_bytes': total_size,
            'total_size_mb': total_size / (1024 * 1024),
            'tier_statistics': tier_stats,
            'compression_statistics': compression_stats
        }
    
    def run_maintenance(self):
        """Run storage maintenance tasks"""
        if not self.hierarchical_storage:
            return
        
        with self.lock:
            migrations = 0
            for object_id, metadata in list(self.metadata.items()):
                target_tier = self.hierarchical_storage.should_migrate(metadata)
                if target_tier:
                    try:
                        updated_metadata = self.hierarchical_storage.migrate_data(
                            object_id, metadata, target_tier
                        )
                        self.metadata[object_id] = updated_metadata
                        migrations += 1
                    except Exception as e:
                        self.logger.error(f"Failed to migrate {object_id}: {e}")
            
            if migrations > 0:
                self._save_metadata()
                self.logger.info(f"Completed maintenance: {migrations} objects migrated")
    
    def _load_metadata(self) -> Dict[str, StorageMetadata]:
        """Load metadata from disk"""
        if not self.metadata_file.exists():
            return {}
        
        try:
            with open(self.metadata_file, 'r') as f:
                metadata_dict = json.load(f)
            
            metadata = {}
            for object_id, meta_data in metadata_dict.items():
                metadata[object_id] = StorageMetadata.from_dict(meta_data)
            
            return metadata
        except Exception as e:
            self.logger.error(f"Failed to load metadata: {e}")
            return {}
    
    def _save_metadata(self):
        """Save metadata to disk"""
        try:
            metadata_dict = {}
            for object_id, metadata in self.metadata.items():
                metadata_dict[object_id] = metadata.to_dict()
            
            with open(self.metadata_file, 'w') as f:
                json.dump(metadata_dict, f, indent=2)
        except Exception as e:
            self.logger.error(f"Failed to save metadata: {e}")
    
    def _start_maintenance_thread(self):
        """Start background maintenance thread"""
        def maintenance_worker():
            while True:
                try:
                    time.sleep(3600)  # Run every hour
                    self.run_maintenance()
                except Exception as e:
                    self.logger.error(f"Maintenance thread error: {e}")
        
        maintenance_thread = threading.Thread(target=maintenance_worker, daemon=True)
        maintenance_thread.start()
    
    @contextmanager
    def transaction(self):
        """Context manager for transactional operations"""
        with self.lock:
            # Create backup of metadata
            backup_metadata = dict(self.metadata)
            try:
                yield self
            except Exception:
                # Restore metadata on error
                self.metadata = backup_metadata
                self._save_metadata()
                raise
            else:
                # Save metadata on success
                self._save_metadata()


# Convenience functions
def create_storage_manager(base_path: Union[str, Path], 
                         hierarchical: bool = True,
                         compression: bool = True) -> DataStorageManager:
    """Create data storage manager with default settings"""
    return DataStorageManager(
        base_storage_path=Path(base_path),
        enable_hierarchical=hierarchical,
        enable_compression=compression
    )


def store_trajectory_data(storage: DataStorageManager, trajectory_id: str, 
                         trajectory_data: bytes) -> str:
    """Store trajectory data with optimized settings"""
    return storage.store_data(
        object_id=trajectory_id,
        data=trajectory_data,
        content_type="trajectory",
        storage_tier=StorageTier.WARM,
        compression=CompressionType.LZMA  # Good compression for trajectory data
    )


def store_analysis_results(storage: DataStorageManager, analysis_id: str,
                          results: Dict[str, Any]) -> str:
    """Store analysis results with metadata"""
    return storage.store_data(
        object_id=analysis_id,
        data=results,
        content_type="analysis",
        storage_tier=StorageTier.HOT,
        compression=CompressionType.GZIP
    )