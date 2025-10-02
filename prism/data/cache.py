#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Caching System

Advanced caching system for PRISM molecular dynamics data with multi-level
caching, intelligent eviction policies, and performance optimization.
"""

import os
import sys
import time
import hashlib
import pickle
import threading
import weakref
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Callable, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
from abc import ABC, abstractmethod
import tempfile
import shutil

from ..utils.logging_system import PrismLogger


class CachePolicy(Enum):
    """Cache eviction policies"""
    LRU = "lru"         # Least Recently Used
    LFU = "lfu"         # Least Frequently Used
    FIFO = "fifo"       # First In First Out
    TTL = "ttl"         # Time To Live
    SIZE_BASED = "size_based"


class CacheLevel(Enum):
    """Cache hierarchy levels"""
    L1_MEMORY = "l1_memory"     # In-memory cache (fastest)
    L2_DISK = "l2_disk"         # Local disk cache
    L3_NETWORK = "l3_network"   # Network/distributed cache


@dataclass
class CacheEntry:
    """Cache entry metadata"""
    key: str
    value: Any
    size_bytes: int
    created_time: float
    accessed_time: float
    access_count: int
    ttl: Optional[float] = None
    priority: int = 0
    
    @property
    def is_expired(self) -> bool:
        """Check if cache entry is expired"""
        if self.ttl is None:
            return False
        return time.time() - self.created_time > self.ttl
    
    @property
    def age(self) -> float:
        """Get age of cache entry in seconds"""
        return time.time() - self.created_time
    
    def touch(self):
        """Update access time and count"""
        self.accessed_time = time.time()
        self.access_count += 1


class CacheStats:
    """Cache performance statistics"""
    
    def __init__(self):
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self.size_bytes = 0
        self.entry_count = 0
        self.lock = threading.RLock()
    
    @property
    def hit_ratio(self) -> float:
        """Calculate cache hit ratio"""
        total = self.hits + self.misses
        return self.hits / total if total > 0 else 0.0
    
    def record_hit(self):
        """Record cache hit"""
        with self.lock:
            self.hits += 1
    
    def record_miss(self):
        """Record cache miss"""
        with self.lock:
            self.misses += 1
    
    def record_eviction(self):
        """Record cache eviction"""
        with self.lock:
            self.evictions += 1
    
    def update_size(self, size_change: int, entry_change: int):
        """Update size statistics"""
        with self.lock:
            self.size_bytes += size_change
            self.entry_count += entry_change
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        with self.lock:
            return {
                'hits': self.hits,
                'misses': self.misses,
                'evictions': self.evictions,
                'hit_ratio': self.hit_ratio,
                'size_bytes': self.size_bytes,
                'size_mb': self.size_bytes / (1024 * 1024),
                'entry_count': self.entry_count
            }


class BaseCache(ABC):
    """Abstract base class for caches"""
    
    def __init__(self, name: str, max_size_bytes: int = 100 * 1024 * 1024):
        self.name = name
        self.max_size_bytes = max_size_bytes
        self.logger = PrismLogger(f"cache.{name}")
        self.stats = CacheStats()
        self.lock = threading.RLock()
    
    @abstractmethod
    def get(self, key: str) -> Optional[Any]:
        """Get value from cache"""
        pass
    
    @abstractmethod
    def put(self, key: str, value: Any, ttl: Optional[float] = None) -> bool:
        """Put value in cache"""
        pass
    
    @abstractmethod
    def delete(self, key: str) -> bool:
        """Delete value from cache"""
        pass
    
    @abstractmethod
    def clear(self):
        """Clear all cache entries"""
        pass
    
    @abstractmethod
    def keys(self) -> List[str]:
        """Get all cache keys"""
        pass
    
    def contains(self, key: str) -> bool:
        """Check if key exists in cache"""
        return self.get(key) is not None
    
    def size_bytes(self) -> int:
        """Get total cache size in bytes"""
        return self.stats.size_bytes
    
    def size_mb(self) -> float:
        """Get total cache size in MB"""
        return self.size_bytes() / (1024 * 1024)
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        stats = self.stats.get_stats()
        stats['cache_name'] = self.name
        stats['max_size_mb'] = self.max_size_bytes / (1024 * 1024)
        return stats


class MemoryCache(BaseCache):
    """In-memory cache with configurable eviction policies"""
    
    def __init__(self, name: str = "memory_cache", 
                 max_size_bytes: int = 100 * 1024 * 1024,
                 policy: CachePolicy = CachePolicy.LRU):
        super().__init__(name, max_size_bytes)
        self.policy = policy
        self.entries: Dict[str, CacheEntry] = {}
        
        # Policy-specific data structures
        if policy == CachePolicy.LRU:
            self.access_order = []  # Most recent last
        elif policy == CachePolicy.FIFO:
            self.insertion_order = []
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from memory cache"""
        with self.lock:
            if key not in self.entries:
                self.stats.record_miss()
                return None
            
            entry = self.entries[key]
            
            # Check if expired
            if entry.is_expired:
                del self.entries[key]
                self.stats.record_miss()
                self.stats.update_size(-entry.size_bytes, -1)
                return None
            
            # Update access statistics
            entry.touch()
            self.stats.record_hit()
            
            # Update policy-specific structures
            if self.policy == CachePolicy.LRU:
                if key in self.access_order:
                    self.access_order.remove(key)
                self.access_order.append(key)
            
            return entry.value
    
    def put(self, key: str, value: Any, ttl: Optional[float] = None) -> bool:
        """Put value in memory cache"""
        with self.lock:
            # Calculate size of value
            try:
                size_bytes = sys.getsizeof(value)
                if hasattr(value, '__len__'):
                    # For collections, estimate total size
                    size_bytes = max(size_bytes, len(value) * 100)  # Rough estimate
            except:
                size_bytes = 1024  # Default size estimate
            
            # Check if value is too large
            if size_bytes > self.max_size_bytes:
                self.logger.warning(f"Value too large for cache: {size_bytes} bytes")
                return False
            
            # Remove existing entry if present
            if key in self.entries:
                old_entry = self.entries[key]
                self.stats.update_size(-old_entry.size_bytes, -1)
            
            # Make space if needed
            while (self.stats.size_bytes + size_bytes > self.max_size_bytes and 
                   len(self.entries) > 0):
                self._evict_one()
            
            # Create new entry
            entry = CacheEntry(
                key=key,
                value=value,
                size_bytes=size_bytes,
                created_time=time.time(),
                accessed_time=time.time(),
                access_count=1,
                ttl=ttl
            )
            
            self.entries[key] = entry
            self.stats.update_size(size_bytes, 1)
            
            # Update policy-specific structures
            if self.policy == CachePolicy.LRU:
                if key in self.access_order:
                    self.access_order.remove(key)
                self.access_order.append(key)
            elif self.policy == CachePolicy.FIFO:
                if key not in self.insertion_order:
                    self.insertion_order.append(key)
            
            return True
    
    def delete(self, key: str) -> bool:
        """Delete value from memory cache"""
        with self.lock:
            if key not in self.entries:
                return False
            
            entry = self.entries[key]
            del self.entries[key]
            self.stats.update_size(-entry.size_bytes, -1)
            
            # Update policy-specific structures
            if self.policy == CachePolicy.LRU and key in self.access_order:
                self.access_order.remove(key)
            elif self.policy == CachePolicy.FIFO and key in self.insertion_order:
                self.insertion_order.remove(key)
            
            return True
    
    def clear(self):
        """Clear all cache entries"""
        with self.lock:
            self.entries.clear()
            self.stats.size_bytes = 0
            self.stats.entry_count = 0
            
            if self.policy == CachePolicy.LRU:
                self.access_order.clear()
            elif self.policy == CachePolicy.FIFO:
                self.insertion_order.clear()
    
    def keys(self) -> List[str]:
        """Get all cache keys"""
        with self.lock:
            return list(self.entries.keys())
    
    def _evict_one(self):
        """Evict one entry based on policy"""
        if not self.entries:
            return
        
        evict_key = None
        
        if self.policy == CachePolicy.LRU:
            # Remove least recently used
            evict_key = self.access_order[0] if self.access_order else next(iter(self.entries))
        
        elif self.policy == CachePolicy.LFU:
            # Remove least frequently used
            evict_key = min(self.entries.keys(), 
                          key=lambda k: self.entries[k].access_count)
        
        elif self.policy == CachePolicy.FIFO:
            # Remove first inserted
            evict_key = self.insertion_order[0] if self.insertion_order else next(iter(self.entries))
        
        elif self.policy == CachePolicy.TTL:
            # Remove oldest or expired entry
            oldest_key = min(self.entries.keys(), 
                           key=lambda k: self.entries[k].created_time)
            if self.entries[oldest_key].is_expired:
                evict_key = oldest_key
            else:
                evict_key = oldest_key
        
        elif self.policy == CachePolicy.SIZE_BASED:
            # Remove largest entry
            evict_key = max(self.entries.keys(), 
                          key=lambda k: self.entries[k].size_bytes)
        
        if evict_key:
            self.delete(evict_key)
            self.stats.record_eviction()


class DiskCache(BaseCache):
    """Disk-based cache with file system storage"""
    
    def __init__(self, name: str = "disk_cache", 
                 cache_dir: Optional[Path] = None,
                 max_size_bytes: int = 1024 * 1024 * 1024):  # 1GB default
        super().__init__(name, max_size_bytes)
        
        if cache_dir is None:
            cache_dir = Path(tempfile.gettempdir()) / "prism_cache" / name
        
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Metadata storage
        self.metadata_file = self.cache_dir / ".cache_metadata.pickle"
        self.entries: Dict[str, CacheEntry] = self._load_metadata()
        
        # Validate existing cache files
        self._validate_cache_files()
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from disk cache"""
        with self.lock:
            if key not in self.entries:
                self.stats.record_miss()
                return None
            
            entry = self.entries[key]
            
            # Check if expired
            if entry.is_expired:
                self.delete(key)
                self.stats.record_miss()
                return None
            
            # Load value from disk
            file_path = self._get_file_path(key)
            if not file_path.exists():
                # File missing, remove from metadata
                del self.entries[key]
                self.stats.record_miss()
                self.stats.update_size(-entry.size_bytes, -1)
                self._save_metadata()
                return None
            
            try:
                with open(file_path, 'rb') as f:
                    value = pickle.load(f)
                
                # Update access statistics
                entry.touch()
                self.stats.record_hit()
                self._save_metadata()
                
                return value
                
            except Exception as e:
                self.logger.error(f"Failed to load cache entry {key}: {e}")
                self.delete(key)
                self.stats.record_miss()
                return None
    
    def put(self, key: str, value: Any, ttl: Optional[float] = None) -> bool:
        """Put value in disk cache"""
        with self.lock:
            file_path = self._get_file_path(key)
            
            try:
                # Serialize value to disk
                with open(file_path, 'wb') as f:
                    pickle.dump(value, f)
                
                # Get file size
                size_bytes = file_path.stat().st_size
                
                # Check if too large
                if size_bytes > self.max_size_bytes:
                    file_path.unlink()
                    self.logger.warning(f"Value too large for disk cache: {size_bytes} bytes")
                    return False
                
                # Remove existing entry if present
                if key in self.entries:
                    old_entry = self.entries[key]
                    self.stats.update_size(-old_entry.size_bytes, -1)
                
                # Make space if needed
                while (self.stats.size_bytes + size_bytes > self.max_size_bytes and 
                       len(self.entries) > 0):
                    self._evict_lru()
                
                # Create new entry
                entry = CacheEntry(
                    key=key,
                    value=None,  # Don't store value in memory for disk cache
                    size_bytes=size_bytes,
                    created_time=time.time(),
                    accessed_time=time.time(),
                    access_count=1,
                    ttl=ttl
                )
                
                self.entries[key] = entry
                self.stats.update_size(size_bytes, 1)
                self._save_metadata()
                
                return True
                
            except Exception as e:
                self.logger.error(f"Failed to store cache entry {key}: {e}")
                if file_path.exists():
                    file_path.unlink()
                return False
    
    def delete(self, key: str) -> bool:
        """Delete value from disk cache"""
        with self.lock:
            if key not in self.entries:
                return False
            
            entry = self.entries[key]
            file_path = self._get_file_path(key)
            
            # Remove file
            if file_path.exists():
                file_path.unlink()
            
            # Remove from metadata
            del self.entries[key]
            self.stats.update_size(-entry.size_bytes, -1)
            self._save_metadata()
            
            return True
    
    def clear(self):
        """Clear all cache entries"""
        with self.lock:
            # Remove all files
            for key in list(self.entries.keys()):
                file_path = self._get_file_path(key)
                if file_path.exists():
                    file_path.unlink()
            
            # Clear metadata
            self.entries.clear()
            self.stats.size_bytes = 0
            self.stats.entry_count = 0
            self._save_metadata()
    
    def keys(self) -> List[str]:
        """Get all cache keys"""
        with self.lock:
            return list(self.entries.keys())
    
    def _get_file_path(self, key: str) -> Path:
        """Get file path for cache key"""
        # Use hash of key to avoid filesystem issues
        key_hash = hashlib.md5(key.encode()).hexdigest()
        return self.cache_dir / f"{key_hash}.cache"
    
    def _evict_lru(self):
        """Evict least recently used entry"""
        if not self.entries:
            return
        
        lru_key = min(self.entries.keys(), 
                     key=lambda k: self.entries[k].accessed_time)
        self.delete(lru_key)
        self.stats.record_eviction()
    
    def _load_metadata(self) -> Dict[str, CacheEntry]:
        """Load metadata from disk"""
        if not self.metadata_file.exists():
            return {}
        
        try:
            with open(self.metadata_file, 'rb') as f:
                entries_data = pickle.load(f)
            
            # Convert back to CacheEntry objects
            entries = {}
            for key, entry_data in entries_data.items():
                entry = CacheEntry(**entry_data)
                entries[key] = entry
                # Update stats
                self.stats.update_size(entry.size_bytes, 1)
            
            return entries
            
        except Exception as e:
            self.logger.error(f"Failed to load cache metadata: {e}")
            return {}
    
    def _save_metadata(self):
        """Save metadata to disk"""
        try:
            # Convert CacheEntry objects to serializable data
            entries_data = {}
            for key, entry in self.entries.items():
                entries_data[key] = asdict(entry)
            
            with open(self.metadata_file, 'wb') as f:
                pickle.dump(entries_data, f)
                
        except Exception as e:
            self.logger.error(f"Failed to save cache metadata: {e}")
    
    def _validate_cache_files(self):
        """Validate that cache files exist and clean up orphaned metadata"""
        with self.lock:
            keys_to_remove = []
            
            for key, entry in self.entries.items():
                file_path = self._get_file_path(key)
                if not file_path.exists():
                    keys_to_remove.append(key)
                    self.stats.update_size(-entry.size_bytes, -1)
            
            for key in keys_to_remove:
                del self.entries[key]
            
            if keys_to_remove:
                self.logger.info(f"Cleaned up {len(keys_to_remove)} orphaned cache entries")
                self._save_metadata()


class CacheManager:
    """Multi-level cache manager with intelligent data placement"""
    
    def __init__(self, l1_size_mb: int = 100, l2_size_mb: int = 1000):
        self.logger = PrismLogger("cache_manager")
        
        # Create cache levels
        self.l1_cache = MemoryCache("l1_memory", l1_size_mb * 1024 * 1024)
        self.l2_cache = DiskCache("l2_disk", max_size_bytes=l2_size_mb * 1024 * 1024)
        
        # Cache hierarchy
        self.caches = [self.l1_cache, self.l2_cache]
        
        # Access pattern tracking
        self.access_patterns: Dict[str, Dict[str, Any]] = {}
        self.lock = threading.RLock()
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from multi-level cache"""
        with self.lock:
            # Try each cache level
            for level, cache in enumerate(self.caches):
                value = cache.get(key)
                if value is not None:
                    # Promote to higher level caches if beneficial
                    self._promote_entry(key, value, level)
                    self._track_access(key, level)
                    return value
            
            self._track_access(key, -1)  # Cache miss
            return None
    
    def put(self, key: str, value: Any, ttl: Optional[float] = None, 
            preferred_level: Optional[int] = None) -> bool:
        """Put value in multi-level cache"""
        with self.lock:
            # Choose cache level based on data characteristics
            if preferred_level is not None:
                target_level = min(preferred_level, len(self.caches) - 1)
            else:
                target_level = self._choose_cache_level(key, value)
            
            # Store in target level and potentially higher levels
            success = False
            
            for level in range(target_level, len(self.caches)):
                if self.caches[level].put(key, value, ttl):
                    success = True
                    break
            
            return success
    
    def delete(self, key: str) -> bool:
        """Delete value from all cache levels"""
        with self.lock:
            deleted = False
            for cache in self.caches:
                if cache.delete(key):
                    deleted = True
            
            # Remove from access patterns
            if key in self.access_patterns:
                del self.access_patterns[key]
            
            return deleted
    
    def clear(self):
        """Clear all cache levels"""
        with self.lock:
            for cache in self.caches:
                cache.clear()
            self.access_patterns.clear()
    
    def get_stats(self) -> Dict[str, Any]:
        """Get comprehensive cache statistics"""
        with self.lock:
            stats = {
                'total_hit_ratio': 0.0,
                'levels': {}
            }
            
            total_hits = 0
            total_requests = 0
            
            for level, cache in enumerate(self.caches):
                cache_stats = cache.get_stats()
                stats['levels'][f'L{level+1}'] = cache_stats
                
                total_hits += cache_stats['hits']
                total_requests += cache_stats['hits'] + cache_stats['misses']
            
            stats['total_hit_ratio'] = total_hits / total_requests if total_requests > 0 else 0.0
            stats['access_patterns_tracked'] = len(self.access_patterns)
            
            return stats
    
    def _choose_cache_level(self, key: str, value: Any) -> int:
        """Choose appropriate cache level for data"""
        try:
            # Estimate value size
            size_bytes = sys.getsizeof(value)
            
            # Small values go to L1 (memory)
            if size_bytes < 1024 * 1024:  # < 1MB
                return 0
            
            # Larger values go to L2 (disk)
            return 1
            
        except:
            # Default to L2 for safety
            return 1
    
    def _promote_entry(self, key: str, value: Any, found_level: int):
        """Promote frequently accessed entries to higher cache levels"""
        if found_level == 0:
            return  # Already in highest level
        
        # Check access pattern
        pattern = self.access_patterns.get(key, {})
        access_count = pattern.get('access_count', 0)
        
        # Promote if accessed frequently
        if access_count > 5:
            # Try to promote to L1
            self.l1_cache.put(key, value)
    
    def _track_access(self, key: str, level: int):
        """Track access patterns for intelligent caching"""
        current_time = time.time()
        
        if key not in self.access_patterns:
            self.access_patterns[key] = {
                'access_count': 0,
                'last_access': current_time,
                'hit_levels': [],
                'miss_count': 0
            }
        
        pattern = self.access_patterns[key]
        pattern['last_access'] = current_time
        
        if level >= 0:
            pattern['access_count'] += 1
            pattern['hit_levels'].append(level)
        else:
            pattern['miss_count'] += 1
        
        # Limit history size
        if len(pattern['hit_levels']) > 100:
            pattern['hit_levels'] = pattern['hit_levels'][-50:]


# Convenience functions
def create_memory_cache(max_size_mb: int = 100, 
                       policy: CachePolicy = CachePolicy.LRU) -> MemoryCache:
    """Create memory cache with specified settings"""
    return MemoryCache(max_size_bytes=max_size_mb * 1024 * 1024, policy=policy)


def create_disk_cache(cache_dir: Optional[Path] = None, 
                     max_size_mb: int = 1000) -> DiskCache:
    """Create disk cache with specified settings"""
    return DiskCache(cache_dir=cache_dir, max_size_bytes=max_size_mb * 1024 * 1024)


def create_cache_manager(l1_size_mb: int = 100, l2_size_mb: int = 1000) -> CacheManager:
    """Create multi-level cache manager"""
    return CacheManager(l1_size_mb=l1_size_mb, l2_size_mb=l2_size_mb)


# Global cache instance
_global_cache_manager = None

def get_global_cache() -> CacheManager:
    """Get global cache manager instance"""
    global _global_cache_manager
    if _global_cache_manager is None:
        _global_cache_manager = create_cache_manager()
    return _global_cache_manager


# Cache decorators
def cached(cache_manager: Optional[CacheManager] = None, 
           ttl: Optional[float] = None):
    """Decorator to cache function results"""
    def decorator(func: Callable) -> Callable:
        cache_mgr = cache_manager or get_global_cache()
        
        def wrapper(*args, **kwargs):
            # Generate cache key
            key_data = {
                'func': func.__name__,
                'args': args,
                'kwargs': sorted(kwargs.items())
            }
            cache_key = hashlib.md5(str(key_data).encode()).hexdigest()
            
            # Try to get from cache
            result = cache_mgr.get(cache_key)
            if result is not None:
                return result
            
            # Compute result and cache it
            result = func(*args, **kwargs)
            cache_mgr.put(cache_key, result, ttl=ttl)
            
            return result
        
        return wrapper
    return decorator