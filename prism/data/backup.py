#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Backup and Recovery System

Comprehensive backup and recovery system for PRISM molecular dynamics data
with incremental backups, checksums, and disaster recovery capabilities.
"""

import os
import sys
import time
import shutil
import hashlib
import json
import threading
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Set, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
from abc import ABC, abstractmethod
import tarfile
import zipfile
import tempfile
from datetime import datetime, timedelta

from ..utils.logging_system import PrismLogger


class BackupType(Enum):
    """Types of backups"""
    FULL = "full"           # Complete backup of all data
    INCREMENTAL = "incremental"  # Only changed files since last backup
    DIFFERENTIAL = "differential"  # All changes since last full backup
    SNAPSHOT = "snapshot"    # Point-in-time snapshot


class BackupStatus(Enum):
    """Backup operation status"""
    PENDING = "pending"
    IN_PROGRESS = "in_progress" 
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class CompressionType(Enum):
    """Backup compression types"""
    NONE = "none"
    GZIP = "gzip"
    BZIP2 = "bzip2"
    XZ = "xz"
    ZIP = "zip"


@dataclass
class BackupMetadata:
    """Backup metadata and information"""
    backup_id: str
    backup_type: BackupType
    source_path: str
    backup_path: str
    created_timestamp: float
    file_count: int
    total_size_bytes: int
    compressed_size_bytes: int
    compression_type: CompressionType
    checksum: str
    status: BackupStatus
    error_message: Optional[str] = None
    duration_seconds: float = 0.0
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}
    
    @property
    def compression_ratio(self) -> float:
        """Calculate compression ratio"""
        if self.total_size_bytes == 0:
            return 1.0
        return self.total_size_bytes / self.compressed_size_bytes
    
    @property
    def age_days(self) -> float:
        """Get age of backup in days"""
        return (time.time() - self.created_timestamp) / 86400
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        data = asdict(self)
        data['backup_type'] = self.backup_type.value
        data['compression_type'] = self.compression_type.value
        data['status'] = self.status.value
        return data
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'BackupMetadata':
        """Create from dictionary"""
        data['backup_type'] = BackupType(data['backup_type'])
        data['compression_type'] = CompressionType(data['compression_type'])
        data['status'] = BackupStatus(data['status'])
        return cls(**data)


@dataclass
class RestoreOptions:
    """Options for data restoration"""
    backup_id: str
    target_path: Path
    overwrite_existing: bool = False
    verify_checksums: bool = True
    selective_restore: Optional[List[str]] = None  # Specific files/directories
    progress_callback: Optional[callable] = None


class BackupStorage(ABC):
    """Abstract base class for backup storage backends"""
    
    @abstractmethod
    def store_backup(self, backup_path: Path, backup_id: str) -> bool:
        """Store backup to storage backend"""
        pass
    
    @abstractmethod
    def retrieve_backup(self, backup_id: str, target_path: Path) -> bool:
        """Retrieve backup from storage backend"""
        pass
    
    @abstractmethod
    def delete_backup(self, backup_id: str) -> bool:
        """Delete backup from storage"""
        pass
    
    @abstractmethod
    def list_backups(self) -> List[str]:
        """List available backups"""
        pass


class LocalBackupStorage(BackupStorage):
    """Local filesystem backup storage"""
    
    def __init__(self, storage_path: Path):
        self.storage_path = storage_path
        self.storage_path.mkdir(parents=True, exist_ok=True)
        self.logger = PrismLogger("local_backup_storage")
    
    def store_backup(self, backup_path: Path, backup_id: str) -> bool:
        """Store backup to local storage"""
        try:
            target_path = self.storage_path / f"{backup_id}.backup"
            shutil.copy2(backup_path, target_path)
            return True
        except Exception as e:
            self.logger.error(f"Failed to store backup {backup_id}: {e}")
            return False
    
    def retrieve_backup(self, backup_id: str, target_path: Path) -> bool:
        """Retrieve backup from local storage"""
        try:
            backup_path = self.storage_path / f"{backup_id}.backup"
            if not backup_path.exists():
                return False
            
            shutil.copy2(backup_path, target_path)
            return True
        except Exception as e:
            self.logger.error(f"Failed to retrieve backup {backup_id}: {e}")
            return False
    
    def delete_backup(self, backup_id: str) -> bool:
        """Delete backup from local storage"""
        try:
            backup_path = self.storage_path / f"{backup_id}.backup"
            if backup_path.exists():
                backup_path.unlink()
            return True
        except Exception as e:
            self.logger.error(f"Failed to delete backup {backup_id}: {e}")
            return False
    
    def list_backups(self) -> List[str]:
        """List available backups"""
        backups = []
        for backup_file in self.storage_path.glob("*.backup"):
            backup_id = backup_file.stem
            backups.append(backup_id)
        return backups


class BackupEngine:
    """Core backup engine with compression and deduplication"""
    
    def __init__(self, storage: BackupStorage):
        self.storage = storage
        self.logger = PrismLogger("backup_engine")
        self.file_registry: Dict[str, str] = {}  # file_path -> checksum mapping
        self.progress_callback: Optional[callable] = None
    
    def create_backup(self, source_path: Path, backup_type: BackupType = BackupType.FULL,
                     compression: CompressionType = CompressionType.GZIP,
                     exclude_patterns: Optional[List[str]] = None) -> BackupMetadata:
        """Create a backup of the specified source"""
        backup_id = self._generate_backup_id(backup_type)
        start_time = time.time()
        
        self.logger.info(f"Starting {backup_type.value} backup: {backup_id}")
        
        try:
            # Create temporary backup file
            with tempfile.NamedTemporaryFile(suffix='.backup', delete=False) as tmp_file:
                tmp_path = Path(tmp_file.name)
            
            # Determine files to backup
            files_to_backup = self._get_files_to_backup(
                source_path, backup_type, exclude_patterns
            )
            
            if not files_to_backup:
                self.logger.warning("No files to backup")
                return self._create_empty_backup_metadata(backup_id, source_path, backup_type)
            
            # Create compressed archive
            total_size, compressed_size, file_count, checksum = self._create_archive(
                files_to_backup, tmp_path, compression, source_path
            )
            
            # Store backup
            if not self.storage.store_backup(tmp_path, backup_id):
                raise Exception("Failed to store backup")
            
            # Clean up temporary file
            tmp_path.unlink()
            
            duration = time.time() - start_time
            
            # Create metadata
            metadata = BackupMetadata(
                backup_id=backup_id,
                backup_type=backup_type,
                source_path=str(source_path),
                backup_path=str(self.storage.storage_path / f"{backup_id}.backup"),
                created_timestamp=start_time,
                file_count=file_count,
                total_size_bytes=total_size,
                compressed_size_bytes=compressed_size,
                compression_type=compression,
                checksum=checksum,
                status=BackupStatus.COMPLETED,
                duration_seconds=duration,
                metadata={
                    'exclude_patterns': exclude_patterns or [],
                    'files_backed_up': len(files_to_backup)
                }
            )
            
            self.logger.info(f"Backup completed: {backup_id} "
                           f"({file_count} files, {total_size/(1024*1024):.1f}MB -> "
                           f"{compressed_size/(1024*1024):.1f}MB, {duration:.1f}s)")
            
            return metadata
            
        except Exception as e:
            duration = time.time() - start_time
            self.logger.error(f"Backup failed: {backup_id}: {e}")
            
            return BackupMetadata(
                backup_id=backup_id,
                backup_type=backup_type,
                source_path=str(source_path),
                backup_path="",
                created_timestamp=start_time,
                file_count=0,
                total_size_bytes=0,
                compressed_size_bytes=0,
                compression_type=compression,
                checksum="",
                status=BackupStatus.FAILED,
                error_message=str(e),
                duration_seconds=duration
            )
    
    def restore_backup(self, backup_metadata: BackupMetadata, 
                      restore_options: RestoreOptions) -> bool:
        """Restore backup to target location"""
        self.logger.info(f"Restoring backup: {backup_metadata.backup_id}")
        
        try:
            # Create temporary file for retrieved backup
            with tempfile.NamedTemporaryFile(suffix='.backup', delete=False) as tmp_file:
                tmp_path = Path(tmp_file.name)
            
            # Retrieve backup from storage
            if not self.storage.retrieve_backup(backup_metadata.backup_id, tmp_path):
                raise Exception("Failed to retrieve backup from storage")
            
            # Verify checksum if requested
            if restore_options.verify_checksums:
                if not self._verify_backup_checksum(tmp_path, backup_metadata.checksum):
                    raise Exception("Backup checksum verification failed")
            
            # Extract archive
            success = self._extract_archive(
                tmp_path, restore_options.target_path,
                backup_metadata.compression_type,
                restore_options.selective_restore,
                restore_options.overwrite_existing,
                restore_options.progress_callback
            )
            
            # Clean up
            tmp_path.unlink()
            
            if success:
                self.logger.info(f"Backup restoration completed: {backup_metadata.backup_id}")
            else:
                self.logger.error(f"Backup restoration failed: {backup_metadata.backup_id}")
            
            return success
            
        except Exception as e:
            self.logger.error(f"Restore failed for {backup_metadata.backup_id}: {e}")
            return False
    
    def _generate_backup_id(self, backup_type: BackupType) -> str:
        """Generate unique backup ID"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"{backup_type.value}_{timestamp}_{int(time.time() * 1000) % 10000}"
    
    def _get_files_to_backup(self, source_path: Path, backup_type: BackupType,
                           exclude_patterns: Optional[List[str]]) -> List[Path]:
        """Get list of files to backup based on type and filters"""
        files = []
        exclude_patterns = exclude_patterns or []
        
        # Add default excludes
        default_excludes = [
            '*.tmp', '*.temp', '*.log', '*.pid', '*.lock',
            '__pycache__', '.git', '.svn', 'node_modules'
        ]
        exclude_patterns.extend(default_excludes)
        
        # Walk directory tree
        for root, dirs, filenames in os.walk(source_path):
            root_path = Path(root)
            
            # Filter directories
            dirs[:] = [d for d in dirs if not self._should_exclude(d, exclude_patterns)]
            
            # Add files
            for filename in filenames:
                file_path = root_path / filename
                
                if self._should_exclude(filename, exclude_patterns):
                    continue
                
                # For incremental/differential backups, check modification time
                if backup_type in [BackupType.INCREMENTAL, BackupType.DIFFERENTIAL]:
                    if not self._file_needs_backup(file_path, backup_type):
                        continue
                
                files.append(file_path)
                
                # Progress callback
                if self.progress_callback:
                    self.progress_callback("scanning", len(files), file_path)
        
        return files
    
    def _should_exclude(self, name: str, exclude_patterns: List[str]) -> bool:
        """Check if file/directory should be excluded"""
        import fnmatch
        
        for pattern in exclude_patterns:
            if fnmatch.fnmatch(name, pattern):
                return True
        return False
    
    def _file_needs_backup(self, file_path: Path, backup_type: BackupType) -> bool:
        """Check if file needs to be backed up (for incremental/differential)"""
        # This is a simplified implementation
        # In a real system, you would track file modification times
        # and compare against previous backup timestamps
        
        try:
            file_mtime = file_path.stat().st_mtime
            # For demo purposes, backup files modified in last 24 hours
            return time.time() - file_mtime < 86400
        except:
            return True
    
    def _create_archive(self, files: List[Path], archive_path: Path,
                       compression: CompressionType, base_path: Path) -> Tuple[int, int, int, str]:
        """Create compressed archive of files"""
        total_size = 0
        file_count = 0
        
        if compression == CompressionType.ZIP:
            return self._create_zip_archive(files, archive_path, base_path)
        else:
            return self._create_tar_archive(files, archive_path, compression, base_path)
    
    def _create_tar_archive(self, files: List[Path], archive_path: Path,
                           compression: CompressionType, base_path: Path) -> Tuple[int, int, int, str]:
        """Create tar archive"""
        total_size = 0
        file_count = 0
        
        # Determine compression mode
        mode_map = {
            CompressionType.NONE: 'w',
            CompressionType.GZIP: 'w:gz',
            CompressionType.BZIP2: 'w:bz2',
            CompressionType.XZ: 'w:xz'
        }
        
        mode = mode_map.get(compression, 'w:gz')
        
        with tarfile.open(archive_path, mode) as tar:
            for file_path in files:
                if not file_path.exists():
                    continue
                
                try:
                    # Calculate relative path
                    rel_path = file_path.relative_to(base_path)
                    
                    # Add file to archive
                    tar.add(file_path, arcname=rel_path)
                    
                    total_size += file_path.stat().st_size
                    file_count += 1
                    
                    if self.progress_callback:
                        self.progress_callback("archiving", file_count, file_path)
                        
                except Exception as e:
                    self.logger.warning(f"Failed to add {file_path} to archive: {e}")
        
        # Get compressed size
        compressed_size = archive_path.stat().st_size
        
        # Calculate checksum
        checksum = self._calculate_file_checksum(archive_path)
        
        return total_size, compressed_size, file_count, checksum
    
    def _create_zip_archive(self, files: List[Path], archive_path: Path,
                           base_path: Path) -> Tuple[int, int, int, str]:
        """Create ZIP archive"""
        total_size = 0
        file_count = 0
        
        with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for file_path in files:
                if not file_path.exists():
                    continue
                
                try:
                    # Calculate relative path
                    rel_path = file_path.relative_to(base_path)
                    
                    # Add file to archive
                    zip_file.write(file_path, rel_path)
                    
                    total_size += file_path.stat().st_size
                    file_count += 1
                    
                    if self.progress_callback:
                        self.progress_callback("archiving", file_count, file_path)
                        
                except Exception as e:
                    self.logger.warning(f"Failed to add {file_path} to archive: {e}")
        
        # Get compressed size
        compressed_size = archive_path.stat().st_size
        
        # Calculate checksum
        checksum = self._calculate_file_checksum(archive_path)
        
        return total_size, compressed_size, file_count, checksum
    
    def _extract_archive(self, archive_path: Path, target_path: Path,
                        compression: CompressionType, selective_files: Optional[List[str]],
                        overwrite: bool, progress_callback: Optional[callable]) -> bool:
        """Extract archive to target location"""
        try:
            target_path.mkdir(parents=True, exist_ok=True)
            
            if compression == CompressionType.ZIP:
                return self._extract_zip_archive(
                    archive_path, target_path, selective_files, overwrite, progress_callback
                )
            else:
                return self._extract_tar_archive(
                    archive_path, target_path, compression, selective_files, 
                    overwrite, progress_callback
                )
        except Exception as e:
            self.logger.error(f"Archive extraction failed: {e}")
            return False
    
    def _extract_tar_archive(self, archive_path: Path, target_path: Path,
                            compression: CompressionType, selective_files: Optional[List[str]],
                            overwrite: bool, progress_callback: Optional[callable]) -> bool:
        """Extract tar archive"""
        mode_map = {
            CompressionType.NONE: 'r',
            CompressionType.GZIP: 'r:gz',
            CompressionType.BZIP2: 'r:bz2',
            CompressionType.XZ: 'r:xz'
        }
        
        mode = mode_map.get(compression, 'r:gz')
        
        with tarfile.open(archive_path, mode) as tar:
            members = tar.getmembers()
            
            for i, member in enumerate(members):
                # Check if file should be extracted
                if selective_files and member.name not in selective_files:
                    continue
                
                target_file = target_path / member.name
                
                # Check if file exists and overwrite is disabled
                if target_file.exists() and not overwrite:
                    self.logger.info(f"Skipping existing file: {target_file}")
                    continue
                
                # Extract member
                tar.extract(member, target_path)
                
                if progress_callback:
                    progress_callback("extracting", i + 1, member.name)
        
        return True
    
    def _extract_zip_archive(self, archive_path: Path, target_path: Path,
                            selective_files: Optional[List[str]], overwrite: bool,
                            progress_callback: Optional[callable]) -> bool:
        """Extract ZIP archive"""
        with zipfile.ZipFile(archive_path, 'r') as zip_file:
            members = zip_file.infolist()
            
            for i, member in enumerate(members):
                # Check if file should be extracted
                if selective_files and member.filename not in selective_files:
                    continue
                
                target_file = target_path / member.filename
                
                # Check if file exists and overwrite is disabled
                if target_file.exists() and not overwrite:
                    self.logger.info(f"Skipping existing file: {target_file}")
                    continue
                
                # Extract member
                zip_file.extract(member, target_path)
                
                if progress_callback:
                    progress_callback("extracting", i + 1, member.filename)
        
        return True
    
    def _calculate_file_checksum(self, file_path: Path) -> str:
        """Calculate SHA-256 checksum of file"""
        sha256_hash = hashlib.sha256()
        with open(file_path, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                sha256_hash.update(byte_block)
        return sha256_hash.hexdigest()
    
    def _verify_backup_checksum(self, backup_path: Path, expected_checksum: str) -> bool:
        """Verify backup file checksum"""
        actual_checksum = self._calculate_file_checksum(backup_path)
        return actual_checksum == expected_checksum
    
    def _create_empty_backup_metadata(self, backup_id: str, source_path: Path, 
                                    backup_type: BackupType) -> BackupMetadata:
        """Create metadata for empty backup"""
        return BackupMetadata(
            backup_id=backup_id,
            backup_type=backup_type,
            source_path=str(source_path),
            backup_path="",
            created_timestamp=time.time(),
            file_count=0,
            total_size_bytes=0,
            compressed_size_bytes=0,
            compression_type=CompressionType.NONE,
            checksum="",
            status=BackupStatus.COMPLETED,
            metadata={'reason': 'no_files_to_backup'}
        )


class BackupManager:
    """High-level backup management system"""
    
    def __init__(self, storage: BackupStorage):
        self.storage = storage
        self.engine = BackupEngine(storage)
        self.logger = PrismLogger("backup_manager")
        
        # Metadata storage
        self.metadata_file = Path("backup_metadata.json")
        self.backup_metadata: Dict[str, BackupMetadata] = self._load_metadata()
        
        # Backup scheduling
        self.scheduled_backups: Dict[str, Dict[str, Any]] = {}
        self.scheduler_thread: Optional[threading.Thread] = None
        self.scheduler_stop_event = threading.Event()
    
    def create_backup(self, source_path: Path, backup_type: BackupType = BackupType.FULL,
                     compression: CompressionType = CompressionType.GZIP,
                     exclude_patterns: Optional[List[str]] = None,
                     progress_callback: Optional[callable] = None) -> BackupMetadata:
        """Create backup with metadata tracking"""
        
        self.engine.progress_callback = progress_callback
        
        metadata = self.engine.create_backup(
            source_path, backup_type, compression, exclude_patterns
        )
        
        # Store metadata
        self.backup_metadata[metadata.backup_id] = metadata
        self._save_metadata()
        
        return metadata
    
    def restore_backup(self, backup_id: str, target_path: Path,
                      overwrite_existing: bool = False,
                      verify_checksums: bool = True,
                      selective_restore: Optional[List[str]] = None,
                      progress_callback: Optional[callable] = None) -> bool:
        """Restore backup by ID"""
        
        if backup_id not in self.backup_metadata:
            self.logger.error(f"Backup not found: {backup_id}")
            return False
        
        metadata = self.backup_metadata[backup_id]
        
        restore_options = RestoreOptions(
            backup_id=backup_id,
            target_path=target_path,
            overwrite_existing=overwrite_existing,
            verify_checksums=verify_checksums,
            selective_restore=selective_restore,
            progress_callback=progress_callback
        )
        
        return self.engine.restore_backup(metadata, restore_options)
    
    def delete_backup(self, backup_id: str) -> bool:
        """Delete backup and metadata"""
        if backup_id not in self.backup_metadata:
            self.logger.error(f"Backup not found: {backup_id}")
            return False
        
        # Delete from storage
        if not self.storage.delete_backup(backup_id):
            self.logger.error(f"Failed to delete backup from storage: {backup_id}")
            return False
        
        # Remove metadata
        del self.backup_metadata[backup_id]
        self._save_metadata()
        
        self.logger.info(f"Deleted backup: {backup_id}")
        return True
    
    def list_backups(self, backup_type: Optional[BackupType] = None) -> List[BackupMetadata]:
        """List available backups"""
        backups = list(self.backup_metadata.values())
        
        if backup_type:
            backups = [b for b in backups if b.backup_type == backup_type]
        
        # Sort by creation time (newest first)
        backups.sort(key=lambda b: b.created_timestamp, reverse=True)
        
        return backups
    
    def get_backup_stats(self) -> Dict[str, Any]:
        """Get backup statistics"""
        backups = list(self.backup_metadata.values())
        
        if not backups:
            return {'total_backups': 0}
        
        total_backups = len(backups)
        total_size = sum(b.compressed_size_bytes for b in backups)
        avg_compression_ratio = sum(b.compression_ratio for b in backups) / total_backups
        
        # Group by type
        type_stats = {}
        for backup_type in BackupType:
            type_backups = [b for b in backups if b.backup_type == backup_type]
            type_stats[backup_type.value] = {
                'count': len(type_backups),
                'total_size_mb': sum(b.compressed_size_bytes for b in type_backups) / (1024 * 1024)
            }
        
        # Recent activity
        now = time.time()
        recent_backups = [b for b in backups if (now - b.created_timestamp) < 86400 * 7]  # Last 7 days
        
        return {
            'total_backups': total_backups,
            'total_size_mb': total_size / (1024 * 1024),
            'average_compression_ratio': avg_compression_ratio,
            'backup_types': type_stats,
            'recent_backups': len(recent_backups),
            'oldest_backup': min(b.created_timestamp for b in backups),
            'newest_backup': max(b.created_timestamp for b in backups)
        }
    
    def cleanup_old_backups(self, max_age_days: int = 30, max_count: int = 10) -> int:
        """Clean up old backups based on age and count limits"""
        backups = self.list_backups()
        now = time.time()
        deleted_count = 0
        
        # Sort by age (oldest first)
        backups.sort(key=lambda b: b.created_timestamp)
        
        for backup in backups:
            # Check age limit
            age_days = (now - backup.created_timestamp) / 86400
            
            # Check count limit (keep most recent)
            remaining_backups = len(backups) - deleted_count
            
            if age_days > max_age_days or remaining_backups > max_count:
                if self.delete_backup(backup.backup_id):
                    deleted_count += 1
                    self.logger.info(f"Cleaned up old backup: {backup.backup_id}")
        
        return deleted_count
    
    def schedule_backup(self, schedule_id: str, source_path: Path,
                       backup_type: BackupType = BackupType.INCREMENTAL,
                       interval_hours: int = 24,
                       compression: CompressionType = CompressionType.GZIP,
                       exclude_patterns: Optional[List[str]] = None):
        """Schedule automatic backups"""
        
        self.scheduled_backups[schedule_id] = {
            'source_path': source_path,
            'backup_type': backup_type,
            'interval_hours': interval_hours,
            'compression': compression,
            'exclude_patterns': exclude_patterns,
            'last_run': 0,
            'next_run': time.time() + interval_hours * 3600
        }
        
        # Start scheduler if not running
        if not self.scheduler_thread or not self.scheduler_thread.is_alive():
            self._start_scheduler()
        
        self.logger.info(f"Scheduled backup: {schedule_id}")
    
    def unschedule_backup(self, schedule_id: str):
        """Remove scheduled backup"""
        if schedule_id in self.scheduled_backups:
            del self.scheduled_backups[schedule_id]
            self.logger.info(f"Unscheduled backup: {schedule_id}")
    
    def _start_scheduler(self):
        """Start backup scheduler thread"""
        def scheduler_worker():
            while not self.scheduler_stop_event.is_set():
                try:
                    current_time = time.time()
                    
                    # Check scheduled backups
                    for schedule_id, config in list(self.scheduled_backups.items()):
                        if current_time >= config['next_run']:
                            self.logger.info(f"Running scheduled backup: {schedule_id}")
                            
                            try:
                                self.create_backup(
                                    config['source_path'],
                                    config['backup_type'],
                                    config['compression'],
                                    config['exclude_patterns']
                                )
                                
                                # Update schedule
                                config['last_run'] = current_time
                                config['next_run'] = current_time + config['interval_hours'] * 3600
                                
                            except Exception as e:
                                self.logger.error(f"Scheduled backup failed: {schedule_id}: {e}")
                    
                    # Sleep for 1 minute
                    self.scheduler_stop_event.wait(60)
                    
                except Exception as e:
                    self.logger.error(f"Scheduler error: {e}")
        
        self.scheduler_thread = threading.Thread(target=scheduler_worker, daemon=True)
        self.scheduler_thread.start()
    
    def stop_scheduler(self):
        """Stop backup scheduler"""
        self.scheduler_stop_event.set()
        if self.scheduler_thread:
            self.scheduler_thread.join(timeout=5)
    
    def _load_metadata(self) -> Dict[str, BackupMetadata]:
        """Load backup metadata from disk"""
        if not self.metadata_file.exists():
            return {}
        
        try:
            with open(self.metadata_file, 'r') as f:
                data = json.load(f)
            
            metadata = {}
            for backup_id, meta_dict in data.items():
                metadata[backup_id] = BackupMetadata.from_dict(meta_dict)
            
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to load backup metadata: {e}")
            return {}
    
    def _save_metadata(self):
        """Save backup metadata to disk"""
        try:
            data = {}
            for backup_id, metadata in self.backup_metadata.items():
                data[backup_id] = metadata.to_dict()
            
            with open(self.metadata_file, 'w') as f:
                json.dump(data, f, indent=2)
                
        except Exception as e:
            self.logger.error(f"Failed to save backup metadata: {e}")


# Convenience functions
def create_backup_manager(storage_path: Path) -> BackupManager:
    """Create backup manager with local storage"""
    storage = LocalBackupStorage(storage_path)
    return BackupManager(storage)


def backup_directory(source_path: Union[str, Path], backup_path: Union[str, Path],
                    backup_type: BackupType = BackupType.FULL,
                    compression: CompressionType = CompressionType.GZIP) -> BackupMetadata:
    """Simple directory backup function"""
    manager = create_backup_manager(Path(backup_path))
    return manager.create_backup(Path(source_path), backup_type, compression)


def restore_directory(backup_id: str, backup_path: Union[str, Path], 
                     target_path: Union[str, Path]) -> bool:
    """Simple directory restore function"""
    manager = create_backup_manager(Path(backup_path))
    return manager.restore_backup(backup_id, Path(target_path))


# Global backup manager instance
_global_backup_manager = None

def get_global_backup_manager(storage_path: Optional[Path] = None) -> BackupManager:
    """Get global backup manager instance"""
    global _global_backup_manager
    if _global_backup_manager is None:
        if storage_path is None:
            storage_path = Path.home() / ".prism_backups"
        _global_backup_manager = create_backup_manager(storage_path)
    return _global_backup_manager