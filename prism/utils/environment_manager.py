#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Environment Manager

Advanced environment configuration and management system for PRISM.
Handles configuration files, environment variables, system detection,
and automatic setup.
"""

import os
import sys
import json
import yaml
import shutil
import platform
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
import tempfile

try:
    from .logging_system import PrismLogger
except ImportError:
    import logging
    PrismLogger = logging.getLogger


class EnvironmentType(Enum):
    """Environment types"""
    DEVELOPMENT = "development"
    TESTING = "testing"
    PRODUCTION = "production"
    RESEARCH = "research"


class ConfigurationLevel(Enum):
    """Configuration precedence levels"""
    SYSTEM = "system"      # System-wide defaults
    USER = "user"          # User-specific settings
    PROJECT = "project"    # Project-specific settings
    RUNTIME = "runtime"    # Runtime overrides


@dataclass
class SystemInfo:
    """System information structure"""
    platform: str
    architecture: str
    python_version: str
    cpu_count: int
    memory_gb: float
    disk_free_gb: float
    environment_type: EnvironmentType
    home_directory: str
    working_directory: str


@dataclass
class EnvironmentConfig:
    """Environment configuration structure"""
    system_info: SystemInfo
    paths: Dict[str, str]
    dependencies: Dict[str, bool]
    settings: Dict[str, Any]
    optimizations: Dict[str, Any]


class PrismEnvironmentManager:
    """Advanced environment configuration and management"""
    
    def __init__(self, environment_type: EnvironmentType = EnvironmentType.RESEARCH,
                 logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("environment_manager")
        self.environment_type = environment_type
        self.config = {}
        self.system_info = None
        
        # Configuration paths
        self.config_paths = self._initialize_config_paths()
        
        # Load configuration hierarchy
        self._load_configuration_hierarchy()
    
    def _initialize_config_paths(self) -> Dict[str, Path]:
        """Initialize configuration file paths"""
        home = Path.home()
        cwd = Path.cwd()
        
        paths = {
            # System-level configs
            "system_config": Path("/etc/prism/config.yaml") if os.name != 'nt' else Path("C:/ProgramData/PRISM/config.yaml"),
            
            # User-level configs
            "user_config": home / ".prism" / "config.yaml",
            "user_settings": home / ".prism" / "settings.json",
            "user_environment": home / ".prism" / "environment.yaml",
            
            # Project-level configs
            "project_config": cwd / ".prism" / "config.yaml",
            "project_settings": cwd / ".prism" / "settings.json",
            "project_environment": cwd / "prism_config.yaml",
            
            # Runtime configs
            "runtime_config": Path(tempfile.gettempdir()) / "prism_runtime.yaml"
        }
        
        # Ensure user config directory exists
        user_config_dir = home / ".prism"
        user_config_dir.mkdir(parents=True, exist_ok=True)
        
        return paths
    
    def _load_configuration_hierarchy(self):
        """Load configuration in hierarchical order"""
        self.config = {}
        
        # Load configurations in order of precedence
        config_levels = [
            ("system", self.config_paths["system_config"]),
            ("user", self.config_paths["user_config"]),
            ("project", self.config_paths["project_config"]),
            ("runtime", self.config_paths["runtime_config"])
        ]
        
        for level, config_path in config_levels:
            if config_path.exists():
                try:
                    config_data = self._load_config_file(config_path)
                    if config_data:
                        self.config = self._merge_configs(self.config, config_data)
                        self.logger.debug(f"Loaded {level} config from {config_path}")
                except Exception as e:
                    self.logger.warning(f"Failed to load {level} config: {e}")
        
        # Apply environment-specific overrides
        self._apply_environment_overrides()
    
    def _load_config_file(self, file_path: Path) -> Optional[Dict[str, Any]]:
        """Load configuration file (YAML or JSON)"""
        try:
            with open(file_path, 'r') as f:
                if file_path.suffix.lower() in ['.yaml', '.yml']:
                    return yaml.safe_load(f)
                elif file_path.suffix.lower() == '.json':
                    return json.load(f)
                else:
                    self.logger.warning(f"Unknown config file format: {file_path}")
                    return None
        except Exception as e:
            self.logger.error(f"Error loading config file {file_path}: {e}")
            return None
    
    def _merge_configs(self, base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
        """Recursively merge configuration dictionaries"""
        result = base.copy()
        
        for key, value in override.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = self._merge_configs(result[key], value)
            else:
                result[key] = value
        
        return result
    
    def _apply_environment_overrides(self):
        """Apply environment-specific configuration overrides"""
        env_overrides = {
            EnvironmentType.DEVELOPMENT: {
                "logging": {"level": "DEBUG", "console": True},
                "optimization": {"cache_enabled": True, "parallel_processing": False},
                "validation": {"strict_mode": False}
            },
            EnvironmentType.TESTING: {
                "logging": {"level": "INFO", "file": False},
                "optimization": {"cache_enabled": False, "parallel_processing": True},
                "validation": {"strict_mode": True}
            },
            EnvironmentType.PRODUCTION: {
                "logging": {"level": "WARNING", "console": False, "file": True},
                "optimization": {"cache_enabled": True, "parallel_processing": True},
                "validation": {"strict_mode": True}
            },
            EnvironmentType.RESEARCH: {
                "logging": {"level": "INFO", "console": True, "file": True},
                "optimization": {"cache_enabled": True, "parallel_processing": True},
                "validation": {"strict_mode": False}
            }
        }
        
        if self.environment_type in env_overrides:
            overrides = env_overrides[self.environment_type]
            self.config = self._merge_configs(self.config, overrides)
    
    def detect_system_info(self) -> SystemInfo:
        """Detect comprehensive system information"""
        # Basic system info
        platform_name = platform.system()
        architecture = platform.machine()
        python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
        
        # CPU information
        cpu_count = os.cpu_count() or 1
        
        # Memory information
        try:
            import psutil
            memory_gb = psutil.virtual_memory().total / (1024**3)
        except ImportError:
            # Fallback estimation
            if platform_name == "Linux":
                try:
                    with open("/proc/meminfo", "r") as f:
                        for line in f:
                            if line.startswith("MemTotal:"):
                                memory_kb = int(line.split()[1])
                                memory_gb = memory_kb / (1024**2)
                                break
                    else:
                        memory_gb = 8.0  # Default assumption
                except:
                    memory_gb = 8.0
            else:
                memory_gb = 8.0  # Default assumption
        
        # Disk space
        try:
            disk_usage = shutil.disk_usage(Path.cwd())
            disk_free_gb = disk_usage.free / (1024**3)
        except:
            disk_free_gb = 10.0  # Default assumption
        
        # Directories
        home_directory = str(Path.home())
        working_directory = str(Path.cwd())
        
        system_info = SystemInfo(
            platform=platform_name,
            architecture=architecture,
            python_version=python_version,
            cpu_count=cpu_count,
            memory_gb=memory_gb,
            disk_free_gb=disk_free_gb,
            environment_type=self.environment_type,
            home_directory=home_directory,
            working_directory=working_directory
        )
        
        self.system_info = system_info
        return system_info
    
    def detect_dependencies(self) -> Dict[str, bool]:
        """Detect available dependencies"""
        dependencies = {}
        
        # Python packages
        python_deps = [
            "numpy", "scipy", "matplotlib", "pandas", "yaml",
            "mdtraj", "openmm", "rdkit", "openff", "pdbfixer"
        ]
        
        for dep in python_deps:
            try:
                __import__(dep)
                dependencies[dep] = True
            except ImportError:
                dependencies[dep] = False
        
        # System dependencies
        system_deps = {
            "gromacs": ["gmx", "--version"],
            "amber": ["antechamber", "-h"],
            "vmd": ["vmd", "-h"],
            "pymol": ["pymol", "-h"]
        }
        
        for dep, command in system_deps.items():
            try:
                result = subprocess.run(command, capture_output=True, timeout=5)
                dependencies[dep] = result.returncode == 0
            except:
                dependencies[dep] = False
        
        return dependencies
    
    def get_optimal_settings(self, system_info: Optional[SystemInfo] = None) -> Dict[str, Any]:
        """Generate optimal settings based on system capabilities"""
        if not system_info:
            system_info = self.system_info or self.detect_system_info()
        
        settings = {}
        
        # CPU optimization
        max_threads = max(1, system_info.cpu_count - 1)  # Leave one core free
        settings["cpu"] = {
            "max_threads": max_threads,
            "parallel_processing": system_info.cpu_count > 2,
            "thread_pool_size": min(8, max_threads)
        }
        
        # Memory optimization
        if system_info.memory_gb >= 16:
            memory_mode = "high"
            cache_size_mb = 2048
        elif system_info.memory_gb >= 8:
            memory_mode = "medium"
            cache_size_mb = 1024
        else:
            memory_mode = "low"
            cache_size_mb = 512
        
        settings["memory"] = {
            "mode": memory_mode,
            "cache_size_mb": cache_size_mb,
            "enable_swapping": system_info.memory_gb < 8
        }
        
        # Storage optimization
        if system_info.disk_free_gb >= 100:
            storage_mode = "generous"
            temp_dir_size_gb = 10
        elif system_info.disk_free_gb >= 20:
            storage_mode = "moderate"
            temp_dir_size_gb = 5
        else:
            storage_mode = "conservative"
            temp_dir_size_gb = 2
        
        settings["storage"] = {
            "mode": storage_mode,
            "temp_dir_size_gb": temp_dir_size_gb,
            "compression_enabled": system_info.disk_free_gb < 50
        }
        
        # Platform-specific optimizations
        if system_info.platform == "Linux":
            settings["platform"] = {
                "use_native_math": True,
                "enable_numa": system_info.cpu_count > 8,
                "io_scheduler": "mq-deadline"
            }
        elif system_info.platform == "Darwin":  # macOS
            settings["platform"] = {
                "use_accelerate": True,
                "metal_performance": True
            }
        elif system_info.platform == "Windows":
            settings["platform"] = {
                "use_mkl": True,
                "high_dpi_aware": True
            }
        
        return settings
    
    def setup_environment_directories(self) -> Dict[str, Path]:
        """Setup all necessary environment directories"""
        base_dirs = {
            "config": Path.home() / ".prism",
            "data": Path.cwd() / "prism_data",
            "logs": Path.home() / ".prism" / "logs",
            "cache": Path.home() / ".prism" / "cache",
            "temp": Path.home() / ".prism" / "temp",
            "plugins": Path.home() / ".prism" / "plugins",
            "profiles": Path.home() / ".prism" / "profiles",
            "backups": Path.home() / ".prism" / "backups"
        }
        
        # Create additional subdirectories
        subdirs = {
            "pmf_results": base_dirs["data"] / "pmf_results",
            "trajectories": base_dirs["data"] / "trajectories",
            "structures": base_dirs["data"] / "structures",
            "analysis": base_dirs["data"] / "analysis",
            "system_logs": base_dirs["logs"] / "system",
            "calculation_logs": base_dirs["logs"] / "calculations",
            "error_logs": base_dirs["logs"] / "errors",
            "performance_cache": base_dirs["cache"] / "performance",
            "data_cache": base_dirs["cache"] / "data",
            "temp_calculations": base_dirs["temp"] / "calculations",
            "temp_processing": base_dirs["temp"] / "processing"
        }
        
        all_dirs = {**base_dirs, **subdirs}
        
        # Create directories
        created_dirs = {}
        for name, path in all_dirs.items():
            try:
                path.mkdir(parents=True, exist_ok=True)
                created_dirs[name] = path
                self.logger.debug(f"Created/verified directory: {name} -> {path}")
            except Exception as e:
                self.logger.error(f"Failed to create directory {name} ({path}): {e}")
        
        return created_dirs
    
    def generate_environment_config(self) -> EnvironmentConfig:
        """Generate comprehensive environment configuration"""
        system_info = self.detect_system_info()
        dependencies = self.detect_dependencies()
        optimal_settings = self.get_optimal_settings(system_info)
        directories = self.setup_environment_directories()
        
        # Convert Path objects to strings for serialization
        paths = {name: str(path) for name, path in directories.items()}
        
        config = EnvironmentConfig(
            system_info=system_info,
            paths=paths,
            dependencies=dependencies,
            settings=optimal_settings,
            optimizations=self._generate_optimization_config(system_info)
        )
        
        return config
    
    def _generate_optimization_config(self, system_info: SystemInfo) -> Dict[str, Any]:
        """Generate optimization configuration"""
        optimizations = {}
        
        # PMF calculation optimizations
        optimizations["pmf"] = {
            "default_num_windows": 30 if system_info.memory_gb >= 8 else 20,
            "batch_processing": system_info.cpu_count > 4,
            "parallel_windows": min(4, system_info.cpu_count // 2),
            "memory_per_window_mb": max(512, int(system_info.memory_gb * 1024 / 40))
        }
        
        # I/O optimizations
        optimizations["io"] = {
            "buffer_size_mb": 64 if system_info.memory_gb >= 8 else 32,
            "async_io": True,
            "compression_level": 6 if system_info.cpu_count > 4 else 3,
            "parallel_io": system_info.cpu_count > 2
        }
        
        # Analysis optimizations
        optimizations["analysis"] = {
            "chunk_size": 10000 if system_info.memory_gb >= 16 else 5000,
            "use_gpu": False,  # Would need GPU detection
            "precision": "double" if system_info.memory_gb >= 16 else "single"
        }
        
        return optimizations
    
    def save_environment_config(self, config: EnvironmentConfig, 
                               config_level: ConfigurationLevel = ConfigurationLevel.USER) -> bool:
        """Save environment configuration to appropriate file"""
        config_dict = asdict(config)
        
        # Choose config file based on level
        if config_level == ConfigurationLevel.USER:
            config_file = self.config_paths["user_environment"]
        elif config_level == ConfigurationLevel.PROJECT:
            config_file = self.config_paths["project_environment"]
        elif config_level == ConfigurationLevel.RUNTIME:
            config_file = self.config_paths["runtime_config"]
        else:
            config_file = self.config_paths["user_environment"]
        
        try:
            # Ensure parent directory exists
            config_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Save as YAML
            with open(config_file, 'w') as f:
                yaml.dump(config_dict, f, default_flow_style=False, indent=2)
            
            self.logger.info(f"Environment configuration saved to {config_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to save environment config: {e}")
            return False
    
    def load_environment_config(self, config_level: ConfigurationLevel = ConfigurationLevel.USER) -> Optional[EnvironmentConfig]:
        """Load environment configuration from file"""
        if config_level == ConfigurationLevel.USER:
            config_file = self.config_paths["user_environment"]
        elif config_level == ConfigurationLevel.PROJECT:
            config_file = self.config_paths["project_environment"]
        elif config_level == ConfigurationLevel.RUNTIME:
            config_file = self.config_paths["runtime_config"]
        else:
            config_file = self.config_paths["user_environment"]
        
        if not config_file.exists():
            return None
        
        try:
            with open(config_file, 'r') as f:
                config_dict = yaml.safe_load(f)
            
            # Convert back to EnvironmentConfig
            system_info = SystemInfo(**config_dict['system_info'])
            config = EnvironmentConfig(
                system_info=system_info,
                paths=config_dict['paths'],
                dependencies=config_dict['dependencies'],
                settings=config_dict['settings'],
                optimizations=config_dict['optimizations']
            )
            
            return config
            
        except Exception as e:
            self.logger.error(f"Failed to load environment config: {e}")
            return None
    
    def validate_environment(self) -> Dict[str, Any]:
        """Validate current environment setup"""
        validation_results = {
            "overall_status": "unknown",
            "system_requirements": {},
            "dependencies": {},
            "configuration": {},
            "directories": {},
            "recommendations": []
        }
        
        # Check system requirements
        system_info = self.detect_system_info()
        
        validation_results["system_requirements"] = {
            "python_version_ok": sys.version_info >= (3, 7),
            "memory_sufficient": system_info.memory_gb >= 4,
            "disk_space_ok": system_info.disk_free_gb >= 5,
            "cpu_count": system_info.cpu_count
        }
        
        # Check dependencies
        dependencies = self.detect_dependencies()
        required_deps = ["numpy", "scipy", "matplotlib", "pandas", "yaml"]
        
        validation_results["dependencies"] = {
            "required_available": all(dependencies.get(dep, False) for dep in required_deps),
            "optional_count": sum(1 for k, v in dependencies.items() if v and k not in required_deps),
            "total_available": sum(dependencies.values()),
            "details": dependencies
        }
        
        # Check configuration files
        config_status = {}
        for name, path in self.config_paths.items():
            config_status[name] = {
                "exists": path.exists(),
                "readable": path.exists() and os.access(path, os.R_OK),
                "writable": path.exists() and os.access(path, os.W_OK) or 
                          os.access(path.parent, os.W_OK)
            }
        
        validation_results["configuration"] = config_status
        
        # Check directories
        directories = self.setup_environment_directories()
        validation_results["directories"] = {
            "created_count": len(directories),
            "all_accessible": all(d.exists() and os.access(d, os.R_OK | os.W_OK) 
                                for d in directories.values())
        }
        
        # Generate recommendations
        recommendations = []
        
        if not validation_results["system_requirements"]["python_version_ok"]:
            recommendations.append("Upgrade Python to version 3.7 or higher")
        
        if not validation_results["system_requirements"]["memory_sufficient"]:
            recommendations.append("Consider adding more RAM (minimum 4 GB recommended)")
        
        if not validation_results["system_requirements"]["disk_space_ok"]:
            recommendations.append("Free up disk space (minimum 5 GB recommended)")
        
        if not validation_results["dependencies"]["required_available"]:
            missing = [dep for dep in required_deps if not dependencies.get(dep, False)]
            recommendations.append(f"Install missing required dependencies: {', '.join(missing)}")
        
        if validation_results["dependencies"]["optional_count"] < 3:
            recommendations.append("Consider installing optional dependencies for enhanced functionality")
        
        validation_results["recommendations"] = recommendations
        
        # Overall status
        if (validation_results["system_requirements"]["python_version_ok"] and
            validation_results["system_requirements"]["memory_sufficient"] and
            validation_results["dependencies"]["required_available"] and
            validation_results["directories"]["all_accessible"]):
            validation_results["overall_status"] = "good"
        elif len(recommendations) <= 2:
            validation_results["overall_status"] = "fair"
        else:
            validation_results["overall_status"] = "poor"
        
        return validation_results
    
    def setup_development_environment(self) -> bool:
        """Setup complete development environment"""
        self.logger.info("Setting up PRISM development environment...")
        
        try:
            # Generate and save environment configuration
            env_config = self.generate_environment_config()
            
            # Save to both user and project levels
            self.save_environment_config(env_config, ConfigurationLevel.USER)
            self.save_environment_config(env_config, ConfigurationLevel.PROJECT)
            
            # Create development-specific configuration
            dev_config = {
                "environment": "development",
                "debug": True,
                "testing": {
                    "enabled": True,
                    "auto_run": False,
                    "coverage": True
                },
                "logging": {
                    "level": "DEBUG",
                    "console": True,
                    "file": True,
                    "detailed": True
                },
                "optimization": {
                    "cache_enabled": True,
                    "parallel_processing": False,
                    "validation_strict": False
                }
            }
            
            # Save development configuration
            dev_config_file = self.config_paths["user_config"]
            with open(dev_config_file, 'w') as f:
                yaml.dump(dev_config, f, default_flow_style=False, indent=2)
            
            self.logger.info("Development environment setup completed")
            return True
            
        except Exception as e:
            self.logger.error(f"Development environment setup failed: {e}")
            return False
    
    def get_configuration_summary(self) -> Dict[str, Any]:
        """Get comprehensive configuration summary"""
        return {
            "environment_type": self.environment_type.value,
            "system_info": asdict(self.system_info) if self.system_info else None,
            "config_paths": {name: str(path) for name, path in self.config_paths.items()},
            "loaded_configs": list(self.config.keys()),
            "validation": self.validate_environment()
        }


def main():
    """Main function for standalone usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="PRISM Environment Manager")
    parser.add_argument("--setup", action="store_true", help="Setup environment")
    parser.add_argument("--validate", action="store_true", help="Validate environment")
    parser.add_argument("--info", action="store_true", help="Show system information")
    parser.add_argument("--config", action="store_true", help="Show configuration")
    parser.add_argument("--env", choices=["development", "testing", "production", "research"],
                       default="research", help="Environment type")
    
    args = parser.parse_args()
    
    env_type = EnvironmentType(args.env)
    manager = PrismEnvironmentManager(env_type)
    
    if args.setup:
        print("🔧 Setting up PRISM environment...")
        success = manager.setup_development_environment()
        print(f"Setup {'completed' if success else 'failed'}")
    
    if args.validate:
        print("🧪 Validating environment...")
        validation = manager.validate_environment()
        print(f"Overall status: {validation['overall_status']}")
        
        if validation['recommendations']:
            print("\\nRecommendations:")
            for rec in validation['recommendations']:
                print(f"  • {rec}")
    
    if args.info:
        print("🖥️  System Information:")
        system_info = manager.detect_system_info()
        for field, value in asdict(system_info).items():
            print(f"  {field}: {value}")
    
    if args.config:
        print("⚙️  Configuration Summary:")
        summary = manager.get_configuration_summary()
        print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()