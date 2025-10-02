#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Plugin System and Extension Framework

Provides a flexible plugin architecture for extending PRISM functionality
with custom modules, analysis tools, and computational methods.
"""

import os
import sys
import importlib
import importlib.util
import inspect
from pathlib import Path
from typing import Dict, Any, List, Optional, Type, Callable, Union
from dataclasses import dataclass
from abc import ABC, abstractmethod
from enum import Enum
import json
import yaml

from .logging_system import PrismLogger


class PluginType(Enum):
    """Types of plugins supported by PRISM"""
    FORCE_FIELD = "force_field"
    ANALYSIS = "analysis"
    VISUALIZATION = "visualization"
    CALCULATION = "calculation"
    IO_HANDLER = "io_handler"
    VALIDATOR = "validator"
    OPTIMIZER = "optimizer"
    WORKFLOW = "workflow"


class PluginStatus(Enum):
    """Plugin status"""
    LOADED = "loaded"
    ACTIVE = "active"
    INACTIVE = "inactive"
    ERROR = "error"
    INCOMPATIBLE = "incompatible"


@dataclass
class PluginMetadata:
    """Metadata for a plugin"""
    name: str
    version: str
    description: str
    author: str
    plugin_type: PluginType
    requires_prism_version: str
    dependencies: List[str]
    api_version: str
    license: Optional[str] = None
    homepage: Optional[str] = None
    keywords: Optional[List[str]] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary"""
        return {
            'name': self.name,
            'version': self.version,
            'description': self.description,
            'author': self.author,
            'plugin_type': self.plugin_type.value,
            'requires_prism_version': self.requires_prism_version,
            'dependencies': self.dependencies,
            'api_version': self.api_version,
            'license': self.license,
            'homepage': self.homepage,
            'keywords': self.keywords or []
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PluginMetadata':
        """Create from dictionary"""
        return cls(
            name=data['name'],
            version=data['version'],
            description=data['description'],
            author=data['author'],
            plugin_type=PluginType(data['plugin_type']),
            requires_prism_version=data['requires_prism_version'],
            dependencies=data.get('dependencies', []),
            api_version=data['api_version'],
            license=data.get('license'),
            homepage=data.get('homepage'),
            keywords=data.get('keywords', [])
        )


class PluginInterface(ABC):
    """Base interface for all PRISM plugins"""
    
    @abstractmethod
    def get_metadata(self) -> PluginMetadata:
        """Get plugin metadata"""
        pass
    
    @abstractmethod
    def initialize(self, config: Optional[Dict[str, Any]] = None) -> bool:
        """Initialize the plugin"""
        pass
    
    @abstractmethod
    def shutdown(self) -> None:
        """Shutdown the plugin"""
        pass
    
    def get_capabilities(self) -> List[str]:
        """Get list of capabilities provided by this plugin"""
        return []
    
    def validate_environment(self) -> Dict[str, bool]:
        """Validate plugin environment and dependencies"""
        return {'valid': True}


class ForceFieldPlugin(PluginInterface):
    """Interface for force field plugins"""
    
    @abstractmethod
    def generate_parameters(self, molecule_data: Dict[str, Any]) -> Dict[str, Any]:
        """Generate force field parameters for a molecule"""
        pass
    
    @abstractmethod
    def supported_molecules(self) -> List[str]:
        """Get list of supported molecule types"""
        pass


class AnalysisPlugin(PluginInterface):
    """Interface for analysis plugins"""
    
    @abstractmethod
    def analyze(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Perform analysis on input data"""
        pass
    
    @abstractmethod
    def get_analysis_types(self) -> List[str]:
        """Get list of supported analysis types"""
        pass


class VisualizationPlugin(PluginInterface):
    """Interface for visualization plugins"""
    
    @abstractmethod
    def create_visualization(self, data: Dict[str, Any], 
                           output_path: Path) -> Path:
        """Create visualization from data"""
        pass
    
    @abstractmethod
    def supported_formats(self) -> List[str]:
        """Get supported output formats"""
        pass


class CalculationPlugin(PluginInterface):
    """Interface for calculation plugins"""
    
    @abstractmethod
    def calculate(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Perform calculation"""
        pass
    
    @abstractmethod
    def get_calculation_methods(self) -> List[str]:
        """Get available calculation methods"""
        pass


class PluginRegistry:
    """Registry for managing plugins"""
    
    def __init__(self):
        self.plugins: Dict[str, PluginInterface] = {}
        self.plugin_metadata: Dict[str, PluginMetadata] = {}
        self.plugin_status: Dict[str, PluginStatus] = {}
        self.plugin_paths: Dict[str, Path] = {}
        self.type_registry: Dict[PluginType, List[str]] = {ptype: [] for ptype in PluginType}
        self.logger = PrismLogger("prism.plugin_registry")
    
    def register_plugin(self, plugin: PluginInterface, 
                       plugin_path: Optional[Path] = None) -> bool:
        """Register a plugin instance"""
        try:
            metadata = plugin.get_metadata()
            
            # Validate metadata
            if not self._validate_metadata(metadata):
                self.logger.error(f"Invalid metadata for plugin {metadata.name}")
                return False
            
            # Check for name conflicts
            if metadata.name in self.plugins:
                self.logger.warning(f"Plugin {metadata.name} already registered, replacing")
            
            # Initialize plugin
            if not plugin.initialize():
                self.logger.error(f"Failed to initialize plugin {metadata.name}")
                self.plugin_status[metadata.name] = PluginStatus.ERROR
                return False
            
            # Register plugin
            self.plugins[metadata.name] = plugin
            self.plugin_metadata[metadata.name] = metadata
            self.plugin_status[metadata.name] = PluginStatus.LOADED
            if plugin_path:
                self.plugin_paths[metadata.name] = plugin_path
            
            # Add to type registry
            if metadata.plugin_type not in self.type_registry:
                self.type_registry[metadata.plugin_type] = []
            self.type_registry[metadata.plugin_type].append(metadata.name)
            
            self.logger.info(f"Registered plugin: {metadata.name} v{metadata.version}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error registering plugin: {e}")
            return False
    
    def unregister_plugin(self, plugin_name: str) -> bool:
        """Unregister a plugin"""
        if plugin_name not in self.plugins:
            self.logger.warning(f"Plugin {plugin_name} not found")
            return False
        
        try:
            # Shutdown plugin
            plugin = self.plugins[plugin_name]
            plugin.shutdown()
            
            # Remove from registries
            metadata = self.plugin_metadata[plugin_name]
            self.type_registry[metadata.plugin_type].remove(plugin_name)
            
            del self.plugins[plugin_name]
            del self.plugin_metadata[plugin_name]
            del self.plugin_status[plugin_name]
            if plugin_name in self.plugin_paths:
                del self.plugin_paths[plugin_name]
            
            self.logger.info(f"Unregistered plugin: {plugin_name}")
            return True
            
        except Exception as e:
            self.logger.error(f"Error unregistering plugin {plugin_name}: {e}")
            return False
    
    def get_plugin(self, plugin_name: str) -> Optional[PluginInterface]:
        """Get plugin by name"""
        return self.plugins.get(plugin_name)
    
    def get_plugins_by_type(self, plugin_type: PluginType) -> List[PluginInterface]:
        """Get all plugins of a specific type"""
        plugin_names = self.type_registry.get(plugin_type, [])
        return [self.plugins[name] for name in plugin_names if name in self.plugins]
    
    def activate_plugin(self, plugin_name: str, config: Optional[Dict[str, Any]] = None) -> bool:
        """Activate a plugin"""
        if plugin_name not in self.plugins:
            self.logger.error(f"Plugin {plugin_name} not found")
            return False
        
        try:
            plugin = self.plugins[plugin_name]
            if plugin.initialize(config):
                self.plugin_status[plugin_name] = PluginStatus.ACTIVE
                self.logger.info(f"Activated plugin: {plugin_name}")
                return True
            else:
                self.plugin_status[plugin_name] = PluginStatus.ERROR
                return False
        except Exception as e:
            self.logger.error(f"Error activating plugin {plugin_name}: {e}")
            self.plugin_status[plugin_name] = PluginStatus.ERROR
            return False
    
    def deactivate_plugin(self, plugin_name: str) -> bool:
        """Deactivate a plugin"""
        if plugin_name not in self.plugins:
            self.logger.error(f"Plugin {plugin_name} not found")
            return False
        
        try:
            plugin = self.plugins[plugin_name]
            plugin.shutdown()
            self.plugin_status[plugin_name] = PluginStatus.INACTIVE
            self.logger.info(f"Deactivated plugin: {plugin_name}")
            return True
        except Exception as e:
            self.logger.error(f"Error deactivating plugin {plugin_name}: {e}")
            return False
    
    def _validate_metadata(self, metadata: PluginMetadata) -> bool:
        """Validate plugin metadata"""
        required_fields = ['name', 'version', 'description', 'author', 'api_version']
        for field in required_fields:
            if not getattr(metadata, field):
                return False
        return True
    
    def get_registry_info(self) -> Dict[str, Any]:
        """Get complete registry information"""
        return {
            'total_plugins': len(self.plugins),
            'plugins_by_type': {ptype.value: len(plugins) 
                               for ptype, plugins in self.type_registry.items()},
            'plugin_status': {name: status.value 
                             for name, status in self.plugin_status.items()},
            'plugin_metadata': {name: metadata.to_dict() 
                               for name, metadata in self.plugin_metadata.items()}
        }


class PluginLoader:
    """Loads plugins from various sources"""
    
    def __init__(self, registry: PluginRegistry):
        self.registry = registry
        self.logger = PrismLogger("prism.plugin_loader")
    
    def load_from_directory(self, plugin_dir: Path, 
                           recursive: bool = True) -> List[str]:
        """Load plugins from directory"""
        loaded_plugins = []
        
        if not plugin_dir.exists():
            self.logger.warning(f"Plugin directory not found: {plugin_dir}")
            return loaded_plugins
        
        # Find plugin files
        patterns = ["*.py", "plugin.yaml", "plugin.json"]
        plugin_files = []
        
        for pattern in patterns:
            if recursive:
                plugin_files.extend(plugin_dir.rglob(pattern))
            else:
                plugin_files.extend(plugin_dir.glob(pattern))
        
        # Load each plugin
        for plugin_file in plugin_files:
            plugin_name = self._load_plugin_file(plugin_file)
            if plugin_name:
                loaded_plugins.append(plugin_name)
        
        self.logger.info(f"Loaded {len(loaded_plugins)} plugins from {plugin_dir}")
        return loaded_plugins
    
    def load_from_file(self, plugin_file: Path) -> Optional[str]:
        """Load plugin from single file"""
        return self._load_plugin_file(plugin_file)
    
    def _load_plugin_file(self, plugin_file: Path) -> Optional[str]:
        """Load plugin from file"""
        try:
            if plugin_file.suffix == '.py':
                return self._load_python_plugin(plugin_file)
            elif plugin_file.name in ['plugin.yaml', 'plugin.yml']:
                return self._load_yaml_plugin(plugin_file)
            elif plugin_file.name == 'plugin.json':
                return self._load_json_plugin(plugin_file)
            else:
                self.logger.debug(f"Skipping non-plugin file: {plugin_file}")
                return None
        except Exception as e:
            self.logger.error(f"Error loading plugin from {plugin_file}: {e}")
            return None
    
    def _load_python_plugin(self, plugin_file: Path) -> Optional[str]:
        """Load Python plugin"""
        try:
            # Load module
            spec = importlib.util.spec_from_file_location("plugin_module", plugin_file)
            if spec is None or spec.loader is None:
                return None
            
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            # Find plugin classes
            plugin_classes = []
            for name, obj in inspect.getmembers(module):
                if (inspect.isclass(obj) and 
                    issubclass(obj, PluginInterface) and 
                    obj != PluginInterface):
                    plugin_classes.append(obj)
            
            if not plugin_classes:
                self.logger.warning(f"No plugin classes found in {plugin_file}")
                return None
            
            # Instantiate and register first plugin class
            plugin_class = plugin_classes[0]
            plugin_instance = plugin_class()
            
            if self.registry.register_plugin(plugin_instance, plugin_file):
                return plugin_instance.get_metadata().name
            
        except Exception as e:
            self.logger.error(f"Error loading Python plugin {plugin_file}: {e}")
        
        return None
    
    def _load_yaml_plugin(self, plugin_file: Path) -> Optional[str]:
        """Load plugin defined in YAML"""
        try:
            with open(plugin_file, 'r') as f:
                plugin_config = yaml.safe_load(f)
            
            return self._load_config_plugin(plugin_config, plugin_file)
        except Exception as e:
            self.logger.error(f"Error loading YAML plugin {plugin_file}: {e}")
            return None
    
    def _load_json_plugin(self, plugin_file: Path) -> Optional[str]:
        """Load plugin defined in JSON"""
        try:
            with open(plugin_file, 'r') as f:
                plugin_config = json.load(f)
            
            return self._load_config_plugin(plugin_config, plugin_file)
        except Exception as e:
            self.logger.error(f"Error loading JSON plugin {plugin_file}: {e}")
            return None
    
    def _load_config_plugin(self, config: Dict[str, Any], 
                           plugin_file: Path) -> Optional[str]:
        """Load plugin from configuration dict"""
        # This would implement loading plugins defined in configuration files
        # For now, just log that we found a config-based plugin
        self.logger.info(f"Found config-based plugin definition: {plugin_file}")
        return None


class PluginManager:
    """High-level plugin management interface"""
    
    def __init__(self, plugin_dirs: Optional[List[Path]] = None):
        self.registry = PluginRegistry()
        self.loader = PluginLoader(self.registry)
        self.logger = PrismLogger("prism.plugin_manager")
        
        # Load plugins from specified directories
        if plugin_dirs:
            for plugin_dir in plugin_dirs:
                self.load_plugins_from_directory(plugin_dir)
    
    def load_plugins_from_directory(self, plugin_dir: Union[str, Path]) -> List[str]:
        """Load all plugins from directory"""
        plugin_dir = Path(plugin_dir)
        return self.loader.load_from_directory(plugin_dir)
    
    def load_plugin_from_file(self, plugin_file: Union[str, Path]) -> Optional[str]:
        """Load single plugin from file"""
        plugin_file = Path(plugin_file)
        return self.loader.load_from_file(plugin_file)
    
    def get_available_plugins(self, plugin_type: Optional[PluginType] = None) -> List[Dict[str, Any]]:
        """Get list of available plugins"""
        if plugin_type:
            plugin_names = self.registry.type_registry.get(plugin_type, [])
        else:
            plugin_names = list(self.registry.plugins.keys())
        
        return [self.registry.plugin_metadata[name].to_dict() 
                for name in plugin_names if name in self.registry.plugin_metadata]
    
    def activate_plugin(self, plugin_name: str, 
                       config: Optional[Dict[str, Any]] = None) -> bool:
        """Activate a plugin"""
        return self.registry.activate_plugin(plugin_name, config)
    
    def deactivate_plugin(self, plugin_name: str) -> bool:
        """Deactivate a plugin"""
        return self.registry.deactivate_plugin(plugin_name)
    
    def get_plugin(self, plugin_name: str) -> Optional[PluginInterface]:
        """Get plugin instance"""
        return self.registry.get_plugin(plugin_name)
    
    def get_plugins_by_type(self, plugin_type: PluginType) -> List[PluginInterface]:
        """Get plugins by type"""
        return self.registry.get_plugins_by_type(plugin_type)
    
    def execute_plugin_method(self, plugin_name: str, method_name: str, 
                             *args, **kwargs) -> Any:
        """Execute method on a plugin"""
        plugin = self.registry.get_plugin(plugin_name)
        if not plugin:
            raise ValueError(f"Plugin {plugin_name} not found")
        
        if not hasattr(plugin, method_name):
            raise ValueError(f"Method {method_name} not found in plugin {plugin_name}")
        
        method = getattr(plugin, method_name)
        return method(*args, **kwargs)
    
    def get_manager_info(self) -> Dict[str, Any]:
        """Get plugin manager information"""
        return {
            'registry_info': self.registry.get_registry_info(),
            'loaded_plugins': list(self.registry.plugins.keys()),
            'active_plugins': [name for name, status in self.registry.plugin_status.items() 
                              if status == PluginStatus.ACTIVE],
            'plugin_directories': []  # Could track loaded directories
        }


# Convenience functions
def create_plugin_manager(plugin_dirs: Optional[List[Union[str, Path]]] = None) -> PluginManager:
    """Create plugin manager with standard plugin directories"""
    if plugin_dirs is None:
        # Default plugin directories
        plugin_dirs = [
            Path.home() / ".prism" / "plugins",
            Path(__file__).parent.parent / "plugins",
            Path.cwd() / "plugins"
        ]
    
    plugin_dirs = [Path(d) for d in plugin_dirs]
    return PluginManager(plugin_dirs)


def discover_plugins(directory: Union[str, Path]) -> List[Dict[str, Any]]:
    """Discover plugins in a directory without loading them"""
    manager = PluginManager()
    loaded = manager.load_plugins_from_directory(directory)
    return manager.get_available_plugins()