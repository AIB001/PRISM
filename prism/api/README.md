# PRISM API Documentation

The PRISM API provides comprehensive programmatic access to PRISM PMF (Potential of Mean Force) calculation functionality through Python library interfaces, internal module communication, and data exchange capabilities.

## Overview

The API system consists of four main components:

1. **Python Library API** - Direct programmatic access through Python classes
2. **Internal Module Communication** - Standardized interfaces for PRISM component integration  
3. **Command-Line API** - Programmatic CLI access for external applications
4. **Data Exchange API** - Format conversion and interoperability with other MD software

## Installation

```python
from prism.api import PMFCalculator, WorkflowManager, DataAnalyzer
from prism.api import PrismAPIClient, FormatConverter, DataExporter
```

## Quick Start

### Basic PMF Calculation

```python
from prism.api import PMFCalculator

# Create calculator
calculator = PMFCalculator()

# Run PMF calculation
result = calculator.run_pmf_calculation(
    system_file="system.gro",
    ligand="LIG",
    protein="Protein",
    method="umbrella_sampling",
    num_windows=30,
    force_constant=1000.0
)

print(f"Binding energy: {result.binding_energy:.2f} kJ/mol")
```

### Using the API Client

```python
from prism.api import create_client

# Create client
with create_client() as client:
    # Run calculation
    result = client.run_pmf_calculation(
        system_file="system.gro",
        ligand="LIG", 
        protein="Protein"
    )
    
    # Export results
    client.export_results(result, "output/", formats=["json", "csv", "xvg"])
```

## Core API Classes

### PMFCalculator

Main class for PMF calculations with support for multiple methods.

```python
from prism.api import PMFCalculator, PMFCalculationConfig

calculator = PMFCalculator()

# Create configuration
config = PMFCalculationConfig(
    system_file="system.gro",
    ligand="LIG",
    protein="Protein", 
    method="umbrella_sampling",
    num_windows=25,
    window_spacing=0.1,
    force_constant=1200.0,
    equilibration_time=2000,  # ps
    sampling_time=15000,      # ps
    temperature=310.15,       # K
    output_directory="pmf_output"
)

# Create calculation
calc_id = calculator.create_calculation(config)

# Submit for execution
result = calculator.submit_calculation(calc_id, blocking=True)

# Check results
if result.status == "completed":
    print(f"PMF calculation completed successfully")
    print(f"Binding energy: {result.binding_energy:.2f} kJ/mol")
    print(f"Convergence: {result.convergence_achieved}")
```

### WorkflowManager  

Orchestrate complex PMF workflows with multiple tasks.

```python
from prism.api import WorkflowManager

manager = WorkflowManager()

# Create workflow
workflow_config = {
    "name": "Protein-Ligand PMF Study",
    "system_file": "complex.gro",
    "ligand": "LIG",
    "protein": "Protein",
    "methods": ["umbrella_sampling"],
    "analysis": ["convergence", "statistics"]
}

workflow_id = manager.create_workflow(workflow_config, priority="high")

# Submit workflow
manager.submit_workflow(workflow_id)

# Monitor progress
status = manager.get_status(workflow_id)
print(f"Workflow {workflow_id}: {status['status']} ({status['progress']:.1f}%)")
```

### DataAnalyzer

Analyze PMF results and generate statistics.

```python
from prism.api import DataAnalyzer

analyzer = DataAnalyzer()

# Analyze PMF profile
pmf_data = [
    {"reaction_coordinate": 0.0, "pmf_energy": 0.0, "error_estimate": 0.3},
    {"reaction_coordinate": 0.1, "pmf_energy": 2.5, "error_estimate": 0.4},
    # ... more data points
]

analysis = analyzer.analyze_pmf_profile(pmf_data)
print(f"Binding energy: {analysis['binding_energy']:.2f} kJ/mol")
print(f"Energy barrier: {analysis['energy_barrier']:.2f} kJ/mol")

# Compare multiple profiles
profiles = [pmf_data_1, pmf_data_2, pmf_data_3]
labels = ["Run 1", "Run 2", "Run 3"]

comparison = analyzer.compare_profiles(profiles, labels)
print(f"Average binding energy: {comparison['statistical_summary']['mean']:.2f}")
```

## Interface APIs

### PMF Interface

Specialized interface for PMF calculations.

```python
from prism.api.interfaces import PMFInterface

# Initialize interface
pmf_interface = PMFInterface()
pmf_interface.initialize({
    'method': 'umbrella_sampling',
    'num_windows': 30,
    'force_constant': 1000.0
})

# Setup umbrella sampling
windows = [i * 0.1 for i in range(30)]  # 0.0 to 2.9 nm
pmf_interface.setup_umbrella_sampling(windows)

# Calculate PMF profile from histogram data  
histogram_data = {
    0.0: [1000, 950, 1050, ...],  # histogram counts
    0.1: [980, 1020, 990, ...],
    # ... more windows
}

profile = pmf_interface.calculate_pmf_profile(histogram_data)
```

### MD Interface

Integration with molecular dynamics simulations.

```python
from prism.api.interfaces import MDInterface, SystemConfig, SimulationParameters

# Initialize MD interface
md_interface = MDInterface()
md_interface.initialize({
    'md_engine': 'gromacs',
    'working_directory': '/path/to/work'
})

# Setup system
system_config = SystemConfig(
    topology_file="topol.top",
    coordinate_file="system.gro", 
    parameter_files=["forcefield.itp"],
    system_name="Protein-Ligand Complex",
    temperature=310.15,
    pressure=1.0
)

md_interface.setup_system(system_config)

# Setup simulation
sim_params = SimulationParameters(
    timestep=0.002,
    total_time=10000,  # ps
    output_frequency=100
)

md_interface.setup_simulation(sim_params)

# Run simulations
md_interface.run_energy_minimization()
md_interface.run_equilibration([
    {"name": "NVT", "time": 1000, "target_temperature": 310.15},
    {"name": "NPT", "time": 1000, "target_pressure": 1.0}
])
md_interface.run_production_simulation("production")
```

## Command-Line API

### Using PrismAPIClient for CLI integration

```python
from prism.api import PrismAPIClient

client = PrismAPIClient()

# Execute CLI commands programmatically
result = client.call_prism_cli([
    'workflow', 'create', 'system.gro',
    '--ligand', 'LIG',
    '--protein', 'Protein',
    '--method', 'umbrella_sampling'
])

if result['success']:
    print("Workflow created successfully")
    print(result['stdout'])

# Get system status via CLI
status = client.get_status_cli()
print(f"System status: {status}")

# Submit workflow via CLI  
workflow_result = client.submit_workflow_cli("workflow_config.yaml", async_mode=True)
```

### Batch Operations

```python
# Run multiple PMF calculations
calculation_configs = [
    {
        "system_file": "system1.gro",
        "ligand": "LIG1", 
        "protein": "Protein",
        "method": "umbrella_sampling"
    },
    {
        "system_file": "system2.gro",
        "ligand": "LIG2",
        "protein": "Protein", 
        "method": "smd"
    },
    # ... more calculations
]

results = client.batch_pmf_calculations(calculation_configs, max_concurrent=3)

for i, result in enumerate(results):
    if result.get('status') == 'completed':
        print(f"Calculation {i+1}: Success - Binding energy = {result.get('binding_energy'):.2f} kJ/mol")
    else:
        print(f"Calculation {i+1}: Failed - {result.get('error')}")
```

## Data Exchange API

### Format Conversion

```python
from prism.api.data import FormatConverter, DataFormat, PMFProfileData

converter = FormatConverter()

# Create PMF data
pmf_data = PMFProfileData(
    reaction_coordinate=[0.0, 0.1, 0.2, 0.3],
    pmf_energy=[0.0, 2.5, 8.1, 15.2],
    error_estimate=[0.3, 0.4, 0.6, 0.8]
)

# Convert to different formats
gromacs_xvg = converter.convert_pmf_profile(pmf_data, DataFormat.GROMACS, "pmf.xvg")
amber_dat = converter.convert_pmf_profile(pmf_data, DataFormat.AMBER, "pmf.dat") 
json_data = converter.convert_pmf_profile(pmf_data, DataFormat.JSON)

print("Converted PMF to multiple formats")
```

### Data Export

```python
from prism.api.data import DataExporter

exporter = DataExporter()

# Export PMF profile
exporter.export_pmf_profile(
    pmf_data=pmf_results,
    output_file="pmf_profile.xvg",
    format=DataFormat.XVG,
    metadata={"method": "umbrella_sampling", "temperature": 310.15}
)

# Export complete calculation results
calculation_results = {
    "pmf_profile": pmf_data,
    "energy_data": energy_data,
    "metadata": {"calculation_id": "pmf_001", "runtime": 7200}
}

exported_files = exporter.export_calculation_results(
    results=calculation_results,
    output_directory="results/",
    formats=[DataFormat.XVG, DataFormat.JSON, DataFormat.CSV]
)

print(f"Exported files: {exported_files}")
```

### Data Import

```python
from prism.api.data import DataImporter

importer = DataImporter()

# Import PMF data from external files
pmf_data = importer.import_pmf_profile("external_pmf.xvg", DataFormat.XVG)
energy_data = importer.import_energy_data("energy.csv", DataFormat.CSV)

# Use imported data
analyzer = DataAnalyzer()
analysis = analyzer.analyze_pmf_profile(pmf_data)
```

## Error Handling

The API uses custom exception classes for different error types:

```python
from prism.api.exceptions import (
    APIError, ValidationError, CalculationError, 
    WorkflowError, DataError, TimeoutError
)

try:
    result = calculator.run_pmf_calculation(
        system_file="missing_file.gro",
        ligand="LIG",
        protein="Protein"
    )
except ValidationError as e:
    print(f"Validation error: {e.message}")
    print(f"Field: {e.field}, Value: {e.value}")
except CalculationError as e:
    print(f"Calculation failed: {e.message}")
    print(f"Stage: {e.stage}, Exit code: {e.exit_code}")
except TimeoutError as e:
    print(f"Operation timed out: {e.message}")
    print(f"Timeout: {e.timeout_seconds}s")
```

## Configuration

### API Client Configuration

```python
config = {
    "timeout": 7200,  # 2 hours
    "working_directory": "/path/to/workdir",
    "temp_directory": "/tmp/prism",
    "max_concurrent": 5,
    "log_level": "INFO"
}

client = PrismAPIClient(config=config)
```

### Logging Configuration

```python
from prism.utils.logging_system import PrismLogger

# Create custom logger
logger = PrismLogger("my_pmf_script", level="DEBUG")

# Use with API components
calculator = PMFCalculator(logger=logger)
```

## Advanced Usage

### Custom Callback Functions

```python
def progress_callback(result):
    print(f"Calculation {result.calculation_id}: {result.status}")

def completion_callback(result):
    if result.status == "completed":
        print(f"✅ Calculation completed: {result.binding_energy:.2f} kJ/mol")
    else:
        print(f"❌ Calculation failed: {result.error_message}")

# Add callbacks
calc_id = calculator.create_calculation(config)
calculator.add_callback(calc_id, completion_callback)

# Submit calculation
calculator.submit_calculation(calc_id, blocking=False)
```

### Context Manager Usage

```python
# Automatic resource cleanup
with create_client() as client:
    # All operations within this block
    result = client.run_pmf_calculation(...)
    client.export_results(result, "output/")
# Resources automatically cleaned up here
```

### Custom Data Converters

```python
def my_custom_converter(data):
    # Custom conversion logic
    converted = f"# Custom Format\n"
    for i, (x, energy) in enumerate(zip(data.reaction_coordinate, data.pmf_energy)):
        converted += f"{x:.3f}\t{energy:.6f}\n"
    return converted

# Register custom converter
converter = FormatConverter()
converter.register_converter("pmf_to_custom", my_custom_converter)

# Use custom converter
result = converter._registered_converters["pmf_to_custom"](pmf_data)
```

## Best Practices

1. **Use context managers** for automatic resource cleanup
2. **Handle exceptions** appropriately for robust error handling  
3. **Set appropriate timeouts** for long-running calculations
4. **Use batch operations** for multiple similar calculations
5. **Export results** in multiple formats for interoperability
6. **Add callbacks** for monitoring long-running processes
7. **Configure logging** for debugging and monitoring

## Examples Directory

See the `examples/api/` directory for complete working examples:

- `basic_pmf_calculation.py` - Simple PMF calculation
- `workflow_management.py` - Complex workflow orchestration
- `data_analysis.py` - PMF data analysis and statistics
- `format_conversion.py` - Data format conversion examples
- `batch_processing.py` - Batch calculation processing
- `error_handling.py` - Comprehensive error handling
- `cli_integration.py` - CLI integration examples

## API Reference

For detailed API reference documentation, see the individual module docstrings:

- `prism.api.core` - Core calculation classes
- `prism.api.interfaces` - Module communication interfaces  
- `prism.api.client` - API client for external access
- `prism.api.data` - Data handling and format conversion
- `prism.api.exceptions` - Exception classes