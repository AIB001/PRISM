# PRISM-FEP Analysis Framework

Complete analysis toolkit for Free Energy Perturbation (FEP) calculations performed with GROMACS.

## Features

- **Multiple Estimators**: Support for BAR, MBAR, and TI methods
- **Comprehensive Reports**: Publication-ready HTML reports with interactive plots
- **Convergence Analysis**: Time-convergence and hysteresis diagnostics
- **Energy Decomposition**: Separate electrostatic and van der Waals contributions
- **CLI Interface**: Easy-to-use command-line tool
- **Python API**: Programmatic access for custom workflows

## Installation

```bash
# Install PRISM with FEP analysis dependencies
pip install -e .[fep_analysis]

# Or install dependencies manually
pip install alchemlyb pandas numpy plotly
```

## Quick Start

### Command Line Interface

```bash
# Basic analysis
prism --fep-analyze \
    --bound-dir fep_project/bound \
    --unbound-dir fep_project/unbound \
    --output report.html

# With custom estimator and temperature
prism --fep-analyze \
    --bound-dir fep_project/bound \
    --unbound-dir fep_project/unbound \
    --estimator MBAR \
    --temperature 300 \
    --output report.html

# Save both HTML and JSON results
prism --fep-analyze \
    --bound-dir fep_project/bound \
    --unbound-dir fep_project/unbound \
    --output report.html \
    --json results.json
```

### Python API

```python
from prism.fep.analysis import FEPAnalyzer

# Create analyzer
analyzer = FEPAnalyzer(
    bound_dir='fep_project/bound',
    unbound_dir='fep_project/unbound',
    temperature=310.0,
    estimator='MBAR'
)

# Run analysis
results = analyzer.analyze()

# Access results
print(f"ΔG = {results.delta_g:.2f} ± {results.delta_g_error:.2f} kcal/mol")
print(f"Bound leg:   {results.delta_g_bound:.2f} kcal/mol")
print(f"Unbound leg: {results.delta_g_unbound:.2f} kcal/mol")

# Generate HTML report
report_path = analyzer.generate_html_report('report.html')

# Save JSON results
analyzer.save_results('results.json')
```

## Directory Structure

The analysis expects the following directory structure:

```
fep_project/
├── bound/
│   ├── window_00/
│   │   ├── dhdl.xvg          # dH/dλ data
│   │   ├── prod.xvg          # Production run data
│   │   └── ...
│   ├── window_01/
│   └── ...
└── unbound/
    ├── window_00/
    │   ├── dhdl.xvg
    │   └── ...
    └── ...
```

## Output Files

### HTML Report (`report.html`)

Comprehensive HTML report containing:

- **Summary**: Quick overview of binding free energy and key parameters
- **Results Table**: Detailed breakdown of all computed quantities
- **Free Energy Profile**: Interactive plot of ΔG vs λ
- **Convergence Analysis**: Time-convergence diagnostics
- **Energy Decomposition**: Contributions from different energy components
- **Methods**: Description of estimators and simulation parameters

### JSON Results (`results.json`)

Machine-readable JSON file with all analysis results:

```json
{
  "delta_g": -2.34,
  "delta_g_error": 0.15,
  "delta_g_bound": 5.67,
  "delta_g_unbound": 3.33,
  "delta_g_components": {
    "elec": -1.23,
    "vdw": -1.11
  },
  "convergence": {
    "n_samples": 100000,
    "time_convergence": "converged"
  },
  "metadata": {
    "temperature": 310.0,
    "estimator": "MBAR",
    "n_lambda_windows": 21
  }
}
```

## Estimators

### MBAR (Multistate Bennett Acceptance Ratio)
- **Recommended**: Best statistical efficiency
- **Uses**: All λ windows simultaneously
- **Pros**: Optimal use of data, robust
- **Cons**: Higher memory usage for many windows

### BAR (Bennett Acceptance Ratio)
- **Uses**: Adjacent λ window pairs
- **Pros**: Fast, low memory
- **Cons**: Lower statistical efficiency than MBAR

### TI (Thermodynamic Integration)
- **Uses**: Numerical integration of dH/dλ
- **Pros**: Simple, intuitive
- **Cons**: Requires careful λ spacing, less robust

## Advanced Usage

### Custom Energy Components

```python
analyzer = FEPAnalyzer(
    bound_dir='bound',
    unbound_dir='unbound',
    energy_components=['elec', 'vdw', 'restraint']
)
```

### Convergence Analysis

```python
results = analyzer.analyze()

# Access convergence metrics
print(f"Total samples: {results.convergence['n_samples']}")
print(f"Time convergence: {results.convergence['time_convergence']}")
```

### Custom Processing

```python
# Parse data without running full analysis
from prism.fep.analysis.parsers import GROMACSParser

parser = GROMACSParser()
bound_data = parser.parse_directory('bound')
unbound_data = parser.parse_directory('unbound')

# Custom analysis...
```

## Troubleshooting

### Missing dhdl.xvg files

**Problem**: `FileNotFoundError: No valid dhdl.xvg files found`

**Solution**: Ensure your GROMACS MDP files include:
```
free-energy = yes
calc-lambda-neighbors = -1   # Output dH/dλ for all windows
```

### Alchemlyb not installed

**Problem**: `ImportError: alchemlyb is required for FEP analysis`

**Solution**: Install alchemlyb:
```bash
pip install alchemlyb
```

### Insufficient convergence

**Problem**: Large error estimates or hysteresis

**Solutions**:
- Increase simulation time per λ window
- Use more λ windows (especially in regions of rapid change)
- Check for adequate sampling in each window
- Consider soft-core potentials for large perturbations

## Performance Tips

1. **Use MBAR** for better statistical efficiency
2. **Subsample data** if you have more than 10,000 samples per window
3. **Parallelize** analysis of bound and unbound legs
4. **Cache parsed data** if running multiple analyses

## Comparison with NAMD FEBuilder

This framework is the GROMACS equivalent of the NAMD FEBuilder analysis module:

| Feature | NAMD FEBuilder | PRISM-FEP |
|---------|----------------|-----------|
| Input files | fepout files | dhdl.xvg files |
| Estimators | BAR, MBAR, TI | BAR, MBAR, TI |
| Convergence analysis | ✓ | ✓ |
| Bootstrap error analysis | ✓ | Planned |
| HTML reports | ✓ | ✓ |
| Energy decomposition | ✓ | Partial |
| Platform | NAMD | GROMACS |

## Citation

If you use PRISM-FEP in your research, please cite:

```
PRISM: Protein Receptor Interaction Simulation Modeler
https://github.com/your-org/prism
```

And also cite alchemlyb:

```
alchemlyb: A Python library for free energy calculations
https://github.com/alchemistry/alchemlyb
```

## License

MIT License - see LICENSE file for details

## Contact

For issues, questions, or contributions:
- GitHub: https://github.com/your-org/prism
- Email: your-email@institution.edu

## Version History

### 0.1.0 (2025-03-19)
- Initial release
- Support for BAR, MBAR, TI estimators
- HTML report generation
- CLI interface
- Basic convergence analysis
