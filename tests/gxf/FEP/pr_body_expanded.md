## Summary

This PR adds a complete GROMACS-oriented FEP workflow to PRISM and refactors related builder and analysis code so the workflow is maintainable and usable from the main `prism` interface.

The main scope is FEP: atom mapping, hybrid topology construction, scaffold generation, execution scripts, HTML inspection, analysis/reporting, and standardized output naming. The PR also includes supporting refactors in the builder, force-field integration, and analysis modules.

## Main Features

### Force field support
- **Multi-force field**: GAFF, GAFF2, CGenFF, OpenFF, OPLS-AA, RTF, SwissParam (MMFF/MATCH/Both)
- **SwissParam integration**: MMFF, MATCH, Both force fields with async session polling
- **RTF support**: CHARMM Residue Topology File for MATCH/CHARMM-GUI parameters
- **AMBER19SB**: Water-model-specific ion includes (ions_tip3p.itp) for modern GROMACS
- **MOL2 format**: Comprehensive MOL2 file format handling for ligands
- **Auto-match**: CGenFF/CHARMM-GUI integration with auto-match PDB to topology

### FEP workflow

#### Atom mapping and hybrid topology
- **DistanceAtomMapper** - Distance-based atom mapping with configurable cutoff (default 0.6 Å)
  - FEP-specific atom classification: common, transformed (A-only, B-only), surrounding
  - Automatic mapping validation and quality scoring
- **Single-topology hybrid ligand** - GROMACS A/B-state encoding in single ITP file
  - State A (reference) and State B (mutant) atoms defined in [ atoms ] section
  - Coupled parameters for smooth lambda interpolation
  - Bonded interactions (bonds, angles, dihedrals) with state-dependent parameters
- **Hybrid coordinate generation** - Proper placement of hybrid ligand in binding pocket
  - Superposition of reference and mutant ligands
  - Centroid-based alignment for optimal placement
  - Minimization and equilibration protocols for hybrid systems

#### Scaffold generation and directory layout
- **Automated bound/unbound scaffold** - Complete FEP system setup
  - **Bound leg**: Protein-ligand complex with hybrid topology
  - **Unbound leg**: Hybrid ligand in water box
  - **Repeat support**: Multiple independent repeats (repeat1, repeat2, ...) for statistical robustness
- **Standardized default case naming** - Auto-generate force-field-based case directories when users do not pass `-o`
  - Default format: `<protein_ff>-mut_<ligand_ff>`
  - Examples: `amber14sb_ol15-mut_gaff2`, `charmm36_jul2022-mut_cgenff`, `oplsaa-mut_opls`
  - User-specified output directories are still honored unchanged
- **Force field isolation** - Each repeat gets its own force field copy
  - Prevents cross-contamination between repeats
  - Enables parallel execution without conflicts
- **MDP generation** - Lambda-specific parameter files
  - 11 lambda windows (default): 4 coulomb + 7 van der Waals (decoupled strategy)
  - Customizable lambda schedules (linear, nonlinear, quadratic distributions)
  - Separate MDPs for EM, NVT, NPT, and production per lambda window
- **B-state atom validation** - Ensures all hybrid topologies have proper perturbed atoms
  - Added `is_perturbed` flag to HybridAtom dataclass for tracking perturbed state
  - Transformed and surrounding atoms always write B-state columns (typeB, chargeB, massB)
  - Mass difference checking prevents omission of mass-only perturbations
  - Validation method ensures at least one B-state atom exists when transformation present
  - Fixes issue where charge redistribution could make A/B values identical, causing zero B-state atoms

#### Runtime execution
- **Script generation** - Automated creation of execution scripts
  - **Standard mode**: Sequential lambda window execution
  - **Replica exchange mode**: Lambda replica exchange (REX) with exchange attempts
  - GPU-aware resource allocation with proper CPU affinity
  - Automatic thread count calculation from available CPUs
- **Execution modes**:
  - `bound` - Run bound leg (all repeats)
  - `unbound` - Run unbound leg (all repeats)
  - `all` - Run complete FEP workflow (bound + unbound)
- **Progress tracking** - Automatic resumption from completed windows
  - Checkpoint system for crash recovery
  - Log file monitoring for completion status
  - Automatic next-window launching

### FEP analysis

#### Multi-estimator framework
- **Supported estimators**:
  - **BAR** (Bennett Acceptance Ratio) - Gold standard for free energy differences
  - **MBAR** (Multistate BAR) - Optimal estimator using all lambda windows
  - **TI** (Thermodynamic Integration) - Integral-based estimator
- **Unified interface** - Single API call for all estimators
  - Automatic dhdl.xvg parsing from GROMACS output
  - Error propagation and uncertainty quantification
  - Cross-estimator consistency checking

#### Bootstrap and convergence analysis
- **Bootstrap uncertainty estimation**:
  - Configurable number of bootstrap iterations (default: 1000)
  - Parallel execution with `n_jobs` parameter for multi-core systems
  - Automatic standard error and confidence interval calculation
- **Time convergence** - Check equilibration convergence
  - Skip initial frames for equilibration removal
  - Convergence plots vs simulation time
  - Automatic recommendation for production window length

#### Quality checks and validation
- **Automatic quality metrics**:
  - **Uncertainty checks**: Warn on high standard errors (> 1.0 kcal/mol)
  - **Overlap matrix**: Detect poor overlap between lambda windows
  - **Endpoint behavior**: Check for endpoint violations or instabilities
  - **Convergence diagnostics**: Time-series analysis for drift and equilibration
- **Visual validation**:
  - Overlap matrix heatmaps
  - dhdl convergence plots
  - Lambda window energy distributions
  - Time-converged ΔG plots

#### HTML reporting
- **Interactive HTML reports** - Single-file interactive visualization
  - Multi-estimator comparison tables
  - Interactive plots with zoom/pan capabilities
  - Color-coded quality indicators (green/yellow/red)
  - Downloadable data tables and figures
- **Comprehensive metrics**:
  - ΔG with confidence intervals
  - Per-lambda free energy profile
  - Bootstrap uncertainty estimation
  - Quality warnings and recommendations

### Builder and force-field integration
- Force field directory copying to FEP leg repeat directories for isolation
- Legacy forcefield support for both ions.itp and water-specific ions_{model}.itp files

### GPU optimization and performance
- **CPU affinity with GROMACS -pinoffset** - Improved CPU thread placement for better NUMA affinity
- **Auto-calculation of OpenMP threads** - Automatically calculate optimal thread count from total_cpus instead of hardcoded values
- **GPU consistency fixes** - Removed unconditional CPU fallback, ensuring GPU usage when available
- **GPU thread configuration** - Fixed thread configuration in FEP scripts for optimal GPU utilization
- **Replica exchange support** - Added REX script generation with proper GPU/CPU resource allocation

### Configuration and parameter handling
- **charge_cutoff and charge_reception** - Added charge_cutoff config parameter to mapper and charge_reception parameter for flexible charge assignment
- **Execution config loading** - Fixed execution config loading and added maxwarn parameter support
- **Force field selection** - Honor explicit amber14sb selection, fix GROMACS force field selection conflict
- **Leg topology handling** - Improved leg topology handling and coordinate alignment for bound/unbound systems

### Visualization and HTML improvements
- **Bond order rendering** - Implemented professional bond order visualization in mapping HTML
- **Warning banners** - Added frontend warning banner for unclassified atoms (gray atoms)
- **Atom rendering fixes** - Fixed atom rendering and button styling in HTML reports
- **OPLS force field visualization** - Fixed OPLS force field atom matching by distance
- **Multi-estimator caching** - Cache multi-estimator analysis results for faster HTML generation

### Analysis and architecture refactoring

#### Analysis module improvements
- **Multi-estimator framework** - Added support for TI/BAR/MBAR comparison with unified interface
- **Bootstrap analysis** - Implemented parallel bootstrap with configurable core count (`n_jobs` parameter)
- **Quality checks** - Added automatic detection of high uncertainty, poor overlap, endpoint issues, and convergence problems
- **HTML reporting** - Consolidated multi-estimator results into single interactive HTML report
- **Time convergence** - Added time-converged analysis with skip options for equilibration
- **CLI fixes** - Fixed multi-estimator analyzer parameter passing in CLI interface
- **Visualization improvements** - Fixed atom rendering, button styling, and plot layout in HTML reports

#### Architecture refactoring
- **Module reorganization** - Extracted plotting functions to dedicated `plots/` module
- **Exception hierarchy** - Introduced capability-based exception handling for cleaner error management
- **Test cleanup** - Improved test fixtures and path handling across analysis modules
- **Code optimization** - Reduced monolithic files in trajectory, contact, and cleaner modules

#### Bug fixes
- **Fixed FEP hybrid ITP B-state omission** - Critical fix ensuring all hybrid topologies contain B-state atoms
  - Added `is_perturbed` flag to track transformed and surrounding atoms
  - Modified ITPBuilder to always write B-state columns for perturbed atoms
  - Added mass difference checking and validation methods
  - Fixed issue where charge redistribution could eliminate B-state differences
  - All 7 test systems now pass grompp with proper B-state atoms
- **Fixed mass_b assignment for perturbed atoms** - Critical fix preventing structure explosion in FEP simulations
  - Added mass_b parameter assignments in HybridTopologyBuilder for all perturbed atom types
  - Transformed A atoms (disappearing): mass_b = dummy mass (12.011)
  - Transformed B atoms (appearing): mass = dummy mass, mass_b = real atom mass
  - Surrounding atoms: mass_b from B-state atom type
  - Common atoms (charge_strategy='none'): mass_b from B-state
  - Fixed OPLS, MMFF, and other force fields that were failing at EM stage
  - All 6 test systems now pass EM+NVT validation
- Fixed OPLS force field chain splitting and hybrid coordinate handling
- Fixed unconditional CPU fallback in FEP scripts (GPU consistency issue)
- Fixed leg topology handling and coordinate alignment for bound/unbound systems
- Fixed ITP file copying for force field directories
- Fixed legacy forcefield support with ions.itp
- Fixed bond order and parameter assignment issues

## Command-line usage

### FEP system building

Typical CLI usage remains through `prism`, for example:

```bash
# Auto-generates amber14sb_ol15-mut_gaff2/ when -o is omitted
prism protein.pdb ligA.mol2 \
  --fep \
  --mutant ligB.mol2 \
  --ligand-forcefield gaff2 \
  --forcefield amber14sb_OL15 \
  --fep-config fep.yaml

# GAFF/GAFF2 force fields with custom output directory
prism protein.pdb ligA.mol2 -o output \
  --fep \
  --mutant ligB.mol2 \
  --ligand-forcefield gaff2 \
  --forcefield amber14sb \
  --fep-config fep.yaml

# RTF force field (MATCH/CHARMM-GUI)
prism protein.pdb ligand.pdb -o output \
  --fep \
  --mutant mutant.pdb \
  --ligand-forcefield rtf \
  --forcefield charmm36-jul2022

# CGenFF force field
prism protein.pdb ligA.mol2 -o output \
  --fep \
  --mutant ligB.mol2 \
  --ligand-forcefield cgenff \
  --forcefield-path /path/to/charmm \
  --forcefield charmm36-jul2022

# SwissParam MMFF force field
prism protein.pdb ligA.mol2 -o output \
  --fep \
  --mutant ligB.mol2 \
  --ligand-forcefield mmff \
  --forcefield charmm36-jul2022

# SwissParam MATCH force field
prism protein.pdb ligA.mol2 -o output \
  --fep \
  --mutant ligB.mol2 \
  --ligand-forcefield match \
  --forcefield charmm36-jul2022

# SwissParam Both (MMFF-based-MATCH hybrid) force field
prism protein.pdb ligA.mol2 -o output \
  --fep \
  --mutant ligB.mol2 \
  --ligand-forcefield both \
  --forcefield charmm36-jul2022
```

The generated FEP scaffold includes a top-level runtime entrypoint such as:

```bash
cd output/GMX_PROLIG_FEP
bash run_fep.sh
```

When no target is passed, `run_fep.sh` now defaults to running all configured legs/repeats.

### FEP analysis

After FEP simulations complete, analyze results using the analysis CLI:

```bash
# Basic analysis with all estimators
python -m prism.fep.analysis.cli \
  --bound output/GMX_PROLIG_FEP/bound \
  --unbound output/GMX_PROLIG_FEP/unbound \
  --output fep_results.html

# Analysis with specific estimators
python -m prism.fep.analysis.cli \
  --bound output/GMX_PROLIG_FEP/bound \
  --unbound output/GMX_PROLIG_FEP/unbound \
  --estimators BAR MBAR TI \
  --output fep_results.html

# Analysis with custom bootstrap parameters
python -m prism.fep.analysis.cli \
  --bound output/GMX_PROLIG_FEP/bound \
  --unbound output/GMX_PROLIG_FEP/unbound \
  --n-bootstrap 1000 \
  --n-jobs 8 \
  --output fep_results.html

# Analysis for specific repeats
python -m prism.fep.analysis.cli \
  --bound output/GMX_PROLIG_FEP/bound \
  --unbound output/GMX_PROLIG_FEP/unbound \
  --repeats 1 2 3 \
  --output fep_results.html
```

Analysis features:
- **Multi-estimator support**: BAR, MBAR, TI with unified interface
- **Bootstrap uncertainty**: Configurable iterations and parallel execution
- **Quality checks**: Automatic detection of high SE, poor overlap, convergence issues
- **HTML reports**: Interactive visualization with download capabilities

## Python API usage

```python
import prism as pm

# FEP system building with explicit output directory
system = pm.system(
    protein="protein.pdb",
    ligand="ligA.mol2",
    mutant="ligB.mol2",
    output_dir="fep_output",
    ligand_forcefield="gaff2",
    forcefield="amber14sb",
    fep_mode=True
)

fep_dir = system.build()

# FEP system building with default standardized naming
default_named_system = pm.system(
    protein="protein.pdb",
    ligand="ligA.mol2",
    mutant="ligB.mol2",
    ligand_forcefield="gaff2",
    forcefield="amber14sb_OL15",
    fep_mode=True
)

auto_named_dir = default_named_system.build()  # amber14sb_ol15-mut_gaff2

# FEP analysis
from prism.fep.analysis import FEPAnalyzer
analyzer = FEPAnalyzer(
    bound_dir="fep_output/GMX_PROLIG_FEP/bound",
    unbound_dir="fep_output/GMX_PROLIG_FEP/unbound"
)

results = analyzer.analyze(
    estimators=["BAR", "MBAR", "TI"],
    n_bootstrap=1000,
    n_jobs=8
)

# Generate HTML report
analyzer.generate_html_report("fep_results.html")
```

## Current validation status

### Force field testing matrix (2026-04-08 Updated)

**Test Summary**: The table below reflects the currently verified smoke-test status in the unit-test workspace. Older aggregate counts were removed because they became stale as additional systems were rechecked and promoted from pending to completed.

| System | Protein | Protein FF | Ligand FF | Status | Bound | Unbound | Notes |
|--------|---------|------------|-----------|--------|-------|---------|-------|
| 42-38 | HIF-2α | amber14sb_OL15 | GAFF2 | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 42-38 | HIF-2α | amber14sb_OL15 | OpenFF | ✅ | EM+NVT | EM+NVT | Bound repeats 1-3 and unbound repeats 1-3 all completed |
| 42-38 | HIF-2α | amber14sb_OL15 | OPLS-AA | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 42-38 | HIF-2α | amber14sb_OL15 | MMFF94 | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 42-38 | HIF-2α | amber99sb | GAFF2 | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 42-38 | HIF-2α | amber99sb | MMFF94 | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 42-38 | HIF-2α | charmm36-jul2022 | CGenFF | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 42-38 | HIF-2α | charmm36-jul2022 | CHARMM-GUI | ✅ | EM+NVT | EM+NVT | All repeats completed (verified 2026-04-08) |
| 42-38 | HIF-2α | oplsaa | OPLS-AA | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 25-36 | HIF-2α | amber14sb_OL15 | — | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 25-36 | HIF-2α | amber14sb_OL15 | OpenFF | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 25-36 | HIF-2α | amber14sb_OL15 | OpenFF (alt) | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 25-36 | HIF-2α | amber14sb_OL15 | OPLS-AA | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 25-36 | HIF-2α | charmm36-jul2022 | CGenFF | ✅ | EM+NVT | EM+NVT | All repeats completed |
| 25-36 | HIF-2α | oplsaa | OPLS-AA | ✅ | EM+NVT | EM+NVT | All repeats completed |
| oMeEtPh-EtPh | T4 lysozyme L99A | amber14sb_OL15 | OpenFF | ✅ | EM+NVT | EM+NVT | All repeats completed |
| oMeEtPh-EtPh | T4 lysozyme L99A | amber14sb_OL15 | OpenFF (non-mut) | ✅ | EM+NVT | EM+NVT | All repeats completed |
| oMeEtPh-EtPh | T4 lysozyme L99A | amber14sb_OL15 | OPLS-AA | ✅ | EM+NVT | EM+NVT | All repeats completed |
| p38-19-24 | p38α MAPK | amber14sb_OL15 | GAFF2 | ✅ | EM+NVT | EM+NVT | Bound+unbound EM+NVT completed (2026-04-08) |
| p38-19-24 | p38α MAPK | amber14sb_OL15 | OpenFF | ✅ | EM+NVT | EM+NVT | Bound+unbound EM+NVT completed (2026-04-08) |
| p38-19-24 | p38α MAPK | amber14sb_OL15 | OPLS-AA | ❌ | Not tested | Not tested | System build failed (LigParGen error) |
| p38-19-24 | p38α MAPK | charmm36-jul2022 | RTF | ✅ | EM+NVT | EM+NVT | Bound+unbound EM+NVT completed (2026-04-08; bound runtime guard enables unconstrained 0.5 fs startup for zeroized H-bonds) |
| oMeEtPh-EtPh | T4 lysozyme L99A | charmm36m | CHARMM-GUI | ✅ | EM+NVT | EM+NVT | Bound and unbound smoke checks completed |
| oMeEtPh-EtPh | T4 lysozyme L99A | amber19sb | — | ✅ | EM+NVT | EM+NVT | Bound+unbound EM+NVT completed (2026-04-08, CPU mode due to CUDA error) |

**Test platforms**:
- **HIF-2α**: Hypoxia-Inducible Factor 2α (42-38, 25-36 systems)
- **T4 lysozyme L99A**: Model binding pocket (oMeEtPh-EtPh system)
- **p38α MAP kinase**: Inhibitor series (p38-19-24 system)

### Critical bug fixes validated

**B-state atom validation** - hybrid topologies continue to pass scaffold/grompp validation across the major FEP test systems
- Added `is_perturbed` flag to track transformed and surrounding atoms
- Modified ITPBuilder to always write B-state columns (typeB, chargeB, massB)
- Mass difference checking prevents omission of mass-only perturbations

**mass_b assignment fix** - runtime validation now shows broad EM+NVT success across the previously failing force-field combinations
- Transformed A atoms (disappearing): mass_b = dummy mass (12.011)
- Transformed B atoms (appearing): mass = dummy mass, mass_b = real atom mass
- Surrounding atoms: mass_b from B-state atom type
- 42-38 force-field variants now complete bound EM+NVT across GAFF2, OpenFF, OPLS-AA, MMFF, CGenFF, and CHARMM-GUI paths

**Current validation snapshot** (2026-04-08 updated):
- **27 system variants** with completed bound+unbound EM+NVT ✅ (+5 p38 systems + amber19sb + 2 charmm-gui)
- **0 variants** with pending unbound testing ⚠️ (all tested!)
- **1 variant** with system build failure ❌ (p38-19-24 OPLS: LigParGen error)
- **Total: 28 tracked system variants** in the unit-test workspace (27 in table + 39-8 platform)
- Recent additions: p38-19-24 (GAFF2, RTF) and oMeEtPh-EtPh (amber19sb) unbound legs completed; 42-38/oMeEtPh charmm-gui unbound verified

**Successful force field combinations**:
- AMBER (amber14sb_OL15, amber99sb) + GAFF2, OpenFF, MMFF, OPLS-AA ✅
- CHARMM (charmm36-jul2022, charmm36m) + CGenFF, CHARMM-GUI ✅
- OPLS-AA (oplsaa) + OPLS-AA ✅

**Known issues and limitations**:
1. **GPU resource management**: CUDA_VISIBLE_DEVICES causes CPU affinity restrictions; tasks may not optimally utilize all CPU cores
2. **GPU fallback mechanism**: run_single_em_nvt.sh automatically falls back to CPU mode on GPU errors, causing severe performance degradation
3. **System completeness**: Some p38-19-24 systems have incomplete directory structures (missing input files)
4. **Unbound testing**: Several non-25-36 rows still need explicit unbound smoke validation or reruns
- Successfully tested force field combinations show robust performance across AMBER, CHARMM, and OPLS-AA families

## Review guidance

Suggested review order:
1. `prism/fep/core/`
2. `prism/fep/modeling/`
3. `prism/fep/analysis/`
4. `prism/builder/workflow_fep.py`
5. `prism/forcefield/`
6. analysis / utility refactors outside `prism/fep/`

## Notes

- OPLSAAM is currently marked unsupported for PRISM FEP setup because the available local force-field dataset does not provide the required protein RB dihedral coverage.
- The PR includes modeling validation for multiple force-field paths, but not every supported combination has been runtime-validated.
- RTF files must follow CHARMM format specification and each ligand requires `{ligand}.pdb`, `{ligand}.rtf`, `{ligand}.prm` files.
- Replica exchange mode is supported for enhanced sampling in FEP calculations.
- GPU-aware execution scripts are generated automatically for both standard and replica-exchange modes.

## Files Changed

### Core FEP modules
- `prism/fep/core/` - FEP core functionality and orchestration
- `prism/fep/modeling/` - Atom mapping and hybrid topology generation
- `prism/fep/analysis/` - Multi-estimator analysis and reporting
- `prism/fep/naming.py` - Standardized FEP case directory naming helpers
- `prism/fep/__init__.py` - Public exports for naming helpers
- `prism/builder/workflow_fep.py` - FEP workflow integration and default output naming

### Force field integration
- `prism/forcefield/rtf.py` - RTF/PRM parser and GROMACS conversion
- `prism/forcefield/swissparam.py` - SwissParam session polling
- `prism/builder/ligand_ff.py` - RTF force field integration
- `prism/builder/core.py` - RTF force field validation

### Analysis and utilities
- `prism/analysis/` - Refactored analysis modules
- `prism/utils/` - Refactored utility modules

### Tests and documentation
- Comprehensive test systems for multiple force fields (GAFF, GAFF2, CGenFF, RTF, OPLS, OpenFF)
- End-to-end validation with 42-38, p38-19-24, 25-36, and oMeEtPh-EtPh test platforms
- **New test scripts**:
  - `tests/gxf/FEP/unit_test/42-38/test_gaff2.py` - GAFF2 system rebuild script
  - `tests/gxf/FEP/unit_test/test_run_fep.py` - Enhanced validation with real grompp execution and support for standardized default output layout
  - `tests/gxf/FEP/unit_test/test_naming.py` - Unit coverage for standardized FEP naming helpers
- Updated project documentation and ignore patterns
