# PRISM-FEP GROMACS Integration Module

This file provides guidance for working with GROMACS integration and MDP templates for FEP calculations.

## Directory Naming Rule (CRITICAL)

- Use simple case names for force-field combinations: `<protein_ff>-mut_<ligand_ff>` (example: `charmm36m-mut_mmff`).
- Never create nested duplicate system directories like `.../GMX_PROLIG_FEP/GMX_PROLIG_FEP/`.
- Avoid ad-hoc suffixes (`_pkgfix*`, `_final*`, `_new*`) in directory names.
- Default output directory for FEP cases should follow `<protein_ff>-mut_<ligand_ff>` when users do not pass an explicit output path.
- For 42-38 maintenance, move legacy non-canonical outputs into `tests/gxf/FEP/unit_test/42-38/Archive/`.

## Module Overview

The GROMACS integration module handles:
- **MDP template generation**: For equilibration and production runs
- **ITP file building**: Hybrid topology file construction
- **Free energy parameters**: Lambda scheduling and soft-core potentials

## Key Components

### MDP Templates
**Location**: `prism/fep/gromacs/mdp_templates.py`

**Templates provided**:
- **EM**: Energy minimization (steepest descent)
- **NVT**: Constant volume equilibration
- **NPT**: Constant pressure equilibration
- **Production**: Free energy calculation with lambda coupling

**Free energy settings**:
```mdp
free-energy = yes
init-lambda-state = 0
calc-lambda-neighbors = -1
delta-lambda = 0
couple-moltype = LIG
couple-lambda0 = vdw-q  # A state: vdw + Coulomb
couple-lambda1 = vdw-q  # B state: vdw + Coulomb
couple-intramol = no

# Soft-core potentials (recommended)
sc-alpha = 0.5
sc-power = 1
sc-sigma = 0.3
```

### ITP Builder
**Location**: `prism/fep/gromacs/itp_builder.py`

**Purpose**: Build hybrid topology ITP files from mapping results

**Key responsibilities**:
- Create [moleculetype] section with hybrid molecule definition
- Generate [atoms] section with A/B state parameters
- Handle [pairs] and [exclusions] for transformed atoms
- Create dummy atom types (DUM) for non-interacting states

**Output structure**:
```itp
[ moleculetype ]
; Name      nrexcl
LIG          3

[ atoms ]
;   nr  type  resnr  residue  atom   cgnr    charge       mass       typeB    chargeB      massB
    1    ca     1     LIG      C1      1     -0.123       12.01      ca       -0.123      12.01
    2    ha     1     LIG      H2      2      0.078        1.008     DUM       0.0        1.008
    3   DUM     1     LIG      C3      3      0.0          12.01      c3        0.123      12.01
```

**Note**: TypeB and chargeB are used for the B state. DUM atoms have zero charge in one state.

## Lambda Scheduling

### Default Lambda Windows
**Standard setup**: 11 windows (λ = 0.0, 0.1, ..., 1.0)

**Recommended**: 21-32 windows for better convergence
```yaml
fep:
  lambda_windows: 32
```

### Lambda State Control
**In MDP files**:
```mdp
init-lambda-state = 0  # Set by script for each window
```

**In run scripts**:
```bash
for window in $(seq 0 31); do
    gmx grompp -f prod.mdp -c npt.gro -p topol.top -o tpr_${window}.tpr \
        -pp topol_processing.mod -defindex ${window}
done
```

## Soft-Core Potentials

**Purpose**: Avoid singularities when atoms appear/disappear

**Parameters** (GROMACS default):
```mdp
sc-alpha = 0.5       # Alpha parameter
sc-power = 1         # Power for lambda dependence
sc-sigma = 0.3       # Soft-core sigma (nm)
```

**Recommended for large perturbations**:
```mdp
sc-alpha = 0.5
sc-power = 2
sc-sigma = 0.3
```

## Force Field Compatibility

### GAFF/GAFF2
**Atom types**: ca, c3, ha, hc, na, oa, etc.

**Nonbonded parameters**:
- [atomtypes] section in hybrid.itp
- [defaults] section (if using explicit combination rules)

### CGenFF
**Atom types**: CA, CB, CG, HA, HB, etc.

**Special handling**:
- Auto-detect OPLS-style types for [defaults] block
- Extract [defaults] to separate file

### OpenFF
**Atom types**: output_0, output_1, output_2, ...

**Mapping**: DistanceAtomMapper skips type checking for OpenFF

### OPLS-AA
**Atom types**: opls_800, opls_801, opls_802, ...

**Special handling**:
- Auto-detect and generate [defaults] block
- DistanceAtomMapper skips type checking

## File Structure

### Hybrid Topology Files
```
common/hybrid/
├── hybrid.itp              # Main hybrid topology
├── atomtypes_hybrid.itp    # Atom type definitions
├── defaults_hybrid.itp     # Nonbonded parameters (optional)
└── hybrid.gro              # Coordinates
```

### System Topology
```
bound/
├── topol.top               # System topology (includes hybrid.itp)
├── topol_processing.mod    # Processing modifiers (lambda states)
└── posre.itp               # Position restraints (if needed)
```

## Common Issues

### Missing [defaults] Block
**Problem**: GROMACS complains about missing combination rules

**Solution**: Extract [defaults] from hybrid.itp to defaults_hybrid.itp

**Code**: `modeling/hybrid_package.py::_extract_defaults_block()`

### Dummy Atom Types Not Defined
**Problem**: "Atom type DUM not defined"

**Solution**: Add DUM type definition to atomtypes_hybrid.itp
```itp
[ atomtypes ]
; name  bond_type  mass  charge  ptype  sigma  epsilon
DUM     DUM         1.0   0.0     A     0.0    0.0
```

### Lambda State Not Applied
**Problem**: All windows produce same results

**Solution**: Check `init-lambda-state` in MDP and `-defindex` in grompp command

## MDP Template Customization

### Temperature
**Default**: 310 K (body temperature)
```mdp
ref-t = 310
gen-temp = 310
```

**Customization**: Set via FEPConfig or YAML
```yaml
fep:
  temperature: 300
```

### Pressure
**Default**: 1 bar
```mdp
ref-p = 1.0
compressibility = 4.5e-5
```

### Time Step
**Default**: 2 fs
```mdp
dt = 0.002
```

**Constraints**: h-bonds (LINCS)
```mdp
constraints = h-bonds
constraint-algorithm = lincs
```

## Code Locations

- **MDP templates**: `prism/fep/gromacs/mdp_templates.py`
- **ITP builder**: `prism/fep/gromacs/itp_builder.py`

## Related Documentation

- `CLAUDE.md` (parent) - FEP module overview
- `../modeling/CLAUDE.md` - System building workflow
- GROMACS Free Energy Documentation: https://manual.gromacs.org/current/fep.html

## GROMACS Version Compatibility

**Tested with**: GROMACS 2020, 2021, 2022, 2023

**Free energy syntax**: Stable across versions (no breaking changes)

**Recommended**: GROMACS 2022 or later for best performance
