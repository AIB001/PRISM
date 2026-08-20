# PRISM-FEP Example Inputs

These version-controlled inputs support complete PRISM-FEP scaffold builds.
They serve as both user examples and developer integration cases; generated
`GMX_PROLIG_FEP/` directories are intentionally not version controlled.

| Case | Protein system | Transformation |
| --- | --- | --- |
| `42-38` | HIF-2alpha | 42 to 38 |
| `p38-19-24` | p38alpha MAPK | 19 to 24 |
| `oMeEtPh-EtPh` | T4 lysozyme L99A | oMeEtPh to EtPh |

Each case contains `input/` (one protein PDB and two aligned MOL2 ligands)
and `configs/`. `fep.yaml` is the recommended GAFF2 profile where present.
Other YAML files preserve the represented force-field-specific settings and
mapping variants.

## Force-Field Pairing

The ligand backend is selected by `--ligand-forcefield`; the matching YAML
records its intended profile. FEP also requires a compatible protein force
field, because the ligand and protein must share non-bonded conventions.

| Ligand backend | Required protein FF family | Example configurations |
| --- | --- | --- |
| GAFF/GAFF2 | AMBER | `fep.yaml`, `config_gaff2.yaml` |
| OpenFF | AMBER | `config_openff.yaml` |
| OPLS-AA | OPLS | `config_opls.yaml`, `case_opls.yaml` |
| MMFF/MATCH/Hybrid (SwissParam) | CHARMM36 | `config_mmff.yaml`, `case_mmff.yaml` |
| CGenFF website | CHARMM36 | `fep_charmm.yaml`, `config_charmm.yaml` |
| RTF/CHARMM | CHARMM36 | `p38-19-24/configs/fep_rtf.yaml` |

PRISM rejects incompatible pairings before it builds a system. In particular,
a successful `grompp` run is not evidence that an Amber+OPLS or Amber+SwissParam
system is physically consistent.

## Build a Scaffold

Run the public CLI from a case directory and write output outside this input
tree. This GAFF2 example uses the AMBER family on both sides of the system:

```bash
cd tests/prism_fep/unit_test/42-38
prism input/receptor.pdb input/42.mol2 \
  --fep --mutant input/38.mol2 \
  --forcefield amber14sb_OL15 \
  --ligand-forcefield gaff2 \
  --fep-config configs/fep.yaml \
  --output /path/to/42-38-gaff2
```

For p38-19-24, use `input/protein.pdb input/19.mol2 --mutant input/24.mol2`.
For oMeEtPh-EtPh, use `input/receptor.pdb input/oMeEtPh.mol2 --mutant
input/EtPh.mol2`. OpenFF uses an AMBER protein force field; choose an OPLS
protein force field for OPLS-AA and CHARMM36 for CGenFF, RTF, or SwissParam.

## CGenFF Website Input

`42-38/input/42_cgenff/` and `42-38/input/38_cgenff/` provide the paired
CGenFF website topology and coordinate files plus the required CHARMM
parameter includes. Supply one `--forcefield-path` for each ligand:

```bash
cd tests/prism_fep/unit_test/42-38
prism input/receptor.pdb input/42.mol2 \
  --fep --mutant input/38.mol2 \
  --forcefield charmm36-jul2022 \
  --ligand-forcefield cgenff \
  --forcefield-path input/42_cgenff \
  --forcefield-path input/38_cgenff \
  --fep-config configs/fep_charmm.yaml \
  --output /path/to/42-38-cgenff
```

Construction writes `GMX_PROLIG_FEP/common/hybrid/mapping.html`. Inspect that
report before production: it displays the common, surrounding, and transformed
atoms used by the hybrid topology.
