# FEP Unit Test Cases

Representative PRISM-FEP test systems with input structures and config files.

## Structure

```
unit_test/
├── 42-38/          # HIF-2α, ligand pair 42→38
├── p38-19-24/      # p38α MAPK, ligand pair 19→24
└── oMeEtPh-EtPh/   # T4 lysozyme L99A, ligand pair oMeEtPh→EtPh
```

Each case contains:
- `input/` — protein PDB + ligand MOL2 files (required to run)
- `configs/` — YAML config files for various force-field combinations

## Usage

```bash
cd tests/prism_fep/unit_test/42-38
python -m prism.fep \
  --protein input/receptor.pdb \
  --ligand-a input/42.mol2 \
  --ligand-b input/38.mol2 \
  -c configs/config_gaff2.yaml \
  -o output/
```

### CGenFF example (42-38 only)

The `42-38/input/42_cgenff/` directory contains a pre-generated CGenFF
toppar for ligand 42, demonstrating manual `--forcefield-path` usage:

```bash
python -m prism.fep \
  --protein input/receptor.pdb \
  --ligand-a input/42.mol2 \
  --ligand-b input/38.mol2 \
  --forcefield-path input/42_cgenff/charmm36.ff \
  -c configs/config_charmm.yaml \
  -o output/
```

## Notes

- Simulation output directories are NOT included (gitignored).
- Only files needed to build and run an FEP system are present.
- Full test suite lives in `tests/prism_fep/` (see sibling test_*.py files).