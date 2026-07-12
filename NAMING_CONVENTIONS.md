# PRISM Naming Conventions

Standardized naming rules for directories, files, and generated artifacts. All
new code MUST follow these. This document is the canonical reference cited by
`CLAUDE.md`.

## 1. System directories

Always use the `GMX_PROLIG_*` prefix for built simulation systems:

| Build mode | Directory |
|---|---|
| Standard MD | `GMX_PROLIG_MD/` |
| PMF (steered MD + umbrella) | `GMX_PROLIG_PMF/` |
| REST2 | `GMX_PROLIG_REST2/` |
| MM/PBSA | `GMX_PROLIG_MMPBSA/` |
| FEP | `GMX_PROLIG_FEP/` |
| Membrane | `GMX_PROLIG_MEMB/` |

❌ Do not invent names like `FEP_SYSTEM/`, `output/`, `my_system/`.

## 2. Force-field directories

Use the generator's method, never a hardcoded literal:

```python
ff_dir = Path(ffgen.output_dir) / ffgen.get_output_dir_name()   # ✅
ff_dir = output_dir / "LIG.amb2gmx"                              # ❌ hardcoded
```

Canonical ligand FF output subdirs: `LIG.amb2gmx` (GAFF/GAFF2), `LIG.openff2gmx`,
`LIG.cgenff2gmx`, `LIG.opls2gmx`, `LIG.sp2gmx` (SwissParam), `LIG.rtf2gmx`.

## 3. Never hardcode paths

Resolve artifacts through helpers (e.g. `prism.forcefield.registry`), not string
literals. No absolute, machine-specific, or user-specific paths in the package.

## 4. No nested duplicate system directories

✅ `.../oMeEtPh-EtPh/charmm36m-mut_mmff/GMX_PROLIG_FEP/`
❌ `.../charmm36m-mut_mmff_pkgfix2/GMX_PROLIG_FEP/GMX_PROLIG_FEP/`

## 5. Case-directory naming stays simple and readable

✅ `<protein_ff>-mut_<ligand_ff>` (e.g. `charmm36m-mut_mmff`)
❌ ad-hoc suffixes: `*_pkgfix*`, `*_final*`, `*_new*`, `*_v2*`, `*_copy*`, `*_tmp*`.

## 6. Test directories

Use descriptive, force-field-specific names (`gaff_test/`, `openff_test/`,
`test_42_38/`), not `*_e2e_test/` / `*_output_final/`.

## 7. No machine-specific or non-software files in the package

Hardware values (core/GPU counts, ids), absolute host paths, benchmark data,
model checkpoints, and planning documents must NOT live inside `PRISM-main/`.
They belong in the surrounding project folder.
