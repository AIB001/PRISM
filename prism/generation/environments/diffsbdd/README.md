# DiffSBDD runtime

DiffSBDD is installed outside PRISM and invoked through the Slurm/Conda
contract in `prism/configs/generation.diffsbdd.slurm.example.yaml`.

## Pinned upstream contract

- Repository: `https://github.com/arneschneuing/DiffSBDD`
- Checkpoint: `crossdocked_fullatom_cond.ckpt` from the upstream release
- Checkpoint SHA256:
  `07f86764bf569aafbc40a9c15fc02de8e2550437dd0f17f657eab3abe66c372c`

The upstream checkout used on the integration cluster was copied without its
`.git` directory, so no source commit is recorded here and `model_commit` is
empty in every DiffSBDD config. That is provenance only -- the run manifest
records what it was given and nothing validates it -- but it means a fresh
clone tracks upstream `main`, so pin a commit locally before you rely on a run
being reproducible.

The source repository is MIT (`Copyright (c) 2022 Arne Schneuing, Yuanqi Du,
Charles Harris`). An earlier note here said the checkout carried no `LICENSE`;
that was the incomplete copy, not upstream. The checkpoint is governed
separately by Zenodo record 8183747 under CC BY 4.0, which is what PRISM
mirrors and what `prism weights list` reports.

## Environment

Reuses the MolCRAFT environment verbatim (`prism-gen-molcraft`: Python 3.9,
PyTorch 1.12.0, CUDA 11.6, PyTorch Geometric 2.1, RDKit). Two upstream imports
are stubbed by the wrapper because they exist only for training-time logging:
`wandb` and `imageio`. Nothing else is installed for DiffSBDD.

CPU inference is supported (`--device cpu`); CUDA works with the MolCRAFT
environment's torch build.

## Input contract

Reference-guided: a full protein PDB and a reference ligand SDF are both
required. The wrapper calls the official CrossDocked full-atom model with
geometry-only sampling, without upstream sanitization or force-field relax
(the PRISM quality-control layer does both jobs, and doing them inside the
wrapper would mix upstream behavior into the QC numbers).

The adapter exposes `diffsbdd_reference_heavy_atoms` (the reference ligand's
heavy-atom count) so a configuration may pin `--num-nodes-lig` and sample the
size the pocket was seen with, instead of the learned size distribution.

## Known behavior

DiffSBDD is the only integrated model whose output occasionally carries
valence errors (~3% of molecules in the multisite benchmark). Those molecules
are quarantined by the QC layer with a written repair proposal, never silently
fixed.
