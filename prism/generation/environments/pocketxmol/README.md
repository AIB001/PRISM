# PocketXMol runtime

PocketXMol is installed outside the PRISM environment and invoked through the
Slurm/Conda contract in
`prism/configs/generation.pocketxmol.slurm.example.yaml`.

## Pinned upstream contract

- Repository: `https://github.com/pengxingang/PocketXMol`
- Source commit: `65488cf635c856101dbe703ac97e2f10f58e005c`
- Source license: MIT
- Model release: Zenodo record `17801271`, version 1.0.0
- Model artifact license: CC BY 4.0
- Official checkpoint: `data/trained_models/pxm/checkpoints/pocketxmol.ckpt`
- Checkpoint SHA256:
  `76d918122e598a3c53c77727898659f059d217e680834b7d223b0205c2715c3a`
- Training config SHA256:
  `7225ef261d7bf20afeec6b41760a7dda8608110a37cb3361395ab49c49a0963d`

The official environment specifies Python 3.8, PyTorch 2.0.0, CUDA 11.7,
PyTorch Geometric 2.3, RDKit 2023.09.3, BioPython 1.83, OpenBabel 3.1.1,
Pandas 1.5.2, SciPy, scikit-learn, NetworkX 2.8, and PeptideBuilder 1.1.0.
Use those versions as the baseline for a fresh environment or image.

## Input contract

PRISM delegates to the official `scripts/sample_use.py` SBDD path. It requires
a full protein PDB and a pocket center. A reference SDF can additionally be
passed to define pocket extraction; MOL and MOL2 inputs still supply a center
but are not silently presented to upstream as SDF. Cropped pocket PDB inputs
are supported for their calculated center while the complete protein remains
the structure passed to PocketXMol. Residue-list pockets are rejected because
the upstream custom-inference path has no residue-list interface.

The adapter uses a 10 Å extraction radius for reference-ligand inputs and
15 Å for center or structure inputs unless the pocket specification supplies
an explicit radius. These values follow the ranges documented by upstream.

The wrapper writes a per-run task configuration, invokes the official sampler,
and collects only rows without an upstream failure tag. PocketXMol's
`cfd_traj`, `cfd_pos`, `cfd_node`, and `cfd_edge` values are retained as SDF
properties. PRISM exposes `cfd_traj` as the model-native score but does not
compare it directly with scores from other generators.

## Inference-only bridge

The official sampling script imports `DataModule` from its PyTorch Lightning
training script, although inference only calls `get_featurizers()` and
`get_in_dims()`. The PRISM wrapper installs a small runtime module containing
those two methods before running the unmodified sampler. Model construction,
checkpoint loading, denoising, molecule reconstruction, and confidence
calculation remain upstream code. This removes the training-only Lightning
dependency from the runtime and reduces the size of offline deployments.

## Security

PyTorch checkpoints use pickle internally. PRISM verifies the configured
artifact SHA256 before execution. The wrapper then loads the published archive
with an allowlist limited to `argparse.Namespace`, `OrderedDict`, `EasyDict`,
the required PyTorch storages, and tensor reconstruction symbols. It validates
the `state_dict` mapping and writes a trusted copy inside the model run before
the official sampler reads it.

## Cluster validation

The Slurm path was accepted end to end on one GPU with Python 3.8, PyTorch
1.13.1, CUDA 11.6, PyTorch Geometric 2.2.0, RDKit 2022.03.2, BioPython 1.79,
NetworkX 2.8, and PeptideBuilder 1.1.0. The official checkpoint loaded all
1,128 model tensors with no missing or unexpected keys. A
`1 pocket x 5 molecules` test completed all 100 denoising steps and produced
five unique molecules; every output passed RDKit sanitization, contained finite
3D coordinates, and retained its `cfd_traj` score.

This compatible stack is useful for the tested offline cluster. The official
versions above remain the preferred reproducibility baseline for new runtime
images, because compatibility outside the tested SBDD path is not implied.
