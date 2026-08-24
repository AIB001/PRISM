# MolCRAFT runtime

MolCRAFT is installed outside PRISM and invoked with the Slurm/Conda contract
in `prism/configs/generation.molcraft.slurm.example.yaml`.

## Pinned upstream contract

- Repository: `https://github.com/GenSI-THUAIR/MolCRAFT`
- Source commit: `743dfb5198c44ca0cec180822bd54258f9d386f1`
- Source subdirectory: `MolCRAFT`
- Source license: CC BY-NC-SA 4.0
- Official checkpoint Google Drive ID:
  `1TcUQM7Lw1klH2wOVBu20cTsvBTcC1WKu`
- Checkpoint SHA256:
  `0a54f9bfb847ea87e7aec538df90a009c8c49d12a7b5d09947d04b7f9444e024`

The NonCommercial and ShareAlike terms are material. PRISM does not bundle or
redistribute the MolCRAFT source or checkpoint; users must review the upstream
license and install those artifacts separately.

## Environment and input contract

The official baseline uses Python 3.9, PyTorch 1.12.0, CUDA 11.6, PyTorch
Geometric 2.1, PyTorch Lightning 2.0.8, and RDKit. Its custom-pocket sampler
requires a PDB structure and a reference ligand SDF, then selects protein atoms
within the model's 10 Å neighborhood.

For reference-ligand pocket input, PRISM passes the complete protein. For a
cropped PDB, center, or residue-list input, PRISM first materializes the
requested pocket as PDB and passes that prepared structure. The reference SDF
is still mandatory because MolCRAFT uses it both to localize the pocket and to
define its ligand-size prior.

## Reproducibility and security

The published checkpoint embeds the training configuration under
`hyper_parameters`; the repository's convenience script instead expects a
relative `checkpoints/config.yaml`. The PRISM wrapper restricted-loads the
checkpoint, verifies the expected Lightning mappings, writes the embedded
configuration into an isolated per-run directory, and rewrites a trusted
checkpoint before the official sampler reads it. The upstream source is not
patched.

The published checkpoint records `time_emb_dim: 1`, while its ligand embedding
weight has input width 13, equal to the 13 ligand classes rather than the 14
inputs that the current source would construct after concatenating time. The
wrapper detects this from the trusted tensor shape and sets the per-run
`time_emb_dim` to zero. Any other config/weight mismatch is rejected instead of
being loaded non-strictly.

The wrapper supplies inert import bridges for the callback's unused Vina,
PoseCheck, W&B, and visualization helpers. It retains the upstream molecule
reconstruction function but skips optional post-generation metrics. Model
construction, sampling, reconstruction, and candidate SDF writing remain
upstream MolCRAFT code.

## Validated smoke run

On 2026-08-14 the pinned source and checkpoint completed a real
`1 molecule x 10 steps` run on a 24,576 MiB RTX 3090 (seed 1234), using the
relocated portable environment described below rather than another model's
environment. The generated SDF sanitized successfully in RDKit and contained
finite 3D coordinates. Its canonical SMILES was
`CC1(NC(=O)c2cc(Cl)c(O)c3ccccc23)CCCNC1`; the validated SDF SHA256 was
`a5ea651f3653039c6cecfcbd81cd69e42bc4b3d99b3bd6d0ca66d13583830843`.

The cluster exposes a system CUDA 11.4 path before environment libraries. To
avoid loading the system `libcublasLt.so.11` alongside the environment's CUDA
11.6 PyTorch build, the wrapper prepends its own `sys.prefix/lib` to
`LD_LIBRARY_PATH` for the isolated upstream sampler process.

The local portable environment artifact
`.prism-build/molcraft/prism-gen-molcraft.tar.gz` is 3,114,745,133 bytes with
SHA256
`7ddf5461497a9e9c0907e247366a8c19586c41d52279af53feb942d78057cfcd`.
Source and checkpoint remain separate because of the upstream license terms.
