# FLOWR runtime

FLOWR is installed outside PRISM and invoked with the Slurm/Conda contract in
`prism/configs/generation.flowr.slurm.example.yaml`.

## Pinned upstream contract

- Repository: `https://github.com/jule-c/flowr`
- Source commit: `bbb8bfb4e8d8dd16bff1d56a7f96246920e643ba`
- Checkpoint record: Zenodo `15737419`
- Selected artifact: `flowr_noHs.ckpt` (629,067,043 bytes)
- Published MD5: `e0d49765629c4b9a6847c4dcdbf678ad`
- Checkpoint SHA256:
  `1e26e06181fa0d6850a159a12acac36847853f979f76b8fdaa148ed3b502a077`

The repository README declares MIT, but the pinned checkout does not contain
the linked `LICENSE` file. The Zenodo checkpoint record also does not publish a
machine-readable artifact license. PRISM therefore does not bundle or
redistribute either source or weights; resolve the upstream licensing gap
before redistribution or commercial use.

## Environment and input contract

The official environment uses Python 3.11, PyTorch 2.5.1, CUDA toolkit 11.8,
Lightning 2.5.0, RDKit 2023.9.6, OpenBabel, Biotite, PyMOL, OpenMM, and the
OpenFF toolkit. Upstream recommends a CUDA GPU with roughly 40 GB VRAM for
generation.

FLOWR accepts PDB/CIF structures and SDF/PDB reference ligands. A full-protein
reference-ligand request is cropped at 6 Å by default. For PRISM structure,
center, and residue-list pockets, the adapter first prepares a PDB pocket and
the wrapper disables upstream re-cropping. A reference ligand remains required
for this initial integration; PRISM does not create a placeholder ligand or
guess a heavy-atom count.

## Narrow generation path

The wrapper calls the official `flowr.gen.generate_from_pdb` module with rigid
holo-pocket generation, sampled molecule sizes, uniqueness/validity filtering,
and the official uniform categorical strategy. It deliberately omits optional
interaction recovery, protonation, docking, and complex optimization so those
evaluation dependencies do not change generation behavior or become runtime
requirements.

The isolated launcher supplies inert bridges for those unused evaluation
interfaces, PyMOL, OpenMM optimization, and PDBInf bond-order optimization.
This also keeps the narrow generation path compatible with glibc 2.17 compute
nodes. The wrapper removes stale PDB `CONECT` records from a per-run copy before
upstream preprocessing; Biotite then reconstructs standard residue bonds. CIF
input is not rewritten.

The checkpoint is first loaded with a strict allowlist limited to ordered
mappings and PyTorch tensor/storage reconstruction, validated for Lightning
`state_dict` and `hyper_parameters` mappings, and rewritten inside the run
directory before upstream loads it.

## Portable environment and validation

The validated local bundle is
`.prism-build/flowr/prism-gen-flowr-v2.tar.gz` (3,316,110,580 bytes), SHA256
`4a21c4d77691480396c987d3228b7e7fd82292c83a4acb34c49772c8b3c32f94`.
After extraction, invoke `conda-unpack` explicitly with the environment's
Python, because some cluster login shells still resolve `python` to Python 2:

```bash
PREFIX="$HOME/conda-envs/prism-gen-flowr"
mkdir -p "$PREFIX"
tar -xzf prism-gen-flowr-v2.tar.gz -C "$PREFIX"
"$PREFIX/bin/python" "$PREFIX/bin/conda-unpack"
```

On 2026-08-14 this exact relocated bundle completed a real
`1 molecule x 10 integration steps` run on a 24,576 MiB RTX 3090 (seed 1234)
without OOM. Upstream sampling took 4.34 seconds. The SDF sanitized in RDKit,
had 27 atoms and finite 3D coordinates, and had SHA256
`1c1be77967b0a73124135c3357810ec88b8e2cb8d7fe91b77af537c12b723104`.
This smoke test establishes operability on 24 GB VRAM, not capacity for the
default 100-step or high-throughput workloads.
