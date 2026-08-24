# Pocket2Mol runtime

Pocket2Mol is intentionally installed outside the PRISM environment. Pin the
official source commit, keep its CUDA/PyTorch/PyG dependencies isolated, and
configure PRISM with
`prism/configs/generation.pocket2mol.slurm.example.yaml`.

## Pinned upstream contract

- Repository: `https://github.com/pengxingang/Pocket2Mol`
- Source commit: `836a0c4ce487297ad24bc54ac2ebd163de13242c`
- Official checkpoint name: `pretrained_Pocket2Mol.pt`
- Official checkpoint SHA256:
  `e9d8152dfdf034dabf48e566904e9d7c6cc1b96843bdc4d0026d891678913911`
- Official tested stack: Python 3.8, PyTorch 1.10.1, CUDA 11.3, PyG 2.0,
  RDKit 2022.03, and BioPython 1.79

The upstream `sample_for_pdb.py` interface accepts a full protein PDB, a pocket
center, and a scalar cubic box size. PRISM does not silently convert CIF/mmCIF
proteins for this adapter. A reference SDF/MOL/MOL2 or cropped pocket PDB can be
used to calculate the center; center YAML is accepted directly. Residue-list
pockets are rejected because the official entry point has no residue-list
contract.

## Security

PyTorch checkpoints use pickle internally. The PRISM wrapper first reads the
official checkpoint with an allowlist containing only `OrderedDict`,
`EasyDict`, and the required PyTorch tensor reconstruction symbols. It verifies
that the result contains model/config mappings and then writes a trusted copy
inside the run directory. Configure the original artifact SHA256 as shown in
the example YAML; the executor verifies it before the wrapper starts.

## Cluster validation

The adapter was accepted end to end on Slurm with PyTorch 1.13.1, CUDA 11.6,
PyG 2.2.0, RDKit 2022.03.2, and BioPython 1.79. A `1 pocket x 5 molecules`
test completed on one GPU; all five outputs were unique, had one 3D conformer,
and passed RDKit sanitization. This compatible stack is useful when an existing
offline cluster environment is available, but the official versions above
remain the reproducibility baseline for a fresh image build.

For a quick smoke test, reduce `--beam-size` in the configured wrapper command
to 50 and request one molecule. For production, retain the upstream default of
300 unless GPU/runtime benchmarking justifies another value.
