# TargetDiff runtime contract

PRISM delegates to the official TargetDiff `scripts/sample_for_pocket.py` entry
point. The integration was verified against upstream commit
`142f1eb7178480d435fe0b8cb95a99beb48997c7`.

The official runtime uses Python 3.8, PyTorch 1.13.1, CUDA 11.6, PyTorch
Geometric 2.2.0, RDKit 2022.03.2, and OpenBabel. Build the image from a checkout
fixed to that commit, copy the PRISM TargetDiff wrapper into the image, and set
the image's working tree to `/opt/targetdiff`.

A Dockerfile is included in this directory. Build it from the PRISM repository
root so its `COPY` path is available:

```bash
docker build \
  --file prism/generation/environments/targetdiff/Dockerfile \
  --tag local/prism-targetdiff:142f1eb \
  .
```

For production, pass a digest-pinned `BASE_IMAGE`, publish or otherwise retain
the built image, and obtain its immutable repository digest with
`docker image inspect`. Put `name@sha256:<digest>` in generation configuration.
PRISM rejects a mutable image tag unless the configuration explicitly opts out.

## Slurm + Conda

On clusters where Docker is unavailable to users, use the bundled
`generation.targetdiff.slurm.example.yaml`. The Slurm backend executes a blocking
chain of `srun -> conda run -> targetdiff wrapper`; model logs and status are
therefore handled exactly like local backends. Configure an absolute Conda
environment path when login and compute nodes do not share named-environment
discovery.

The executor accepts explicit Slurm partition, account, optional QoS, walltime,
GPU count, CPU count, and memory. It deliberately does not accept arbitrary
`srun` arguments. The smoke-test profile requests one GPU in `debug`, four CPU
cores, 24 GB RAM, and 45 minutes. The model sees its allocated device as
`cuda:0`, while Slurm controls the actual GPU assignment through
`CUDA_VISIBLE_DEVICES`.

The checkpoint is intentionally not copied into the image. Configure it as a
PRISM artifact with its exact SHA256 and mount it read-only at runtime. The
wrapper writes upstream outputs below `<output>/sdf/*.sdf`, which is the pattern
collected by `TargetDiffAdapter`.

The upstream checkpoint is a PyTorch pickle. Treat it as executable input and
verify the recorded SHA256 before every run. PRISM performs that verification
before starting the model process. The wrapper then loads the archive through a
strict pickle allowlist (EasyDict, OrderedDict, and the two PyTorch tensor
reconstruction globals used by the published checkpoint), validates the
checkpoint shape, and rewrites a trusted per-run copy before invoking upstream
TargetDiff. Any other pickle global is rejected.

The upstream README still points to a Google Drive folder, but that folder was
returning HTTP 404 during integration in August 2026. The cluster smoke test
therefore uses the Zenodo model record `10.5281/zenodo.14041881`, whose archive
describes itself as adopted from the official TargetDiff repository. This is a
third-party preservation record rather than an author-authenticated release.
Record both the archive provenance and the extracted checkpoint SHA256 in the
run manifest; replace it with an author-authenticated checkpoint if one becomes
available.

Minimal configured command:

```yaml
command:
  - python
  - /opt/prism_wrappers/targetdiff.py
  - --targetdiff-root
  - /opt/targetdiff
  - --checkpoint
  - "{checkpoint}"
  - --pocket
  - "{pocket}"
  - --output
  - "{output}"
  - --num-samples
  - "{num_samples}"
  - --seed
  - "{seed}"
  - --device
  - "{device}"
```

TargetDiff expects an already cropped pocket PDB. PRISM passes a structure
pocket unchanged; for center, reference-ligand, or residue specifications it
extracts whole protein residues into `models/targetdiff/prepared/pocket.pdb`.
The default extraction radius is 10 Å when the pocket specification does not
provide one.
