# PRISM-Generation Module Instructions

Guidance for working on structure-based ligand generation in PRISM.

This file records the decisions that are **not** recoverable from the code, the
measurements behind them, and the traps that have already been paid for. Read
the invariants section before changing anything in this module.

Last substantive update: 2026-08-24.

## Module Overview

Orchestration layer for six structure-based generative models. The model
environments, source trees and checkpoints deliberately live **outside** the
PRISM Python environment and are reached through configured Conda, Docker or
Slurm backends.

| Model | Mode | Conditioning ligand | Notes |
|---|---|---|---|
| Pocket2Mol | direct | optional (locates pocket) | autoregressive, see the sampling trap below |
| TargetDiff | direct | optional | drops its own reconstruction failures before writing SDF |
| PocketXMol | direct | optional | most stable delivery rate |
| MolCRAFT | reference_guided | required (.sdf) | |
| FLOWR | reference_guided | required (.sdf/.pdb) | CUDA-only, ~24 GB GPU |
| DiffSBDD | reference_guided | required (.sdf) | only model producing valence errors (~3%) |

**Architecture**

- `types.py` — serializable `GenerationRequest` / `PocketSpec` / `ModelRunResult` /
  `CandidateRecord` / `QuarantineRecord`
- `registry.py` — model IDs, `all` expansion, capability metadata
- `adapters/` — per-model capability contracts and command placeholders
- `wrappers/` — the actual call into each upstream entry point
- `execution.py` — command rendering, artifact verification, process control
- `weights.py` — checkpoint inventory, discovery, preflight, download
- `pocket.py` / `protein.py` — pocket normalization, pocket PDB extraction
- `postprocess.py` — record splitting, QC invocation, quarantine, reports
- `quality.py` — the quality-control engine
- `runner.py` — orchestration, resume/overwrite, manifest
- `handoff.py` — selects MD-ready ligands from a finished run (see below)

**Entry points**

```bash
prism generate --model diffsbdd --protein p.pdb --pocket ref.sdf -o out \
  --generation-config cfg.yaml --num-samples 100 --qc full
prism --generation true --model ...          # compatibility form
python -m prism.generation ...               # no console script needed

prism prepare-md out -o md_inputs            # export MD-ready ligands
prism prepare-md out -o md_inputs --only-pass --limit 10
```

MCP tools: `list_generation_models`, `generate_ligands`, `check_generated_ligands`,
`prepare_generated_for_md`.

---

## Invariants — do not undo these without re-measuring

### 1. Detect, propose, never silently repair

`quality.py` splits transformations into *apply* and *detect*:

- **apply**: only things that cannot change the molecular graph or the pose —
  SDF layout normalization and hydrogen addition (written to a separate file).
- **detect**: everything else, **including speculative valence repair**.

When a molecule fails sanitization, PRISM searches for minimal edits, writes
every equally-ranked variant to `quarantine/`, and applies none of them. A bond
order that sanitizes is not evidence that it is the bond order the model
intended: one real case had three equally minimal tautomers and a 3D SDF cannot
distinguish them.

`RepairProposal.to_dict()` carries `"applied": False` so every consumer —
manifest, CSV, MCP response — sees it.

### 2. Thresholds are calibrated against experimentally determined ligands

Every threshold in `DEFAULT_THRESHOLDS` is a claim about where correct
chemistry ends. `tests/test_generation_calibration.py` holds the answer key:
crystal ligands must pass with zero flags, marketed drugs must not be called
reactive.

Measured on 30 organic crystal ligands: worst bond angle deviation median
14.3°, max 26.1°; aromatic planarity RMS < 0.10 Å. The 35° warning threshold
sits above the observed maximum on purpose.

**The calibration set must cover the upper tail.** A first version used only
well-behaved ligands (max 11°) and did not notice the angle threshold being
tightened to 20°. 2XNB (26.1°) and 1D3P (0.035 Å) were added for this reason.

### 3. Reactive-group SMARTS must not flag marketed drugs

The original patterns flagged 5 of 20 approved drugs: `[CX3;!R]=[NX2;!R]` read
amidines and guanidines as Schiff bases (benzamidine, metformin), `[OX2][NX3]`
read hydroxamic acids as labile N–O (vorinostat), `[NX3][NX3]` read hydrazides
as free hydrazines (isoniazid). False-positive rate went 25% → 5% after
tightening; the residual is hydralazine, a genuine free arylhydrazine.

These are **advisory**. Whether a reactive group disqualifies a candidate is a
project policy decision, not a correctness question, so they raise
`WARN_PLAUSIBILITY` and never withhold a molecule.

### 4. Checkpoint allowlisting cannot use `torch.load(pickle_module=...)`

torch < 2.0 discards the argument:

```python
if weights_only:
    ...
else:
    pickle_module = pickle      # torch 1.13.1 serialization.py
```

`weights_only` defaults to False, so on torch 1.13 — five of the six model
environments — the custom `find_class` is **never called**. Verified: 0 calls,
and a poisoned checkpoint executed its payload. `weights_only=True` is not an
alternative: every published checkpoint here carries a non-tensor config object
(EasyDict, argparse.Namespace) and refuses to load.

`wrappers/_checkpoint.py` therefore walks the pickle opcode stream with
`pickletools.genops` **before torch touches the file**. Nothing executes until
the allowlist passes. Real checkpoints reference only 3–4 benign globals.

### 5. Wrappers must run as standalone scripts

Every generation config invokes a wrapper by path
(`python /path/wrappers/xxx.py`) inside the *model's* environment, which does
not have PRISM installed. A bare `from . import _device` breaks all six models
at once, and package-import unit tests cannot see it.

Both shared modules are imported with a fallback:

```python
try:
    from . import _checkpoint
    from . import _device as _device_support
except ImportError:                      # standalone script invocation
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import _checkpoint
    import _device as _device_support
```

`tests/test_generation.py::test_every_wrapper_runs_as_a_standalone_script`
executes each wrapper the way production does. Keep it.

### 6. `--device auto` must probe the host

`auto` previously returned `cuda:0` unconditionally, which fails on a CPU-only
node and violates the project rule against hardcoding GPU ids. `wrappers/_device.py`
now probes with `nvidia-smi -L`, honors `CUDA_VISIBLE_DEVICES` (including `-1`
and the empty string), and treats a failed probe as "no GPU".

An **explicit** `cuda` / `cuda:N` request is passed through unchecked: some CUDA
containers ship without `nvidia-smi`, and refusing a deliberate request on a
failed probe would be worse.

Measured CPU support: PocketXMol ✅, DiffSBDD ✅, TargetDiff ✅,
Pocket2Mol ❌ (`models/maskfill.py` hardcodes `.to('cuda')`),
MolCRAFT ❌ (Lightning Trainer built without an accelerator setting),
FLOWR ❌ (by design). The three without a CPU path refuse in argparse with the
reason, rather than failing minutes later with a CUDA traceback.

### 7. Timeouts must kill the process group

The direct child is only a launcher (`conda run`, `srun`, `docker run`); the
model runs as a grandchild. `subprocess.run(timeout=...)` kills the launcher
only — an abandoned AM1-BCC job was observed holding a core for nearly two
hours after its task had already been recorded as timed out.

`execution.py` uses `Popen(start_new_session=True)` and `os.killpg` on timeout
**and on any other exception**, so a cancelled run does not leak either.

Known gap: the Docker backend still cannot stop a container this way (it runs
under the daemon, not as a child). Would need a deterministic `--name` plus
`docker kill`.

### 8. Natural sort, everywhere output files are numbered

Pocket2Mol writes `0.sdf` … `N.sdf` without zero padding. Lexicographic order
puts 10 before 2, and truncating that order at `--num-samples` silently drops a
contiguous middle block. Fixed in `wrappers/pocket2mol.py` and in
`postprocess.collect_candidates`.

### 9. The MD hand-off verifies hydrogens on the file, not in the manifest

`handoff.py` is the only supported path from a finished run to system building.
It re-reads every hydrogenated file it is about to hand over and counts the
hydrogens there, rather than trusting the manifest's record that hydrogenation
succeeded. A run directory can be copied, pruned, or partially regenerated
between the two steps, and a missing hydrogen is invisible downstream — see the
MD handoff measurements below.

Two deliberate constraints:

- **It does not import the builder.** It writes the ligands and emits the
  `prism -pf ... -lf ...` command that consumes them. Building is long-running
  with its own failure modes; invoking it here would make a parameterization
  error surface as a generation error, and would couple the two modules.
- **It reports every exclusion.** `skipped.csv` carries a reason code per
  withheld candidate (`FAILED_QC`, `NO_EXPLICIT_HYDROGENS`, `NO_HYDROGEN_FILE`,
  `NOT_PASS`, `OVER_LIMIT`, `MODEL_NOT_SELECTED`, `UNREADABLE_LIGAND`). A
  hand-off that quietly delivered fewer ligands than the run produced would
  read as "these are the candidates" rather than "these are the survivors".

`--limit` truncates a QC-ranked order, not collection order: grade predicts
build success, so an arbitrary slice would discard good molecules first.

### 10. Record splitting must not strip the terminating blank line

`_sdf_records` used `fragment.rstrip("\r\n")`, which deleted the blank line that
terminates a trailing SD data field. That manufactured a whole class of
"SDF format defect" that was not in the source files — an earlier report
claiming FLOWR emits malformed SDF was wrong for this reason. Guarded by
`test_record_splitting_preserves_a_well_formed_trailing_data_field`.

---

## The Pocket2Mol sampling trap

Upstream's loop is `while len(pool.finished) < num_samples`, and `global_step`
counts atoms added. A molecule finishing at step k has roughly k heavy atoms,
so **stopping when enough molecules have accumulated stops at small sizes**.

Measured on 1IEP (reference ligand 37 heavy atoms):

| | request 10 | request 100 | run to max_steps |
|---|---:|---:|---:|
| heavy-atom median | 7.5 | 18 | 31 |
| within ±40% of reference | 0% | — | **84%** |

Completion order vs size: Pearson r = 0.956, monotonically non-decreasing
across all 100 molecules. The run reached step 22 of 50 and the largest
molecule had exactly 22 heavy atoms.

**Fix in place**: the wrapper passes `_UNBOUNDED_SAMPLES` so the stopping
condition never fires, lets `max_steps` bound the size, and selects
`--num-samples` afterwards via `--select {stratified,largest,target,completion}`.
`sdf_all/` keeps everything; `selection.json` records the full distribution.

**Semantic difference to remember**: for Pocket2Mol, `--num-samples` selects
from a completed run. For the other five it means "generate this many".

**Upstream bug this exposed**: running to `max_steps` drains the beam queue,
and `np.random.choice` then raises `probabilities do not sum to 1`. Upstream
only guards its loop against `KeyboardInterrupt`, so the process dies before
its own SDF writer and every finished molecule is lost. The wrapper recovers
from the per-step `samples_N.pt` checkpoint (needs the Pocket2Mol root on
`sys.path`, or unpickling fails with `No module named 'utils'`).

---

## Measured findings — do not re-derive these

Scale: 3204 molecules, 5 binding sites (1IEP Abl kinase, 1HPV HIV-1 protease,
1OYT thrombin, 3ERT estrogen receptor, 4EY7 acetylcholinesterase), 6 models,
100 samples each, single fixed seed.

**Defect rates**

| Defect | Rate | Repairable |
|---|---:|---|
| reactive groups | 23.3% | no — human judgment |
| bond angle > 35° | 18.0% | yes — restrained minimization → 1.7% |
| unassignable stereocentre | 14.5% | yes → 0.1% |
| pocket centre offset | 9.2% | no — re-dock |
| protein hard clash | 7.1% | **only with the protein present** |
| aromatic non-planar | 3.2% | yes → 0% |
| bond length outlier | 1.2% | yes → 0% |
| valence failure | 0.5% | proposal only |
| integrity defect | **0%** | — |

**Hydrogens**: no model emits explicit hydrogens (0/3204). 3188/3188 hydrogenate
successfully with **heavy-atom displacement exactly 0.000000 Å**.

**Restrained MMFF94s minimization** (heavy atoms tethered, k = 5.0) fixes the
intramolecular geometry at a median 0.24 Å cost. **Free minimization is
catastrophic**: protein clashes go 17.5% → 73.9% because the force field does
not know the protein exists. Ligand-only restrained minimization also makes
clashes *worse* (17.5% → 21.0%); those need the complex.

**MD handoff** (82 builds, gaff2):

| Input | build returns success | topology hydrogens correct |
|---|---:|---:|
| raw model output | 19/41 | **0/41** |
| hydrogenated | 35/41 | **35/41 (85%)** |

The raw arm's 19 "successes" produced topologies with **zero hydrogens** —
correct file layout, wrong chemistry, no warning. This included a crystal
ligand that QC rated PASS.

Closed by `handoff.py`, which exports only hydrogen-verified files. The same
accounting is available to the builder as `require_explicit_hydrogens()` and
`verify_topology_hydrogens()` in `prism/forcefield/base.py`; those are additions
to the shared base class and nothing calls them yet, so whether GAFF should
refuse such a ligand outright stays the builder's decision. Verified on 36 real
candidates that the two agree: exported files report 0 missing hydrogens, raw
model output 25–36 missing each.

QC grade predicts build success: PASS → 100%, angle < 35° → 100%,
angle > 35° → 73%, quarantined originals → 0%, **repair variants → 100%**.

**Reactive groups do NOT break GAFF2 parameterization** (12/12 succeeded,
allene and Schiff bases included). An earlier contrary impression came from
concurrency-induced timeouts, not chemistry.

---

## Cluster deployment

Host details (address, port, account) are deployment-specific and belong in your
own notes or SSH config, not in this repository. What matters here is the
layout, which the wrappers and the weights manifest both assume.

Clusters typically have **no outbound internet**: run `prism weights download`
on a machine that does, then transfer `$PRISM_MODELS_DIR` across, or serve it
with `PRISM_MODELS_MIRROR=file:///path/to/mirror`.

| Item | Path (under the deployment root) |
|---|---|
| Model sources + checkpoints | `<root>/prism-models/<model>/` — point `PRISM_MODELS_DIR` here |
| Conda environments | `<root>/conda-envs/prism-gen-<model>/` |
| DiffSBDD | reuses `prism-gen-molcraft` (zero extra install; stubs `wandb`, `imageio`) |
| Benchmark configs | `<root>/prism-benchmarks/<run-id>/config/` |

**Scheduling notes**

- `--cpus-per-task=1..2`. GPUs are often free while CPUs are exhausted; asking
  for 4–8 leaves jobs pending on `(Resources)` indefinitely.
- Request a walltime close to the real runtime so backfill can schedule it.
- `PYTORCH_CUDA_ALLOC_CONF=expandable_segments` is **only** valid for the FLOWR
  environment (torch 2.5.1); the torch 1.13 environments abort on it.
- FLOWR OOMs at `--batch-cost 100` on a 24 GB card; use 20.

Approximate runtimes for 100 samples: PocketXMol ~2 min, DiffSBDD ~3.5 min,
MolCRAFT ~5 min, FLOWR ~6 min, TargetDiff ~25 min, Pocket2Mol ~4 min at
beam 50 (2h46m at beam 300 running to completion).

---

## Weights distribution

The six checkpoints total 966 MiB. They are not in the package, not in the
repository, and not in any wheel — `prism/data/model_weights/manifest.json` is
the only thing that ships, and it is the single source of truth for what a
model needs.

| Concern | Where it lives |
|---|---|
| What each model needs | `manifest.json` — relpath, sha256, size, licence, upstream URL |
| Where it goes | `PRISM_MODELS_DIR` → `$XDG_CACHE_HOME/prism/models` → `~/.cache/prism/models` |
| Where it comes from | the artifact's one acquisition route (below) |
| User-facing commands | `prism weights {list,path,download,verify}` |

This mirrors `prism/membrane/patch.py`, which solved the same problem for the
bilayer patches: env-var override, package-relative default, and a descriptive
not-available error instead of a `FileNotFoundError`.

### Three routes, one per artifact

`manifest.json` is `"version": 2`. Every artifact declares **exactly one** of:

| Route | Field(s) | Used for |
|---|---|---|
| PRISM mirror | `mirror_asset` | TargetDiff, Pocket2Mol, MolCRAFT |
| publisher | `download_url` | FLOWR, DiffSBDD |
| publisher bundle | `archive` + `archive_member` | PocketXMol (both artifacts) |

No route at all means the user places the file by hand, and `upstream_url` /
`upstream_instructions` say how. Two routes on one artifact raises at manifest
load: it is not a richer manifest, it is an unanswered question about which one
is authoritative.

**The split is about the host, not the licence.** All six licences permit
redistribution. TargetDiff, Pocket2Mol and MolCRAFT publish through Google
Drive, which cannot be scripted — virus-scan interstitials on large files, a
per-file quota that locks out for a day once tripped, and no access at all from
some networks. Those three are mirrored. The other three are on Zenodo, which
serves stable direct URLs, so mirroring them would only add a copy to keep in
step with upstream. Mirroring is the fallback, not the default; if a publisher
moves to a scriptable host, its asset should leave the mirror.

An earlier version mirrored all seven artifacts (966 MiB). The hybrid cuts the
upload to 117 MiB — the three unscriptable checkpoints happen to be the three
smallest. The measured counter-argument, recorded because it is a real cost:
on the development network Zenodo delivered 11.5–16.3 KB/s and dropped
connections mid-transfer, against 28.4–47.8 KB/s from the GitHub release CDN.
Where bandwidth to Zenodo is poor, `PRISM_MODELS_MIRROR` pointed at a local or
institutional copy is the answer, not widening the mirror.

**PocketXMol costs more to fetch than it occupies.** Zenodo publishes every
trained model in one 611 MiB bundle; PRISM wants 232 MiB of it. Both members
are taken in a single pass, each verified against its own SHA256, and the
bundle is deleted once nothing needs it. The CLI states the bundle size before
starting, because watching 611 MiB arrive for a "232 MiB missing" artifact
otherwise reads as a bug.

### Rules that are enforced, not documented

1. *Only redistributable weights may carry a mirror asset.* Adding a
   `mirror_asset` to a model whose `redistributable` is false raises at
   **manifest load**, so the mistake surfaces on any command that touches
   weights rather than at upload time. A non-redistributable model must instead
   carry an `upstream_url` or `upstream_instructions`. Asserted in
   `tests/test_generation_weights.py` against a synthetic manifest, so the rule
   keeps being tested after the licence facts change.

   Mirroring is redistribution, not relicensing — each artifact keeps its
   upstream terms, and `attribution` is printed on download:

   | Weights | Licence | Note |
   |---|---|---|
   | TargetDiff, Pocket2Mol | MIT | TargetDiff's licence file is spelled `LICIENCE` upstream, so GitHub's detector reports "no licence" — it is MIT (Copyright (c) 2023 Jiaqi Guan) |
   | PocketXMol, FLOWR, DiffSBDD | CC BY 4.0 | all three repos are MIT, but the **weights** are published on Zenodo under CC BY 4.0 (records 17801271, 15737419, 8183747); the Zenodo licence governs the checkpoint |
   | MolCRAFT | CC BY-NC-SA 4.0 | NonCommercial and ShareAlike — a PRISM mirror of it may not be used commercially |

   The code/weights split matters: a model's repository licence does not
   automatically cover a checkpoint hosted elsewhere. Check the host of the
   file you are actually redistributing.
2. *A download is installed only after it verifies.* Bytes stream into
   `<target>.part` and are moved into place only when both the SHA256 and the
   byte count match. Any failure unlinks the partial file, so a tampered mirror
   cannot replace a checkpoint that is already good. Download URLs — mirror,
   publisher and bundle alike — are restricted to `https` and `file` (plus
   `http` on loopback, for tests), and that check runs at manifest load as well
   as at download time.

   The digest is computed by **re-reading the finished file**, not by hashing
   blocks as they stream. Across a resumed transfer an incremental digest and
   the bytes on disk can disagree, and the digest is the one thing here that has
   to be trustworthy.
3. *A dropped transfer resumes; it does not restart.* Verified against the real
   thing: a DiffSBDD download died after 2.8 of 17.9 MB at 14.5 KB/s, and the
   hash check correctly discarded the truncated file — but the next attempt
   started from zero, which on a 600 MiB checkpoint means never finishing. Each
   attempt now sends a `Range` header for what is already on disk, with a budget
   of `_DOWNLOAD_ATTEMPTS` tries; a server that ignores the range forces a clean
   restart rather than splicing two copies together.

   `http.client.IncompleteRead` is caught explicitly. A server that announces a
   Content-Length and then hangs up raises it, it is **not** an `OSError`, and
   it is exactly the failure the retry loop exists to survive — before it was
   named, that case escaped as a traceback.

   Two constants make the difference between "slow" and "hung", and both were
   found the hard way against Zenodo at ~20 KB/s:

   - `_STREAM_BLOCK` is 256 KiB, not the 4 MiB used for reading local files.
     `response.read(n)` blocks until it has `n` bytes, so a 4 MiB socket block
     on a 20 KB/s link means three minutes between writes: the progress bar
     sits still, the `.part` file stays at zero, and the download looks dead
     when it is working perfectly. A dropped connection also loses the whole
     unwritten block.
   - `_SOCKET_TIMEOUT` is 60 s per read. `urlopen` has **no** timeout by
     default, so a connection that opens and then delivers nothing blocks
     forever — no retry, no error, no output. It is a per-read timeout, so a
     slow-but-progressing transfer is never cut off.

**Preflight covers the dry run too.** `execution._artifacts()` verifies
checkpoints, but it is never reached with `--dry-run`, so a validation pass
used to report PLANNED for a model whose weights did not exist. `runner.run()`
now calls `weights.preflight()` on both paths, and a missing checkpoint is
reported as `MODEL_WEIGHTS_MISSING` in `status.json`.

**Config placeholders are deliberately narrow.** `load_generation_config()`
expands exactly two names — `${PRISM_MODELS_DIR}` and `${PRISM_WRAPPERS_DIR}` —
and not `os.path.expandvars`. Every string in a generation config can end up as
a command-line argument, so general environment expansion would let a config
file exfiltrate the environment it runs in.

`prism/configs/generation.local.example.yaml` is the portable six-model config
built on those placeholders; a test asserts its paths and hashes agree with the
manifest, so the two cannot drift.

`PRISM_MODELS_DIR` pointed at an existing `prism-models/<model>/` tree — the
cluster layout above — resolves all seven artifacts with no file moves. That
was a layout constraint on the manifest, not a coincidence.

**Where the mirror lives.** `AIB002/PRISM_Model_Weights`, as GitHub **release
assets** rather than files in a git tree. Release assets allow 2 GiB each and
are served over plain HTTPS with no authentication; a git tree rejects any file
over 100 MiB on push, and a repository that once held a checkpoint stays large
forever. This is why `mirror_asset` is a flat file name rather than a path — a
release is a flat namespace, which is also why the names are prefixed with the
model (`molcraft-molcraft.ckpt`) and why a test asserts they are unique.

Git LFS was rejected: the free tier is 1 GiB of storage and 1 GiB of bandwidth
per month, so even the 117 MiB mirror would allow roughly nine full downloads
per month, and the 966 MiB one would have exhausted storage on upload.

Uploading is a REST call, so an SSH key is not enough — `scripts/upload_release.sh`
in that repository needs a token with `contents: write`. It verifies every file
against `SHA256SUMS` before uploading anything: a wrong byte in the mirror
becomes a failed hash check on every user's machine, and the mirror is the one
place where that cannot be caught later.

---

## Testing

```bash
pytest tests/test_generation.py            # orchestration, wrappers, execution
pytest tests/test_generation_weights.py    # manifest, discovery, preflight, download
pytest tests/test_generation_quality.py    # QC engine, quarantine, repair proposals
pytest tests/test_generation_calibration.py  # threshold answer keys
pytest tests/test_mcp_generation.py        # MCP tool contracts
pytest tests/test_generation_handoff.py    # MD hand-off selection and hydrogen guarantee
pytest tests/test_ligand_hydrogens.py      # ligand hydrogen validation
```

Several guards deliberately exercise the **real** execution path because the
corresponding bugs were invisible to package-level unit tests: standalone
wrapper execution, real poisoned checkpoints, real grandchild processes, and
crystal ligands plus marketed drugs as calibration inputs. If one of them looks
slow or awkward, that is why.

Pre-existing failures unrelated to this module: `tests/gxf/FEP/unit_test/42-38/test_gaff2.py`
and `tests/test_rtf_forcefield.py` (its input `24.rtf` is not tracked in git).
`tests/test_javascript_code.py` was fixed upstream in 5967cd5.

---

## Open items

1. **Restrained geometry minimization is measured but not implemented.** Should
   follow the existing convention: explicit `--repair-geometry` switch, off by
   default, separate output file, original untouched.
2. **A size check relative to the reference ligand.** Pocket2Mol's original
   output scored full marks on every chemical metric while being unusable
   (benzene against a 37-heavy-atom reference). Would be `WARN_SIZE`, advisory.
3. **TargetDiff returns rc=0 when zero molecules survive reconstruction.**
   `collect_candidates` catches it via `NO_VALID_MOLECULES`, but the wrapper's
   own exit code does not reflect the outcome.
4. **Docker backend cannot stop its container on timeout** (see invariant 7).
5. **README generation section** does not yet cover DiffSBDD, the Pocket2Mol
   selection semantics, the QC layer, or `prism prepare-md`.
6. **MolCRAFT CPU support may be recoverable** by injecting
   `accelerator="cpu"` into its Lightning Trainer from the launcher. Untested.
7. **DiffSBDD's `model_commit` is unpinned** — empty in the configs, `null` in
   the manifest. The other five are pinned to the commits the wrappers were
   written against.
8. **FLOWR's upstream repo has no committed LICENSE file.** Its README states
   MIT, and the checkpoint on Zenodo is explicitly CC BY 4.0 (which is what
   PRISM mirrors under), but the source tree itself is only licensed by that
   README sentence.

---

## Artifacts

Experiment outputs live under `.prism-build/` (gitignored, outside the package):

| Directory | Contents |
|---|---|
| `multisite-20260818/` | 5-site × 6-model run, `records.json`, `defect_cases.png` (14 annotated real defects), analysis scripts |
| `five-model-audit-1iep-20260817/` | 1IEP run, Pocket2Mol before/after comparison |
| `five-model-audit-8c7y-20260816/` | first audit; the only site so far with PocketXMol valence failures |
| `validate_thresholds/` | 33 crystal ligands used to calibrate thresholds |
| `agent-md-handoff/` | 82 MD builds, per-molecule logs |
| `agent-security/`, `agent-cpu-inference/`, `agent-pipeline/` | robustness test reports |
