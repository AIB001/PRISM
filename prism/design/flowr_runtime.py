"""Standalone bridge executed by FLOWR.root's isolated Python environment.

Uses the upstream loader, interpolant and sampler without modifying their code.
Writes every raw decoded slot before any validity/uniqueness/property filter.
M0 is standalone. Optional M1 imports only PRISM design modules from the source
checkout; the model environment does not need the MD application's dependencies.
"""

import argparse
from contextlib import nullcontext
import hashlib
import json
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace


def write_json(path, value):
    path = Path(path)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def check_file(ref):
    path = Path(ref["path"])
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != ref["sha256"]:
        raise RuntimeError(f"Input changed after compilation: {path}")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compiled", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20261013)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--batch-cost", type=int, default=1)
    parser.add_argument("--grow-size", type=int, default=2)
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--probe-only", action="store_true")
    parser.add_argument("--m1-checkpoint", type=Path)
    parser.add_argument("--m1-instruction-file", type=Path)
    parser.add_argument("--m1-checkpoint-sha256")
    parser.add_argument("--m1-instruction-sha256")
    parser.add_argument("--base-sha256")
    parser.add_argument("--source-commit")
    cli = parser.parse_args()
    compiled = json.loads(cli.compiled.read_text())
    if cli.m1_checkpoint and (compiled["task"] != "fragment_replacement" or not cli.m1_instruction_file):
        raise RuntimeError("Current M1 weights support explicit-text fragment replacement only")
    if cli.m1_checkpoint:
        check_file({"path": str(cli.m1_checkpoint), "sha256": cli.m1_checkpoint_sha256})
        check_file({"path": str(cli.m1_instruction_file), "sha256": cli.m1_instruction_sha256})
        check_file({"path": str(cli.checkpoint), "sha256": cli.base_sha256})
    cli.output.mkdir(parents=True, exist_ok=True)
    reference_path = check_file(compiled["inputs"]["reference_ligand"])
    pocket_path = check_file(compiled["inputs"]["backend_pocket"])

    import numpy as np
    import torch
    from rdkit import Chem, rdBase
    import flowr
    from flowr.data.interpolate import GeometricInterpolant
    import flowr.gen.generate_from_pdb as api

    if cli.source_root and cli.source_root.resolve() not in Path(flowr.__file__).resolve().parents:
        raise RuntimeError("Imported FLOWR module is outside the pinned source checkout")

    reference = Chem.MolFromMolFile(str(reference_path), removeHs=True)
    mode = "substructure_inpainting" if compiled["task"] == "fragment_replacement" else "fragment_growing"
    flags = {
        name: False
        for name in (
            "scaffold_hopping",
            "scaffold_elaboration",
            "linker_inpainting",
            "core_growing",
            "fragment_growing",
            "fragment_inpainting",
            "substructure_inpainting",
            "graph_inpainting",
            "mixed_uncond_inpaint",
        )
    }
    flags.update(
        {mode: True, "inference": True, "substructure": compiled["editable_indices_0based"], "grow_size": cli.grow_size}
    )
    actual_mode, mask, is_local = GeometricInterpolant._determine_modes_and_extract_masks(
        SimpleNamespace(**flags), [SimpleNamespace(seq_length=reference.GetNumAtoms())], [reference], None
    )[0]
    actual_fixed = torch.where(mask)[0].tolist()
    if actual_fixed != sorted(compiled["fixed_indices_0based"]) or actual_mode != mode:
        raise RuntimeError("Native conditional mask does not preserve the requested side")
    unsafe = torch.serialization.get_unsafe_globals_in_checkpoint(str(cli.checkpoint))
    if unsafe:
        raise RuntimeError("Checkpoint contains unreviewed pickle globals: " + repr(unsafe))
    checkpoint = torch.load(cli.checkpoint, map_location="cpu", weights_only=True)
    hparams = checkpoint["hyper_parameters"]
    if not hparams.get("fragment_inpainting") or not hparams.get("remove_hs"):
        raise RuntimeError("Checkpoint does not declare the required heavy-atom fragment training")
    del checkpoint
    probe = {
        "status": "PASS",
        "native_mode": actual_mode,
        "is_local": is_local,
        "fixed_indices_0based": actual_fixed,
        "editable_indices_0based": compiled["editable_indices_0based"],
        "source_module": str(Path(flowr.__file__).resolve()),
        "torch_version": torch.__version__,
        "rdkit_version": rdBase.rdkitVersion,
        "cuda_available": torch.cuda.is_available(),
        "checkpoint_fragment_inpainting": bool(hparams["fragment_inpainting"]),
        "native_final_inpaint": True,
        "valence_repair": False,
        "exact_attachment_enforcement": "PRISM independent verifier; growth location is a soft native prior",
    }
    write_json(cli.output / "capability_probe.json", probe)
    if cli.probe_only:
        print(json.dumps(probe))
        return 0
    if cli.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")

    argv = [
        "flowr-root",
        "--pdb_file",
        str(pocket_path),
        "--ligand_file",
        str(reference_path),
        "--arch",
        "pocket",
        "--pocket_type",
        "holo",
        "--ckpt_path",
        str(cli.checkpoint),
        "--save_dir",
        str(cli.output),
        "--gpus",
        "1" if cli.device == "cuda" else "0",
        "--num_workers",
        "0",
        "--batch_cost",
        str(cli.batch_cost),
        "--max_sample_iter",
        "0",
        "--sample_n_molecules_per_target",
        str(cli.samples),
        "--seed",
        str(cli.seed),
        "--integration_steps",
        str(cli.steps),
        "--final_inpaint",
        "--no_ligand_valence_repair",
    ]
    if compiled["task"] == "fragment_replacement":
        argv += ["--substructure_inpainting", "--substructure"] + list(map(str, compiled["editable_indices_0based"]))
    else:
        lo, hi = compiled["new_fragment_heavy_atoms_range"]
        if not lo <= cli.grow_size <= hi:
            raise RuntimeError("Requested grow size exceeds DesignSpec range")
        prior = check_file(compiled["growth_prior"]["file"])
        argv += ["--fragment_growing", "--grow_size", str(cli.grow_size), "--prior_center_file", str(prior)]
    original_argv = sys.argv
    sys.argv = argv
    try:
        args = api.get_args()
    finally:
        sys.argv = original_argv
    write_json(cli.output / "upstream_arguments.json", vars(args))
    # This is the narrow inference path from upstream evaluate(): no repair,
    # minimization, final sanitization filter, reference-PDB rendering or retry.
    torch.set_float32_matmul_precision("high")
    loaded = api.load_model(args)
    model, hparams, vocab, charges, hybridization, aromatic, _, _ = loaded
    model = model.to(cli.device).eval()
    transform, interpolant = api.util.load_util(args, hparams, vocab, charges, hybridization, aromatic)
    text_session = None
    if cli.m1_checkpoint:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from prism.design.m1.inference import TextConditionedSampling

        text_session = TextConditionedSampling(
            model,
            cli.m1_checkpoint,
            cli.m1_instruction_file.read_text(encoding="utf-8").strip(),
            cli.device,
            cli.base_sha256,
            cli.source_commit,
            interpolant,
            compiled["attachment_index_0based"],
            len(compiled["reference_atom_map_to_backend_index_0based"]),
        )
        write_json(cli.output / "m1_conditioning.json", text_session.metadata)
    system = api.util.load_data_from_pdb(
        args,
        remove_hs=hparams["remove_hs"],
        remove_aromaticity=hparams["remove_aromaticity"],
        ligand_idx=0,
        chain_id=None,
        canonicalize_conformer=False,
    )
    # Verify index and frame preservation through actual upstream preprocessing.
    system_mol = system.ligand.to_rdkit(vocab, charges, hybridization, sanitise=True)
    if text_session and system_mol is not None:
        Chem.AssignStereochemistryFrom3D(system_mol, replaceExistingTags=True)
    if system_mol is None or [a.GetAtomicNum() for a in system_mol.GetAtoms()] != [
        a.GetAtomicNum() for a in reference.GetAtoms()
    ]:
        raise RuntimeError("Upstream preprocessing changed the reference atom order/elements")
    if not np.allclose(system_mol.GetConformer().GetPositions(), reference.GetConformer().GetPositions(), atol=1e-3):
        raise RuntimeError("Upstream preprocessing changed reference coordinates/order")
    if Chem.MolToSmiles(system_mol) != Chem.MolToSmiles(reference):
        raise RuntimeError("Upstream preprocessing changed reference chemistry")
    probe["real_preprocessing_atom_order_and_frame"] = "PASS"
    write_json(cli.output / "capability_probe.json", probe)
    dataset = api.get_dataset(system, transform, vocab, interpolant, args, hparams)
    loader = api.util.get_dataloader(args, dataset, interpolant, iter=0)
    guidance = api.util.get_guidance_params(args)
    raw_dir = cli.output / "raw_slots"
    raw_dir.mkdir(exist_ok=True)
    events = []
    start = time.monotonic()
    try:
        with (cli.output / "attempts.jsonl").open("w", encoding="utf-8") as ledger:
            for batch_index, batch in enumerate(loader):
                prior, data, _, _ = batch
                with text_session.condition() if text_session else nullcontext():
                    sampled = api.generate_ligands_per_target(
                        args,
                        model,
                        prior=prior,
                        posterior=data,
                        pocket_noise="fix",
                        save_traj=False,
                        iter=f"0_{batch_index}",
                        guidance_params=guidance,
                        device=cli.device,
                    )
                for mol in sampled:
                    slot = len(events) + 1
                    event = {
                        "slot": slot,
                        "batch_index": batch_index,
                        "status": "DECODE_FAILED" if mol is None else "DECODED",
                        "sdf": None,
                        "error": None,
                    }
                    if mol is not None:
                        target = raw_dir / f"slot_{slot:04d}.sdf"
                        try:
                            writer = Chem.SDWriter(str(target))
                            writer.SetKekulize(False)
                            writer.write(mol)
                            writer.close()
                            event["sdf"] = str(target)
                        except Exception as exc:
                            event["status"] = "SERIALIZE_FAILED"
                            event["error"] = str(exc)
                    events.append(event)
                    ledger.write(json.dumps(event) + "\n")
                    ledger.flush()
                print(f"PRISM raw slots completed: {len(events)}/{cli.samples}", flush=True)
        if len(events) != cli.samples:
            raise RuntimeError(f"Raw slot accounting mismatch: {len(events)} != {cli.samples}")
    finally:
        write_json(
            cli.output / "runtime_summary.json",
            {
                "requested": cli.samples,
                "raw_slots": len(events),
                "decoded": sum(e["status"] == "DECODED" for e in events),
                "decode_failed": sum(e["status"] == "DECODE_FAILED" for e in events),
                "serialize_failed": sum(e["status"] == "SERIALIZE_FAILED" for e in events),
                "wall_seconds_sampling": time.monotonic() - start,
                "retries": 0,
                "internal_validity_or_uniqueness_filters": False,
                "chemical_repair": False,
                "native_final_inpaint": True,
                "repair_stats": api.util.repair_stats_summary(model),
                "peak_cuda_allocated_bytes": torch.cuda.max_memory_allocated() if cli.device == "cuda" else None,
            },
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
