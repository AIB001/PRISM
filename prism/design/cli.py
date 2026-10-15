"""Command line entry point for reviewed, local molecular edits."""

import argparse
import json
from pathlib import Path

from .errors import DesignError
from .io import write_json


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="prism design", description="Compile, generate and independently verify local molecular edits"
    )
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("plan-prepare", help="Build a hash-bound molecular evidence packet for an LLM")
    prepare.add_argument("--spec", required=True, type=Path)
    prepare.add_argument("--input-root", type=Path)
    prepare.add_argument("--output", required=True, type=Path)
    goals = prepare.add_mutually_exclusive_group()
    goals.add_argument("--goal", help="Design goal; defaults to the existing DesignSpec prompt")
    goals.add_argument("--goal-file", type=Path, help="UTF-8 design goal")
    plan = commands.add_parser("plan", help="Call a configured LLM or validate an imported edit proposal")
    plan.add_argument("--request", required=True, type=Path)
    plan.add_argument("--output", required=True, type=Path)
    source = plan.add_mutually_exclusive_group(required=True)
    source.add_argument("--llm-config", type=Path)
    source.add_argument("--response", type=Path, help="Offline JSON response; never treated as a live API call")
    check = commands.add_parser("plan-check", help="Replay edits and optionally verify a 3D candidate matches the plan")
    check.add_argument("--plan", required=True, type=Path)
    check.add_argument("--candidate", type=Path)
    check.add_argument("--output", type=Path)
    for name in ("validate", "compile", "verify", "run"):
        sub = commands.add_parser(name)
        sub.add_argument("--spec", required=True, type=Path)
        sub.add_argument("--input-root", type=Path)
        sub.add_argument("--output", type=Path, required=name in {"compile", "run"})
        if name == "validate":
            sub.add_argument("--for-execution", action="store_true")
        if name == "verify":
            sub.add_argument("--candidate", required=True, type=Path, help="Exactly one raw SDF molecule")
        if name == "run":
            sub.add_argument("--backend-config", required=True, type=Path)
            sub.add_argument("--num-samples", type=int, default=16)
            sub.add_argument("--seed", type=int, default=20261013)
    backend = commands.add_parser("backend-check")
    backend.add_argument("--backend-config", required=True, type=Path)
    backend.add_argument("--output", type=Path)
    pilot = commands.add_parser("pilot")
    pilot.add_argument("--queue", required=True, type=Path)
    pilot.add_argument("--backend-config", required=True, type=Path)
    pilot.add_argument("--output", required=True, type=Path)
    pilot.add_argument("--seed", type=int, default=20261013)
    dataset = commands.add_parser("dataset-check", help="Validate portable M1 records and create grouped splits")
    dataset.add_argument("--manifest", required=True, type=Path)
    dataset.add_argument("--split-output", required=True, type=Path)
    dataset.add_argument("--seed", type=int, default=20261013)
    for name in ("train", "preflight"):
        sub = commands.add_parser(name, help="Run M1 in the configured isolated FLOWR environment")
        sub.add_argument("--config", required=True, type=Path)
        if name == "train":
            sub.add_argument("--resume", type=Path)
            sub.add_argument("--stop-after-step", type=int)
        else:
            sub.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "plan-prepare":
            from .planner import prepare_planning

            goal = args.goal_file.read_text(encoding="utf-8-sig") if args.goal_file else args.goal
            result = prepare_planning(args.spec, args.output, goal, args.input_root)
        elif args.command == "plan":
            from .planner import plan_design

            result = plan_design(args.request, args.output, args.llm_config, args.response)
        elif args.command == "plan-check":
            from .planner import check_plan

            result = check_plan(args.plan, args.candidate)
        elif args.command in {"train", "preflight"}:
            import os
            import subprocess
            from .backends import load_backend
            from .io import read_json, resolve_path

            config = read_json(args.config)
            backend = load_backend(resolve_path(args.config.resolve().parent, config["backend_config"]))
            command = [
                backend["python"],
                "-m",
                "prism.design.m1." + args.command,
                "--config",
                str(args.config.resolve()),
            ]
            if args.command == "train":
                if args.resume:
                    command += ["--resume", str(args.resume.resolve())]
                if args.stop_after_step is not None:
                    command += ["--stop-after-step", str(args.stop_after_step)]
            else:
                command += ["--output", str(args.output.resolve())]
            env = os.environ.copy()
            env.pop("PYTHONHOME", None)
            env["PYTHONPATH"] = str(Path(__file__).resolve().parents[2])
            return subprocess.run(command, env=env, check=False).returncode
        elif args.command == "dataset-check":
            from .m1.data import prepare_manifest

            split = prepare_manifest(args.manifest, args.split_output, args.seed)
            result = {k: v for k, v in split.items() if k != "assignments"}
        elif args.command == "validate":
            from .validate import validate_design

            result = validate_design(args.spec, args.input_root, args.for_execution)
        elif args.command == "compile":
            from .compile import compile_design
            from .schema import load_design
            from .validate import ground_design

            compiled = compile_design(
                ground_design(load_design(args.spec, args.input_root), for_execution=True), args.output
            )
            result = {
                "status": "COMPILED",
                "compiled_design": str(args.output / "compiled_design.json"),
                "task": compiled["task"],
            }
        elif args.command == "verify":
            from rdkit import Chem
            from .errors import require
            from .schema import load_design
            from .validate import ground_design
            from .verify import verify_candidate

            supplier = Chem.SDMolSupplier(str(args.candidate), sanitize=False, removeHs=False)
            require(len(supplier) == 1, "--candidate must contain exactly one raw molecule", "CANDIDATE_COUNT")
            result = verify_candidate(ground_design(load_design(args.spec, args.input_root)), supplier[0])
        elif args.command == "backend-check":
            from .backends import check_backend, load_backend

            result = check_backend(load_backend(args.backend_config))
        elif args.command == "run":
            from .runner import run_design

            result = run_design(
                args.spec, args.backend_config, args.output, args.num_samples, args.seed, args.input_root
            )
        else:
            from .runner import run_pilot

            result = run_pilot(args.queue, args.backend_config, args.output, args.seed)
        if args.command in {"validate", "verify", "backend-check", "plan-check"} and args.output:
            write_json(args.output, result)
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
        return 2 if result.get("status") in {"FAILED", "REJECTED"} else 0
    except (DesignError, OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        parser.exit(2, f"prism design: {getattr(exc, 'code', 'INVALID_INPUT')}: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
