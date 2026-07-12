#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GROMACS simulator for PRISM.

Generates and runs a stage-aware, hardware-adaptive workflow (EM -> NVT -> NPT ->
production). It honors the requested `stages` subset and auto-detects the thread
count and GPU availability at runtime — nothing about the host machine is
hardcoded. Each stage resumes from a checkpoint if one is present.
"""

import os
import subprocess


# Per-stage grompp/mdrun recipe. {flags} is the hardware-adaptive mdrun flag set.
_STAGE_RECIPE = {
    "em": dict(prev="solv_ions.gro", mdp="em.mdp", restr="solv_ions.gro", deffnm="em/em", use_t=False),
    "nvt": dict(prev="em/em.gro", mdp="nvt.mdp", restr="em/em.gro", deffnm="nvt/nvt", use_t=False),
    "npt": dict(prev="nvt/nvt.gro", mdp="npt.mdp", restr="nvt/nvt.gro", deffnm="npt/npt", use_t="nvt/nvt.cpt"),
    "prod": dict(prev="npt/npt.gro", mdp="md.mdp", restr="npt/npt.gro", deffnm="prod/md", use_t=False),
}


class GMXSimulator:
    """Run a stage-aware, hardware-adaptive GROMACS MD workflow."""

    def __init__(self, gmx_dir, system_files, ff_dir):
        self.gmx_dir = gmx_dir
        self.system_files = system_files
        self.ff_dir = ff_dir
        self.gmx = self._check_gromacs()

    def _check_gromacs(self):
        for cmd in ("gmx", "gmx_mpi"):
            try:
                r = subprocess.run([cmd, "--version"], capture_output=True, text=True, check=False, timeout=30)
                if r.returncode == 0:
                    for line in r.stdout.split("\n"):
                        if "GROMACS version" in line:
                            print(f"Found {line.strip()}")
                            break
                    return cmd
            except (FileNotFoundError, subprocess.TimeoutExpired):
                continue
        raise RuntimeError("GROMACS not found on PATH. Please install GROMACS (provides `gmx`).")

    # ------------------------------------------------------------------ #
    @staticmethod
    def _detect_threads():
        """A sensible default thread count for a single job (all logical cores)."""
        return max(1, os.cpu_count() or 1)

    @staticmethod
    def _gpu_available():
        try:
            r = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, check=False, timeout=10)
            return r.returncode == 0 and "GPU 0" in r.stdout
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    def _mdrun_flags(self, gpu_id, ntomp, ntmpi):
        flags = f"-ntmpi {ntmpi} -ntomp {ntomp}"
        if gpu_id is not None:
            flags += f" -nb gpu -bonded gpu -pme gpu -gpu_id {gpu_id}"
        return flags + " -v"

    def run(self, stages=None, gpu_id="auto", ntomp=None, ntmpi=1, continue_from=None, maxwarn=10, **kwargs):
        """Run the GROMACS workflow.

        Parameters
        ----------
        stages : list[str]    subset of ['em','nvt','npt','prod'] (default: all)
        gpu_id : int|str|None 'auto' detects a GPU; an int forces a device id; None forces CPU-only
        ntomp : int|None      OpenMP threads; None auto-detects (all logical cores)
        ntmpi : int           thread-MPI ranks (default 1)
        """
        stages = stages or ["em", "nvt", "npt", "prod"]
        # validate/order stages
        order = ["em", "nvt", "npt", "prod"]
        stages = [s for s in order if s in stages]
        if not stages:
            raise ValueError("No valid stages requested. Choose from em/nvt/npt/prod.")

        if ntomp is None:
            ntomp = self._detect_threads()
        if gpu_id == "auto":
            gpu_id = 0 if self._gpu_available() else None
        if kwargs:
            print(f"Ignoring extra parameters: {list(kwargs)}")

        flags = self._mdrun_flags(gpu_id, ntomp, ntmpi)
        print(f"GROMACS run: stages={stages}, threads={ntomp}, "
              f"{'GPU '+str(gpu_id) if gpu_id is not None else 'CPU-only'}")

        script_path = os.path.join(self.gmx_dir, "localrun.sh")
        self._write_script(script_path, stages, flags, maxwarn)
        os.chmod(script_path, 0o755)

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(ntomp)
        original = os.getcwd()
        os.chdir(self.gmx_dir)
        try:
            subprocess.run(["bash", "localrun.sh"], check=True, text=True, env=env)
            return self._collect_outputs(stages)
        except subprocess.CalledProcessError as e:
            print(f"Error running GROMACS simulation: {e}")
            raise
        finally:
            os.chdir(original)

    def _write_script(self, path, stages, flags, maxwarn):
        blocks = ["#!/bin/bash", "set -e", "", "# Hardware-adaptive, stage-aware GROMACS workflow (PRISM)", ""]
        for stage in stages:
            r = _STAGE_RECIPE[stage]
            sd = os.path.dirname(r["deffnm"])
            t_flag = f" -t {r['use_t']}" if r["use_t"] else ""
            blocks.append(f"""# {stage.upper()}
mkdir -p {sd}
if [ -f ./{r['deffnm']}.gro ]; then
    echo "{stage.upper()} already completed, skipping..."
elif [ -f ./{r['deffnm']}.tpr ]; then
    echo "{stage.upper()} tpr found, continuing from checkpoint..."
    {self.gmx} mdrun -s ./{r['deffnm']}.tpr -deffnm ./{r['deffnm']} {flags} $([ -f ./{r['deffnm']}.cpt ] && echo "-cpi ./{r['deffnm']}.cpt")
else
    echo "Starting {stage.upper()}..."
    {self.gmx} grompp -f ../mdps/{r['mdp']} -c ./{r['prev']} -r ./{r['restr']}{t_flag} -p topol.top -o ./{r['deffnm']}.tpr -maxwarn {maxwarn}
    {self.gmx} mdrun -s ./{r['deffnm']}.tpr -deffnm ./{r['deffnm']} {flags}
fi
""")
        with open(path, "w") as f:
            f.write("\n".join(blocks))

    def _collect_outputs(self, stages):
        outputs = {}
        ext_map = {"em": ("gro", "edr", "log", "tpr"),
                   "nvt": ("gro", "edr", "log", "tpr", "cpt"),
                   "npt": ("gro", "edr", "log", "tpr", "cpt"),
                   "prod": ("gro", "xtc", "edr", "log", "tpr", "cpt")}
        for stage in stages:
            deffnm = _STAGE_RECIPE[stage]["deffnm"]
            d = {}
            for ext in ext_map[stage]:
                p = os.path.join(self.gmx_dir, f"{deffnm}.{ext}")
                if os.path.exists(p):
                    d[ext] = p
            if d:
                outputs[stage] = d
        return outputs
