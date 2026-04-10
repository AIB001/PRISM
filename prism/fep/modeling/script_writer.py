#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FEP script generation utilities.

This module provides functions for generating shell scripts and SLURM submission
scripts for FEP calculations.

IMPORTANT: PRISM_MDRUN_NSTEPS Environment Variable
-------------------------------------------------
The PRISM_MDRUN_NSTEPS environment variable is intended for TESTING ONLY:
- Applied to: per-window NPT short equilibration and production runs
- NOT applied to: EM, NVT, NPT (these must complete fully for proper system equilibration)

IMPORTANT: PRISM_EQ_NSTEPS Environment Variable
----------------------------------------------
The PRISM_EQ_NSTEPS environment variable is intended for fast infrastructure
validation only:
- Applied to: EM, NVT, NPT via temporary test-mode MDP files
- NOT valid for any scientific result or publication

Usage for testing:
    export PRISM_MDRUN_NSTEPS=100  # Run only 100 steps for quick validation
    cd GMX_PROLIG_FEP && bash run_fep.sh bound

For production FEP calculations:
    unset PRISM_MDRUN_NSTEPS  # Ensure it's not set
    cd GMX_PROLIG_FEP && bash run_fep.sh bound

WARNING: Results obtained with PRISM_MDRUN_NSTEPS set are NOT valid for publication!
"""

import os
import textwrap
from pathlib import Path
from typing import Optional


def _detect_num_gpus(exec_config: dict) -> int:
    """Resolve the total GPU count for generated scripts."""
    num_gpus_config = exec_config.get("num_gpus")
    if num_gpus_config is not None:
        return int(num_gpus_config)

    import os

    cuda_devices = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if cuda_devices:
        return len([gpu for gpu in cuda_devices.split(",") if gpu.strip()])
    return 1


def _get_forcefield_name(config: Optional[dict]) -> str:
    """Extract the configured protein force field name as a normalized string."""
    config = config or {}
    protein_ff_info = config.get("forcefield", "")
    if isinstance(protein_ff_info, dict) and "name" in protein_ff_info:
        protein_ff = protein_ff_info["name"]
    elif isinstance(protein_ff_info, str):
        protein_ff = protein_ff_info
    else:
        protein_ff = str(protein_ff_info)
    return protein_ff.lower() if isinstance(protein_ff, str) else ""


def _forcefield_uses_cmap(config: Optional[dict]) -> bool:
    """Return True for force fields that need conservative CMAP runtime defaults."""
    protein_ff = _get_forcefield_name(config)
    return "charmm" in protein_ff or protein_ff == "amber19sb"


def _get_execution_settings(config: Optional[dict]) -> dict:
    """Normalize execution-layer settings used by script generation.

    Resource allocation strategy:
    1. If execution.total_cpus is set: omp_threads = total_cpus // num_gpus
    2. If execution.omp_threads is set: use that value directly (manual override)
    3. Otherwise: auto-detect from system or use default
    """
    config = config or {}
    exec_config = dict(config.get("execution", {}))
    fep_config = dict(config.get("fep", {}))

    num_gpus = _detect_num_gpus(exec_config)
    mode = str(exec_config.get("mode", "standard")).strip().lower()
    if mode not in {"standard", "repex"}:
        raise ValueError(f"Unsupported execution.mode: {mode}")

    parallel_windows = exec_config.get("parallel_windows", num_gpus)
    if parallel_windows is None:
        parallel_windows = num_gpus
    parallel_windows = max(1, int(parallel_windows))

    # Calculate omp_threads from total_cpus if available.
    # Standard mode runs one window per active GPU/process, so divide by the
    # effective concurrent window count. Replica exchange uses one shared job
    # across the configured GPUs, so divide by the GPU count.
    total_cpus = exec_config.get("total_cpus")
    if total_cpus is not None:
        if mode == "standard":
            divisor = max(1, min(int(num_gpus), int(parallel_windows)))
        else:
            divisor = max(1, int(num_gpus))
        omp_threads = max(1, int(total_cpus) // divisor)
    elif exec_config.get("omp_threads"):
        # Manual override
        omp_threads = int(exec_config.get("omp_threads"))
    else:
        # Default fallback
        omp_threads = 8

    # CMAP force fields are less robust with GPU-side update/constraints during
    # early equilibration on hybrid topologies, so prefer CPU updates by default.
    protein_ff = _get_forcefield_name(config)
    is_charmm = "charmm" in protein_ff
    uses_cmap_runtime_guard = _forcefield_uses_cmap(config)

    update_mode = str(exec_config.get("mdrun_update_mode", "auto")).strip().lower()
    if update_mode not in {"auto", "gpu", "cpu", "none"}:
        raise ValueError(f"Invalid execution.mdrun_update_mode: {update_mode}")

    if update_mode == "auto":
        # Keep the run GPU-resident, but move update/constraint work to CPU for
        # CMAP systems where GPU update groups are more failure-prone.
        update_mode = (
            "cpu" if uses_cmap_runtime_guard else ("gpu" if bool(exec_config.get("use_gpu_update", False)) else "none")
        )

    update_flag = "" if update_mode == "none" else f"-update {update_mode}"
    use_gpu_update = update_mode == "gpu"
    min_gromacs_version, min_gromacs_reason = _infer_min_gromacs_version(config)

    return {
        "mode": mode,
        "use_gpu_pme": bool(exec_config.get("use_gpu_pme", True)),
        "use_gpu_update": use_gpu_update,
        "mdrun_update_mode": update_mode,
        "update_flag": update_flag,
        "num_gpus": max(1, int(num_gpus)),
        "parallel_windows": parallel_windows,
        "omp_threads": omp_threads,
        "total_cpus": total_cpus,  # For logging/validation
        "replicas": int(config.get("replicas", fep_config.get("replicas", 3))),
        "is_charmm": is_charmm,  # Pass for logging
        "min_gromacs_version": min_gromacs_version,
        "min_gromacs_reason": min_gromacs_reason,
    }


def _infer_min_gromacs_version(config: Optional[dict]) -> tuple[Optional[str], Optional[str]]:
    """Return an optional minimum GROMACS version required by this scaffold."""
    config = config or {}
    exec_config = dict(config.get("execution", {}))
    explicit = exec_config.get("min_gromacs_version")
    if explicit:
        return str(explicit), "configured via execution.min_gromacs_version"

    protein_ff = _get_forcefield_name(config)

    if protein_ff.strip() == "amber19sb":
        return "2026.0", "amber19sb CMAP wildcard support requires GROMACS 2026.0+"
    return None, None


def _write_leg_execution_scripts(
    leg_dir: Path,
    execution_mode: str,  # noqa: ARG001
    *,
    use_gpu_pme: bool,
    update_mode: str,
    num_gpus: int,
    parallel_windows: int,
    omp_threads: int,
) -> None:
    """Generate per-leg production helper scripts for standard and repex modes."""
    pme_gpu_flag = "-pme gpu" if use_gpu_pme else ""
    standard_script = leg_dir / "run_prod_standard.sh"
    repex_script = leg_dir / "run_prod_repex.sh"
    runtime_resource_block = f"""\
    DEFAULT_NUM_GPUS={num_gpus}
    DEFAULT_PARALLEL_WINDOWS={parallel_windows}
    DEFAULT_OMP_THREADS={omp_threads}

    get_positive_int() {{
        local value="$1"
        local name="$2"
        if ! [[ "${{value}}" =~ ^[0-9]+$ ]] || [ "${{value}}" -lt 1 ]; then
            echo "Error: ${{name}} must be a positive integer (got '${{value}}')"
            exit 1
        fi
    }}

    NUM_GPUS="${{PRISM_NUM_GPUS:-${{DEFAULT_NUM_GPUS}}}}"
    PARALLEL_WINDOWS="${{PRISM_PARALLEL_WINDOWS:-${{DEFAULT_PARALLEL_WINDOWS}}}}"
    OMP_THREADS="${{PRISM_OMP_THREADS:-${{OMP_NUM_THREADS:-${{DEFAULT_OMP_THREADS}}}}}}"

    get_positive_int "${{NUM_GPUS}}" "PRISM_NUM_GPUS"
    get_positive_int "${{PARALLEL_WINDOWS}}" "PRISM_PARALLEL_WINDOWS"
    get_positive_int "${{OMP_THREADS}}" "PRISM_OMP_THREADS/OMP_NUM_THREADS"

    export OMP_NUM_THREADS="${{OMP_THREADS}}"
    """

    standard_content = textwrap.dedent(f"""\
    #!/usr/bin/env bash
    set -euo pipefail
    
    LEG_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    MDP_DIR="${{LEG_DIR}}/mdps"
    BUILD_DIR="${{LEG_DIR}}/build"
    
    if [ -n "${{PRISM_MDRUN_NSTEPS:-}}" ]; then
        MDRUN_NSTEPS_ARG="-nsteps ${{PRISM_MDRUN_NSTEPS}}"
    else
        MDRUN_NSTEPS_ARG=""
    fi

{runtime_resource_block}
    
    UPDATE_MODE="${{PRISM_MDRUN_UPDATE_MODE:-{update_mode}}}"
    case "${{UPDATE_MODE}}" in
        gpu|cpu)
            UPDATE_FLAG="-update ${{UPDATE_MODE}}"
            ;;
        none)
            UPDATE_FLAG=""
            ;;
        *)
            echo "Error: unsupported PRISM_MDRUN_UPDATE_MODE '${{UPDATE_MODE}}'"
            echo "Supported values: gpu, cpu, none"
            exit 1
            ;;
    esac
    
    MDRUN_GPU_ARGS="-nb gpu -bonded gpu ${{UPDATE_FLAG}} {pme_gpu_flag} -pin on"
    
    check_em_short() {{
        local window_dir="$1"
        local lambda_name="$2"
        local em_log="${{window_dir}}/em_short.log"
    
        if [ ! -f "${{em_log}}" ]; then
            echo "Error: ${{lambda_name}} missing em_short.log"
            return 1
        fi
    
        local max_force
        max_force=$(python - "${{em_log}}" <<'PY'
    import re
    import sys
    from pathlib import Path
    
    text = Path(sys.argv[1]).read_text(errors="ignore")
    match = re.search(r"Maximum force\\s*=\\s*([0-9.eE+-]+)", text)
    print(match.group(1) if match else "")
    PY
    )
    
        if [ -z "${{max_force}}" ]; then
            echo "Error: ${{lambda_name}} em_short maximum force not found"
            return 1
        fi
    
        python - "${{max_force}}" <<'PY'
    import sys
    sys.exit(0 if float(sys.argv[1]) <= 1000.0 else 1)
    PY
        if [ $? -ne 0 ]; then
            echo "Error: ${{lambda_name}} em_short did not converge (Fmax=${{max_force}})"
            return 1
        fi
    }}
    
    
    echo "Production mode: standard"
    echo "Concurrent windows: ${{PARALLEL_WINDOWS}}"
    echo "Configured GPUs: ${{NUM_GPUS}}"
    echo "OpenMP threads per GPU: ${{OMP_THREADS}}"
    echo "CPU affinity: enabled (each GPU gets dedicated CPU cores)"
    
    JOB_COUNT=0
    PIDS=()
    for lambda_mdp in "${{MDP_DIR}}"/prod_*.mdp; do
        lambda_name=$(basename "${{lambda_mdp}}" .mdp)
        lambda_idx="${{lambda_name#prod_}}"
        window_dir="${{LEG_DIR}}/window_${{lambda_idx}}"
        mkdir -p "${{window_dir}}"
    
        if [ -f "${{window_dir}}/prod.gro" ]; then
            echo "  ${{lambda_name}}: already completed"
            continue
        fi
    
        if [ "${{PARALLEL_WINDOWS}}" -gt 1 ]; then
            running_jobs=$(find "${{LEG_DIR}}" -maxdepth 1 -name ".run_*" -type f 2>/dev/null | wc -l)
            while [ "${{running_jobs}}" -ge "${{PARALLEL_WINDOWS}}" ]; do
                sleep 1
                running_jobs=$(find "${{LEG_DIR}}" -maxdepth 1 -name ".run_*" -type f 2>/dev/null | wc -l)
            done
        fi
    
        if [ "${{NUM_GPUS}}" -gt 1 ]; then
            gpu_id=$((JOB_COUNT % NUM_GPUS))
        else
            gpu_id=0
        fi
    
        run_window() {{
            local lambda_idx="$1"
            local lambda_name="$2"
            local lambda_mdp="$3"
            local window_dir="$4"
            local gpu_id="$5"
    
            export CUDA_VISIBLE_DEVICES="${{gpu_id}}"
            local cpu_offset=$((gpu_id * OMP_THREADS))

            echo "$$" > "${{LEG_DIR}}/.run_${{lambda_idx}}"
            trap 'rm -f "${{LEG_DIR}}/.run_${{lambda_idx}}"' EXIT
    
            if [ -f "${{window_dir}}/npt_short.gro" ]; then
                echo "  ${{lambda_name}}: NPT short already completed"
            elif [ -f "${{window_dir}}/em_short.gro" ]; then
                echo "  ${{lambda_name}}: lambda-specific EM already completed"
                if [ -f "${{window_dir}}/npt_short.tpr" ] && [ -f "${{window_dir}}/npt_short.cpt" ]; then
                    echo "  ${{lambda_name}}: resuming NPT short on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                    gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -cpi "${{window_dir}}/npt_short.cpt" -v ${{MDRUN_NSTEPS_ARG}}
                else
                    echo "  ${{lambda_name}}: starting NPT short on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                    gmx grompp -f "${{MDP_DIR}}/npt_short_${{lambda_idx}}.mdp" -c "${{window_dir}}/em_short.gro" -r "${{window_dir}}/em_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/npt_short.tpr" -maxwarn 20
                    gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -v ${{MDRUN_NSTEPS_ARG}}
                fi
            elif [ -f "${{window_dir}}/em_short.tpr" ]; then
                echo "  ${{lambda_name}}: resuming lambda-specific EM"
                gmx mdrun -deffnm "${{window_dir}}/em_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" -v
                check_em_short "${{window_dir}}" "${{lambda_name}}"
                echo "  ${{lambda_name}}: starting NPT short on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                gmx grompp -f "${{MDP_DIR}}/npt_short_${{lambda_idx}}.mdp" -c "${{window_dir}}/em_short.gro" -r "${{window_dir}}/em_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/npt_short.tpr" -maxwarn 20
                gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -v ${{MDRUN_NSTEPS_ARG}}
            elif [ -f "${{window_dir}}/npt_short.tpr" ] && [ -f "${{window_dir}}/npt_short.cpt" ]; then
                echo "  ${{lambda_name}}: resuming NPT short on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -cpi "${{window_dir}}/npt_short.cpt" -v ${{MDRUN_NSTEPS_ARG}}
            else
                echo "  ${{lambda_name}}: starting lambda-specific EM"
                gmx grompp -f "${{MDP_DIR}}/em_short_${{lambda_idx}}.mdp" -c "${{BUILD_DIR}}/npt.gro" -r "${{BUILD_DIR}}/npt.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/em_short.tpr" -maxwarn 20
                gmx mdrun -deffnm "${{window_dir}}/em_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" -v
                check_em_short "${{window_dir}}" "${{lambda_name}}"
                echo "  ${{lambda_name}}: starting NPT short on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                gmx grompp -f "${{MDP_DIR}}/npt_short_${{lambda_idx}}.mdp" -c "${{window_dir}}/em_short.gro" -r "${{window_dir}}/em_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/npt_short.tpr" -maxwarn 20
                gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -v ${{MDRUN_NSTEPS_ARG}}
            fi
    
            if [ -f "${{window_dir}}/prod.gro" ]; then
                echo "  ${{lambda_name}}: production already completed"
            elif [ -f "${{window_dir}}/prod.tpr" ] && [ -f "${{window_dir}}/prod.cpt" ]; then
                echo "  ${{lambda_name}}: resuming production on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                gmx mdrun -deffnm "${{window_dir}}/prod" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -cpi "${{window_dir}}/prod.cpt" -v ${{MDRUN_NSTEPS_ARG}}
            else
                echo "  ${{lambda_name}}: starting production on GPU ${{gpu_id}} (CPUs ${{cpu_offset}}-$((cpu_offset + OMP_THREADS - 1)))"
                if [ ! -f "${{window_dir}}/prod.tpr" ]; then
                    gmx grompp -f "${{lambda_mdp}}" -c "${{window_dir}}/npt_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/prod.tpr" -maxwarn 20
                fi
                gmx mdrun -deffnm "${{window_dir}}/prod" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -pinoffset "${{cpu_offset}}" -pinstride 1 -v ${{MDRUN_NSTEPS_ARG}}
            fi
    
            rm -f "${{LEG_DIR}}/.run_${{lambda_idx}}"
            trap - EXIT
            echo "  ${{lambda_name}}: ✓ Complete"
        }}
    
        if [ "${{PARALLEL_WINDOWS}}" -gt 1 ]; then
            run_window "${{lambda_idx}}" "${{lambda_name}}" "${{lambda_mdp}}" "${{window_dir}}" "${{gpu_id}}" &
            PIDS+=($!)
            JOB_COUNT=$((JOB_COUNT + 1))
        else
            run_window "${{lambda_idx}}" "${{lambda_name}}" "${{lambda_mdp}}" "${{window_dir}}" "${{gpu_id}}"
        fi
    done
    
    if [ "${{PARALLEL_WINDOWS}}" -gt 1 ]; then
        echo "Waiting for all windows to complete..."
        failure=0
        for job_pid in "${{PIDS[@]}}"; do
            if ! wait "${{job_pid}}"; then
                failure=1
            fi
        done
        if [ "${{failure}}" -ne 0 ]; then
            echo "Error: one or more lambda windows failed."
            exit 1
        fi
    fi
    """)

    repex_content = textwrap.dedent(f"""\
    #!/usr/bin/env bash
    set -euo pipefail
    
    LEG_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
    MDP_DIR="${{LEG_DIR}}/mdps"
    BUILD_DIR="${{LEG_DIR}}/build"
    
    if [ -n "${{PRISM_MDRUN_NSTEPS:-}}" ]; then
        MDRUN_NSTEPS_ARG="-nsteps ${{PRISM_MDRUN_NSTEPS}}"
    else
        MDRUN_NSTEPS_ARG=""
    fi

{runtime_resource_block}
    
    UPDATE_MODE="${{PRISM_MDRUN_UPDATE_MODE:-{update_mode}}}"
    case "${{UPDATE_MODE}}" in
        gpu|cpu)
            UPDATE_FLAG="-update ${{UPDATE_MODE}}"
            ;;
        none)
            UPDATE_FLAG=""
            ;;
        *)
            echo "Error: unsupported PRISM_MDRUN_UPDATE_MODE '${{UPDATE_MODE}}'"
            echo "Supported values: gpu, cpu, none"
            exit 1
            ;;
    esac
    
    MDRUN_GPU_ARGS="-nb gpu -bonded gpu ${{UPDATE_FLAG}} {pme_gpu_flag} -pin on"
    
    check_em_short() {{
        local window_dir="$1"
        local lambda_name="$2"
        local em_log="${{window_dir}}/em_short.log"
    
        if [ ! -f "${{em_log}}" ]; then
            echo "Error: ${{lambda_name}} missing em_short.log"
            return 1
        fi
    
        local max_force
        max_force=$(python - "${{em_log}}" <<'PY'
    import re
    import sys
    from pathlib import Path
    
    text = Path(sys.argv[1]).read_text(errors="ignore")
    match = re.search(r"Maximum force\\s*=\\s*([0-9.eE+-]+)", text)
    print(match.group(1) if match else "")
    PY
    )
    
        if [ -z "${{max_force}}" ]; then
            echo "Error: ${{lambda_name}} em_short maximum force not found"
            return 1
        fi
    
        python - "${{max_force}}" <<'PY'
    import sys
    sys.exit(0 if float(sys.argv[1]) <= 1000.0 else 1)
    PY
        if [ $? -ne 0 ]; then
            echo "Error: ${{lambda_name}} em_short did not converge (Fmax=${{max_force}})"
            return 1
        fi
    }}
    
    window_dirs=()
    for lambda_mdp in "${{MDP_DIR}}"/prod_*.mdp; do
        lambda_name=$(basename "${{lambda_mdp}}" .mdp)
        lambda_idx="${{lambda_name#prod_}}"
        window_dir="${{LEG_DIR}}/window_${{lambda_idx}}"
        window_dirs+=("${{window_dir}}")
        mkdir -p "${{window_dir}}"
    
        if [ ! -f "${{window_dir}}/npt_short.gro" ]; then
            if [ ! -f "${{window_dir}}/em_short.gro" ]; then
                if [ ! -f "${{window_dir}}/em_short.tpr" ]; then
                    gmx grompp -f "${{MDP_DIR}}/em_short_${{lambda_idx}}.mdp" -c "${{BUILD_DIR}}/npt.gro" -r "${{BUILD_DIR}}/npt.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/em_short.tpr" -maxwarn 20
                fi
                gmx mdrun -deffnm "${{window_dir}}/em_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" -v
                check_em_short "${{window_dir}}" "${{lambda_name}}"
            fi
            if [ ! -f "${{window_dir}}/npt_short.tpr" ]; then
                gmx grompp -f "${{MDP_DIR}}/npt_short_${{lambda_idx}}.mdp" -c "${{window_dir}}/em_short.gro" -r "${{window_dir}}/em_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/npt_short.tpr" -maxwarn 20
            fi
            gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp "${{OMP_THREADS}}" ${{MDRUN_GPU_ARGS}} -v ${{MDRUN_NSTEPS_ARG}}
        fi
    
        if [ ! -f "${{window_dir}}/prod.tpr" ]; then
            gmx grompp -f "${{lambda_mdp}}" -c "${{window_dir}}/npt_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/prod.tpr" -maxwarn 20
        fi
    done
    
    num_windows="${{#window_dirs[@]}}"
    if [ "${{num_windows}}" -eq 0 ]; then
        echo "Error: no prod_*.mdp files found in $MDP_DIR"
        exit 1
    fi
    
    gpu_ids=$(seq 0 $((NUM_GPUS - 1)) | awk '{{printf "%s", $0}}')
    window_cpi="${{window_dirs[0]}}/prod.cpt"
    cpi_args=()
    if [ -f "${{window_cpi}}" ]; then
        cpi_args=(-cpi prod.cpt)
    fi
    
    echo "Production mode: repex"
    echo "Replica-exchange windows: ${{num_windows}}"
    echo "Configured GPUs: $NUM_GPUS"
    echo "OpenMP threads: $OMP_THREADS"
    echo "GPU id string: $gpu_ids"
    
    mpirun -oversubscribe -np "${{num_windows}}" \\
        gmx_mpi mdrun -v -deffnm prod ${{MDRUN_GPU_ARGS}} -replex 1000 \\
        -multidir "${{window_dirs[@]}}" -gpu_id "${{gpu_ids}}" -pin on -dhdl dhdl \\
        "${{cpi_args[@]}}" ${{MDRUN_NSTEPS_ARG}}
    """)

    for path, content in ((standard_script, standard_content), (repex_script, repex_content)):
        with open(path, "w") as handle:
            handle.write(content)
        path.chmod(0o755)


def write_root_scripts(layout, config: Optional[dict] = None) -> None:
    """
    Write root-level run_fep.sh master script.

    Parameters
    ----------
    layout : FEPScaffoldLayout
        The scaffold layout containing paths to bound/unbound legs
    config : dict, optional
        Configuration dictionary with GPU and parallel execution settings

    Note
    ----
    This generates a single master script (run_fep.sh) that can run bound, unbound,
    or both legs via command-line arguments. Legacy run_all.sh and submit_all.slurm.sh
    scripts are no longer generated as they are redundant with run_fep.sh.

    Configuration Options
    ---------------------
    execution.total_cpus : int
        Total number of CPU cores available. When set, OpenMP threads per GPU
        is automatically calculated as: total_cpus // num_gpus
        Example: 40 total CPUs with 4 GPUs → 10 threads per GPU
    execution.num_gpus : int
        Number of GPUs available (default: 1, or from CUDA_VISIBLE_DEVICES)
    execution.omp_threads : int
        OpenMP threads per GPU (manual override). If total_cpus is set, this
        is ignored and auto-calculated instead. If neither is set, defaults to 8.
    execution.parallel_windows : int
        Number of concurrent lambda windows in ``standard`` mode (default: num_gpus)
    execution.use_gpu_pme : bool
        Enable GPU PME (default: True)
    execution.mdrun_update_mode : str
        Update/constraint placement for GPU runs: ``auto``, ``gpu``, ``cpu``,
        or ``none``. ``auto`` maps CHARMM to ``cpu`` and other force fields to
        ``none`` by default for non-CHARMM FEP scaffolds unless legacy ``use_gpu_update: true`` enables GPU-side update explicitly.
    execution.mode : str
        Production execution mode: ``standard`` (independent windows) or
        ``repex`` (lambda replica exchange via ``gmx_mpi -multidir``)
    fep.replicas : int
        Number of replica runs for error estimation (default: 1)
        When >1, creates bound1, bound2, ... and unbound1, unbound2, ... directories

    Examples
    --------
    Example 1: Auto-calculate threads from total CPUs
    >>> config = {
    ...     "execution": {
    ...         "total_cpus": 40,  # Total CPU cores
    ...         "num_gpus": 4,     # 4 GPUs
    ...         # omp_threads will be auto-calculated as 40//4 = 10
    ...     }
    ... }

    Example 2: Manual thread specification
    >>> config = {
    ...     "execution": {
    ...         "omp_threads": 12,  # Explicit value
    ...         "num_gpus": 4,
    ...     }
    ... }

    Example 3: Default configuration
    >>> config = {
    ...     "execution": {
    ...         "num_gpus": 4,
    ...         # Uses default omp_threads=8
    ...     }
    ... }
    """
    settings = _get_execution_settings(config)
    for leg_dir in (layout.bound_dir, layout.unbound_dir):
        _write_leg_execution_scripts(
            leg_dir,
            settings["mode"],
            use_gpu_pme=settings["use_gpu_pme"],
            update_mode=settings["mdrun_update_mode"],
            num_gpus=settings["num_gpus"],
            parallel_windows=settings["parallel_windows"],
            omp_threads=settings["omp_threads"],
        )
    _propagate_leg_scripts_to_repeats(layout)
    write_fep_master_script(layout.root, config)


def _propagate_leg_scripts_to_repeats(layout) -> None:
    """Make per-repeat production scripts available in repeat2..N.

    repeat1 owns the real generated scripts. Other repeats can safely reuse them
    because the scripts resolve their working directory from the repeat they are
    executed in, while the static script content is identical across repeats.
    """
    script_names = ("run_prod_standard.sh", "run_prod_repex.sh")
    for leg_root in (layout.root / "bound", layout.root / "unbound"):
        repeat1 = leg_root / "repeat1"
        if not repeat1.exists():
            continue
        for repeat_dir in sorted(leg_root.glob("repeat*")):
            if repeat_dir.name == "repeat1" or not repeat_dir.is_dir():
                continue
            for script_name in script_names:
                source = repeat1 / script_name
                if not source.exists():
                    continue
                target = repeat_dir / script_name
                if target.exists() or target.is_symlink():
                    continue
                target.symlink_to(os.path.relpath(source, repeat_dir))


def write_fep_master_script(fep_dir: Path, config: Optional[dict] = None) -> None:
    """
    Generate master run_fep.sh script at FEP root level.
    This single script can run bound, unbound, or both legs via command-line arguments.

    Parameters
    ----------
    fep_dir : Path
        Path to the FEP root directory (GMX_PROLIG_FEP/)
    config : dict, optional
        Configuration dictionary with GPU settings
    """
    settings = _get_execution_settings(config)
    execution_mode = settings["mode"]
    use_gpu_pme = settings["use_gpu_pme"]
    update_flag = settings["update_flag"]
    num_gpus = settings["num_gpus"]
    parallel_windows = settings["parallel_windows"]
    omp_threads = settings["omp_threads"]
    replicas = settings["replicas"]
    update_mode = settings["mdrun_update_mode"]
    pme_gpu_flag = "-pme gpu" if use_gpu_pme else ""

    script_path = fep_dir / "run_fep.sh"
    script_content = textwrap.dedent(f"""\
#!/usr/bin/env bash
# PRISM-FEP Master Script
# GPU PME: {"enabled" if use_gpu_pme else "disabled (CPU PME for small systems)"}
set -euo pipefail

# Test mode: limit production run steps (DO NOT use for production FEP calculations)
if [ -n "${{PRISM_MDRUN_NSTEPS:-}}" ]; then
    echo "⚠️  WARNING: PRISM_MDRUN_NSTEPS is set to ${{PRISM_MDRUN_NSTEPS}}"
    echo "⚠️  This is for TESTING ONLY - results will NOT be valid for publication!"
    MDRUN_NSTEPS_ARG="-nsteps ${{PRISM_MDRUN_NSTEPS}}"
else
    MDRUN_NSTEPS_ARG=""
fi

UPDATE_MODE="${{PRISM_MDRUN_UPDATE_MODE:-{update_mode}}}"
case "${{UPDATE_MODE}}" in
    gpu|cpu)
        UPDATE_FLAG="-update ${{UPDATE_MODE}}"
        ;;
    none)
        UPDATE_FLAG=""
        ;;
    *)
        echo "Error: unsupported PRISM_MDRUN_UPDATE_MODE '${{UPDATE_MODE}}'"
        echo "Supported values: gpu, cpu, none"
        exit 1
        ;;
esac

MDRUN_GPU_ARGS="-nb gpu -bonded gpu ${{UPDATE_FLAG}} {pme_gpu_flag} -pin on"

EXECUTION_MODE="${{PRISM_FEP_MODE:-{execution_mode}}}"
if [[ "${{EXECUTION_MODE}}" != "standard" && "${{EXECUTION_MODE}}" != "repex" ]]; then
    echo "Error: unsupported execution mode '${{EXECUTION_MODE}}'"
    echo "Supported modes: standard, repex"
    exit 1
fi

    DEFAULT_NUM_GPUS={num_gpus}
    DEFAULT_PARALLEL_WINDOWS={parallel_windows}
    DEFAULT_OMP_THREADS={omp_threads}

    get_positive_int() {{
        local value="$1"
        local name="$2"
        if ! [[ "${{value}}" =~ ^[0-9]+$ ]] || [ "${{value}}" -lt 1 ]; then
            echo "Error: ${{name}} must be a positive integer (got '${{value}}')"
            exit 1
        fi
    }}

    GPU_COUNT="${{PRISM_NUM_GPUS:-${{DEFAULT_NUM_GPUS}}}}"
    PARALLEL_WINDOWS="${{PRISM_PARALLEL_WINDOWS:-${{DEFAULT_PARALLEL_WINDOWS}}}}"

    get_positive_int "${{GPU_COUNT}}" "PRISM_NUM_GPUS"
    get_positive_int "${{PARALLEL_WINDOWS}}" "PRISM_PARALLEL_WINDOWS"

    if [ -n "${{PRISM_OMP_THREADS:-}}" ]; then
        export OMP_NUM_THREADS="${{PRISM_OMP_THREADS}}"
    elif [ -n "${{OMP_NUM_THREADS:-}}" ]; then
        export OMP_NUM_THREADS="${{OMP_NUM_THREADS}}"
    elif [ -n "${{PRISM_TOTAL_CPUS:-}}" ]; then
        get_positive_int "${{PRISM_TOTAL_CPUS}}" "PRISM_TOTAL_CPUS"
        if [ "${{EXECUTION_MODE}}" = "standard" ]; then
            ACTIVE_WORKERS="${{GPU_COUNT}}"
            if [ "${{PARALLEL_WINDOWS}}" -lt "${{ACTIVE_WORKERS}}" ]; then
                ACTIVE_WORKERS="${{PARALLEL_WINDOWS}}"
            fi
        else
            ACTIVE_WORKERS="${{GPU_COUNT}}"
        fi
        export OMP_NUM_THREADS="$((PRISM_TOTAL_CPUS / ACTIVE_WORKERS))"
        if [ "${{OMP_NUM_THREADS}}" -lt 1 ]; then
            export OMP_NUM_THREADS=1
        fi
    else
        export OMP_NUM_THREADS="{omp_threads}"
    fi

    export OMP_PROC_BIND="${{OMP_PROC_BIND:-close}}"
    export OMP_PLACES="${{OMP_PLACES:-cores}}"

resolve_gpu_id() {{
    if [ -n "${{PRISM_GPU_ID:-}}" ]; then
        printf '%s\n' "${{PRISM_GPU_ID}}"
        return 0
    fi

    if [ -n "${{GMX_GPU_ID:-}}" ]; then
        printf '%s\n' "${{GMX_GPU_ID}}"
        return 0
    fi

    if [ -n "${{CUDA_VISIBLE_DEVICES:-}}" ] && [[ "${{CUDA_VISIBLE_DEVICES}}" != *,* ]]; then
        printf '%s\n' "${{CUDA_VISIBLE_DEVICES}}"
        return 0
    fi

    printf '%s\n' 0
}}

# Resolve a run target to its on-disk leg directory.
resolve_leg_dir() {{
    local target="$1"

    if [[ "${{target}}" =~ ^(bound|unbound)([0-9]+)$ ]]; then
        local leg="${{BASH_REMATCH[1]}}"
        local replica_idx="${{BASH_REMATCH[2]}}"
        local nested_dir="${{FEP_DIR}}/${{leg}}/repeat${{replica_idx}}"
        if [ -d "${{nested_dir}}" ]; then
            printf '%s\\n' "${{nested_dir}}"
            return 0
        fi
    fi

    if [[ "${{target}}" =~ ^(bound|unbound)$ ]]; then
        local nested_dir="${{FEP_DIR}}/${{target}}/repeat1"
        if [ -d "${{nested_dir}}" ]; then
            printf '%s\\n' "${{nested_dir}}"
            return 0
        fi
    fi

    if [ -d "${{FEP_DIR}}/${{target}}" ]; then
        printf '%s\\n' "${{FEP_DIR}}/${{target}}"
        return 0
    fi

    return 1
}}

# Function to run a single leg / repeat target
run_leg() {{
    local target_name="$1"
    local leg_dir
    leg_dir=$(resolve_leg_dir "${{target_name}}") || {{
        echo "Error: Leg directory not found for target: ${{target_name}}"
        return 1
    }}
    local gpu_id
    gpu_id=$(resolve_gpu_id)
    local cpu_offset="${{PRISM_CPU_OFFSET:-$((gpu_id * OMP_NUM_THREADS))}}"
    export CUDA_VISIBLE_DEVICES="${{gpu_id}}"

    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                    PRISM-FEP ${{target_name}} LEG                                   ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "GPU assignment: physical GPU ${{gpu_id}} (visible as local GPU 0)"
    echo "OpenMP threads: ${{OMP_NUM_THREADS}}"
    echo "Configured GPU pool: ${{GPU_COUNT}}"
    echo "Parallel windows: ${{PARALLEL_WINDOWS}}"
    echo "CPU pin offset: ${{cpu_offset}}"

    cd "${{leg_dir}}"

    local MDP_DIR="${{leg_dir}}/mdps"
    local INPUT_DIR="${{leg_dir}}/input"
    local BUILD_DIR="${{leg_dir}}/build"

    # Energy minimization
    mkdir -p ${{BUILD_DIR}}

    if [ -f ${{BUILD_DIR}}/em.gro ]; then
        echo "✓ EM already completed, skipping..."
    elif [ -f ${{BUILD_DIR}}/em.tpr ]; then
        echo "EM checkpoint found, resuming..."
        gmx mdrun -deffnm ${{BUILD_DIR}}/em -cpi ${{BUILD_DIR}}/em.cpt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} -pin on -pinoffset "${{cpu_offset}}" -pinstride 1 -v
    else
        echo "Running energy minimization..."
        gmx grompp -f ${{MDP_DIR}}/em.mdp -c ${{INPUT_DIR}}/conf.gro -p topol.top -o ${{BUILD_DIR}}/em.tpr -maxwarn 20
        gmx mdrun -deffnm ${{BUILD_DIR}}/em -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} -pin on -pinoffset "${{cpu_offset}}" -pinstride 1 -v
    fi

    # Check EM convergence
    if [ -f ${{BUILD_DIR}}/em.log ]; then
        max_force=$(grep -E "Maximum force|Fmax=" ${{BUILD_DIR}}/em.log | tail -1 | grep -oE "[0-9]+\\\\.?[0-9]*e[+-][0-9]+" | tail -1 || true)
        if [ -n "$max_force" ]; then
            echo "EM final max force: $max_force kJ/mol/nm"
            # Warning if force is still high (> 1000 kJ/mol/nm)
            # Use awk to handle scientific notation comparison
            if [ "$(echo "$max_force" | awk '{{print ($1 > 1000)}}')" -eq 1 ]; then
                echo "⚠️  WARNING: EM converged with high force. System may need more minimization."
            fi
        fi
    fi

    # NVT equilibration (must run to completion)
    if [ -f ${{BUILD_DIR}}/nvt.gro ]; then
        echo "✓ NVT already completed, skipping..."
    elif [ -f ${{BUILD_DIR}}/nvt.tpr ]; then
        if [ -f ${{BUILD_DIR}}/nvt.cpt ]; then
            echo "NVT checkpoint found, resuming on GPU..."
            (
                cd ${{BUILD_DIR}}
                gmx mdrun -deffnm nvt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} ${{MDRUN_GPU_ARGS}} -gpu_id 0 -pinoffset "${{cpu_offset}}" -pinstride 1 -cpi nvt.cpt -v
            )
        else
            echo "⚠️  Previous NVT run failed (no checkpoint), cleaning up and restarting..."
            rm -f ${{BUILD_DIR}}/nvt.* 2>/dev/null || true
            echo "Running NVT equilibration..."
            gmx grompp -f ${{MDP_DIR}}/nvt.mdp -c ${{BUILD_DIR}}/em.gro -r ${{BUILD_DIR}}/em.gro -p topol.top -o ${{BUILD_DIR}}/nvt.tpr -maxwarn 20
            gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} ${{MDRUN_GPU_ARGS}} -gpu_id 0 -pinoffset "${{cpu_offset}}" -pinstride 1 -v
        fi
    else
        echo "Running NVT equilibration..."
        gmx grompp -f ${{MDP_DIR}}/nvt.mdp -c ${{BUILD_DIR}}/em.gro -r ${{BUILD_DIR}}/em.gro -p topol.top -o ${{BUILD_DIR}}/nvt.tpr -maxwarn 20
        gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} ${{MDRUN_GPU_ARGS}} -gpu_id 0 -pinoffset "${{cpu_offset}}" -pinstride 1 -v
    fi

    # NPT equilibration (must run to completion)
    if [ -f ${{BUILD_DIR}}/npt.gro ]; then
        echo "✓ NPT already completed, skipping..."
    elif [ -f ${{BUILD_DIR}}/npt.tpr ]; then
        # Check if NPT checkpoint is valid (has corresponding cpt file)
        if [ -f ${{BUILD_DIR}}/npt.cpt ]; then
            echo "NPT checkpoint found, resuming on GPU..."
            (
                cd ${{BUILD_DIR}}
                gmx mdrun -deffnm npt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} ${{MDRUN_GPU_ARGS}} -gpu_id 0 -pinoffset "${{cpu_offset}}" -pinstride 1 -cpi npt.cpt -v
            )
        else
            # NPT TPR exists but no checkpoint - previous run failed, clean up and restart
            echo "⚠️  Previous NPT run failed (no checkpoint), cleaning up and restarting..."
            rm -f ${{BUILD_DIR}}/npt.* 2>/dev/null || true
            echo "Running NPT equilibration..."
            gmx grompp -f ${{MDP_DIR}}/npt.mdp -c ${{BUILD_DIR}}/nvt.gro -r ${{BUILD_DIR}}/nvt.gro -t ${{BUILD_DIR}}/nvt.cpt -p topol.top -o ${{BUILD_DIR}}/npt.tpr -maxwarn 20
            gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} ${{MDRUN_GPU_ARGS}} -gpu_id 0 -pinoffset "${{cpu_offset}}" -pinstride 1 -v
        fi
    else
        echo "Running NPT equilibration..."
        gmx grompp -f ${{MDP_DIR}}/npt.mdp -c ${{BUILD_DIR}}/nvt.gro -r ${{BUILD_DIR}}/nvt.gro -t ${{BUILD_DIR}}/nvt.cpt -p topol.top -o ${{BUILD_DIR}}/npt.tpr -maxwarn 20
        gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp ${{OMP_NUM_THREADS}} ${{MDRUN_GPU_ARGS}} -gpu_id 0 -pinoffset "${{cpu_offset}}" -pinstride 1 -v
    fi

    echo "Running FEP production for all lambda windows..."
    echo "Execution mode: ${{EXECUTION_MODE}}"
    if [ "${{EXECUTION_MODE}}" = "repex" ]; then
        ./run_prod_repex.sh
    else
        ./run_prod_standard.sh
    fi

    # Clean up intermediate step PDB files
    echo "Cleaning up intermediate PDB files..."
    rm -f step*.pdb 2>/dev/null || true

    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                    ${{target_name}} LEG COMPLETE                                   ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"

    cd "${{FEP_DIR}}"
}}

# Main script logic
FEP_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"

# Number of replica runs (configured at build time)
REPLICAS={replicas}

# Parse target specification (e.g., bound1, bound1-3, unbound1-2,unbound3, etc.)
parse_targets() {{
    local spec=$1
    local targets=()

    # Check if it's a simple leg name (bound, unbound, bound, unbound with replicas)
    if [[ "$spec" =~ ^(bound|unbound)$ ]]; then
        if [ ${{REPLICAS}} -gt 1 ]; then
            # Expand to all replicas: bound -> bound1 bound2 bound3
            for i in $(seq 1 ${{REPLICAS}}); do
                targets+=("${{spec}}${{i}}")
            done
        else
            # No replicas, just the leg name
            targets+=("$spec")
        fi
    elif [[ "$spec" =~ ^(bound|unbound)([0-9]+)$ ]]; then
        # Single replica: bound1, unbound2, etc.
        targets+=("$spec")
    elif [[ "$spec" =~ ^(bound|unbound)([0-9]+)-([0-9]+)$ ]]; then
        # Range of replicas: bound1-3, unbound1-2, etc.
        local leg="${{BASH_REMATCH[1]}}"
        local start="${{BASH_REMATCH[2]}}"
        local end="${{BASH_REMATCH[3]}}"
        for i in $(seq $start $end); do
            targets+=("${{leg}}${{i}}")
        done
    else
        echo "Error: Invalid target specification: $spec"
        return 1
    fi

    # Output targets as space-separated list
    echo "${{targets[@]}}"
}}

# Default behavior with no arguments: run all configured legs/repeats.
if [ $# -eq 0 ]; then
    echo "No target specified; defaulting to 'all'."
    set -- all
fi

# Special case: 'all' means all bound and unbound replicas
if [ "$1" = "all" ]; then
    all_targets=()
    if [ ${{REPLICAS}} -gt 1 ]; then
        for i in $(seq 1 ${{REPLICAS}}); do
            all_targets+=("bound${{i}}")
            all_targets+=("unbound${{i}}")
        done
    else
        all_targets=("bound" "unbound")
    fi
else
    # Parse all target specifications
    all_targets=()
    for spec in "$@"; do
        targets=$(parse_targets "$spec")
        if [ $? -ne 0 ]; then
            exit 1
        fi
        for target in $targets; do
            all_targets+=("$target")
        done
    done
fi

# Run all requested targets
for target in "${{all_targets[@]}}"; do
    if resolve_leg_dir "${{target}}" >/dev/null; then
        run_leg "${{target}}"
    else
        echo "⚠️  Warning: Directory ${{target}} not found, skipping..."
    fi
done

echo ""
echo "✓ All requested FEP calculations completed."
""")

    script_path.write_text(script_content)
    script_path.chmod(0o755)


def write_fep_slurm_script(
    leg_dir: Path,
    leg_name: str,
    job_name: str = "fep",
    partition: str = "gpu",
    time: str = "48:00:00",
    nodes: int = 1,
    ntasks: int = 1,
    cpus_per_task: int = 16,
    gpus: int = 1,
) -> None:
    """
    Generate SLURM submission script for FEP calculations.

    Parameters
    ----------
    leg_dir : Path
        Path to the leg directory (bound or unbound)
    leg_name : str
        Name of the leg ('bound' or 'unbound')
    job_name : str
        SLURM job name
    partition : str
        SLURM partition name
    time : str
        Wall time limit
    nodes : int
        Number of nodes
    ntasks : int
        Number of tasks
    cpus_per_task : int
        CPUs per task
    gpus : int
        Number of GPUs
    """
    script_path = leg_dir / "submit.slurm.sh"

    script_content = f"""#!/bin/bash
#SBATCH -J {job_name}_{leg_name}
#SBATCH -p {partition}
#SBATCH --time={time}
#SBATCH --nodes={nodes}
#SBATCH --ntasks={ntasks}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --gres=gpu:{gpus}

# FEP {leg_name} leg SLURM submission script

# Load GROMACS module (adjust for your cluster)
# module load gromacs/2023.2

# Set OpenMP threads to match SLURM allocation
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Set number of lambda windows
NUM_WINDOWS=21

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Energy Minimization
echo 'Starting energy minimization...'
for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do
    if [ ! -f lambda_$i/em.gro ]; then
        gmx grompp -f lambda_$i/em.mdp -c conf.gro -r conf.gro -p topol.top -o lambda_$i/em.tpr -maxwarn 3
        gmx mdrun -s lambda_$i/em.tpr -deffnm lambda_$i/em -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK
    fi
done
echo 'Energy minimization completed.'

# NVT Equilibration
echo 'Starting NVT equilibration...'
for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do
    if [ ! -f lambda_$i/nvt.gro ]; then
        gmx grompp -f lambda_$i/nvt.mdp -c lambda_$i/em.gro -r lambda_$i/em.gro -p topol.top -o lambda_$i/nvt.tpr -maxwarn 2
        gmx mdrun -s lambda_$i/nvt.tpr -deffnm lambda_$i/nvt -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0
    fi
done
echo 'NVT equilibration completed.'

# NPT Equilibration (short)
echo 'Starting NPT equilibration (short)...'
for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do
    if [ -f lambda_$i/npt_short.mdp ]; then
        if [ ! -f lambda_$i/npt_short.gro ]; then
            gmx grompp -f lambda_$i/npt_short.mdp -c lambda_$i/nvt.gro -r lambda_$i/nvt.gro -p topol.top -o lambda_$i/npt_short.tpr -maxwarn 2
            gmx mdrun -s lambda_$i/npt_short.tpr -deffnm lambda_$i/npt_short -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0
        fi
    fi
done
echo 'NPT equilibration (short) completed.'

# NPT Equilibration
echo 'Starting NPT equilibration...'
for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do
    if [ ! -f lambda_$i/npt.gro ]; then
        # Use npt_short.gro if available, otherwise nvt.gro
        if [ -f lambda_$i/npt_short.gro ]; then
            INPUT_GRO=lambda_$i/npt_short.gro
        else
            INPUT_GRO=lambda_$i/nvt.gro
        fi
        gmx grompp -f lambda_$i/npt.mdp -c $INPUT_GRO -r $INPUT_GRO -p topol.top -o lambda_$i/npt.tpr -maxwarn 2
        gmx mdrun -s lambda_$i/npt.tpr -deffnm lambda_$i/npt -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0
    fi
done
echo 'NPT equilibration completed.'

# Production MD
echo 'Starting production MD...'
for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do
    if [ ! -f lambda_$i/prod.gro ]; then
        gmx grompp -f lambda_$i/prod.mdp -c lambda_$i/npt.gro -r lambda_$i/npt.gro -p topol.top -o lambda_$i/prod.tpr -maxwarn 2
        gmx mdrun -s lambda_$i/prod.tpr -deffnm lambda_$i/prod -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0
    fi
done
echo 'Production MD completed.'

echo "End time: $(date)"
echo 'All FEP calculations for {leg_name} leg completed.'
"""
    script_path.write_text(script_content)
    script_path.chmod(0o755)
