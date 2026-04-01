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

Usage for testing:
    export PRISM_MDRUN_NSTEPS=100  # Run only 100 steps for quick validation
    cd bound && bash localrun.sh

For production FEP calculations:
    unset PRISM_MDRUN_NSTEPS  # Ensure it's not set
    cd bound && bash localrun.sh

WARNING: Results obtained with PRISM_MDRUN_NSTEPS set are NOT valid for publication!
"""

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


def _get_execution_settings(config: Optional[dict]) -> dict:
    """Normalize execution-layer settings used by script generation."""
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

    return {
        "mode": mode,
        "use_gpu_pme": bool(exec_config.get("use_gpu_pme", True)),
        "num_gpus": max(1, int(num_gpus)),
        "parallel_windows": parallel_windows,
        "omp_threads": int(exec_config.get("omp_threads", 8)),  # GPU optimization: -ntmpi 1 -ntomp 8
        "replicas": int(fep_config.get("replicas", 1)),
    }


def _write_leg_execution_scripts(
    leg_dir: Path,
    execution_mode: str,
    *,
    use_gpu_pme: bool,
    num_gpus: int,
    parallel_windows: int,
    omp_threads: int,
) -> None:
    """Generate per-leg production helper scripts for standard and repex modes."""
    pme_gpu_flag = "-pme gpu" if use_gpu_pme else ""
    standard_script = leg_dir / "run_prod_standard.sh"
    repex_script = leg_dir / "run_prod_repex.sh"

    standard_content = f"""#!/usr/bin/env bash
set -euo pipefail

LEG_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
MDP_DIR="${{LEG_DIR}}/mdps"
BUILD_DIR="${{LEG_DIR}}/build"

if [ -n "${{PRISM_MDRUN_NSTEPS:-}}" ]; then
    MDRUN_NSTEPS_ARG="-nsteps ${{PRISM_MDRUN_NSTEPS}}"
else
    MDRUN_NSTEPS_ARG=""
fi

echo "Production mode: standard"
echo "Concurrent windows: {parallel_windows}"
echo "Configured GPUs: {num_gpus}"

JOB_COUNT=0
for lambda_mdp in "${{MDP_DIR}}"/prod_*.mdp; do
    lambda_name=$(basename "${{lambda_mdp}}" .mdp)
    lambda_idx="${{lambda_name#prod_}}"
    window_dir="${{LEG_DIR}}/window_${{lambda_idx}}"
    mkdir -p "${{window_dir}}"

    if [ -f "${{window_dir}}/prod.gro" ]; then
        echo "  ${{lambda_name}}: already completed"
        continue
    fi

    if [ {parallel_windows} -gt 1 ]; then
        running_jobs=$(find "${{LEG_DIR}}" -maxdepth 1 -name ".run_*" -type f 2>/dev/null | wc -l)
        while [ "${{running_jobs}}" -ge {parallel_windows} ]; do
            sleep 1
            running_jobs=$(find "${{LEG_DIR}}" -maxdepth 1 -name ".run_*" -type f 2>/dev/null | wc -l)
        done
    fi

    if [ {num_gpus} -gt 1 ]; then
        gpu_id=$((JOB_COUNT % {num_gpus}))
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

        echo "$$" > "${{LEG_DIR}}/.run_${{lambda_idx}}"
        trap 'rm -f "${{LEG_DIR}}/.run_${{lambda_idx}}"' EXIT

        if [ -f "${{window_dir}}/npt_short.gro" ]; then
            echo "  ${{lambda_name}}: NPT short already completed"
        elif [ -f "${{window_dir}}/npt_short.tpr" ] && [ -f "${{window_dir}}/npt_short.cpt" ]; then
            echo "  ${{lambda_name}}: resuming NPT short"
            gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp {omp_threads} -nb gpu -bonded gpu {pme_gpu_flag} -cpi "${{window_dir}}/npt_short.cpt" -v ${{MDRUN_NSTEPS_ARG}} || \
                gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp 10 -cpi "${{window_dir}}/npt_short.cpt" -v ${{MDRUN_NSTEPS_ARG}}
        else
            echo "  ${{lambda_name}}: starting NPT short on GPU ${{gpu_id}}"
            gmx grompp -f "${{MDP_DIR}}/npt_short_${{lambda_idx}}.mdp" -c "${{BUILD_DIR}}/npt.gro" -r "${{BUILD_DIR}}/npt.gro" -t "${{BUILD_DIR}}/npt.cpt" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/npt_short.tpr" -maxwarn 20
            gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp {omp_threads} -nb gpu -bonded gpu {pme_gpu_flag} -v ${{MDRUN_NSTEPS_ARG}} || \\
                gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp 10 -v ${{MDRUN_NSTEPS_ARG}}
        fi

        if [ -f "${{window_dir}}/prod.gro" ]; then
            echo "  ${{lambda_name}}: production already completed"
        elif [ -f "${{window_dir}}/prod.tpr" ] && [ -f "${{window_dir}}/prod.cpt" ]; then
            echo "  ${{lambda_name}}: resuming production"
            gmx mdrun -deffnm "${{window_dir}}/prod" -cpi "${{window_dir}}/prod.cpt" -v ${{MDRUN_NSTEPS_ARG}}
        else
            echo "  ${{lambda_name}}: starting production on GPU ${{gpu_id}}"
            if [ ! -f "${{window_dir}}/prod.tpr" ]; then
                gmx grompp -f "${{lambda_mdp}}" -c "${{window_dir}}/npt_short.gro" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/prod.tpr" -maxwarn 20
            fi
            gmx mdrun -deffnm "${{window_dir}}/prod" -ntmpi 1 -ntomp {omp_threads} -nb gpu -bonded gpu {pme_gpu_flag} -v ${{MDRUN_NSTEPS_ARG}} || \\
                gmx mdrun -deffnm "${{window_dir}}/prod" -ntmpi 1 -ntomp 10 -v ${{MDRUN_NSTEPS_ARG}}
        fi

        rm -f "${{LEG_DIR}}/.run_${{lambda_idx}}"
        trap - EXIT
        echo "  ${{lambda_name}}: ✓ Complete"
    }}

    if [ {parallel_windows} -gt 1 ]; then
        run_window "${{lambda_idx}}" "${{lambda_name}}" "${{lambda_mdp}}" "${{window_dir}}" "${{gpu_id}}" &
        JOB_COUNT=$((JOB_COUNT + 1))
    else
        run_window "${{lambda_idx}}" "${{lambda_name}}" "${{lambda_mdp}}" "${{window_dir}}" "${{gpu_id}}"
    fi
done

if [ {parallel_windows} -gt 1 ]; then
    echo "Waiting for all windows to complete..."
    wait
fi
"""

    repex_content = f"""#!/usr/bin/env bash
set -euo pipefail

LEG_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
MDP_DIR="${{LEG_DIR}}/mdps"
BUILD_DIR="${{LEG_DIR}}/build"

if [ -n "${{PRISM_MDRUN_NSTEPS:-}}" ]; then
    MDRUN_NSTEPS_ARG="-nsteps ${{PRISM_MDRUN_NSTEPS}}"
else
    MDRUN_NSTEPS_ARG=""
fi

NUM_GPUS={num_gpus}

window_dirs=()
for lambda_mdp in "${{MDP_DIR}}"/prod_*.mdp; do
    lambda_name=$(basename "${{lambda_mdp}}" .mdp)
    lambda_idx="${{lambda_name#prod_}}"
    window_dir="${{LEG_DIR}}/window_${{lambda_idx}}"
    window_dirs+=("${{window_dir}}")
    mkdir -p "${{window_dir}}"

    if [ ! -f "${{window_dir}}/npt_short.gro" ]; then
        if [ ! -f "${{window_dir}}/npt_short.tpr" ]; then
            gmx grompp -f "${{MDP_DIR}}/npt_short_${{lambda_idx}}.mdp" -c "${{BUILD_DIR}}/npt.gro" -r "${{BUILD_DIR}}/npt.gro" -t "${{BUILD_DIR}}/npt.cpt" -p "${{LEG_DIR}}/topol.top" -o "${{window_dir}}/npt_short.tpr" -maxwarn 20
        fi
        gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp {omp_threads} -nb gpu -bonded gpu {pme_gpu_flag} -v ${{MDRUN_NSTEPS_ARG}} || \\
            gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp 10 -v ${{MDRUN_NSTEPS_ARG}}
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
echo "GPU id string: $gpu_ids"

mpirun -oversubscribe -np "${{num_windows}}" \\
    gmx_mpi mdrun -v -deffnm prod -nb gpu -bonded gpu {pme_gpu_flag} -replex 1000 \\
    -multidir "${{window_dirs[@]}}" -gpu_id "${{gpu_ids}}" -pin on -dhdl dhdl \\
    "${{cpi_args[@]}}" ${{MDRUN_NSTEPS_ARG}}
"""

    for path, content in ((standard_script, standard_content), (repex_script, repex_content)):
        with open(path, "w") as handle:
            handle.write(content)
        path.chmod(0o755)

    mode_file = leg_dir / "run_prod.sh"
    mode_file.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f'exec "$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)/run_prod_{execution_mode}.sh" "$@"\n'
    )
    mode_file.chmod(0o755)


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
    execution.use_gpu_pme : bool
        Enable GPU PME (default: True)
    execution.mode : str
        Production execution mode: ``standard`` (independent windows) or
        ``repex`` (lambda replica exchange via ``gmx_mpi -multidir``)
    execution.num_gpus : int
        Number of GPUs available (default: 1)
    execution.parallel_windows : int
        Number of concurrent lambda windows in ``standard`` mode. In ``repex``
        mode, all configured GPUs are assigned to one multidir production job.
    execution.omp_threads : int
        OpenMP threads per GPU (default: 14)
    fep.replicas : int
        Number of replica runs for error estimation (default: 1)
        When >1, creates bound1, bound2, ... and unbound1, unbound2, ... directories
    """
    settings = _get_execution_settings(config)
    for leg_dir in (layout.bound_dir, layout.unbound_dir):
        _write_leg_execution_scripts(
            leg_dir,
            settings["mode"],
            use_gpu_pme=settings["use_gpu_pme"],
            num_gpus=settings["num_gpus"],
            parallel_windows=settings["parallel_windows"],
            omp_threads=settings["omp_threads"],
        )
    write_fep_master_script(layout.root, config)


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
    num_gpus = settings["num_gpus"]
    parallel_windows = settings["parallel_windows"]
    omp_threads = settings["omp_threads"]
    replicas = settings["replicas"]
    pme_gpu_flag = "-pme gpu" if use_gpu_pme else ""

    script_path = fep_dir / "run_fep.sh"

    script_content = f"""#!/usr/bin/env bash
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

EXECUTION_MODE="${{PRISM_FEP_MODE:-{execution_mode}}}"
if [[ "${{EXECUTION_MODE}}" != "standard" && "${{EXECUTION_MODE}}" != "repex" ]]; then
    echo "Error: unsupported execution mode '${{EXECUTION_MODE}}'"
    echo "Supported modes: standard, repex"
    exit 1
fi

# Function to run a single leg
run_leg() {{
    local leg_name=$1
    local leg_dir="${{FEP_DIR}}/${{leg_name}}"

    if [ ! -d "${{leg_dir}}" ]; then
        echo "Error: Leg directory not found: ${{leg_dir}}"
        return 1
    fi

    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                    PRISM-FEP ${{leg_name}} LEG                                   ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
    echo ""

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
        gmx mdrun -deffnm ${{BUILD_DIR}}/em -cpi ${{BUILD_DIR}}/em.cpt -v
    else
        echo "Running energy minimization..."
        gmx grompp -f ${{MDP_DIR}}/em.mdp -c ${{INPUT_DIR}}/conf.gro -p topol.top -o ${{BUILD_DIR}}/em.tpr -maxwarn 20
        gmx mdrun -deffnm ${{BUILD_DIR}}/em -ntmpi 1 -ntomp 10 -v
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
        echo "NVT checkpoint found, resuming..."
        (
            cd ${{BUILD_DIR}}
            gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -cpi nvt.cpt -v || gmx mdrun -deffnm nvt -ntmpi 1 -ntomp 10 -cpi nvt.cpt -v
        )
    else
        echo "Running NVT equilibration..."
        gmx grompp -f ${{MDP_DIR}}/nvt.mdp -c ${{BUILD_DIR}}/em.gro -r ${{BUILD_DIR}}/em.gro -p topol.top -o ${{BUILD_DIR}}/nvt.tpr -maxwarn 20
        gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v || gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -ntmpi 1 -ntomp 10 -v
    fi

    # NPT equilibration (must run to completion)
    if [ -f ${{BUILD_DIR}}/npt.gro ]; then
        echo "✓ NPT already completed, skipping..."
    elif [ -f ${{BUILD_DIR}}/npt.tpr ]; then
        # Check if NPT checkpoint is valid (has corresponding cpt file)
        if [ -f ${{BUILD_DIR}}/npt.cpt ]; then
            echo "NPT checkpoint found, resuming..."
            (
                cd ${{BUILD_DIR}}
                gmx mdrun -deffnm npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -cpi npt.cpt -v || gmx mdrun -deffnm npt -ntmpi 1 -ntomp 10 -cpi npt.cpt -v
            )
        else
            # NPT TPR exists but no checkpoint - previous run failed, clean up and restart
            echo "⚠️  Previous NPT run failed (no checkpoint), cleaning up and restarting..."
            rm -f ${{BUILD_DIR}}/npt.* 2>/dev/null || true
            echo "Running NPT equilibration..."
            gmx grompp -f ${{MDP_DIR}}/npt.mdp -c ${{BUILD_DIR}}/nvt.gro -r ${{BUILD_DIR}}/nvt.gro -t ${{BUILD_DIR}}/nvt.cpt -p topol.top -o ${{BUILD_DIR}}/npt.tpr -maxwarn 20
            gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v || gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 10 -v
        fi
    else
        echo "Running NPT equilibration..."
        gmx grompp -f ${{MDP_DIR}}/npt.mdp -c ${{BUILD_DIR}}/nvt.gro -r ${{BUILD_DIR}}/nvt.gro -t ${{BUILD_DIR}}/nvt.cpt -p topol.top -o ${{BUILD_DIR}}/npt.tpr -maxwarn 20
        gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v || gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 10 -v
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
    echo "║                    ${{leg_name}} LEG COMPLETE                                   ║"
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

# Show usage if no argument or invalid argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <target> [target ...]"
    echo ""
    echo "Target specifications:"
    echo "  bound           - Run all bound replicas (bound1, bound2, ...)"
    echo "  unbound         - Run all unbound replicas (unbound1, unbound2, ...)"
    echo "  all             - Run all bound and unbound replicas"
    echo "  bound1          - Run specific replica"
    echo "  bound1-3        - Run range of replicas (bound1, bound2, bound3)"
    echo "  unbound1-2      - Run range of replicas"
    echo "  bound1,unbound2 - Mix different legs and replicas"
    echo "  bound1-3,unbound1-2 - Mix ranges"
    echo ""
    echo "Configuration:"
    echo "  Replicas configured: ${{REPLICAS}}"
    if [ ${{REPLICAS}} -gt 1 ]; then
        echo "  Available: bound1..bound${{REPLICAS}}, unbound1..unbound${{REPLICAS}}"
    else
        echo "  Available: bound, unbound"
    fi
    echo ""
    echo "Environment variables:"
    echo "  PRISM_MDRUN_NSTEPS - Limit production steps (testing only, DO NOT use for production)"
    echo "  PRISM_FEP_MODE     - Override execution mode (standard or repex)"
    exit 1
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
    if [ -d "${{FEP_DIR}}/${{target}}" ]; then
        run_leg "${{target}}"
    else
        echo "⚠️  Warning: Directory ${{target}} not found, skipping..."
    fi
done

echo ""
echo "✓ All requested FEP calculations completed."
"""

    with open(script_path, "w") as f:
        f.write(script_content)

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

    with open(script_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH -J {job_name}_{leg_name}\n")
        f.write(f"#SBATCH -p {partition}\n")
        f.write(f"#SBATCH --time={time}\n")
        f.write(f"#SBATCH --nodes={nodes}\n")
        f.write(f"#SBATCH --ntasks={ntasks}\n")
        f.write(f"#SBATCH --cpus-per-task={cpus_per_task}\n")
        f.write(f"#SBATCH --gres=gpu:{gpus}\n\n")

        f.write(f"# FEP {leg_name} leg SLURM submission script\n\n")
        f.write("# Load GROMACS module (adjust for your cluster)\n")
        f.write("# module load gromacs/2023.2\n\n")

        f.write("# Set OpenMP threads to match SLURM allocation\n")
        f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n\n")

        f.write("# Set number of lambda windows\n")
        f.write("NUM_WINDOWS=21\n\n")

        f.write('echo "Start time: $(date)"\n')
        f.write('echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"\n')
        f.write('echo "hostname: $(hostname)"\n')
        f.write('echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"\n')
        f.write('echo "Job directory: $(pwd)"\n\n')

        f.write("# Energy Minimization\n")
        f.write("echo 'Starting energy minimization...'\n")
        f.write("for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do\n")
        f.write("    if [ ! -f lambda_$i/em.gro ]; then\n")
        f.write(
            "        gmx grompp -f lambda_$i/em.mdp -c conf.gro -r conf.gro -p topol.top -o lambda_$i/em.tpr -maxwarn 3\n"
        )
        f.write("        gmx mdrun -s lambda_$i/em.tpr -deffnm lambda_$i/em -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK\n")
        f.write("    fi\n")
        f.write("done\n")
        f.write("echo 'Energy minimization completed.'\n\n")

        f.write("# NVT Equilibration\n")
        f.write("echo 'Starting NVT equilibration...'\n")
        f.write("for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do\n")
        f.write("    if [ ! -f lambda_$i/nvt.gro ]; then\n")
        f.write(
            "        gmx grompp -f lambda_$i/nvt.mdp -c lambda_$i/em.gro -r lambda_$i/em.gro -p topol.top -o lambda_$i/nvt.tpr -maxwarn 2\n"
        )
        f.write(
            "        gmx mdrun -s lambda_$i/nvt.tpr -deffnm lambda_$i/nvt -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0\n"
        )
        f.write("    fi\n")
        f.write("done\n")
        f.write("echo 'NVT equilibration completed.'\n\n")

        f.write("# NPT Equilibration (short)\n")
        f.write("echo 'Starting NPT equilibration (short)...'\n")
        f.write("for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do\n")
        f.write("    if [ -f lambda_$i/npt_short.mdp ]; then\n")
        f.write("        if [ ! -f lambda_$i/npt_short.gro ]; then\n")
        f.write(
            "            gmx grompp -f lambda_$i/npt_short.mdp -c lambda_$i/nvt.gro -r lambda_$i/nvt.gro -p topol.top -o lambda_$i/npt_short.tpr -maxwarn 2\n"
        )
        f.write(
            "            gmx mdrun -s lambda_$i/npt_short.tpr -deffnm lambda_$i/npt_short -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0\n"
        )
        f.write("        fi\n")
        f.write("    fi\n")
        f.write("done\n")
        f.write("echo 'NPT equilibration (short) completed.'\n\n")

        f.write("# NPT Equilibration\n")
        f.write("echo 'Starting NPT equilibration...'\n")
        f.write("for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do\n")
        f.write("    if [ ! -f lambda_$i/npt.gro ]; then\n")
        f.write("        # Use npt_short.gro if available, otherwise nvt.gro\n")
        f.write("        if [ -f lambda_$i/npt_short.gro ]; then\n")
        f.write("            INPUT_GRO=lambda_$i/npt_short.gro\n")
        f.write("        else\n")
        f.write("            INPUT_GRO=lambda_$i/nvt.gro\n")
        f.write("        fi\n")
        f.write(
            "        gmx grompp -f lambda_$i/npt.mdp -c $INPUT_GRO -r $INPUT_GRO -p topol.top -o lambda_$i/npt.tpr -maxwarn 2\n"
        )
        f.write(
            "        gmx mdrun -s lambda_$i/npt.tpr -deffnm lambda_$i/npt -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0\n"
        )
        f.write("    fi\n")
        f.write("done\n")
        f.write("echo 'NPT equilibration completed.'\n\n")

        f.write("# Production MD\n")
        f.write("echo 'Starting production MD...'\n")
        f.write("for i in $(seq -f '%02g' 0 $((NUM_WINDOWS-1))); do\n")
        f.write("    if [ ! -f lambda_$i/prod.gro ]; then\n")
        f.write(
            "        gmx grompp -f lambda_$i/prod.mdp -c lambda_$i/npt.gro -r lambda_$i/npt.gro -p topol.top -o lambda_$i/prod.tpr -maxwarn 2\n"
        )
        f.write(
            "        gmx mdrun -s lambda_$i/prod.tpr -deffnm lambda_$i/prod -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu -bonded gpu -pme gpu -gpu_id 0\n"
        )
        f.write("    fi\n")
        f.write("done\n")
        f.write("echo 'Production MD completed.'\n\n")

        f.write('echo "End time: $(date)"\n')
        f.write(f"echo 'All FEP calculations for {leg_name} leg completed.'\n")

    script_path.chmod(0o755)
