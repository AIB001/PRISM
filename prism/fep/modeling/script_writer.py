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


def write_root_scripts(layout) -> None:
    """
    Write root-level run scripts (run_all.sh, submit_all.slurm.sh).

    Parameters
    ----------
    layout : FEPScaffoldLayout
        The scaffold layout containing paths to bound/unbound legs
    """
    # Write run_all.sh
    run_all_script = layout.root / "run_all.sh"
    with open(run_all_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Run all FEP calculations\n\n")
        f.write("echo 'Starting bound leg...'\n")
        f.write("cd bound && bash localrun.sh\n")
        f.write("cd ..\n\n")
        f.write("echo 'Starting unbound leg...'\n")
        f.write("cd unbound && bash localrun.sh\n")
        f.write("cd ..\n\n")
        f.write("echo 'All FEP calculations completed.'\n")
    run_all_script.chmod(0o755)

    # Write submit_all.slurm.sh
    submit_all_script = layout.root / "submit_all.slurm.sh"
    with open(submit_all_script, "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Submit all FEP calculations to SLURM\n\n")
        f.write("echo 'Submitting bound leg...'\n")
        f.write("cd bound && sbatch submit.slurm.sh\n")
        f.write("cd ..\n\n")
        f.write("echo 'Submitting unbound leg...'\n")
        f.write("cd unbound && sbatch submit.slurm.sh\n")
        f.write("cd ..\n\n")
        f.write("echo 'All jobs submitted.'\n")
    submit_all_script.chmod(0o755)


def write_fep_run_script(leg_dir: Path, leg_name: str, config: Optional[dict] = None) -> None:
    """
    Generate localrun.sh script for FEP calculations.

    Parameters
    ----------
    leg_dir : Path
        Path to the leg directory (bound or unbound)
    leg_name : str
        Name of the leg ('bound' or 'unbound')
    config : dict, optional
        Configuration dictionary with GPU settings
    """
    # Get GPU settings from config (default to GPU enabled for backward compatibility)
    config = config or {}
    exec_config = config.get("execution", {})
    use_gpu_pme = exec_config.get("use_gpu_pme", True)

    # Build GPU PME flag for mdrun commands
    pme_gpu_flag = "-pme gpu" if use_gpu_pme else ""

    script_path = leg_dir / "localrun.sh"

    with open(script_path, "w") as f:
        # Script header with environment variable support
        f.write("#!/usr/bin/env bash\n")
        f.write("# GPU PME: " + ("enabled" if use_gpu_pme else "disabled (CPU PME for small systems)") + "\n")
        f.write("set -euo pipefail\n\n")
        f.write('LEG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"\n')
        f.write('MDP_DIR="${LEG_DIR}/mdps"\n')
        f.write('INPUT_DIR="${LEG_DIR}/input"\n')
        f.write('BUILD_DIR="${LEG_DIR}/build"\n')
        # WARNING: PRISM_MDRUN_NSTEPS should ONLY be used for quick testing of production runs
        # It is intentionally NOT applied to EM/NVT/NPT equilibration - these must complete fully
        f.write("# Test mode: limit production run steps (DO NOT use for production FEP calculations)\n")
        f.write('if [ -n "${PRISM_MDRUN_NSTEPS:-}" ]; then\n')
        f.write('    echo "⚠️  WARNING: PRISM_MDRUN_NSTEPS is set to ${PRISM_MDRUN_NSTEPS}"\n')
        f.write('    echo "⚠️  This is for TESTING ONLY - results will NOT be valid for publication!"\n')
        f.write('    MDRUN_NSTEPS_ARG="-nsteps ${PRISM_MDRUN_NSTEPS}"\n')
        f.write("else\n")
        f.write('    MDRUN_NSTEPS_ARG=""\n')
        f.write("fi\n\n")

        # Banner
        f.write('echo "╔═══════════════════════════════════════════════════════════════════════════╗"\n')
        f.write(f'echo "║                    PRISM-FEP {leg_name.upper()} LEG                                    ║"\n')
        f.write('echo "╚═══════════════════════════════════════════════════════════════════════════╝"\n')
        f.write('echo ""\n\n')

        f.write('cd "${LEG_DIR}"\n\n')

        # Energy Minimization
        f.write("# Energy minimization\n")
        f.write("mkdir -p ${BUILD_DIR}\n\n")
        f.write("if [ -f ${BUILD_DIR}/em.gro ]; then\n")
        f.write('    echo "✓ EM already completed, skipping..."\n')
        f.write("elif [ -f ${BUILD_DIR}/em.tpr ]; then\n")
        f.write('    echo "EM checkpoint found, resuming..."\n')
        f.write("    gmx mdrun -deffnm ${BUILD_DIR}/em -cpi ${BUILD_DIR}/em.cpt -v\n")
        f.write("else\n")
        f.write('    echo "Running energy minimization..."\n')
        f.write(
            "    gmx grompp -f ${MDP_DIR}/em.mdp -c ${INPUT_DIR}/conf.gro -p topol.top -o ${BUILD_DIR}/em.tpr -maxwarn 2\n"
        )
        f.write("    gmx mdrun -deffnm ${BUILD_DIR}/em -ntmpi 1 -ntomp 10 -v\n")
        f.write("fi\n\n")
        # Validate EM convergence
        f.write("# Check EM convergence\n")
        f.write("if [ -f ${BUILD_DIR}/em.log ]; then\n")
        f.write(
            '    max_force=$(grep -E "Maximum force|Fmax=" ${BUILD_DIR}/em.log | tail -1 | grep -oE "[0-9]+\\.?[0-9]*e[+-][0-9]+" | tail -1 || true)\n'
        )
        f.write('    if [ -n "$max_force" ]; then\n')
        f.write('        echo "EM final max force: $max_force kJ/mol/nm"\n')
        f.write("        # Warning if force is still high (> 1000 kJ/mol/nm)\n")
        f.write('        if (( $(echo "$max_force > 1000" | bc -l) )); then\n')
        f.write('            echo "⚠️  WARNING: EM converged with high force. System may need more minimization."\n')
        f.write("        fi\n")
        f.write("    fi\n")
        f.write("fi\n\n")

        # NVT equilibration
        f.write("# NVT equilibration (must run to completion)\n")
        f.write("if [ -f ${BUILD_DIR}/nvt.gro ]; then\n")
        f.write('    echo "✓ NVT already completed, skipping..."\n')
        f.write("elif [ -f ${BUILD_DIR}/nvt.tpr ]; then\n")
        f.write('    echo "NVT checkpoint found, resuming..."\n')
        f.write("    gmx mdrun -deffnm ${BUILD_DIR}/nvt -cpi ${BUILD_DIR}/nvt.cpt -v\n")
        f.write("else\n")
        f.write('    echo "Running NVT equilibration..."\n')
        f.write(
            "    gmx grompp -f ${MDP_DIR}/nvt.mdp -c ${BUILD_DIR}/em.gro -r ${BUILD_DIR}/em.gro -p topol.top -o ${BUILD_DIR}/nvt.tpr -maxwarn 2\n"
        )
        f.write(
            f"    gmx mdrun -deffnm ${{BUILD_DIR}}/nvt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v || "
        )
        f.write("gmx mdrun -deffnm ${BUILD_DIR}/nvt -ntmpi 1 -ntomp 10 -v\n")
        f.write("fi\n\n")

        # NPT equilibration
        f.write("# NPT equilibration (must run to completion)\n")
        f.write("if [ -f ${BUILD_DIR}/npt.gro ]; then\n")
        f.write('    echo "✓ NPT already completed, skipping..."\n')
        f.write("elif [ -f ${BUILD_DIR}/npt.tpr ]; then\n")
        f.write("    # Check if NPT checkpoint is valid (has corresponding cpt file)\n")
        f.write("    if [ -f ${BUILD_DIR}/npt.cpt ]; then\n")
        f.write('        echo "NPT checkpoint found, resuming..."\n')
        f.write("        gmx mdrun -deffnm ${BUILD_DIR}/npt -cpi ${BUILD_DIR}/npt.cpt -v\n")
        f.write("    else\n")
        f.write("        # NPT TPR exists but no checkpoint - previous run failed, clean up and restart\n")
        f.write('        echo "⚠️  Previous NPT run failed (no checkpoint), cleaning up and restarting..."\n')
        f.write("        rm -f ${BUILD_DIR}/npt.* 2>/dev/null || true\n")
        f.write('        echo "Running NPT equilibration..."\n')
        f.write(
            "        gmx grompp -f ${MDP_DIR}/npt.mdp -c ${BUILD_DIR}/nvt.gro -r ${BUILD_DIR}/nvt.gro -t ${BUILD_DIR}/nvt.cpt -p topol.top -o ${BUILD_DIR}/npt.tpr -maxwarn 2\n"
        )
        f.write(
            f"        gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v || "
        )
        f.write("gmx mdrun -deffnm ${BUILD_DIR}/npt -ntmpi 1 -ntomp 10 -v\n")
        f.write("    fi\n")
        f.write("else\n")
        f.write('    echo "Running NPT equilibration..."\n')
        f.write(
            "    gmx grompp -f ${MDP_DIR}/npt.mdp -c ${BUILD_DIR}/nvt.gro -r ${BUILD_DIR}/nvt.gro -t ${BUILD_DIR}/nvt.cpt -p topol.top -o ${BUILD_DIR}/npt.tpr -maxwarn 2\n"
        )
        f.write(
            f"    gmx mdrun -deffnm ${{BUILD_DIR}}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v || "
        )
        f.write("gmx mdrun -deffnm ${BUILD_DIR}/npt -ntmpi 1 -ntomp 10 -v\n")
        f.write("fi\n\n")

        # FEP production runs for each lambda window
        f.write("# FEP production runs for each lambda window\n")
        f.write('echo "Running FEP production for all lambda windows..."\n')
        f.write("for lambda_mdp in ${MDP_DIR}/prod_*.mdp; do\n")
        f.write("    lambda_name=$(basename ${lambda_mdp} .mdp)\n")
        f.write("    lambda_idx=$(echo ${lambda_name} | sed 's/prod_//')\n")
        f.write('    window_dir="window_${lambda_idx}"\n\n')
        f.write('    mkdir -p "${window_dir}"\n\n')
        f.write("    # Per-window NPT short equilibration (quick, can use test mode)\n")
        f.write('    if [ -f "${window_dir}/npt_short.gro" ]; then\n')
        f.write('        echo "  ${lambda_name}: NPT short already completed"\n')
        f.write('    elif [ -f "${window_dir}/npt_short.tpr" ]; then\n')
        f.write('        echo "  ${lambda_name}: NPT short resuming from checkpoint"\n')
        f.write(
            '        gmx mdrun -deffnm "${window_dir}/npt_short" -cpi "${window_dir}/npt_short.cpt" -v ${MDRUN_NSTEPS_ARG}\n'
        )
        f.write("    else\n")
        f.write('        echo "  ${lambda_name}: starting NPT short equilibration"\n')
        f.write(
            '        gmx grompp -f ${MDP_DIR}/npt_short_${lambda_idx}.mdp -c ${BUILD_DIR}/npt.gro -r ${BUILD_DIR}/npt.gro -t ${BUILD_DIR}/nvt.cpt -p topol.top -o "${window_dir}/npt_short.tpr" -maxwarn 2\n'
        )
        f.write(
            f'        gmx mdrun -deffnm "${{window_dir}}/npt_short" -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v ${{MDRUN_NSTEPS_ARG}} || '
        )
        f.write('gmx mdrun -deffnm "${window_dir}/npt_short" -ntmpi 1 -ntomp 10 -v ${MDRUN_NSTEPS_ARG}\n')
        f.write("    fi\n\n")
        f.write("    # Production run (can use test mode for quick validation)\n")
        f.write('    if [ -f "${window_dir}/prod.gro" ]; then\n')
        f.write('        echo "  ${lambda_name}: production already completed"\n')
        f.write('    elif [ -f "${window_dir}/prod.tpr" ]; then\n')
        f.write('        echo "  ${lambda_name}: production resuming from checkpoint"\n')
        f.write('        gmx mdrun -deffnm "${window_dir}/prod" -cpi "${window_dir}/prod.cpt" -v ${MDRUN_NSTEPS_ARG}\n')
        f.write("    else\n")
        f.write('        echo "  ${lambda_name}: starting production"\n')
        f.write(
            '        gmx grompp -f ${lambda_mdp} -c "${window_dir}/npt_short.gro" -p topol.top -o "${window_dir}/prod.tpr" -maxwarn 2\n'
        )
        f.write(
            f'        gmx mdrun -deffnm "${{window_dir}}/prod" -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu {pme_gpu_flag} -v ${{MDRUN_NSTEPS_ARG}} || '
        )
        f.write('gmx mdrun -deffnm "${window_dir}/prod" -ntmpi 1 -ntomp 10 -v ${MDRUN_NSTEPS_ARG}\n')
        f.write("    fi\n")
        f.write("done\n\n")

        # Cleanup and completion message
        f.write("# Clean up intermediate step PDB files\n")
        f.write('echo "Cleaning up intermediate PDB files..."\n')
        f.write("rm -f step*.pdb 2>/dev/null || true\n\n")
        f.write('echo ""\n')
        f.write('echo "╔═══════════════════════════════════════════════════════════════════════════╗"\n')
        f.write(f'echo "║                    {leg_name.upper()} LEG COMPLETE                                   ║"\n')
        f.write('echo "╚═══════════════════════════════════════════════════════════════════════════╝"\n')

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
