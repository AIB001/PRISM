#!/usr/bin/env bash
set -euo pipefail

LEG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDP_DIR="${LEG_DIR}/mdps"
BUILD_DIR="${LEG_DIR}/build"

if [ -n "${PRISM_MDRUN_NSTEPS:-}" ]; then
    MDRUN_NSTEPS_ARG="-nsteps ${PRISM_MDRUN_NSTEPS}"
else
    MDRUN_NSTEPS_ARG=""
fi

echo "Production mode: standard"
echo "Concurrent windows: 4"
echo "Configured GPUs: 4"

JOB_COUNT=0
for lambda_mdp in "${MDP_DIR}"/prod_*.mdp; do
    lambda_name=$(basename "${lambda_mdp}" .mdp)
    lambda_idx="${lambda_name#prod_}"
    window_dir="${LEG_DIR}/window_${lambda_idx}"
    mkdir -p "${window_dir}"

    if [ -f "${window_dir}/prod.gro" ]; then
        echo "  ${lambda_name}: already completed"
        continue
    fi

    if [ 4 -gt 1 ]; then
        running_jobs=$(find "${LEG_DIR}" -maxdepth 1 -name ".run_*" -type f 2>/dev/null | wc -l)
        while [ "${running_jobs}" -ge 4 ]; do
            sleep 1
            running_jobs=$(find "${LEG_DIR}" -maxdepth 1 -name ".run_*" -type f 2>/dev/null | wc -l)
        done
    fi

    if [ 4 -gt 1 ]; then
        gpu_id=$((JOB_COUNT % 4))
    else
        gpu_id=0
    fi

    run_window() {
        local lambda_idx="$1"
        local lambda_name="$2"
        local lambda_mdp="$3"
        local window_dir="$4"
        local gpu_id="$5"

        export CUDA_VISIBLE_DEVICES="${gpu_id}"
        export OMP_NUM_THREADS=2

        echo "$$" > "${LEG_DIR}/.run_${lambda_idx}"
        trap 'rm -f "${LEG_DIR}/.run_${lambda_idx}"' EXIT

        if [ -f "${window_dir}/npt_short.gro" ]; then
            echo "  ${lambda_name}: NPT short already completed"
        elif [ -f "${window_dir}/npt_short.tpr" ] && [ -f "${window_dir}/npt_short.cpt" ]; then
            echo "  ${lambda_name}: resuming NPT short"
            gmx mdrun -deffnm "${window_dir}/npt_short" -cpi "${window_dir}/npt_short.cpt" -v ${MDRUN_NSTEPS_ARG}
        else
            echo "  ${lambda_name}: starting NPT short on GPU ${gpu_id}"
            gmx grompp -f "${MDP_DIR}/npt_short_${lambda_idx}.mdp" -c "${BUILD_DIR}/npt.gro" -r "${BUILD_DIR}/npt.gro" -t "${BUILD_DIR}/nvt.cpt" -p "${LEG_DIR}/topol.top" -o "${window_dir}/npt_short.tpr" -maxwarn 10
            gmx mdrun -deffnm "${window_dir}/npt_short" -ntmpi 1 -ntomp 2 -nb gpu -bonded gpu -pme gpu -v ${MDRUN_NSTEPS_ARG} || \
                gmx mdrun -deffnm "${window_dir}/npt_short" -ntmpi 1 -ntomp 10 -v ${MDRUN_NSTEPS_ARG}
        fi

        if [ -f "${window_dir}/prod.gro" ]; then
            echo "  ${lambda_name}: production already completed"
        elif [ -f "${window_dir}/prod.tpr" ] && [ -f "${window_dir}/prod.cpt" ]; then
            echo "  ${lambda_name}: resuming production"
            gmx mdrun -deffnm "${window_dir}/prod" -cpi "${window_dir}/prod.cpt" -v ${MDRUN_NSTEPS_ARG}
        else
            echo "  ${lambda_name}: starting production on GPU ${gpu_id}"
            if [ ! -f "${window_dir}/prod.tpr" ]; then
                gmx grompp -f "${lambda_mdp}" -c "${window_dir}/npt_short.gro" -p "${LEG_DIR}/topol.top" -o "${window_dir}/prod.tpr" -maxwarn 10
            fi
            gmx mdrun -deffnm "${window_dir}/prod" -ntmpi 1 -ntomp 2 -nb gpu -bonded gpu -pme gpu -v ${MDRUN_NSTEPS_ARG} || \
                gmx mdrun -deffnm "${window_dir}/prod" -ntmpi 1 -ntomp 10 -v ${MDRUN_NSTEPS_ARG}
        fi

        rm -f "${LEG_DIR}/.run_${lambda_idx}"
        trap - EXIT
        echo "  ${lambda_name}: ✓ Complete"
    }

    if [ 4 -gt 1 ]; then
        run_window "${lambda_idx}" "${lambda_name}" "${lambda_mdp}" "${window_dir}" "${gpu_id}" &
        JOB_COUNT=$((JOB_COUNT + 1))
    else
        run_window "${lambda_idx}" "${lambda_name}" "${lambda_mdp}" "${window_dir}" "${gpu_id}"
    fi
done

if [ 4 -gt 1 ]; then
    echo "Waiting for all windows to complete..."
    wait
fi
