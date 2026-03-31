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

NUM_GPUS=4
export OMP_NUM_THREADS=2

window_dirs=()
for lambda_mdp in "${MDP_DIR}"/prod_*.mdp; do
    lambda_name=$(basename "${lambda_mdp}" .mdp)
    lambda_idx="${lambda_name#prod_}"
    window_dir="${LEG_DIR}/window_${lambda_idx}"
    window_dirs+=("${window_dir}")
    mkdir -p "${window_dir}"

    if [ ! -f "${window_dir}/npt_short.gro" ]; then
        if [ ! -f "${window_dir}/npt_short.tpr" ]; then
            gmx grompp -f "${MDP_DIR}/npt_short_${lambda_idx}.mdp" -c "${BUILD_DIR}/npt.gro" -r "${BUILD_DIR}/npt.gro" -t "${BUILD_DIR}/nvt.cpt" -p "${LEG_DIR}/topol.top" -o "${window_dir}/npt_short.tpr" -maxwarn 10
        fi
        gmx mdrun -deffnm "${window_dir}/npt_short" -ntmpi 1 -ntomp 2 -nb gpu -bonded gpu -pme gpu -v ${MDRUN_NSTEPS_ARG} || \
            gmx mdrun -deffnm "${window_dir}/npt_short" -ntmpi 1 -ntomp 10 -v ${MDRUN_NSTEPS_ARG}
    fi

    if [ ! -f "${window_dir}/prod.tpr" ]; then
        gmx grompp -f "${lambda_mdp}" -c "${window_dir}/npt_short.gro" -p "${LEG_DIR}/topol.top" -o "${window_dir}/prod.tpr" -maxwarn 10
    fi
done

num_windows="${#window_dirs[@]}"
if [ "${num_windows}" -eq 0 ]; then
    echo "Error: no prod_*.mdp files found in $MDP_DIR"
    exit 1
fi

gpu_ids=$(seq 0 $((NUM_GPUS - 1)) | awk '{printf "%s", $0}')
window_cpi="${window_dirs[0]}/prod.cpt"
cpi_args=()
if [ -f "${window_cpi}" ]; then
    cpi_args=(-cpi prod.cpt)
fi

echo "Production mode: repex"
echo "Replica-exchange windows: ${num_windows}"
echo "Configured GPUs: $NUM_GPUS"
echo "GPU id string: $gpu_ids"

mpirun -oversubscribe -np "${num_windows}" \
    gmx_mpi mdrun -v -deffnm prod -nb gpu -bonded gpu -pme gpu -replex 1000 \
    -multidir "${window_dirs[@]}" -gpu_id "${gpu_ids}" -pin on -dhdl dhdl \
    "${cpi_args[@]}" ${MDRUN_NSTEPS_ARG}
