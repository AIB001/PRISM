#!/usr/bin/env bash
# PRISM-FEP Master Script
# GPU PME: enabled
set -euo pipefail

# Test mode: limit production run steps (DO NOT use for production FEP calculations)
if [ -n "${PRISM_MDRUN_NSTEPS:-}" ]; then
    echo "⚠️  WARNING: PRISM_MDRUN_NSTEPS is set to ${PRISM_MDRUN_NSTEPS}"
    echo "⚠️  This is for TESTING ONLY - results will NOT be valid for publication!"
    MDRUN_NSTEPS_ARG="-nsteps ${PRISM_MDRUN_NSTEPS}"
else
    MDRUN_NSTEPS_ARG=""
fi

EXECUTION_MODE="${PRISM_FEP_MODE:-standard}"
if [[ "${EXECUTION_MODE}" != "standard" && "${EXECUTION_MODE}" != "repex" ]]; then
    echo "Error: unsupported execution mode '${EXECUTION_MODE}'"
    echo "Supported modes: standard, repex"
    exit 1
fi

# Function to run a single leg
run_leg() {
    local leg_name=$1
    local leg_dir="${FEP_DIR}/${leg_name}"

    if [ ! -d "${leg_dir}" ]; then
        echo "Error: Leg directory not found: ${leg_dir}"
        return 1
    fi

    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                    PRISM-FEP ${leg_name} LEG                                   ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
    echo ""

    cd "${leg_dir}"

    local MDP_DIR="${leg_dir}/mdps"
    local INPUT_DIR="${leg_dir}/input"
    local BUILD_DIR="${leg_dir}/build"

    # Energy minimization
    mkdir -p ${BUILD_DIR}

    if [ -f ${BUILD_DIR}/em.gro ]; then
        echo "✓ EM already completed, skipping..."
    elif [ -f ${BUILD_DIR}/em.tpr ]; then
        echo "EM checkpoint found, resuming..."
        gmx mdrun -deffnm ${BUILD_DIR}/em -cpi ${BUILD_DIR}/em.cpt -v
    else
        echo "Running energy minimization..."
        gmx grompp -f ${MDP_DIR}/em.mdp -c ${INPUT_DIR}/conf.gro -p topol.top -o ${BUILD_DIR}/em.tpr -maxwarn 10
        gmx mdrun -deffnm ${BUILD_DIR}/em -ntmpi 1 -ntomp 10 -v
    fi

    # Check EM convergence
    if [ -f ${BUILD_DIR}/em.log ]; then
        max_force=$(grep -E "Maximum force|Fmax=" ${BUILD_DIR}/em.log | tail -1 | grep -oE "[0-9]+\\.?[0-9]*e[+-][0-9]+" | tail -1 || true)
        if [ -n "$max_force" ]; then
            echo "EM final max force: $max_force kJ/mol/nm"
            # Warning if force is still high (> 1000 kJ/mol/nm)
            # Use awk to handle scientific notation comparison
            if [ "$(echo "$max_force" | awk '{print ($1 > 1000)}')" -eq 1 ]; then
                echo "⚠️  WARNING: EM converged with high force. System may need more minimization."
            fi
        fi
    fi

    # NVT equilibration (must run to completion)
    if [ -f ${BUILD_DIR}/nvt.gro ]; then
        echo "✓ NVT already completed, skipping..."
    elif [ -f ${BUILD_DIR}/nvt.tpr ]; then
        echo "NVT checkpoint found, resuming..."
        gmx mdrun -deffnm ${BUILD_DIR}/nvt -cpi ${BUILD_DIR}/nvt.cpt -v
    else
        echo "Running NVT equilibration..."
        gmx grompp -f ${MDP_DIR}/nvt.mdp -c ${BUILD_DIR}/em.gro -r ${BUILD_DIR}/em.gro -p topol.top -o ${BUILD_DIR}/nvt.tpr -maxwarn 10
        gmx mdrun -deffnm ${BUILD_DIR}/nvt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -v || gmx mdrun -deffnm ${BUILD_DIR}/nvt -ntmpi 1 -ntomp 10 -v
    fi

    # NPT equilibration (must run to completion)
    if [ -f ${BUILD_DIR}/npt.gro ]; then
        echo "✓ NPT already completed, skipping..."
    elif [ -f ${BUILD_DIR}/npt.tpr ]; then
        # Check if NPT checkpoint is valid (has corresponding cpt file)
        if [ -f ${BUILD_DIR}/npt.cpt ]; then
            echo "NPT checkpoint found, resuming..."
            gmx mdrun -deffnm ${BUILD_DIR}/npt -cpi ${BUILD_DIR}/npt.cpt -v
        else
            # NPT TPR exists but no checkpoint - previous run failed, clean up and restart
            echo "⚠️  Previous NPT run failed (no checkpoint), cleaning up and restarting..."
            rm -f ${BUILD_DIR}/npt.* 2>/dev/null || true
            echo "Running NPT equilibration..."
            gmx grompp -f ${MDP_DIR}/npt.mdp -c ${BUILD_DIR}/nvt.gro -r ${BUILD_DIR}/nvt.gro -t ${BUILD_DIR}/nvt.cpt -p topol.top -o ${BUILD_DIR}/npt.tpr -maxwarn 10
            gmx mdrun -deffnm ${BUILD_DIR}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -v || gmx mdrun -deffnm ${BUILD_DIR}/npt -ntmpi 1 -ntomp 10 -v
        fi
    else
        echo "Running NPT equilibration..."
        gmx grompp -f ${MDP_DIR}/npt.mdp -c ${BUILD_DIR}/nvt.gro -r ${BUILD_DIR}/nvt.gro -t ${BUILD_DIR}/nvt.cpt -p topol.top -o ${BUILD_DIR}/npt.tpr -maxwarn 10
        gmx mdrun -deffnm ${BUILD_DIR}/npt -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -v || gmx mdrun -deffnm ${BUILD_DIR}/npt -ntmpi 1 -ntomp 10 -v
    fi

    echo "Running FEP production for all lambda windows..."
    echo "Execution mode: ${EXECUTION_MODE}"
    if [ "${EXECUTION_MODE}" = "repex" ]; then
        ./run_prod_repex.sh
    else
        ./run_prod_standard.sh
    fi

    # Clean up intermediate step PDB files
    echo "Cleaning up intermediate PDB files..."
    rm -f step*.pdb 2>/dev/null || true

    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                    ${leg_name} LEG COMPLETE                                   ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"

    cd "${FEP_DIR}"
}

# Main script logic
FEP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Number of replica runs (configured at build time)
REPLICAS=1

# Parse target specification (e.g., bound1, bound1-3, unbound1-2,unbound3, etc.)
parse_targets() {
    local spec=$1
    local targets=()

    # Check if it's a simple leg name (bound, unbound, bound, unbound with replicas)
    if [[ "$spec" =~ ^(bound|unbound)$ ]]; then
        if [ ${REPLICAS} -gt 1 ]; then
            # Expand to all replicas: bound -> bound1 bound2 bound3
            for i in $(seq 1 ${REPLICAS}); do
                targets+=("${spec}${i}")
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
        local leg="${BASH_REMATCH[1]}"
        local start="${BASH_REMATCH[2]}"
        local end="${BASH_REMATCH[3]}"
        for i in $(seq $start $end); do
            targets+=("${leg}${i}")
        done
    else
        echo "Error: Invalid target specification: $spec"
        return 1
    fi

    # Output targets as space-separated list
    echo "${targets[@]}"
}

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
    echo "  Replicas configured: ${REPLICAS}"
    if [ ${REPLICAS} -gt 1 ]; then
        echo "  Available: bound1..bound${REPLICAS}, unbound1..unbound${REPLICAS}"
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
    if [ ${REPLICAS} -gt 1 ]; then
        for i in $(seq 1 ${REPLICAS}); do
            all_targets+=("bound${i}")
            all_targets+=("unbound${i}")
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
for target in "${all_targets[@]}"; do
    if [ -d "${FEP_DIR}/${target}" ]; then
        run_leg "${target}"
    else
        echo "⚠️  Warning: Directory ${target} not found, skipping..."
    fi
done

echo ""
echo "✓ All requested FEP calculations completed."
