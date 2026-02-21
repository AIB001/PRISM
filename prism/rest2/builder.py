"""
Build the REST2 output directory structure.

Creates scaled topologies, MDP files, and run scripts for REST2
replica exchange enhanced sampling with GROMACS.
"""

import math
import os
import re
import shutil
from pathlib import Path
from typing import List, Tuple

from .partial_tempering import partial_tempering


# ---------------------------------------------------------------------------
# Temperature ladder
# ---------------------------------------------------------------------------

def calculate_temperature_ladder(
    t_ref: float,
    t_max: float,
    n_replicas: int,
) -> Tuple[List[float], List[float]]:
    """Calculate exponential temperature ladder for REST2.

    T[i] = T_ref * (T_max / T_ref) ^ (i / (n_replicas - 1))
    lambda[i] = T_ref / T[i]

    Args:
        t_ref: Reference (thermostat) temperature in K
        t_max: Maximum effective temperature in K
        n_replicas: Number of replicas

    Returns:
        (temperatures, lambdas) - both as lists of floats
    """
    temps = []
    lambdas = []

    for i in range(n_replicas):
        if n_replicas > 1:
            t_i = t_ref * (t_max / t_ref) ** (i / (n_replicas - 1))
        else:
            t_i = t_ref
        lam_i = t_ref / t_i
        temps.append(round(t_i, 2))
        lambdas.append(round(lam_i, 6))

    return temps, lambdas


# ---------------------------------------------------------------------------
# MDP file generators
# ---------------------------------------------------------------------------

def write_em_mdp(path: str) -> None:
    """Write energy minimization MDP file."""
    content = """\
; Energy minimization parameters
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0

; Van der Waals
rvdw        = 1.0
DispCorr    = EnerPres

; Constraints
constraints = h-bonds
"""
    with open(path, 'w') as f:
        f.write(content)


def write_nvt_mdp(path: str, t_ref: float, nsteps: int = 50000) -> None:
    """Write NVT equilibration MDP file."""
    content = f"""\
; NVT equilibration
integrator  = md
dt          = 0.002
nsteps      = {nsteps}
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0

; Van der Waals
rvdw        = 1.0
DispCorr    = EnerPres

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = {t_ref}

; Pressure coupling - off for NVT
pcoupl      = no

; Constraints
constraints = h-bonds

; Velocity generation
gen_vel     = yes
gen_temp    = {t_ref}
gen_seed    = -1

; Position restraints
define      = -DPOSRES
"""
    with open(path, 'w') as f:
        f.write(content)


def write_npt_mdp(path: str, t_ref: float, nsteps: int = 50000) -> None:
    """Write NPT equilibration MDP file."""
    content = f"""\
; NPT equilibration
integrator  = md
dt          = 0.002
nsteps      = {nsteps}
nstxout-compressed = 5000
nstenergy   = 5000
nstlog      = 5000

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0

; Van der Waals
rvdw        = 1.0
DispCorr    = EnerPres

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = {t_ref}

; Pressure coupling (C-rescale for equilibration stability)
pcoupl      = C-rescale
pcoupltype  = isotropic
tau_p       = 5.0
ref_p       = 1.0
compressibility = 4.5e-5
refcoord_scaling = com

; Constraints
constraints = h-bonds

; Velocity generation
gen_vel     = no

; Position restraints
define      = -DPOSRES
"""
    with open(path, 'w') as f:
        f.write(content)


def write_prod_mdp(path: str, t_ref: float, nsteps: int = 5000000) -> None:
    """Write production MD MDP file for REST2."""
    content = f"""\
; REST2 production MD
integrator  = md
dt          = 0.002
nsteps      = {nsteps}
nstxout-compressed = 5000
nstenergy   = 1000
nstlog      = 1000

; Neighbor searching
nstlist     = 10
cutoff-scheme = Verlet

; Electrostatics
coulombtype = PME
rcoulomb    = 1.0

; Van der Waals
rvdw        = 1.0
DispCorr    = EnerPres

; Temperature coupling
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = {t_ref}

; Pressure coupling
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5

; Constraints
constraints = h-bonds

; Velocity generation
gen_vel     = no
"""
    with open(path, 'w') as f:
        f.write(content)


# ---------------------------------------------------------------------------
# Script generator
# ---------------------------------------------------------------------------

def write_rest2_run_script(
    path: str, n_replicas: int, t_ref: float, replex: int = 1000,
) -> None:
    """Write a single all-in-one REST2 run script.

    The script is placed directly inside GMX_PROLIG_REST2/ so BASEDIR is the
    directory containing the script itself (no '..' needed).

    Phases:
      1. Energy minimization  (shared, uses topol_original.top)
      2. NVT equilibration    (shared, uses topol_original.top)
      3. NPT equilibration    (shared, uses topol_original.top)
      4. Generate per-replica TPR files (each replica uses its own scaled topol.top)
      5. Production REST2 replica-exchange MD
    """
    content = f"""\
#!/bin/bash
# =============================================================================
#  REST2 (Replica Exchange with Solute Tempering v2) - Complete Run Script
# =============================================================================
#  Generated by PRISM
#
#  Directory layout expected:
#    GMX_PROLIG_REST2/
#    ├── rest2_run.sh          <-- this script
#    ├── solv_ions.gro
#    ├── topol_original.top    (with #include, used for equilibration)
#    ├── mdps/
#    │   ├── em.mdp, nvt.mdp, npt.mdp, prod.mdp
#    ├── replica_0/topol.top   (scaled topologies for production)
#    └── ...
#
#  Usage:
#    cd GMX_PROLIG_REST2
#    bash rest2_run.sh
# =============================================================================

set -e

GMX=${{GMX:-gmx}}
GMX_MPI=${{GMX_MPI:-gmx_mpi}}

# BASEDIR = directory where this script lives (GMX_PROLIG_REST2/)
BASEDIR=$(cd "$(dirname "$0")" && pwd)
cd "$BASEDIR"

echo "============================================================"
echo "  REST2 Enhanced Sampling - Complete Workflow"
echo "============================================================"
echo "  Base directory : $BASEDIR"
echo "  Replicas       : {n_replicas}"
echo "  T_ref          : {t_ref} K"
echo "  Exchange every  : {replex} steps"
echo ""

# -----------------------------------------------------------------------
#  Phase 1 : Energy Minimization
# -----------------------------------------------------------------------
echo "=== Phase 1/5 : Energy Minimization ==="
if [ -f em.gro ]; then
    echo "  em.gro already exists, skipping EM."
else
    $GMX grompp -f mdps/em.mdp -c solv_ions.gro -p topol_original.top -o em.tpr -maxwarn 2
    $GMX mdrun -deffnm em -v
fi
echo "  EM done."
echo ""

# -----------------------------------------------------------------------
#  Phase 2 : NVT Equilibration (position restraints)
# -----------------------------------------------------------------------
echo "=== Phase 2/5 : NVT Equilibration ==="
if [ -f nvt.gro ]; then
    echo "  nvt.gro already exists, skipping NVT."
else
    $GMX grompp -f mdps/nvt.mdp -c em.gro -r em.gro -p topol_original.top -o nvt.tpr -maxwarn 2
    $GMX mdrun -deffnm nvt -v
fi
echo "  NVT done."
echo ""

# -----------------------------------------------------------------------
#  Phase 3 : NPT Equilibration (position restraints)
# -----------------------------------------------------------------------
echo "=== Phase 3/5 : NPT Equilibration ==="
if [ -f npt.gro ]; then
    echo "  npt.gro already exists, skipping NPT."
else
    $GMX grompp -f mdps/npt.mdp -c nvt.gro -r nvt.gro -p topol_original.top -o npt.tpr -maxwarn 2
    $GMX mdrun -deffnm npt -v
fi
echo "  NPT done."
echo ""

# -----------------------------------------------------------------------
#  Phase 4 : Generate per-replica TPR files
# -----------------------------------------------------------------------
echo "=== Phase 4/5 : Generating per-replica TPR files ==="
for i in $(seq 0 {n_replicas - 1}); do
    echo "  Replica $i ..."
    cd "$BASEDIR/replica_$i"
    $GMX grompp \\
        -f "$BASEDIR/mdps/prod.mdp" \\
        -c "$BASEDIR/npt.gro" \\
        -p topol.top \\
        -o prod.tpr \\
        -maxwarn 2
done
cd "$BASEDIR"
echo "  All TPR files generated."
echo ""

# -----------------------------------------------------------------------
#  Phase 5 : Production REST2 Replica Exchange MD
# -----------------------------------------------------------------------
echo "=== Phase 5/5 : Production REST2 Replica Exchange MD ==="

# Build -multidir argument
MULTIDIR=""
for i in $(seq 0 {n_replicas - 1}); do
    MULTIDIR="$MULTIDIR $BASEDIR/replica_$i"
done

mpirun -np {n_replicas} $GMX_MPI mdrun \\
    -deffnm prod \\
    -multidir $MULTIDIR \\
    -hrex \\
    -replex {replex} \\
    -v

echo ""
echo "============================================================"
echo "  REST2 Production Complete!"
echo "============================================================"
"""
    with open(path, 'w') as f:
        f.write(content)
    os.chmod(path, 0o755)


# ---------------------------------------------------------------------------
# Main build function
# ---------------------------------------------------------------------------

def build_rest2(
    md_dir: str,
    output_dir: str,
    marked_top_path: str,
    original_top_path: str,
    gro_path: str,
    n_replicas: int,
    t_ref: float,
    t_max: float,
    processed_top_path: str,
) -> None:
    """Build the complete REST2 output directory.

    Args:
        md_dir: Source GMX_PROLIG_MD directory
        output_dir: Output GMX_PROLIG_REST2 directory
        marked_top_path: Path to REST2-marked topology
        original_top_path: Path to original topology with #includes
        gro_path: Path to solv_ions.gro
        n_replicas: Number of replicas
        t_ref: Reference temperature (K)
        t_max: Maximum effective temperature (K)
        processed_top_path: Path to processed topology (for reference)
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- Temperature ladder ---
    temps, lambdas = calculate_temperature_ladder(t_ref, t_max, n_replicas)

    print(f"\n  Temperature ladder ({n_replicas} replicas):")
    for i, (t, l) in enumerate(zip(temps, lambdas)):
        print(f"    Replica {i}: T_eff = {t:.2f} K, lambda = {l:.6f}")

    # Write temperature ladder file
    ladder_path = os.path.join(output_dir, 'temperature_ladder.txt')
    with open(ladder_path, 'w') as f:
        f.write(f"# REST2 Temperature Ladder\n")
        f.write(f"# T_ref = {t_ref} K, T_max = {t_max} K, N_replicas = {n_replicas}\n")
        f.write(f"# Thermostat temperature for ALL replicas = {t_ref} K\n")
        f.write(f"# replica  T_effective(K)  lambda\n")
        for i, (t, l) in enumerate(zip(temps, lambdas)):
            f.write(f"  {i:>7d}  {t:>14.2f}  {l:>8.6f}\n")

    # --- Copy GRO file ---
    shutil.copy2(gro_path, os.path.join(output_dir, 'solv_ions.gro'))

    # --- Copy original topology with all dependencies for equilibration ---
    _copy_original_topology(md_dir, original_top_path, output_dir)

    # --- Copy processed and marked topologies for reference ---
    dest_processed = os.path.join(output_dir, 'processed.top')
    dest_marked = os.path.join(output_dir, 'rest2_marked.top')
    if os.path.abspath(processed_top_path) != os.path.abspath(dest_processed):
        shutil.copy2(processed_top_path, dest_processed)
    if os.path.abspath(marked_top_path) != os.path.abspath(dest_marked):
        shutil.copy2(marked_top_path, dest_marked)

    # --- Generate scaled topologies ---
    print(f"\n  Generating scaled topologies...")
    with open(marked_top_path, 'r') as f:
        marked_lines = f.readlines()

    for i in range(n_replicas):
        replica_dir = os.path.join(output_dir, f'replica_{i}')
        os.makedirs(replica_dir, exist_ok=True)
        scaled_lines = partial_tempering(marked_lines, lambdas[i])
        with open(os.path.join(replica_dir, 'topol.top'), 'w') as f:
            f.writelines(scaled_lines)
        print(f"    Replica {i}: lambda = {lambdas[i]:.6f}")

    # --- MDP files ---
    mdp_dir = os.path.join(output_dir, 'mdps')
    os.makedirs(mdp_dir, exist_ok=True)
    write_em_mdp(os.path.join(mdp_dir, 'em.mdp'))
    write_nvt_mdp(os.path.join(mdp_dir, 'nvt.mdp'), t_ref)
    write_npt_mdp(os.path.join(mdp_dir, 'npt.mdp'), t_ref)
    write_prod_mdp(os.path.join(mdp_dir, 'prod.mdp'), t_ref)
    print(f"  MDP files written to: {mdp_dir}")

    # --- Run script ---
    script_path = os.path.join(output_dir, 'rest2_run.sh')
    write_rest2_run_script(script_path, n_replicas, t_ref)
    print(f"  Run script written to: {script_path}")


def _is_forcefield_include(inc_path: str) -> bool:
    """Check if an #include path refers to a GROMACS force field (not a local file)."""
    # Force field includes contain ".ff/" (e.g. "amber14sb.ff/forcefield.itp",
    # "charmm36.ff/tip3p.itp", "oplsaa.ff/ions.itp")
    return '.ff/' in inc_path


def _copy_original_topology(md_dir: str, original_top_path: str, output_dir: str) -> None:
    """Copy the original topology and all its #include dependencies.

    Recursively resolves nested #include directives so that files like
    posre_LIG.itp (included from within LIG.itp) are also copied.
    """
    md_dir = os.path.abspath(md_dir)
    output_top = os.path.join(output_dir, 'topol_original.top')
    _INCLUDE_RE = re.compile(r'#include\s+"([^"]+)"')

    # Track files already copied to avoid duplicates
    copied_basenames: set = set()

    def _copy_and_rewrite(src_path: str, dst_path: str) -> None:
        """Copy an .itp/.top file, rewriting local #include paths to basenames.

        Also recursively copies any local files referenced by nested #includes.
        """
        src_dir = os.path.dirname(os.path.abspath(src_path))

        with open(src_path, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            stripped = line.strip()
            match = _INCLUDE_RE.search(stripped)
            if match and stripped.startswith('#'):
                inc_path = match.group(1)
                if _is_forcefield_include(inc_path):
                    # Keep force field includes as-is
                    new_lines.append(line)
                    continue

                # Resolve the included file relative to the source file's directory
                abs_inc = os.path.normpath(os.path.join(src_dir, inc_path))
                if os.path.exists(abs_inc):
                    basename = os.path.basename(abs_inc)
                    dest_inc = os.path.join(output_dir, basename)

                    # Copy the included file (and its own nested includes) if not done yet
                    if basename not in copied_basenames:
                        copied_basenames.add(basename)
                        if abs_inc.endswith('.itp') or abs_inc.endswith('.top'):
                            _copy_and_rewrite(abs_inc, dest_inc)
                        else:
                            shutil.copy2(abs_inc, dest_inc)

                    # Rewrite the include to use just the basename
                    # Preserve any surrounding preprocessor context (#ifdef etc.)
                    new_lines.append(line[:match.start(1)] + basename + line[match.end(1):])
                    continue

            new_lines.append(line)

        with open(dst_path, 'w') as f:
            f.writelines(new_lines)

    _copy_and_rewrite(original_top_path, output_top)

    # Also copy any posre files from the MD directory that weren't already copied
    # (they might not be directly #included from topol.top)
    for fname in os.listdir(md_dir):
        if fname.startswith('posre_') and fname.endswith('.itp'):
            dst = os.path.join(output_dir, fname)
            if not os.path.exists(dst):
                shutil.copy2(os.path.join(md_dir, fname), dst)
