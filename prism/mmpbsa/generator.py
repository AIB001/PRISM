#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MMPBSA Generator - Generates mmpbsa.in and mmpbsa_run.sh for MM/PBSA calculations.

Supports two backends:
- gmx_MMPBSA (default): Uses gmx_MMPBSA with GROMACS topology directly.
- AMBER MMPBSA.py (--gmx2amber): Uses parmed to convert GROMACS -> AMBER topology,
  then calls AMBER's native MMPBSA.py. Requires AmberTools with parmed.

Each backend supports two sub-modes:
- Single-frame (default): EM -> NVT -> NPT -> MM/PBSA on equilibrated structure
- Trajectory-based: EM -> NVT -> NPT -> Production MD -> MM/PBSA on trajectory
"""

import os
import shutil


class MMPBSAGenerator:
    """Generates MM/PBSA input and run script files.

    Parameters
    ----------
    mmpbsa_dir : str
        Output directory for generated files.
    protein_ff_name : str
        GROMACS protein force field name (e.g., 'amber14sb').
    ligand_forcefield : str
        Ligand force field type (e.g., 'gaff', 'gaff2').
    temperature : float
        Simulation temperature in K.
    traj_ns : float or None
        Production MD length in ns for trajectory mode. None = single-frame mode.
    trajectory_interval_ps : float
        Trajectory save interval in ps (default: 500).
    gmx2amber : bool
        If True, use parmed + AMBER MMPBSA.py instead of gmx_MMPBSA (default: False).
    """

    # GROMACS protein FF -> AMBER leaprc mapping
    PROTEIN_FF_MAP = {
        "amber14sb": "leaprc.protein.ff14SB",
        "amber99sb": "oldff/leaprc.ff99SB",
        "amber99sb-ildn": "oldff/leaprc.ff99SBildn",
        "charmm27": "leaprc.protein.ff14SB",  # fallback for CHARMM
    }

    # Ligand FF -> AMBER leaprc mapping
    LIGAND_FF_MAP = {
        "gaff": "leaprc.gaff",
        "gaff2": "leaprc.gaff2",
    }

    def __init__(
        self,
        mmpbsa_dir,
        protein_ff_name,
        ligand_forcefield,
        temperature,
        traj_ns=None,
        trajectory_interval_ps=500,
        gmx2amber=False,
    ):
        self.mmpbsa_dir = mmpbsa_dir
        self.protein_ff_name = protein_ff_name.lower()
        self.ligand_forcefield = ligand_forcefield.lower()
        self.temperature = temperature
        self.traj_ns = traj_ns
        self.trajectory_interval_ps = trajectory_interval_ps
        self.gmx2amber = gmx2amber

    @property
    def is_trajectory(self):
        return self.traj_ns is not None

    def _resolve_protein_leaprc(self):
        """Resolve GROMACS protein FF name to AMBER leaprc string."""
        leaprc = self.PROTEIN_FF_MAP.get(self.protein_ff_name)
        if not leaprc:
            for key, val in self.PROTEIN_FF_MAP.items():
                if key in self.protein_ff_name:
                    leaprc = val
                    break
        return leaprc or "leaprc.protein.ff14SB"

    def _resolve_ligand_leaprc(self):
        """Resolve ligand FF type to AMBER leaprc string.

        Returns None for force fields without an AMBER leaprc equivalent
        (e.g., OpenFF, CGenFF), signaling that gmx_MMPBSA should use
        parmed-based topology conversion instead of tleap.
        """
        return self.LIGAND_FF_MAP.get(self.ligand_forcefield)

    # ------------------------------------------------------------------
    # mmpbsa.in
    # ------------------------------------------------------------------

    def generate_input(self):
        """Generate mmpbsa.in input file.

        When gmx2amber=False, generates gmx_MMPBSA format (with forcefields, solvated_trajectory).
        When gmx2amber=True, generates AMBER MMPBSA.py format (no gmx_MMPBSA-specific fields).

        Returns
        -------
        str
            Path to the generated mmpbsa.in file.
        """
        if self.is_trajectory:
            interval = max(1, int(1000 / self.trajectory_interval_ps))
            startframe, endframe = 1, 9999999
        else:
            startframe, endframe, interval = 1, 1, 1

        if self.gmx2amber:
            # strip_mask: MMPBSA.py uses this to strip solvent/ions from the
            # solvated trajectory (-sp) to match the unsolvated complex (-cp).
            # Must include GROMACS residue names (SOL, NA, CL) since parmed
            # preserves them when converting to prmtop format.
            content = f"""&general
  startframe = {startframe}
  endframe = {endframe}
  interval = {interval}
  strip_mask = :WAT,SOL,HOH,Cl-,CL,Na+,NA,K+,K,MG,CA,ZN,SOD,CLA
  keep_files = 2
  verbose = 2
/
&pb
  istrng = 0.15
  inp = 2
  radiopt = 1
/
"""
        else:
            protein_leaprc = self._resolve_protein_leaprc()
            ligand_leaprc = self._resolve_ligand_leaprc()

            # For GAFF/GAFF2: specify forcefields so gmx_MMPBSA uses tleap
            # For non-GAFF (OpenFF, CGenFF, etc.): omit forcefields so
            # gmx_MMPBSA falls back to parmed-based topology conversion
            if ligand_leaprc:
                forcefields_line = f'  forcefields = "{protein_leaprc},{ligand_leaprc}"'
            else:
                forcefields_line = f'  forcefields = "{protein_leaprc}"'

            content = f"""&general
  sys_name = "PRISM_MMPBSA"
  startframe = {startframe}
  endframe = {endframe}
  interval = {interval}
{forcefields_line}
  temperature = {self.temperature}
  PBRadii = 4
  solvated_trajectory = 1
  keep_files = 2
  verbose = 2
/
&pb
  istrng = 0.15
  inp = 2
  radiopt = 0
/
"""
        input_path = os.path.join(self.mmpbsa_dir, "mmpbsa.in")
        with open(input_path, "w") as f:
            f.write(content)
        return input_path

    # ------------------------------------------------------------------
    # mmpbsa_run.sh
    # ------------------------------------------------------------------

    def generate_script(self):
        """Generate mmpbsa_run.sh bash script.

        Returns
        -------
        str
            Path to the generated mmpbsa_run.sh file.
        """
        sim_header = self._sim_header_block()

        if self.gmx2amber:
            backend_name = "AMBER MMPBSA.py"
            mmpbsa_tail = self._amber_mmpbsa_tail_block()
            # Copy utils.py to output dir for prmtop periodicity fix
            utils_src = os.path.join(os.path.dirname(__file__), "utils.py")
            utils_dst = os.path.join(self.mmpbsa_dir, "utils.py")
            shutil.copy2(utils_src, utils_dst)
        else:
            backend_name = "gmx_MMPBSA"
            mmpbsa_tail = self._mmpbsa_tail_block()

        if self.is_trajectory:
            workflow_desc = f"EM -> NVT -> NPT -> Production MD ({self.traj_ns} ns) -> {backend_name}"
            middle = self._production_md_block()
            cs_path, ct_path = "./prod/md.tpr", "./prod/md.xtc"
        else:
            workflow_desc = f"EM -> NVT -> NPT -> {backend_name} (single-frame)"
            middle = self._trjconv_block()
            cs_path, ct_path = "./npt/npt.tpr", "single_frame.xtc"

        script_content = (
            sim_header.format(workflow_desc=workflow_desc)
            + middle
            + mmpbsa_tail.format(cs_path=cs_path, ct_path=ct_path)
        )

        script_path = os.path.join(self.mmpbsa_dir, "mmpbsa_run.sh")
        with open(script_path, "w") as f:
            f.write(script_content)
        os.chmod(script_path, 0o755)
        return script_path

    # ------------------------------------------------------------------
    # Shell script building blocks
    # ------------------------------------------------------------------

    @staticmethod
    def _sim_header_block():
        """EM + NVT + NPT equilibration block (shared by both modes)."""
        return """#!/bin/bash

######################################################
# PRISM MM/PBSA Run Script
# Runs: {workflow_desc}
######################################################

GMX=${{GMX:-gmx}}

######################################################
# SIMULATION
######################################################

# Energy Minimization (EM)
mkdir -p em
if [ -f ./em/em.gro ]; then
    echo "EM already completed, skipping..."
elif [ -f ./em/em.tpr ]; then
    echo "EM tpr file found, continuing from checkpoint..."
    $GMX mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v -cpi ./em/em.cpt
else
    echo "Starting EM from scratch..."
    $GMX grompp -f ../mdps/em.mdp -c solv_ions.gro -r solv_ions.gro -p topol.top -o ./em/em.tpr -maxwarn 999
    $GMX mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi

# NVT Equilibration
mkdir -p nvt
if [ -f ./nvt/nvt.gro ]; then
    echo "NVT already completed, skipping..."
elif [ -f ./nvt/nvt.tpr ]; then
    echo "NVT tpr file found, continuing from checkpoint..."
    $GMX mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v -cpi ./nvt/nvt.cpt
else
    echo "Starting NVT from scratch..."
    $GMX grompp -f ../mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -maxwarn 999
    $GMX mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt -v
fi

# NPT Equilibration
mkdir -p npt
if [ -f ./npt/npt.gro ]; then
    echo "NPT already completed, skipping..."
elif [ -f ./npt/npt.tpr ]; then
    echo "NPT tpr file found, continuing from checkpoint..."
    $GMX mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v -cpi ./npt/npt.cpt
else
    echo "Starting NPT from scratch..."
    $GMX grompp -f ../mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr -maxwarn 999
    $GMX mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt -v
fi
"""

    @staticmethod
    def _production_md_block():
        """Production MD block (trajectory mode only)."""
        return """
# Production MD
mkdir -p prod
if [ -f ./prod/md.gro ]; then
    echo "Production MD already completed, skipping..."
elif [ -f ./prod/md.tpr ]; then
    echo "Production MD tpr file found, continuing from checkpoint..."
    $GMX mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v -cpi ./prod/md.cpt
else
    echo "Starting Production MD from scratch..."
    $GMX grompp -f ../mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro -p topol.top -o ./prod/md.tpr -maxwarn 999
    $GMX mdrun -ntmpi 1 -ntomp 10 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./prod/md.tpr -deffnm ./prod/md -v
fi
"""

    @staticmethod
    def _trjconv_block():
        """Convert NPT final structure to single-frame XTC (single-frame mode only)."""
        return """
# Convert NPT final structure to single-frame XTC for MM/PBSA
echo "Converting NPT structure to single-frame trajectory..."
echo "System" | $GMX trjconv -f ./npt/npt.gro -s ./npt/npt.tpr -o single_frame.xtc
"""

    @staticmethod
    def _mmpbsa_tail_block():
        """Index generation + group detection + gmx_MMPBSA invocation block.

        Note: -lm (ligand mol2) is intentionally omitted. When -cp topol.top
        is provided, gmx_MMPBSA extracts all parameters (including ligand)
        from the GROMACS topology via parmed conversion. This works with any
        ligand force field (GAFF, GAFF2, OpenFF, CGenFF, etc.) and avoids
        parmchk2 failures from atom type mismatches in the original mol2.
        """
        return """
######################################################
# MM/PBSA CALCULATION
######################################################

echo ""
echo "======================================================"
echo "Starting MM/PBSA calculation..."
echo "======================================================"

# Generate index file from solvated system
echo "Generating index groups..."
echo "q" | $GMX make_ndx -f solv_ions.gro -o index.ndx

# Find Protein and LIG group numbers from index.ndx
# Groups are numbered starting from 0, in order of appearance
PROTEIN_GROUP=""
LIG_GROUP=""
GROUP_NUM=-1

while IFS= read -r line; do
    if [[ "$line" =~ ^\\[.*\\]$ ]]; then
        GROUP_NUM=$((GROUP_NUM + 1))
        GROUP_NAME=$(echo "$line" | sed 's/\\[//;s/\\]//' | xargs)
        if [ "$GROUP_NAME" = "Protein" ]; then
            PROTEIN_GROUP=$GROUP_NUM
        elif [ "$GROUP_NAME" = "LIG" ]; then
            LIG_GROUP=$GROUP_NUM
        fi
    fi
done < index.ndx

if [ -z "$PROTEIN_GROUP" ] || [ -z "$LIG_GROUP" ]; then
    echo "ERROR: Could not find Protein and/or LIG groups in index.ndx"
    echo "Found: Protein=$PROTEIN_GROUP, LIG=$LIG_GROUP"
    exit 1
fi

echo "Using groups: Protein=$PROTEIN_GROUP, LIG=$LIG_GROUP"

gmx_MMPBSA -O \\
  -i mmpbsa.in \\
  -cs {cs_path} \\
  -ci index.ndx \\
  -cg $PROTEIN_GROUP $LIG_GROUP \\
  -ct {ct_path} \\
  -cp topol.top \\
  -o FINAL_RESULTS_MMPBSA.dat \\
  -eo FINAL_RESULTS_MMPBSA.csv

echo ""
echo "======================================================"
echo "MM/PBSA calculation complete!"
echo "Results: FINAL_RESULTS_MMPBSA.dat"
echo "======================================================"
"""

    @staticmethod
    def _amber_mmpbsa_tail_block():
        """parmed GROMACS->AMBER conversion + MMPBSA.py invocation block."""
        return """
######################################################
# GROMACS -> AMBER TOPOLOGY CONVERSION
######################################################

echo ""
echo "======================================================"
echo "Converting GROMACS topology to AMBER format (parmed)..."
echo "======================================================"

python3 << 'PYEOF'
import parmed as pmd
import sys
try:
    system = pmd.load_file('topol.top', xyz='solv_ions.gro')
    system.save('solvated.prmtop', overwrite=True)
    print("  -> solvated.prmtop")

    # Detect solvent/ion residue names to strip
    known_solvent_ions = {{'SOL','HOH','WAT','NA','CL','K','MG','ZN','CA',
                          'Na+','Cl-','K+','Mg2+','Zn2+','Ca2+','SOD','CLA'}}
    strip_resnames = sorted(r for r in set(res.name for res in system.residues)
                            if r in known_solvent_ions)
    strip_mask = ':' + ','.join(strip_resnames)

    complex_sys = system['!' + strip_mask]
    complex_sys.save('complex.prmtop', overwrite=True)
    print("  -> complex.prmtop")

    receptor = complex_sys['!:LIG']
    receptor.save('receptor.prmtop', overwrite=True)
    print("  -> receptor.prmtop")

    ligand = complex_sys[':LIG']
    ligand.save('ligand.prmtop', overwrite=True)
    print("  -> ligand.prmtop")

    print("Topology conversion complete.")
except Exception as e:
    print(f"ERROR: {{e}}", file=sys.stderr)
    sys.exit(1)
PYEOF

######################################################
# FIX PRMTOP PERIODICITY VALUES
######################################################

echo "Checking prmtop files for invalid periodicity values..."
python3 utils.py solvated.prmtop complex.prmtop receptor.prmtop ligand.prmtop

######################################################
# MM/PBSA CALCULATION (AMBER MMPBSA.py)
######################################################

echo ""
echo "======================================================"
echo "Starting MM/PBSA calculation (AMBER MMPBSA.py)..."
echo "======================================================"

MMPBSA.py -O \\
  -i mmpbsa.in \\
  -sp solvated.prmtop \\
  -cp complex.prmtop \\
  -rp receptor.prmtop \\
  -lp ligand.prmtop \\
  -y {ct_path} \\
  -o FINAL_RESULTS_MMPBSA.dat

echo ""
echo "======================================================"
echo "MM/PBSA calculation complete!"
echo "Results: FINAL_RESULTS_MMPBSA.dat"
echo "======================================================"
"""
