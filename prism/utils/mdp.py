#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM MDP Generator - Generate MDP files for MD simulations
"""

import os

# Import color utilities
try:
    from .colors import print_success, print_info, print_subheader, path, number
except ImportError:
    try:
        from prism.utils.colors import print_success, print_info, print_subheader, path, number
    except ImportError:
        # Fallback if colors not available
        def print_success(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}✓ {x}")

        def print_info(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}ℹ {x}")

        def print_subheader(x, **kwargs):
            prefix = kwargs.get("prefix", "")
            print(f"{prefix}\n=== {x} ===")

        def path(x):
            return x

        def number(x):
            return x


class MDPGenerator:
    """Generate MDP files for MD simulations"""

    def __init__(self, config, output_dir):
        self.config = config
        self.mdp_dir = os.path.join(output_dir, "mdps")
        os.makedirs(self.mdp_dir, exist_ok=True)

    def generate_all(self):
        """Generate all MDP files"""
        print_subheader("Generating MDP Files")

        self._generate_ions_mdp()
        self._generate_em_mdp()
        self._generate_nvt_mdp()
        self._generate_npt_mdp()
        self._generate_production_mdp()

        print_success(f"MDP files generated in {path(self.mdp_dir)}")
        self._print_summary()

    def _generate_ions_mdp(self):
        """Generate ions.mdp file"""
        em_config = self.config["energy_minimization"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]

        content = f"""; ions.mdp - Parameters for adding ions
integrator  = {em_config['integrator']}
emtol       = {em_config['emtol']}
emstep      = {em_config['emstep']}
nsteps      = {em_config['nsteps']}

; Output
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rcoulomb        = {elec_config['rcoulomb']}
rvdw            = {vdw_config['rvdw']}

; Electrostatics
coulombtype     = {elec_config['coulombtype']}
pme_order       = {elec_config['pme_order']}
fourierspacing  = {elec_config['fourierspacing']}

; Temperature and pressure coupling are off
tcoupl      = no
pcoupl      = no

; Periodic boundary conditions
pbc         = xyz

; Constraints
constraints             = {self.config['constraints']['type']}
constraint-algorithm    = {self.config['constraints']['algorithm']}
"""

        mdp_path = os.path.join(self.mdp_dir, "ions.mdp")
        with open(mdp_path, "w") as f:
            f.write(content)

    def _generate_em_mdp(self):
        """Generate energy minimization MDP"""
        em_config = self.config["energy_minimization"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]

        em_mdp = os.path.join(self.mdp_dir, "em.mdp")
        with open(em_mdp, "w") as f:
            f.write(f"""; em.mdp - Energy minimization
title           = Energy minimization
define          = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator      = {em_config['integrator']}
emtol           = {em_config['emtol']}
emstep          = {em_config['emstep']}
nsteps          = {em_config['nsteps']}

nstxout         = 0
nstvout         = 0
nstenergy       = 0
nstlog          = 0
nstxout-compressed = 0

cutoff-scheme   = Verlet
coulombtype     = {elec_config['coulombtype']}
rcoulomb        = {elec_config['rcoulomb']}
rvdw            = {vdw_config['rvdw']}
pbc             = xyz
""")

    def _generate_nvt_mdp(self):
        """Generate NVT equilibration MDP"""
        sim_config = self.config["simulation"]
        const_config = self.config["constraints"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]
        tc_config = self.config["temperature_coupling"]
        output_config = self.config["output"]

        nvt_steps = int(sim_config["equilibration_nvt_time_ps"] / sim_config["dt"])
        energy_interval = int(output_config["energy_interval_ps"] / sim_config["dt"])
        log_interval = int(output_config["log_interval_ps"] / sim_config["dt"])

        tc_grps = " ".join(tc_config["tc_grps"])
        tau_t = "     ".join(map(str, tc_config["tau_t"]))
        ref_t = "     ".join([str(sim_config["temperature"])] * len(tc_config["tc_grps"]))

        nvt_mdp = os.path.join(self.mdp_dir, "nvt.mdp")
        with open(nvt_mdp, "w") as f:
            f.write(f"""; nvt.mdp - NVT equilibration
title               = NVT equilibration with the complex restrained
define              = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator          = md
nsteps              = {nvt_steps}
dt                  = {sim_config['dt']}

nstxout             = 0
nstvout             = 0
nstenergy           = {energy_interval}
nstlog              = {log_interval}

continuation            = no
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = no
pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}

gen_vel                 = yes
gen_temp                = {sim_config['temperature']}
gen_seed                = -1

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {tc_grps}
""")

    def _generate_npt_mdp(self):
        """Generate NPT equilibration MDP"""
        sim_config = self.config["simulation"]
        const_config = self.config["constraints"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]
        tc_config = self.config["temperature_coupling"]
        pc_config = self.config["pressure_coupling"]
        output_config = self.config["output"]

        npt_steps = int(sim_config["equilibration_npt_time_ps"] / sim_config["dt"])
        energy_interval = int(output_config["energy_interval_ps"] / sim_config["dt"])
        log_interval = int(output_config["log_interval_ps"] / sim_config["dt"])

        tc_grps = " ".join(tc_config["tc_grps"])
        tau_t = "     ".join(map(str, tc_config["tau_t"]))
        ref_t = "     ".join([str(sim_config["temperature"])] * len(tc_config["tc_grps"]))

        npt_mdp = os.path.join(self.mdp_dir, "npt.mdp")
        with open(npt_mdp, "w") as f:
            f.write(f"""; npt.mdp - NPT equilibration
title               = NPT equilibration with protein restrained
define              = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator          = md
nsteps              = {npt_steps}
dt                  = {sim_config['dt']}

nstxout             = 0
nstvout             = 0
nstenergy           = {energy_interval}
nstlog              = {log_interval}

continuation            = yes
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = {pc_config['pcoupl']}
pcoupltype              = {pc_config['pcoupltype']}
tau_p                   = {pc_config['tau_p']}
ref_p                   = {sim_config['pressure']}
compressibility         = {pc_config['compressibility']}

pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {tc_grps}
""")

    def _generate_production_mdp(self):
        """Generate production MD MDP"""
        sim_config = self.config["simulation"]
        const_config = self.config["constraints"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]
        tc_config = self.config["temperature_coupling"]
        pc_config = self.config["pressure_coupling"]
        output_config = self.config["output"]

        prod_steps = int(sim_config["production_time_ns"] * 1000 / sim_config["dt"])
        log_interval = int(output_config["log_interval_ps"] / sim_config["dt"])
        traj_interval = int(output_config["trajectory_interval_ps"] / sim_config["dt"])

        tc_grps = " ".join(tc_config["tc_grps"])
        tau_t = "     ".join(map(str, tc_config["tau_t"]))
        ref_t = "     ".join([str(sim_config["temperature"])] * len(tc_config["tc_grps"]))

        output_line = "nstxtcout" if output_config["compressed_trajectory"] else "nstxout"

        md_mdp = os.path.join(self.mdp_dir, "md.mdp")
        with open(md_mdp, "w") as f:
            f.write(f"""; md.mdp - Production MD
title               = Production run ({sim_config['production_time_ns']} ns)
integrator          = md
nsteps              = {prod_steps}
dt                  = {sim_config['dt']}

{output_line}       = {traj_interval}
nstvout             = 0
nstenergy           = 0
nstlog              = {log_interval}

continuation            = yes
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

cutoff-scheme           = Verlet
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = {pc_config['pcoupl']}
pcoupltype              = {pc_config['pcoupltype']}
tau_p                   = {pc_config['tau_p']}
ref_p                   = {sim_config['pressure']}
compressibility         = {pc_config['compressibility']}

pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}
gen_vel                 = no

refcoord_scaling        = com
comm-mode               = Linear
comm-grps               = {tc_grps}
""")

    def generate_smd_mdp(
        self,
        pull_rate=0.01,
        pull_k=1000.0,
        pull_distance=3.6,
        ref_group="Protein",
        pull_group="LIG",
        pbcatom=None,
        freeze_group=None,
    ):
        """
        Generate SMD (Steered Molecular Dynamics) MDP file for pulling simulations.

        Parameters:
        -----------
        pull_rate : float
            Pull rate in nm/ps (default: 0.01 nm/ps)
        pull_k : float
            Pull force constant in kJ/mol/nm^2 (default: 1000.0)
        pull_distance : float
            Total pull distance in nm (default: 5.0 nm)
        ref_group : str
            Name of the reference (anchor) pull group (default: 'Protein')
        pull_group : str
            Name of the pulled group (default: 'LIG')
        pbcatom : int or None
            PBC reference atom index (1-based) for the reference group
        freeze_group : str or None
            Name of the group to freeze during SMD (e.g. 'CA_freeze')
        """
        sim_config = self.config["simulation"]
        const_config = self.config["constraints"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]
        tc_config = self.config["temperature_coupling"]
        pc_config = self.config["pressure_coupling"]
        output_config = self.config["output"]

        # Calculate SMD simulation time based on pull distance and rate
        # Time = distance / rate
        smd_time_ps = pull_distance / pull_rate  # in ps
        smd_steps = int(smd_time_ps / sim_config["dt"])

        # Output settings - more frequent for pulling
        traj_interval = int(10.0 / sim_config["dt"])  # Every 10 ps
        energy_interval = int(1.0 / sim_config["dt"])  # Every 1 ps for pull force
        log_interval = int(10.0 / sim_config["dt"])

        tc_grps = " ".join(tc_config["tc_grps"])
        tau_t = "     ".join(map(str, tc_config["tau_t"]))
        ref_t = "     ".join([str(sim_config["temperature"])] * len(tc_config["tc_grps"]))

        # Build conditional sections
        if freeze_group:
            comm_section = "; COM motion removal - disabled due to frozen atoms\nnstcomm                 = 0\ncomm-mode               = None"
            freeze_section = (
                "\n; ===========================================\n"
                "; Freeze group - keep protein backbone fixed during pulling\n"
                "; ===========================================\n"
                f"freezegrps              = {freeze_group}\n"
                "freezedim               = Y Y Y\n"
            )
        else:
            comm_section = (
                f"; COM motion removal\ncomm-mode               = Linear\ncomm-grps               = {tc_grps}"
            )
            freeze_section = ""

        pbcatom_line = f"pull-group1-pbcatom     = {pbcatom}" if pbcatom else "; pull-group1-pbcatom auto-selected"

        smd_mdp = os.path.join(self.mdp_dir, "smd.mdp")
        with open(smd_mdp, "w") as f:
            f.write(f"""; smd.mdp - Steered Molecular Dynamics (SMD) for PMF calculation
title               = SMD pulling simulation along Z-axis
integrator          = md
nsteps              = {smd_steps}
dt                  = {sim_config['dt']}

; Output control - frequent output for pulling
nstxout-compressed  = {traj_interval}
nstvout             = 0
nstfout             = 0
nstenergy           = {energy_interval}
nstlog              = {log_interval}

; Bond parameters
continuation            = yes
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

; Nonbonded settings
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

; Electrostatics
coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

; Temperature coupling
tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

; Pressure coupling
pcoupl                  = {pc_config['pcoupl']}
pcoupltype              = {pc_config['pcoupltype']}
tau_p                   = {pc_config['tau_p']}
ref_p                   = {sim_config['pressure']}
compressibility         = {pc_config['compressibility']}

; PBC
pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}

; Velocity generation
gen_vel                 = no

{comm_section}

; ===========================================
; Pull code - Steered MD along Z-axis
; ===========================================
pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2

; Pull groups: {ref_group} (reference) and {pull_group} (pulled)
pull-group1-name        = {ref_group}
{pbcatom_line}
pull-group2-name        = {pull_group}

; Pull coordinate 1: Pull ligand away from protein along Z-axis
pull-coord1-type        = umbrella
pull-coord1-geometry    = direction
pull-coord1-dim         = N N Y
pull-coord1-vec         = 0.0 0.0 1.0
pull-coord1-groups      = 1 2
pull-coord1-start       = yes
pull-coord1-rate        = {pull_rate}
pull-coord1-k           = {pull_k}

; Pull output
pull-nstxout            = 500
pull-nstfout            = 500
{freeze_section}""")

        print_success(f"SMD MDP generated: {path(smd_mdp)}")
        print_info(f"  Pull rate: {number(f'{pull_rate} nm/ps')}")
        print_info(f"  Pull constant: {number(f'{pull_k} kJ/mol/nm²')}")
        print_info(f"  Pull distance: {number(f'{pull_distance} nm')}")
        print_info(f"  SMD time: {number(f'{smd_time_ps:.1f} ps')} ({number(smd_steps)} steps)")

        return smd_mdp

    def generate_umbrella_mdp(
        self,
        umbrella_time_ns=10.0,
        pull_k=1000.0,
        ref_group="Protein_near_LIG",
        pull_group="LIG",
        pbcatom=None,
        freeze_group=None,
    ):
        """
        Generate umbrella sampling MDP file for PMF calculations.

        Parameters:
        -----------
        umbrella_time_ns : float
            Simulation time per umbrella window in ns (default: 10.0)
        pull_k : float
            Umbrella restraint force constant in kJ/mol/nm^2 (default: 1000.0)
        ref_group : str
            Name of the reference pull group (default: 'Protein_near_LIG')
        pull_group : str
            Name of the pulled group (default: 'LIG')
        pbcatom : int or None
            PBC reference atom index (1-based) for the reference group
        freeze_group : str or None
            Name of the group to freeze (e.g. 'CA_freeze')
        """
        sim_config = self.config["simulation"]
        const_config = self.config["constraints"]
        elec_config = self.config["electrostatics"]
        vdw_config = self.config["vdw"]
        tc_config = self.config["temperature_coupling"]
        pc_config = self.config["pressure_coupling"]

        # Calculate umbrella simulation steps (ns -> ps -> steps)
        umbrella_steps = int(umbrella_time_ns * 1000 / sim_config["dt"])

        # Output settings - same as SMD
        traj_interval = int(10.0 / sim_config["dt"])  # Every 10 ps
        energy_interval = int(1.0 / sim_config["dt"])  # Every 1 ps for pull data
        log_interval = int(10.0 / sim_config["dt"])

        tc_grps = " ".join(tc_config["tc_grps"])
        tau_t = "     ".join(map(str, tc_config["tau_t"]))
        ref_t = "     ".join([str(sim_config["temperature"])] * len(tc_config["tc_grps"]))

        # Build conditional sections (same logic as SMD)
        if freeze_group:
            comm_section = "; COM motion removal - disabled due to frozen atoms\nnstcomm                 = 0\ncomm-mode               = None"
            freeze_section = (
                "\n; ===========================================\n"
                "; Freeze group - keep protein backbone fixed during umbrella sampling\n"
                "; ===========================================\n"
                f"freezegrps              = {freeze_group}\n"
                "freezedim               = Y Y Y\n"
            )
        else:
            comm_section = (
                f"; COM motion removal\ncomm-mode               = Linear\ncomm-grps               = {tc_grps}"
            )
            freeze_section = ""

        pbcatom_line = f"pull-group1-pbcatom     = {pbcatom}" if pbcatom else "; pull-group1-pbcatom auto-selected"

        umbrella_mdp = os.path.join(self.mdp_dir, "umbrella.mdp")
        with open(umbrella_mdp, "w") as f:
            f.write(f"""; umbrella.mdp - Umbrella sampling window for PMF calculation
title               = Umbrella sampling window ({umbrella_time_ns} ns)
integrator          = md
nsteps              = {umbrella_steps}
dt                  = {sim_config['dt']}

; Output control - frequent output for pull data
nstxout-compressed  = {traj_interval}
nstvout             = 0
nstfout             = 0
nstenergy           = {energy_interval}
nstlog              = {log_interval}

; Bond parameters
continuation            = yes
constraint_algorithm    = {const_config['algorithm']}
constraints             = {const_config['type']}
lincs_iter              = {const_config['lincs_iter']}
lincs_order             = {const_config['lincs_order']}

; Nonbonded settings
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = {elec_config['rcoulomb']}
rvdw                    = {vdw_config['rvdw']}

; Electrostatics
coulombtype             = {elec_config['coulombtype']}
pme_order               = {elec_config['pme_order']}
fourierspacing          = {elec_config['fourierspacing']}

; Temperature coupling
tcoupl                  = {tc_config['tcoupl']}
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

; Pressure coupling
pcoupl                  = {pc_config['pcoupl']}
pcoupltype              = {pc_config['pcoupltype']}
tau_p                   = {pc_config['tau_p']}
ref_p                   = {sim_config['pressure']}
compressibility         = {pc_config['compressibility']}

; PBC
pbc                     = xyz
DispCorr                = {vdw_config['dispcorr']}

; Velocity generation
gen_vel                 = no

{comm_section}

; ===========================================
; Pull code - Umbrella restraint (no pulling)
; ===========================================
pull                    = yes
pull-ncoords            = 1
pull-ngroups            = 2

; Pull groups: {ref_group} (reference) and {pull_group} (restrained)
pull-group1-name        = {ref_group}
{pbcatom_line}
pull-group2-name        = {pull_group}

; Pull coordinate 1: Umbrella restraint at initial distance
pull-coord1-type        = umbrella
pull-coord1-geometry    = direction
pull-coord1-dim         = N N Y
pull-coord1-vec         = 0.0 0.0 1.0
pull-coord1-groups      = 1 2
pull-coord1-start       = yes
pull-coord1-rate        = 0.0
pull-coord1-k           = {pull_k}

; Pull output
pull-nstxout            = 500
pull-nstfout            = 500
{freeze_section}""")

        print_success(f"Umbrella MDP generated: {path(umbrella_mdp)}")
        print_info(f"  Restraint constant: {number(f'{pull_k} kJ/mol/nm²')}")
        print_info(f"  Window time: {number(f'{umbrella_time_ns} ns')} ({number(umbrella_steps)} steps)")

        return umbrella_mdp

    def _print_summary(self):
        """Print summary of generated MDP files"""
        sim_config = self.config["simulation"]
        em_config = self.config["energy_minimization"]

        print_info("Simulation parameters:")
        print(f"  - Energy minimization: {number(em_config['nsteps'])} steps")
        print(f"  - NVT equilibration: {number(sim_config['equilibration_nvt_time_ps'])} ps")
        print(f"  - NPT equilibration: {number(sim_config['equilibration_npt_time_ps'])} ps")
        print(f"  - Production: {number(sim_config['production_time_ns'])} ns")
