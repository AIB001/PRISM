#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM MDP Generator - Generate MDP files for MD simulations
"""

import os


class MDPGenerator:
    """Generate MDP files for MD simulations"""

    def __init__(self, config, output_dir):
        self.config = config
        self.mdp_dir = os.path.join(output_dir, "mdps")
        os.makedirs(self.mdp_dir, exist_ok=True)

    def generate_all(self):
        """Generate all MDP files"""
        print("\n=== Generating MDP Files ===")

        self._generate_ions_mdp()
        self._generate_em_mdp()
        self._generate_nvt_mdp()
        self._generate_npt_mdp()
        self._generate_production_mdp()

        print(f"MDP files generated in {self.mdp_dir}")
        self._print_summary()

    def _generate_ions_mdp(self):
        """Generate ions.mdp file"""
        em_config = self.config['energy_minimization']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']

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
        with open(mdp_path, 'w') as f:
            f.write(content)

    def _generate_em_mdp(self):
        """Generate energy minimization MDP"""
        em_config = self.config['energy_minimization']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']

        em_mdp = os.path.join(self.mdp_dir, "em.mdp")
        with open(em_mdp, 'w') as f:
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
        sim_config = self.config['simulation']
        const_config = self.config['constraints']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']
        tc_config = self.config['temperature_coupling']
        output_config = self.config['output']

        nvt_steps = int(sim_config['equilibration_nvt_time_ps'] / sim_config['dt'])
        energy_interval = int(output_config['energy_interval_ps'] / sim_config['dt'])
        log_interval = int(output_config['log_interval_ps'] / sim_config['dt'])

        tc_grps = ' '.join(tc_config['tc_grps'])
        tau_t = '     '.join(map(str, tc_config['tau_t']))
        ref_t = '     '.join([str(sim_config['temperature'])] * len(tc_config['tc_grps']))

        nvt_mdp = os.path.join(self.mdp_dir, "nvt.mdp")
        with open(nvt_mdp, 'w') as f:
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
        sim_config = self.config['simulation']
        const_config = self.config['constraints']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']
        tc_config = self.config['temperature_coupling']
        pc_config = self.config['pressure_coupling']
        output_config = self.config['output']

        npt_steps = int(sim_config['equilibration_npt_time_ps'] / sim_config['dt'])
        energy_interval = int(output_config['energy_interval_ps'] / sim_config['dt'])
        log_interval = int(output_config['log_interval_ps'] / sim_config['dt'])

        tc_grps = ' '.join(tc_config['tc_grps'])
        tau_t = '     '.join(map(str, tc_config['tau_t']))
        ref_t = '     '.join([str(sim_config['temperature'])] * len(tc_config['tc_grps']))

        npt_mdp = os.path.join(self.mdp_dir, "npt.mdp")
        with open(npt_mdp, 'w') as f:
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
        sim_config = self.config['simulation']
        const_config = self.config['constraints']
        elec_config = self.config['electrostatics']
        vdw_config = self.config['vdw']
        tc_config = self.config['temperature_coupling']
        pc_config = self.config['pressure_coupling']
        output_config = self.config['output']

        prod_steps = int(sim_config['production_time_ns'] * 1000 / sim_config['dt'])
        log_interval = int(output_config['log_interval_ps'] / sim_config['dt'])
        traj_interval = int(output_config['trajectory_interval_ps'] / sim_config['dt'])

        tc_grps = ' '.join(tc_config['tc_grps'])
        tau_t = '     '.join(map(str, tc_config['tau_t']))
        ref_t = '     '.join([str(sim_config['temperature'])] * len(tc_config['tc_grps']))

        output_line = "nstxtcout" if output_config['compressed_trajectory'] else "nstxout"

        md_mdp = os.path.join(self.mdp_dir, "md.mdp")
        with open(md_mdp, 'w') as f:
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

    def _print_summary(self):
        """Print summary of generated MDP files"""
        sim_config = self.config['simulation']
        em_config = self.config['energy_minimization']

        print(f"  - Energy minimization: {em_config['nsteps']} steps")
        print(f"  - NVT equilibration: {sim_config['equilibration_nvt_time_ps']} ps")
        print(f"  - NPT equilibration: {sim_config['equilibration_npt_time_ps']} ps")
        print(f"  - Production: {sim_config['production_time_ns']} ns")