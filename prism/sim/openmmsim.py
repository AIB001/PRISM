#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenMM simulator for PRISM.

Runs simulations from a GROMACS-prepared system (topol.top + solv_ions.gro) with
a broad set of run types — not just plain MD:

  * minimize                energy minimization
  * nvt / npt               equilibration (optionally position-restrained)
  * npt-membrane            NPT with a MonteCarloMembraneBarostat (bilayers)
  * production              NVT or NPT (isotropic / membrane) production
  * anneal                  simulated-annealing temperature ramp

Settings (nonbonded cutoff, constraints, dt, temperature, pressure, ensemble)
are read from the PRISM-generated MDP files rather than hardcoded. Extras:
position restraints, modern LangevinMiddle integrator, optional hydrogen-mass
repartitioning (HMR) for larger time steps, checkpoint/restart resume, and
DCD + state-data + checkpoint reporters.

Optional plugins (installed separately — see setup.py extras `openmm-plugins`):
  * openmm-plumed  -> metadynamics / steered MD via a PLUMED script
  * openmm-ml      -> hybrid ML/MM potentials (NNP) for the ligand region
Both are gated: absent plugins simply disable the corresponding feature.

Nothing here is machine-specific: platform, device, and thread counts are
auto-detected or passed in, never hardcoded.
"""

import os
import subprocess
import tempfile
import re
import time
from datetime import timedelta

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    OPENMM_AVAILABLE = True
except ImportError:
    OPENMM_AVAILABLE = False
    print("Warning: OpenMM not available. Optional install:  conda install -c conda-forge openmm")

# --- optional plugins (gated; absence only disables the feature) ---
try:
    import openmmplumed  # noqa: F401  (provides openmm.app PlumedForce via the plugin)
    from openmmplumed import PlumedForce

    PLUMED_AVAILABLE = True
except Exception:
    PLUMED_AVAILABLE = False

try:
    from openmmml import MLPotential  # noqa: F401

    OPENMMML_AVAILABLE = True
except Exception:
    OPENMMML_AVAILABLE = False

try:
    from tqdm import tqdm

    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False


class ProgressReporter:
    """Custom reporter for displaying simulation progress with tqdm (or text)."""

    def __init__(self, total_steps, update_interval=1000, desc="Simulation"):
        self.total_steps = total_steps
        self.update_interval = update_interval
        self.desc = desc
        self.start_time = None
        self.pbar = None

    def __enter__(self):
        if TQDM_AVAILABLE:
            self.pbar = tqdm(total=self.total_steps, desc=self.desc, unit="steps", unit_scale=True)
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.pbar:
            self.pbar.close()

    def update(self, current_step):
        if self.pbar:
            increment = current_step - self.pbar.n
            if increment > 0:
                self.pbar.update(increment)
        elif self.total_steps and current_step % max(1, self.total_steps // 20) == 0:
            elapsed = time.time() - self.start_time
            progress = current_step / self.total_steps * 100
            eta = elapsed / (current_step / self.total_steps) - elapsed if current_step > 0 else 0
            print(f"{self.desc}: {progress:.1f}% ({current_step}/{self.total_steps}) ETA {timedelta(seconds=int(eta))}")


class TqdmReporter:
    """OpenMM reporter that drives a ProgressReporter."""

    def __init__(self, progress_reporter, reportInterval):
        self.progress_reporter = progress_reporter
        self.reportInterval = reportInterval

    def describeNextReport(self, simulation):
        steps = self.reportInterval - simulation.currentStep % self.reportInterval
        return (steps, False, False, False, False, None)

    def report(self, simulation, state):
        self.progress_reporter.update(simulation.currentStep)


# Default per-stage step counts when an MDP value is missing.
_DEFAULT_NSTEPS = {"nvt": 250000, "npt": 250000, "prod": 5000000, "anneal": 500000}


class OpenMMSimulator:
    """Run flexible OpenMM simulations from a GROMACS-prepared PRISM system."""

    def __init__(self, gmx_dir, system_files, ff_dir):
        if not OPENMM_AVAILABLE:
            raise ImportError(
                "OpenMM is not installed. Optional install: conda install -c conda-forge openmm"
            )
        self.gmx_dir = gmx_dir
        self.system_files = system_files
        self.ff_dir = ff_dir
        self.gmx_data_dir = self._find_gromacs_data_dir()
        self.mdp_params = self._parse_mdp_files()
        self.cuda_info = self._check_cuda_availability()

        print("OpenMM simulator initialized")
        print(f"  Platforms: {self._get_platform_info()}")
        if self.cuda_info["available"]:
            print(f"  CUDA: available (driver {self.cuda_info['driver_version']}, "
                  f"{self.cuda_info['device_count']} device(s))")
        else:
            print(f"  CUDA: not available ({self.cuda_info['reason']})")
        plugins = [n for n, ok in (("openmm-plumed", PLUMED_AVAILABLE), ("openmm-ml", OPENMMML_AVAILABLE)) if ok]
        print(f"  Optional plugins: {', '.join(plugins) if plugins else 'none (metadynamics/ML disabled)'}")

    # ================================================================== #
    #  Topology / environment helpers (load a GROMACS system into OpenMM) #
    # ================================================================== #
    def _check_cuda_availability(self):
        info = {"available": False, "driver_version": None, "device_count": 0, "devices": [], "reason": "not checked"}
        try:
            for i in range(mm.Platform.getNumPlatforms()):
                if mm.Platform.getPlatform(i).getName() == "CUDA":
                    info["available"] = True
                    info["reason"] = "CUDA platform found"
                    try:
                        r = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                                           capture_output=True, text=True, timeout=10)
                        if r.returncode == 0:
                            info["devices"] = [d for d in r.stdout.strip().split("\n") if d]
                            info["device_count"] = len(info["devices"])
                        r = subprocess.run(["nvidia-smi", "--query-gpu=driver_version", "--format=csv,noheader"],
                                           capture_output=True, text=True, timeout=10)
                        if r.returncode == 0:
                            info["driver_version"] = r.stdout.strip().split("\n")[0]
                    except Exception:
                        pass
                    break
            if not info["available"]:
                info["reason"] = "CUDA platform not present in this OpenMM build"
        except Exception as e:
            info["reason"] = str(e)
        return info

    def _find_gromacs_data_dir(self):
        try:
            r = subprocess.run(["which", "gmx"], capture_output=True, text=True)
            if r.returncode == 0:
                gmx_bin = r.stdout.strip()
                bin_dir = os.path.dirname(gmx_bin)
                for d in (
                    os.path.join(os.path.dirname(bin_dir), "share", "gromacs", "top"),
                    os.path.join(os.path.dirname(bin_dir), "share", "top"),
                    "/usr/share/gromacs/top", "/usr/local/share/gromacs/top",
                ):
                    if os.path.isdir(d) and any(f.endswith(".ff") for f in os.listdir(d)):
                        return d
        except Exception:
            pass
        return None

    def _create_modified_topology(self, original_top):
        with open(original_top) as f:
            lines = f.readlines()
        top_dir = os.path.dirname(original_top)
        modified = []
        for line in lines:
            if line.strip().startswith("#include"):
                if line.count('"') == 1:
                    line = line.rstrip() + '"\n'
                m = re.search(r'#include\s+"([^"]+)"', line)
                if m:
                    new = self._resolve_include_path(m.group(1), top_dir)
                    if new != m.group(1):
                        line = f'#include "{new}"\n'
            modified.append(line)
        if modified != lines:
            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".top", delete=False)
            tmp.write("".join(modified))
            tmp.close()
            return tmp.name
        return original_top

    def _resolve_include_path(self, include_file, top_dir):
        if os.path.isabs(include_file):
            return include_file
        if ".ff/" in include_file and self.gmx_data_dir:
            p = os.path.join(self.gmx_data_dir, include_file)
            if os.path.exists(p):
                return p
        candidates = [
            os.path.join(top_dir, include_file),
            os.path.join(self.ff_dir, include_file) if self.ff_dir else None,
            os.path.join(os.path.dirname(top_dir), include_file),
            os.path.abspath(os.path.join(top_dir, include_file)),
        ]
        for c in candidates:
            if c and os.path.exists(c):
                return os.path.abspath(c)
        return include_file

    def _get_platform_info(self):
        return ", ".join(mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms()))

    def _parse_mdp_files(self):
        params = {}
        for stage in ["em", "nvt", "npt", "md"]:
            f = self.system_files.get(f"{stage}_mdp")
            if f and os.path.exists(f):
                params[stage] = self._parse_mdp(f)
        return params

    @staticmethod
    def _parse_mdp(mdp_file):
        params = {}
        with open(mdp_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith(";") and "=" in line:
                    k, v = line.split("=", 1)
                    params[k.strip()] = v.split(";")[0].strip()
        return params

    def _mdp(self, stage, key, default=None):
        return self.mdp_params.get(stage, {}).get(key, default)

    # ================================================================== #
    #  System / integrator / restraint / barostat builders               #
    # ================================================================== #
    def _nonbonded_from_mdp(self):
        """Derive nonbonded cutoff (nm) and constraints from the MDP files."""
        # cutoff: prefer rvdw/rcoulomb from any stage (md > npt > nvt > em)
        cutoff = 1.0
        for stage in ("md", "npt", "nvt", "em"):
            for key in ("rvdw", "rcoulomb", "rlist"):
                v = self._mdp(stage, key)
                if v:
                    try:
                        cutoff = float(v)
                        break
                    except ValueError:
                        pass
            else:
                continue
            break
        cons_str = (self._mdp("md", "constraints") or self._mdp("npt", "constraints")
                    or self._mdp("nvt", "constraints") or "h-bonds").lower()
        cons = {"h-bonds": app.HBonds, "hbonds": app.HBonds, "all-bonds": app.AllBonds,
                "h-angles": app.HAngles, "none": None}.get(cons_str, app.HBonds)
        return cutoff * unit.nanometer, cons

    def _build_system(self, top, box_vectors, hmr=False, hydrogen_mass_amu=1.5, ml_region=None):
        cutoff, constraints = self._nonbonded_from_mdp()
        kwargs = dict(nonbondedMethod=app.PME, nonbondedCutoff=cutoff, constraints=constraints, rigidWater=True)
        if hmr:
            kwargs["hydrogenMass"] = hydrogen_mass_amu * unit.amu
        system = top.createSystem(**kwargs)
        if ml_region is not None:
            system = self._apply_ml_potential(system, top.topology, ml_region)
        return system

    @staticmethod
    def _make_integrator(temperature, dt, kind="langevin-middle", friction=1.0):
        fr = friction / unit.picosecond
        if kind == "langevin":
            return mm.LangevinIntegrator(temperature, fr, dt)
        if kind == "nose-hoover":
            return mm.NoseHooverIntegrator(temperature, fr, dt)
        # default: LangevinMiddle (more accurate sampling than legacy Langevin)
        return mm.LangevinMiddleIntegrator(temperature, fr, dt)

    def _add_position_restraints(self, system, topology, ref_positions, selection="backbone", fc=1000.0):
        """Harmonic position restraints (kJ/mol/nm^2) on a selection of atoms."""
        if not selection or selection == "none":
            return None
        force = mm.CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
        force.addGlobalParameter("k", fc * unit.kilojoule_per_mole / unit.nanometer**2)
        for p in ("x0", "y0", "z0"):
            force.addPerParticleParameter(p)
        bb = {"CA", "C", "N", "O"}
        std_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
                  "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                  "HID", "HIE", "HIP", "CYX", "SEP", "TPO", "PTR"}
        n = 0
        for atom in topology.atoms():
            keep = False
            rn = atom.residue.name
            if selection == "heavy":
                keep = atom.element is not None and atom.element.symbol != "H"
            elif selection == "protein":
                keep = rn in std_aa
            else:  # backbone (default)
                keep = rn in std_aa and atom.name in bb
            if keep:
                pos = ref_positions[atom.index].value_in_unit(unit.nanometer)
                force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
                n += 1
        if n == 0:
            return None
        system.addForce(force)
        print(f"  Position restraints: {n} atoms ({selection}, k={fc} kJ/mol/nm^2)")
        return force

    @staticmethod
    def _add_barostat(system, ensemble, temperature, pressure, surface_tension=0.0, freq=25):
        if ensemble == "npt":
            system.addForce(mm.MonteCarloBarostat(pressure, temperature, freq))
        elif ensemble in ("npt-membrane", "membrane"):
            system.addForce(mm.MonteCarloMembraneBarostat(
                pressure, surface_tension * unit.bar * unit.nanometer, temperature,
                mm.MonteCarloMembraneBarostat.XYIsotropic, mm.MonteCarloMembraneBarostat.ZFree, freq))

    def _apply_plumed(self, system, plumed_input):
        if not PLUMED_AVAILABLE:
            print("  PLUMED requested but openmm-plumed is not installed — skipping (see setup.py extras).")
            return
        script = open(plumed_input).read() if os.path.exists(plumed_input) else plumed_input
        system.addForce(PlumedForce(script))
        print("  PLUMED bias added (metadynamics / steered MD).")

    def _apply_ml_potential(self, system, topology, ml_region):
        if not OPENMMML_AVAILABLE:
            print("  ML potential requested but openmm-ml is not installed — using pure MM (see setup.py extras).")
            return system
        # ml_region: name of the ML model, e.g. 'ani2x' or 'mace-off'; applied to the ligand residue
        lig_atoms = [a.index for a in topology.atoms() if a.residue.name in ("LIG", "MOL", "UNK")]
        if not lig_atoms:
            print("  ML potential: no ligand atoms found; using pure MM.")
            return system
        potential = MLPotential(ml_region if isinstance(ml_region, str) else "ani2x")
        return potential.createMixedSystem(topology, system, lig_atoms)

    def _setup_platform(self, platform_name, device_index, precision="mixed"):
        try:
            if platform_name == "auto":
                platform_name = "CUDA" if self.cuda_info["available"] else (
                    "OpenCL" if "OpenCL" in self._get_platform_info() else "CPU")
            platform = mm.Platform.getPlatformByName(platform_name)
            props = {}
            if platform_name in ("CUDA", "OpenCL"):
                props["DeviceIndex"] = str(device_index)
                if platform_name == "CUDA":
                    props["Precision"] = precision
            print(f"  Platform: {platform_name}" + (f" (device {device_index}, {precision})" if props else ""))
            return platform, props
        except Exception as e:
            print(f"  Platform '{platform_name}' unavailable ({e}); falling back to CPU.")
            return mm.Platform.getPlatformByName("CPU"), {}

    # ================================================================== #
    #  Public run API                                                     #
    # ================================================================== #
    def run(self, stages=None, platform="auto", device_index=0, precision="mixed",
            ensemble="npt", restraints="backbone", restraint_fc=1000.0, equil_restraints=True,
            integrator="langevin-middle", hmr=False, hydrogen_mass_amu=1.5,
            plumed=None, ml_potential=None, resume=True, num_threads=None,
            gmx_data_dir=None, **kwargs):
        """Run an OpenMM workflow.

        Parameters
        ----------
        stages : list[str]   subset of ['em','nvt','npt','prod','anneal'] (default em,nvt,npt,prod)
        ensemble : str       production ensemble: 'nvt' | 'npt' | 'npt-membrane'
        restraints : str     'backbone' | 'heavy' | 'protein' | 'none' (applied during nvt/npt if equil_restraints)
        hmr : bool           hydrogen-mass repartitioning (allows a larger time step)
        integrator : str     'langevin-middle' (default) | 'langevin' | 'nose-hoover'
        plumed : str         path to a PLUMED input (metadynamics / steered MD; needs openmm-plumed)
        ml_potential : str   ML model name for a hybrid ML/MM ligand region (needs openmm-ml)
        resume : bool        resume production from a checkpoint if present
        """
        if gmx_data_dir:
            self.gmx_data_dir = gmx_data_dir
        if stages is None:
            stages = ["em", "nvt", "npt", "prod"]
        if num_threads:
            os.environ.setdefault("OPENMM_CPU_THREADS", str(num_threads))
        if kwargs:
            print(f"Ignoring extra parameters: {list(kwargs)}")

        platform_obj = self._setup_platform(platform, device_index, precision)

        # Load GROMACS system
        gro = app.GromacsGroFile(self.system_files["gro"])
        modified_top = self._create_modified_topology(self.system_files["top"])
        try:
            top = app.GromacsTopFile(modified_top, periodicBoxVectors=gro.getPeriodicBoxVectors())
        except Exception as e:
            if modified_top != self.system_files["top"] and os.path.exists(modified_top):
                os.unlink(modified_top)
            raise RuntimeError(f"Could not load GROMACS topology in OpenMM ({e}). "
                               f"Consider the GROMACS engine: sim.run(engine='gmx').")

        try:
            outputs = {}
            positions = gro.positions
            velocities = None
            box = gro.getPeriodicBoxVectors()

            for stage in stages:
                print(f"\n{'='*60}\n {stage.upper()} stage\n{'='*60}")
                if stage == "em":
                    system = self._build_system(top, box, hmr, hydrogen_mass_amu, ml_potential)
                    if plumed:
                        self._apply_plumed(system, plumed)
                    positions, box = self._run_minimization(system, top.topology, positions, box, platform_obj)
                    outputs["em"] = {"gro": os.path.join(self.gmx_dir, "em_openmm.gro")}
                elif stage in ("nvt", "npt"):
                    system = self._build_system(top, box, hmr, hydrogen_mass_amu, ml_potential)
                    if plumed:
                        self._apply_plumed(system, plumed)
                    sel = restraints if equil_restraints else "none"
                    positions, velocities, box = self._run_md(
                        system, top.topology, positions, velocities, box, platform_obj,
                        stage=stage, ensemble=("npt" if stage == "npt" else "nvt"),
                        integrator=integrator, restraints=sel, restraint_fc=restraint_fc)
                    outputs[stage] = {"gro": os.path.join(self.gmx_dir, f"{stage}_openmm.gro")}
                elif stage in ("prod", "production"):
                    system = self._build_system(top, box, hmr, hydrogen_mass_amu, ml_potential)
                    if plumed:
                        self._apply_plumed(system, plumed)
                    positions, velocities, box = self._run_md(
                        system, top.topology, positions, velocities, box, platform_obj,
                        stage="prod", ensemble=ensemble, integrator=integrator,
                        restraints="none", production=True, resume=resume)
                    outputs["prod"] = {"dcd": os.path.join(self.gmx_dir, "prod_openmm.dcd"),
                                       "log": os.path.join(self.gmx_dir, "prod_openmm.log"),
                                       "gro": os.path.join(self.gmx_dir, "prod_openmm.gro")}
                elif stage == "anneal":
                    system = self._build_system(top, box, hmr, hydrogen_mass_amu, ml_potential)
                    positions, velocities, box = self._run_anneal(
                        system, top.topology, positions, velocities, box, platform_obj, integrator=integrator)
                    outputs["anneal"] = {"gro": os.path.join(self.gmx_dir, "anneal_openmm.gro")}
                else:
                    print(f"  Unknown stage '{stage}', skipping.")
            return outputs
        finally:
            if modified_top != self.system_files["top"] and os.path.exists(modified_top):
                os.unlink(modified_top)

    # ================================================================== #
    #  Stage implementations                                              #
    # ================================================================== #
    def _run_minimization(self, system, topology, positions, box, platform_info):
        platform, props = platform_info
        em = self.mdp_params.get("em", {})
        max_it = int(float(em.get("nsteps", 10000)))
        tol = float(em.get("emtol", 200.0)) * unit.kilojoule_per_mole / unit.nanometer
        integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
        sim = app.Simulation(topology, system, integrator, platform, props)
        sim.context.setPositions(positions)
        sim.context.setPeriodicBoxVectors(*box)
        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
        print(f"  Initial PE: {e0}")
        sim.minimizeEnergy(tolerance=tol, maxIterations=max_it)
        state = sim.context.getState(getPositions=True, getEnergy=True, enforcePeriodicBox=True)
        print(f"  Final PE:   {state.getPotentialEnergy()}")
        out = os.path.join(self.gmx_dir, "em_openmm.gro")
        self._write_gro_file(out, topology, state.getPositions(), state.getPeriodicBoxVectors(), "EM (OpenMM)")
        return state.getPositions(), state.getPeriodicBoxVectors()

    def _run_md(self, system, topology, positions, velocities, box, platform_info,
                stage="nvt", ensemble="nvt", integrator="langevin-middle",
                restraints="none", restraint_fc=1000.0, production=False, resume=False):
        platform, props = platform_info
        mdp_stage = "md" if production else stage
        p = self.mdp_params.get(mdp_stage, {})
        dt = float(p.get("dt", 0.002)) * unit.picoseconds
        nsteps = int(float(p.get("nsteps", _DEFAULT_NSTEPS.get("prod" if production else stage, 250000))))
        temperature = float(str(p.get("ref_t", "310")).split()[0]) * unit.kelvin
        pressure = float(str(p.get("ref_p", "1.0")).split()[0]) * unit.bar

        # restraints must be added before the context is created
        if restraints and restraints != "none":
            self._add_position_restraints(system, topology, positions, restraints, restraint_fc)
        # detect membrane from MDP (semiisotropic) unless explicitly set
        if ensemble == "npt" and "semiisotropic" in str(p.get("pcoupltype", "")).lower():
            ensemble = "npt-membrane"
        if ensemble in ("npt", "npt-membrane", "membrane"):
            self._add_barostat(system, ensemble, temperature, pressure)

        integ = self._make_integrator(temperature, dt, integrator)
        sim = app.Simulation(topology, system, integ, platform, props)

        deffnm = os.path.join(self.gmx_dir, f"{'prod' if production else stage}_openmm")
        chk = deffnm + ".chk"
        start = 0
        if production and resume and os.path.exists(chk):
            try:
                sim.loadCheckpoint(chk)
                start = sim.currentStep
                print(f"  Resuming from checkpoint at step {start}")
            except Exception as e:
                print(f"  Checkpoint load failed ({e}); starting fresh.")
        if start == 0:
            sim.context.setPositions(positions)
            sim.context.setPeriodicBoxVectors(*box)
            if velocities is not None:
                sim.context.setVelocities(velocities)
            else:
                sim.context.setVelocitiesToTemperature(temperature)

        total_time = (nsteps * dt).in_units_of(unit.nanoseconds)
        print(f"  {stage.upper()} {ensemble}: {nsteps:,} steps ({total_time}), T={temperature}"
              + (f", P={pressure}" if ensemble != "nvt" else ""))

        report_every = int(float(p.get("nstxout-compressed", p.get("nstxout", 5000)) or 5000))
        if report_every <= 0:   # a literal "0" in the mdp would give a zero-interval reporter
            report_every = 5000
        sim.reporters.append(app.StateDataReporter(
            deffnm + ".log", max(1000, report_every), step=True, time=True, temperature=True,
            potentialEnergy=True, density=(ensemble != "nvt"), speed=True, totalSteps=nsteps, separator="\t"))
        if production:
            sim.reporters.append(app.DCDReporter(deffnm + ".dcd", report_every, append=(start > 0)))
            sim.reporters.append(app.CheckpointReporter(chk, max(report_every, 50000)))

        with ProgressReporter(nsteps, desc=f"{stage.upper()} {ensemble}") as prog:
            sim.reporters.append(TqdmReporter(prog, 1000))
            remaining = nsteps - start
            if remaining > 0:
                sim.step(remaining)

        state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        self._write_gro_file(deffnm + ".gro", topology, state.getPositions(),
                             state.getPeriodicBoxVectors(), f"{stage} (OpenMM)")
        return state.getPositions(), state.getVelocities(), state.getPeriodicBoxVectors()

    def _run_anneal(self, system, topology, positions, velocities, box, platform_info,
                    integrator="langevin-middle", t_low=100.0, t_high=None, n_cycles=1):
        """Simulated annealing: ramp temperature low->high->low over the MD nsteps."""
        platform, props = platform_info
        p = self.mdp_params.get("md", self.mdp_params.get("npt", {}))
        dt = float(p.get("dt", 0.002)) * unit.picoseconds
        nsteps = int(float(p.get("nsteps", _DEFAULT_NSTEPS["anneal"])))
        t_high = t_high or float(str(p.get("ref_t", "310")).split()[0])
        integ = self._make_integrator(t_low * unit.kelvin, dt, integrator)
        sim = app.Simulation(topology, system, integ, platform, props)
        sim.context.setPositions(positions)
        sim.context.setPeriodicBoxVectors(*box)
        sim.context.setVelocitiesToTemperature(t_low * unit.kelvin)
        n_windows = 20
        steps_per = max(1, nsteps // (n_windows * 2 * n_cycles))
        print(f"  Simulated annealing {t_low}->{t_high}->{t_low} K, {nsteps:,} steps")
        with ProgressReporter(nsteps, desc="Anneal") as prog:
            sim.reporters.append(TqdmReporter(prog, 1000))
            for _ in range(n_cycles):
                for ramp in (list(range(n_windows)), list(range(n_windows, 0, -1))):
                    for w in ramp:
                        T = t_low + (t_high - t_low) * (w / n_windows)
                        integ.setTemperature(T * unit.kelvin)
                        sim.step(steps_per)
        state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
        out = os.path.join(self.gmx_dir, "anneal_openmm.gro")
        self._write_gro_file(out, topology, state.getPositions(), state.getPeriodicBoxVectors(), "Anneal (OpenMM)")
        return state.getPositions(), state.getVelocities(), state.getPeriodicBoxVectors()

    # ================================================================== #
    def _write_gro_file(self, filename, topology, positions, box_vectors, title="OpenMM"):
        try:
            with open(filename, "w") as f:
                f.write(f"{title}\n{topology.getNumAtoms():5d}\n")
                for i, atom in enumerate(topology.atoms()):
                    pos = positions[i].value_in_unit(unit.nanometer)
                    f.write(f"{atom.residue.index + 1:5d}{atom.residue.name[:5]:<5s}{atom.name[:5]:>5s}"
                            f"{i + 1:5d}{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n")
                if box_vectors is not None:
                    box = box_vectors.value_in_unit(unit.nanometer)
                    f.write(f"{box[0][0]:10.5f}{box[1][1]:10.5f}{box[2][2]:10.5f}\n")
        except Exception as e:
            print(f"Warning: failed to write GRO ({e}); writing PDB instead.")
            with open(filename.replace(".gro", ".pdb"), "w") as f:
                app.PDBFile.writeFile(topology, positions, f)
