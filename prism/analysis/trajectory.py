import os
import logging
import MDAnalysis as mda
from typing import Optional, List, Dict
from pathlib import Path

try:
    import mdtraj as md
    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
    md = None

logger = logging.getLogger(__name__)

class TrajectoryManager:
    """Manage trajectory file loading and caching for both MDTraj and MDAnalysis"""

    # Supported trajectory formats
    MDTRAJ_FORMATS = {
        '.dcd': 'DCD (NAMD/CHARMM)',
        '.xtc': 'XTC (GROMACS)',
        '.trr': 'TRR (GROMACS)',
        '.h5': 'HDF5 (MDTraj)',
        '.pdb': 'PDB',
        '.nc': 'NetCDF (Amber)',
        '.crd': 'CRD (CHARMM)',
        '.mdcrd': 'MDCRD (Amber)'
    }

    MDANALYSIS_FORMATS = {
        '.dcd': 'DCD (NAMD/CHARMM)',
        '.xtc': 'XTC (GROMACS)',
        '.trr': 'TRR (GROMACS)',
        '.nc': 'NetCDF (Amber)',
        '.crd': 'CRD (CHARMM)',
        '.pdb': 'PDB'
    }

    def __init__(self):
        self._cached_mdtraj_trajectories = {}
        self._cached_universes = {}

    def detect_format(self, trajectory_file: str) -> str:
        """Detect trajectory file format from extension"""
        ext = Path(trajectory_file).suffix.lower()
        if ext in self.MDTRAJ_FORMATS:
            return self.MDTRAJ_FORMATS[ext]
        elif ext in self.MDANALYSIS_FORMATS:
            return self.MDANALYSIS_FORMATS[ext]
        else:
            raise ValueError(f"Unsupported trajectory format: {ext}")

    def is_supported_format(self, trajectory_file: str) -> bool:
        """Check if trajectory format is supported"""
        ext = Path(trajectory_file).suffix.lower()
        return ext in self.MDTRAJ_FORMATS or ext in self.MDANALYSIS_FORMATS

    def get_supported_formats(self) -> Dict[str, List[str]]:
        """Get list of supported formats for each backend"""
        return {
            'mdtraj': list(self.MDTRAJ_FORMATS.values()),
            'mdanalysis': list(self.MDANALYSIS_FORMATS.values())
        }

    def load_mdtraj_trajectory(self, topology_file: str, trajectory_file: str, cache_key: Optional[str] = None):
        """Load MDTraj trajectory for structural analysis"""
        if not MDTRAJ_AVAILABLE:
            raise ImportError("MDTraj is not available. Install with: pip install mdtraj")

        # Check format support
        if not self.is_supported_format(trajectory_file):
            ext = Path(trajectory_file).suffix.lower()
            raise ValueError(f"Unsupported trajectory format: {ext}")

        if cache_key is None:
            cache_key = f"mdtraj_{topology_file}_{trajectory_file}"

        if cache_key in self._cached_mdtraj_trajectories:
            logger.info(f"Using cached MDTraj trajectory: {cache_key}")
            return self._cached_mdtraj_trajectories[cache_key]

        try:
            if not os.path.exists(topology_file):
                raise FileNotFoundError(f"Topology file not found: {topology_file}")
            if not os.path.exists(trajectory_file):
                raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")

            format_name = self.detect_format(trajectory_file)
            logger.info(f"Loading {format_name} trajectory: {trajectory_file}")

            traj = md.load(trajectory_file, top=topology_file)
            if traj.n_frames == 0:
                raise ValueError("Trajectory has no frames")

            self._cached_mdtraj_trajectories[cache_key] = traj
            logger.info(f"Loaded MDTraj trajectory: {traj.n_frames} frames, {traj.n_atoms} atoms")
            return traj
        except Exception as e:
            logger.error(f"Failed to load MDTraj trajectory: {e}")
            raise
    
    def load_universe(self, topology_file: str, trajectory_file: str, cache_key: Optional[str] = None):
        """Load MDAnalysis Universe for RMSD and hydrogen bond analysis"""
        # Check format support
        if not self.is_supported_format(trajectory_file):
            ext = Path(trajectory_file).suffix.lower()
            raise ValueError(f"Unsupported trajectory format: {ext}")

        if cache_key is None:
            cache_key = f"mda_{topology_file}_{trajectory_file}"

        if cache_key in self._cached_universes:
            logger.info(f"Using cached MDAnalysis Universe: {cache_key}")
            return self._cached_universes[cache_key]

        try:
            if not os.path.exists(topology_file):
                raise FileNotFoundError(f"Topology file not found: {topology_file}")
            if not os.path.exists(trajectory_file):
                raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")

            format_name = self.detect_format(trajectory_file)
            logger.info(f"Loading {format_name} universe: {trajectory_file}")

            universe = mda.Universe(topology_file, trajectory_file)
            self._cached_universes[cache_key] = universe
            logger.info(f"Loaded MDAnalysis Universe: {len(universe.trajectory)} frames, {universe.atoms.n_atoms} atoms")
            return universe
        except Exception as e:
            logger.error(f"Failed to load universe: {e}")
            raise
