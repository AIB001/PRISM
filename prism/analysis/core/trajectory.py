import os
import logging
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
    """Manage trajectory file loading and caching for MDTraj only"""

    # Supported trajectory formats (MDTraj only)
    MDTRAJ_FORMATS = {
        '.dcd': 'DCD (NAMD/CHARMM)',
        '.xtc': 'XTC (GROMACS)',
        '.trr': 'TRR (GROMACS)',
        '.pdb': 'PDB',
        '.nc': 'NetCDF (Amber)',
        '.crd': 'CRD (CHARMM)',
        '.mdcrd': 'MDCRD (Amber)'
    }

    def __init__(self):
        self._cached_mdtraj_trajectories = {}

    def detect_format(self, trajectory_file: str) -> str:
        """Detect trajectory file format from extension"""
        ext = Path(trajectory_file).suffix.lower()
        if ext in self.MDTRAJ_FORMATS:
            return self.MDTRAJ_FORMATS[ext]
        else:
            raise ValueError(f"Unsupported trajectory format: {ext}")

    def is_supported_format(self, trajectory_file: str) -> bool:
        """Check if trajectory format is supported"""
        ext = Path(trajectory_file).suffix.lower()
        return ext in self.MDTRAJ_FORMATS

    def get_supported_formats(self) -> Dict[str, List[str]]:
        """Get list of supported formats"""
        return {
            'mdtraj': list(self.MDTRAJ_FORMATS.values())
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

            trajectory = md.load(trajectory_file, top=topology_file)
            self._cached_mdtraj_trajectories[cache_key] = trajectory
            logger.info(f"Loaded MDTraj trajectory: {trajectory.n_frames} frames, {trajectory.n_atoms} atoms")
            return trajectory
        except Exception as e:
            logger.error(f"Failed to load trajectory: {e}")
            raise