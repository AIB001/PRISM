import os
import logging
import MDAnalysis as mda
from typing import Optional

try:
    import mdtraj as md
    MDTRAJ_AVAILABLE = True
except ImportError:
    MDTRAJ_AVAILABLE = False
    md = None

logger = logging.getLogger(__name__)

class TrajectoryManager:
    """Manage trajectory file loading and caching for both MDTraj and MDAnalysis"""
    
    def __init__(self):
        self._cached_mdtraj_trajectories = {}
        self._cached_universes = {}
    
    def load_mdtraj_trajectory(self, topology_file: str, trajectory_file: str, cache_key: Optional[str] = None):
        """Load MDTraj trajectory for new calculation methods"""
        if not MDTRAJ_AVAILABLE:
            raise ImportError("MDTraj is not available. Install with: pip install mdtraj")

        if cache_key is None:
            cache_key = f"mdtraj_{topology_file}_{trajectory_file}"

        if cache_key in self._cached_mdtraj_trajectories:
            return self._cached_mdtraj_trajectories[cache_key]

        try:
            if not os.path.exists(topology_file):
                raise FileNotFoundError(f"Topology file not found: {topology_file}")
            if not os.path.exists(trajectory_file):
                raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")

            traj = md.load(trajectory_file, top=topology_file)
            if traj.n_frames == 0:
                raise ValueError("Trajectory has no frames")

            self._cached_mdtraj_trajectories[cache_key] = traj
            logger.warning(f"Loaded MDTraj trajectory: {traj.n_frames} frames")
            return traj
        except Exception as e:
            logger.error(f"Failed to load MDTraj trajectory: {e}")
            raise
    
    def load_universe(self, topology_file: str, trajectory_file: str, cache_key: Optional[str] = None):
        """Load MDAnalysis Universe for hydrogen bond analysis"""
        if cache_key is None:
            cache_key = f"mda_{topology_file}_{trajectory_file}"
        
        if cache_key in self._cached_universes:
            return self._cached_universes[cache_key]
        
        try:
            universe = mda.Universe(topology_file, trajectory_file)
            self._cached_universes[cache_key] = universe
            return universe
        except Exception as e:
            logger.error(f"Failed to load universe: {e}")
            raise
