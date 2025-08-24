import os
from dataclasses import dataclass
from typing import Optional

@dataclass
class AnalysisConfig:
    """Configuration parameters for analysis"""
    ligand_resname: str = "LIG"
    contact_cutoff_nm: float = 0.4
    contact_enter_threshold_nm: float = 0.3
    contact_exit_threshold_nm: float = 0.45
    hbond_distance_cutoff_nm: float = 0.35
    hbond_angle_cutoff_deg: float = 120.0
    bond_length_threshold_nm: float = 0.15
    smooth_window: int = 11
    min_frames_for_smoothing: int = 10
    parallel_workers: Optional[int] = None
    distance_analysis: bool = True  
    distance_cutoff_nm: float = 0.5  
    
    def __post_init__(self):
        if self.parallel_workers is None:
            self.parallel_workers = min(os.cpu_count(), 8)