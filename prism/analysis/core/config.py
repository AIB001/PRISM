import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Optional
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="MDAnalysis")

def convert_numpy_types(obj):
    """Recursively convert numpy types to native Python types for JSON serialization"""
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    else:
        return obj

# Removed hardcoded rcParams - global publication style is now applied by individual plotting modules

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
