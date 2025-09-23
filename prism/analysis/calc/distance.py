import logging
import numpy as np
from typing import Dict
from ..config import AnalysisConfig

logger = logging.getLogger(__name__)

class DistanceAnalyzer:
    """Distance statistics analysis module"""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
    
    def calculate_distance_statistics(self, distance_results: Dict[str, np.ndarray]) -> Dict[str, Dict[str, float]]:
        """Calculate distance statistics from pre-calculated results"""
        distance_stats = {}
        
        for residue_key, distances in distance_results.items():
            try:
                distances_angstrom = distances * 10.0
                
                distance_stats[residue_key] = {
                    'mean': float(np.mean(distances_angstrom)),
                    'min': float(np.min(distances_angstrom)),
                    'max': float(np.max(distances_angstrom)),
                    'std': float(np.std(distances_angstrom))
                }
            except Exception as e:
                logger.warning(f"Failed to analyze distance statistics for {residue_key}: {e}")
                continue
        
        return distance_stats
