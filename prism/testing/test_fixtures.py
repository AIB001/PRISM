#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Test Fixtures and Mock Data Generation

Provides test fixtures, mock data generators, and test system builders
for comprehensive testing of PRISM components.
"""

import os
import json
import random
import tempfile
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from dataclasses import dataclass

from ..utils.logging_system import PrismLogger


@dataclass
class MockSystemConfig:
    """Configuration for mock molecular systems"""
    protein_atoms: int = 100
    ligand_atoms: int = 20
    water_molecules: int = 1000
    box_size: Tuple[float, float, float] = (5.0, 5.0, 5.0)
    temperature: float = 300.0
    force_field: str = "amber99sb-ildn"
    water_model: str = "tip3p"


class MockDataGenerator:
    """Generates mock data for testing PRISM components"""
    
    def __init__(self, seed: Optional[int] = None):
        self.logger = PrismLogger("mock_data_generator")
        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)
    
    def generate_coordinates(self, num_atoms: int, 
                           box_size: Tuple[float, float, float] = (5.0, 5.0, 5.0)) -> np.ndarray:
        """Generate random atomic coordinates"""
        coords = np.random.uniform(
            low=[0, 0, 0], 
            high=box_size, 
            size=(num_atoms, 3)
        )
        return coords
    
    def generate_velocities(self, num_atoms: int, temperature: float = 300.0) -> np.ndarray:
        """Generate Maxwell-Boltzmann distributed velocities"""
        # Simple approximation: normally distributed velocities
        # Real implementation would use proper Maxwell-Boltzmann distribution
        sigma = np.sqrt(temperature / 100.0)  # Simplified relation
        velocities = np.random.normal(0, sigma, size=(num_atoms, 3))
        return velocities
    
    def generate_forces(self, num_atoms: int, magnitude_range: Tuple[float, float] = (0.1, 10.0)) -> np.ndarray:
        """Generate random forces on atoms"""
        forces = np.random.uniform(
            low=-magnitude_range[1], 
            high=magnitude_range[1], 
            size=(num_atoms, 3)
        )
        return forces
    
    def generate_energy_trajectory(self, num_frames: int, 
                                 base_energy: float = -50000.0,
                                 fluctuation: float = 1000.0) -> List[float]:
        """Generate realistic energy trajectory"""
        energies = []
        current_energy = base_energy
        
        for i in range(num_frames):
            # Add random walk with mean reversion
            change = np.random.normal(0, fluctuation * 0.1)
            current_energy += change
            
            # Mean reversion towards base energy
            current_energy += (base_energy - current_energy) * 0.01
            
            energies.append(current_energy)
        
        return energies
    
    def generate_pmf_profile(self, reaction_coordinate: np.ndarray,
                           binding_energy: float = -8.5,
                           barrier_height: float = 5.0) -> np.ndarray:
        """Generate realistic PMF profile"""
        # Create a profile with a binding minimum and dissociation barrier
        profile = np.zeros_like(reaction_coordinate)
        
        for i, coord in enumerate(reaction_coordinate):
            if coord < 1.0:  # Bound state
                profile[i] = binding_energy + 2.0 * coord**2
            elif coord < 2.5:  # Barrier region
                profile[i] = binding_energy + barrier_height * np.exp(-(coord - 1.8)**2 / 0.2)
            else:  # Unbound state
                profile[i] = 0.0
        
        # Add some noise
        noise = np.random.normal(0, 0.5, size=len(profile))
        profile += noise
        
        return profile
    
    def generate_pull_data(self, num_steps: int, pull_rate: float = 0.01,
                          dt: float = 0.002) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate SMD pull data (time, distance, force)"""
        time_points = np.arange(num_steps) * dt
        distances = time_points * pull_rate
        
        # Generate realistic pull forces
        forces = np.zeros(num_steps)
        for i, dist in enumerate(distances):
            if dist < 1.0:  # Bound region - low force
                forces[i] = 50 + 20 * np.random.normal()
            elif dist < 2.0:  # Breaking region - high force
                forces[i] = 200 + 100 * np.exp(-(dist - 1.5)**2 / 0.1) + 30 * np.random.normal()
            else:  # Dissociated - low force
                forces[i] = 30 + 15 * np.random.normal()
        
        return time_points, distances, forces
    
    def generate_umbrella_windows(self, num_windows: int, 
                                distance_range: Tuple[float, float] = (0.5, 4.0)) -> List[float]:
        """Generate umbrella sampling window positions"""
        min_dist, max_dist = distance_range
        windows = np.linspace(min_dist, max_dist, num_windows)
        return windows.tolist()


class TestFixtures:
    """Manages test fixtures and temporary test environments"""
    
    def __init__(self, base_dir: Optional[Path] = None):
        self.base_dir = base_dir or Path(tempfile.mkdtemp(prefix="prism_test_fixtures_"))
        self.base_dir.mkdir(exist_ok=True)
        self.logger = PrismLogger("test_fixtures")
        self.mock_generator = MockDataGenerator()
        
        # Track created fixtures for cleanup
        self.created_fixtures: List[Path] = []
    
    def create_mock_system(self, system_name: str, 
                          config: Optional[MockSystemConfig] = None) -> Path:
        """Create a complete mock molecular system"""
        if config is None:
            config = MockSystemConfig()
        
        system_dir = self.base_dir / system_name
        system_dir.mkdir(exist_ok=True)
        self.created_fixtures.append(system_dir)
        
        # Generate coordinates
        protein_coords = self.mock_generator.generate_coordinates(
            config.protein_atoms, config.box_size
        )
        ligand_coords = self.mock_generator.generate_coordinates(
            config.ligand_atoms, config.box_size
        )
        
        # Create GRO file
        self._create_gro_file(system_dir / "system.gro", protein_coords, ligand_coords, config)
        
        # Create topology file
        self._create_topology_file(system_dir / "system.top", config)
        
        # Create index file
        self._create_index_file(system_dir / "index.ndx", config)
        
        # Create MDP files
        self._create_mdp_files(system_dir, config)
        
        # Create mock trajectory
        self._create_mock_trajectory(system_dir, config)
        
        self.logger.info(f"Created mock system: {system_dir}")
        return system_dir
    
    def create_mock_md_results(self, system_name: str) -> Path:
        """Create mock MD simulation results"""
        md_dir = self.base_dir / system_name / "GMX_PROLIG_MD"
        md_dir.mkdir(parents=True, exist_ok=True)
        self.created_fixtures.append(md_dir)
        
        # Create result files
        files_to_create = {
            "production.xtc": self._mock_trajectory_content(),
            "production.gro": self._mock_final_structure(),
            "production.edr": b"mock energy file",  # Binary mock
            "production.log": self._mock_log_content(),
            "system.top": self._mock_topology_content(),
            "index.ndx": self._mock_index_content()
        }
        
        for filename, content in files_to_create.items():
            file_path = md_dir / filename
            if isinstance(content, str):
                with open(file_path, 'w') as f:
                    f.write(content)
            else:
                with open(file_path, 'wb') as f:
                    f.write(content)
        
        return md_dir
    
    def create_mock_pmf_results(self, system_name: str) -> Dict[str, Any]:
        """Create mock PMF calculation results"""
        # Generate mock PMF data
        reaction_coord = np.linspace(0.5, 4.0, 40)
        pmf_profile = self.mock_generator.generate_pmf_profile(reaction_coord)
        
        # Find binding energy (minimum of PMF)
        binding_energy = float(np.min(pmf_profile))
        
        # Calculate error estimates
        error_estimates = np.random.uniform(0.5, 2.0, size=len(pmf_profile))
        
        pmf_results = {
            'binding_energy': binding_energy,
            'pmf_profile': pmf_profile.tolist(),
            'reaction_coordinate': reaction_coord.tolist(),
            'error_estimate': error_estimates.tolist(),
            'convergence': True,
            'method': 'WHAM',
            'temperature': 300.0,
            'units': 'kcal/mol'
        }
        
        # Save to file
        results_dir = self.base_dir / system_name / "pmf_results"
        results_dir.mkdir(parents=True, exist_ok=True)
        self.created_fixtures.append(results_dir)
        
        with open(results_dir / "pmf_analysis.json", 'w') as f:
            json.dump(pmf_results, f, indent=2)
        
        return pmf_results
    
    def create_performance_benchmark_data(self, num_systems: int = 5) -> List[Dict[str, Any]]:
        """Create benchmark data for performance testing"""
        benchmark_data = []
        
        for i in range(num_systems):
            system_size = random.randint(1000, 50000)  # atoms
            simulation_time = random.uniform(1.0, 100.0)  # ns
            
            # Generate realistic performance metrics
            ns_per_day = max(0.1, 50 - system_size / 1000 + random.uniform(-5, 5))
            memory_usage = system_size * 0.1 + random.uniform(0, 500)  # MB
            cpu_efficiency = random.uniform(70, 95)  # percent
            
            benchmark_data.append({
                'system_id': f"benchmark_system_{i+1}",
                'system_size': system_size,
                'simulation_time_ns': simulation_time,
                'ns_per_day': ns_per_day,
                'memory_usage_mb': memory_usage,
                'cpu_efficiency_percent': cpu_efficiency,
                'execution_time_hours': simulation_time / ns_per_day * 24,
                'force_field': random.choice(['amber99sb-ildn', 'charmm36', 'opls-aa'])
            })
        
        return benchmark_data
    
    def _create_gro_file(self, file_path: Path, protein_coords: np.ndarray, 
                        ligand_coords: np.ndarray, config: MockSystemConfig):
        """Create mock GRO file"""
        total_atoms = len(protein_coords) + len(ligand_coords)
        
        content = f"Mock system for testing\n{total_atoms}\n"
        
        # Add protein atoms
        for i, coord in enumerate(protein_coords):
            content += f"{i+1:5d}PRO     CA{i+1:5d}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}\n"
        
        # Add ligand atoms
        for i, coord in enumerate(ligand_coords):
            atom_num = len(protein_coords) + i + 1
            content += f"{atom_num:5d}LIG     C{i+1:6d}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}\n"
        
        # Add box dimensions
        content += f"{config.box_size[0]:10.5f}{config.box_size[1]:10.5f}{config.box_size[2]:10.5f}\n"
        
        with open(file_path, 'w') as f:
            f.write(content)
    
    def _create_topology_file(self, file_path: Path, config: MockSystemConfig):
        """Create mock topology file"""
        content = f"""; Mock topology file
#include "{config.force_field}.ff/forcefield.itp"

[ moleculetype ]
; Name      nrexcl
Protein     3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1   CT    1     PRO    CA    1    0.0   12.01

[ moleculetype ]
; Name      nrexcl  
LIG         3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1   C     1     LIG    C1    1    0.0   12.01

[ system ]
Mock System

[ molecules ]
Protein  1
LIG      1
"""
        
        with open(file_path, 'w') as f:
            f.write(content)
    
    def _create_index_file(self, file_path: Path, config: MockSystemConfig):
        """Create mock index file"""
        protein_indices = " ".join(str(i+1) for i in range(config.protein_atoms))
        ligand_indices = " ".join(str(i+config.protein_atoms+1) for i in range(config.ligand_atoms))
        
        content = f"""[ System ]
{protein_indices} {ligand_indices}

[ Protein ]
{protein_indices}

[ LIG ]
{ligand_indices}

[ Protein_LIG ]
{protein_indices} {ligand_indices}
"""
        
        with open(file_path, 'w') as f:
            f.write(content)
    
    def _create_mdp_files(self, system_dir: Path, config: MockSystemConfig):
        """Create mock MDP files"""
        mdp_content = f"""; Mock MDP file
integrator = md
dt = 0.002
nsteps = 1000
ref_t = {config.temperature}
"""
        
        mdp_files = ["em.mdp", "nvt.mdp", "npt.mdp", "production.mdp"]
        for mdp_file in mdp_files:
            with open(system_dir / mdp_file, 'w') as f:
                f.write(mdp_content)
    
    def _create_mock_trajectory(self, system_dir: Path, config: MockSystemConfig):
        """Create mock trajectory files"""
        # For testing, create simple text files
        # Real implementation would create proper XTC files
        
        with open(system_dir / "trajectory.xtc", 'w') as f:
            f.write("Mock XTC trajectory file")
        
        with open(system_dir / "trajectory.trr", 'w') as f:
            f.write("Mock TRR trajectory file")
    
    def _mock_trajectory_content(self) -> str:
        """Mock trajectory file content"""
        return "Mock trajectory data for testing"
    
    def _mock_final_structure(self) -> str:
        """Mock final structure content"""
        return """Mock final structure
    2
    1PRO     CA    1   1.500   1.500   1.500
    2LIG     C1    2   2.500   2.500   2.500
   5.00000   5.00000   5.00000
"""
    
    def _mock_log_content(self) -> str:
        """Mock log file content"""
        return """GROMACS log file mock
Step    Time    Potential    Kinetic    Total
0       0.000   -50000       5000       -45000
1000    2.000   -49800       4900       -44900
"""
    
    def _mock_topology_content(self) -> str:
        """Mock topology content"""
        return """; Mock topology
[ system ]
Mock System

[ molecules ]
Protein 1
LIG     1
"""
    
    def _mock_index_content(self) -> str:
        """Mock index content"""
        return """[ System ]
1 2

[ Protein ]
1

[ LIG ]
2
"""
    
    def cleanup(self):
        """Cleanup all created fixtures"""
        import shutil
        
        for fixture_path in self.created_fixtures:
            try:
                if fixture_path.exists():
                    if fixture_path.is_dir():
                        shutil.rmtree(fixture_path)
                    else:
                        fixture_path.unlink()
            except Exception as e:
                self.logger.warning(f"Failed to cleanup {fixture_path}: {e}")
        
        self.created_fixtures.clear()
        
        # Cleanup base directory
        try:
            if self.base_dir.exists():
                shutil.rmtree(self.base_dir)
        except Exception as e:
            self.logger.warning(f"Failed to cleanup base directory: {e}")


class SystemFixtures:
    """High-level fixtures for common test scenarios"""
    
    def __init__(self):
        self.fixtures = TestFixtures()
        self.logger = PrismLogger("system_fixtures")
    
    def small_protein_ligand_system(self) -> Path:
        """Create small protein-ligand system for fast tests"""
        config = MockSystemConfig(
            protein_atoms=50,
            ligand_atoms=10,
            water_molecules=100,
            box_size=(3.0, 3.0, 3.0)
        )
        return self.fixtures.create_mock_system("small_system", config)
    
    def medium_protein_ligand_system(self) -> Path:
        """Create medium protein-ligand system for standard tests"""
        config = MockSystemConfig(
            protein_atoms=200,
            ligand_atoms=30,
            water_molecules=1000,
            box_size=(5.0, 5.0, 5.0)
        )
        return self.fixtures.create_mock_system("medium_system", config)
    
    def large_protein_ligand_system(self) -> Path:
        """Create large protein-ligand system for performance tests"""
        config = MockSystemConfig(
            protein_atoms=1000,
            ligand_atoms=100,
            water_molecules=10000,
            box_size=(10.0, 10.0, 10.0)
        )
        return self.fixtures.create_mock_system("large_system", config)
    
    def pmf_test_scenario(self) -> Dict[str, Any]:
        """Create complete PMF test scenario"""
        # Create system
        system_dir = self.medium_protein_ligand_system()
        
        # Create MD results
        md_results = self.fixtures.create_mock_md_results("medium_system")
        
        # Create PMF results
        pmf_results = self.fixtures.create_mock_pmf_results("medium_system")
        
        return {
            'system_dir': system_dir,
            'md_results': md_results,
            'pmf_results': pmf_results,
            'expected_binding_energy': pmf_results['binding_energy']
        }
    
    def cleanup_all(self):
        """Cleanup all fixtures"""
        self.fixtures.cleanup()


# Convenience functions
def create_test_fixtures() -> TestFixtures:
    """Create test fixtures manager"""
    return TestFixtures()


def create_system_fixtures() -> SystemFixtures:
    """Create system fixtures manager"""
    return SystemFixtures()