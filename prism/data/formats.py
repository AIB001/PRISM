#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data Format Handlers

Specialized handlers for molecular dynamics data formats including
trajectories, energy data, PMF results, and analysis outputs.
"""

import os
import struct
import json
import csv
import numpy as np
from pathlib import Path
from typing import Dict, Any, List, Optional, Union, Tuple, Iterator
from dataclasses import dataclass, asdict
from enum import Enum
from abc import ABC, abstractmethod
import re
import time

from ..utils.logging_system import PrismLogger


class DataFormat(Enum):
    """Supported molecular dynamics data formats"""
    XTC = "xtc"           # GROMACS compressed trajectory
    TRR = "trr"           # GROMACS full precision trajectory  
    GRO = "gro"           # GROMACS structure format
    PDB = "pdb"           # Protein Data Bank format
    XVG = "xvg"           # GROMACS plot format
    EDR = "edr"           # GROMACS energy format
    TOP = "top"           # GROMACS topology
    NDX = "ndx"           # GROMACS index format
    PMF = "pmf"           # PMF analysis results
    WHAM = "wham"         # WHAM output format
    JSON = "json"         # JSON structured data
    CSV = "csv"           # Comma-separated values


@dataclass
class FrameData:
    """Single trajectory frame data"""
    frame_number: int
    time_ps: float
    coordinates: np.ndarray  # Shape: (n_atoms, 3)
    velocities: Optional[np.ndarray] = None
    forces: Optional[np.ndarray] = None
    box_vectors: Optional[np.ndarray] = None
    
    def __post_init__(self):
        """Validate frame data"""
        if not isinstance(self.coordinates, np.ndarray):
            self.coordinates = np.array(self.coordinates)
        
        if self.coordinates.ndim != 2 or self.coordinates.shape[1] != 3:
            raise ValueError("Coordinates must be (n_atoms, 3) array")


@dataclass
class EnergyData:
    """Energy data from MD simulations"""
    time_ps: List[float]
    potential_energy: List[float]
    kinetic_energy: List[float]
    total_energy: List[float]
    temperature: List[float]
    pressure: List[float]
    volume: List[float]
    density: List[float]
    custom_terms: Dict[str, List[float]] = None
    
    def __post_init__(self):
        if self.custom_terms is None:
            self.custom_terms = {}
    
    @property
    def n_frames(self) -> int:
        return len(self.time_ps)
    
    def get_energy_term(self, term_name: str) -> Optional[List[float]]:
        """Get specific energy term"""
        standard_terms = {
            'potential': self.potential_energy,
            'kinetic': self.kinetic_energy,
            'total': self.total_energy,
            'temperature': self.temperature,
            'pressure': self.pressure,
            'volume': self.volume,
            'density': self.density
        }
        
        return standard_terms.get(term_name) or self.custom_terms.get(term_name)


@dataclass
class PMFResult:
    """PMF calculation results"""
    reaction_coordinate: List[float]
    pmf_values: List[float]
    error_estimates: List[float]
    binding_energy: float
    convergence_status: bool
    method: str
    temperature: float
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}
    
    @property
    def n_points(self) -> int:
        return len(self.reaction_coordinate)


class BaseFormatHandler(ABC):
    """Abstract base class for data format handlers"""
    
    def __init__(self, format_type: DataFormat):
        self.format_type = format_type
        self.logger = PrismLogger(f"format_handler.{format_type.value}")
    
    @abstractmethod
    def can_handle(self, file_path: Path) -> bool:
        """Check if this handler can process the file"""
        pass
    
    @abstractmethod
    def read(self, file_path: Path, **kwargs) -> Any:
        """Read data from file"""
        pass
    
    @abstractmethod
    def write(self, data: Any, file_path: Path, **kwargs) -> bool:
        """Write data to file"""
        pass
    
    def validate_data(self, data: Any) -> bool:
        """Validate data format"""
        return True


class TrajectoryHandler(BaseFormatHandler):
    """Handler for trajectory files (XTC, TRR, GRO)"""
    
    def __init__(self, format_type: DataFormat = DataFormat.XTC):
        super().__init__(format_type)
    
    def can_handle(self, file_path: Path) -> bool:
        """Check if file is a supported trajectory format"""
        suffix = file_path.suffix.lower()
        return suffix in ['.xtc', '.trr', '.gro', '.pdb']
    
    def read(self, file_path: Path, **kwargs) -> Iterator[FrameData]:
        """Read trajectory frames"""
        suffix = file_path.suffix.lower()
        
        if suffix == '.gro':
            yield from self._read_gro_file(file_path, **kwargs)
        elif suffix == '.pdb':
            yield from self._read_pdb_file(file_path, **kwargs)
        elif suffix in ['.xtc', '.trr']:
            yield from self._read_binary_trajectory(file_path, **kwargs)
        else:
            raise ValueError(f"Unsupported trajectory format: {suffix}")
    
    def write(self, frames: List[FrameData], file_path: Path, **kwargs) -> bool:
        """Write trajectory frames to file"""
        suffix = file_path.suffix.lower()
        
        try:
            if suffix == '.gro':
                return self._write_gro_file(frames, file_path, **kwargs)
            elif suffix == '.pdb':
                return self._write_pdb_file(frames, file_path, **kwargs)
            else:
                self.logger.error(f"Writing not supported for format: {suffix}")
                return False
        except Exception as e:
            self.logger.error(f"Failed to write trajectory: {e}")
            return False
    
    def _read_gro_file(self, file_path: Path, **kwargs) -> Iterator[FrameData]:
        """Read GROMACS GRO format"""
        try:
            with open(file_path, 'r') as f:
                frame_num = 0
                
                while True:
                    # Read header
                    title = f.readline().strip()
                    if not title:
                        break
                    
                    # Read number of atoms
                    n_atoms_line = f.readline().strip()
                    if not n_atoms_line:
                        break
                    
                    n_atoms = int(n_atoms_line)
                    
                    # Read coordinates
                    coordinates = []
                    for _ in range(n_atoms):
                        line = f.readline()
                        if not line:
                            return
                        
                        # Parse GRO line: resid resname atomname atomid x y z [vx vy vz]
                        parts = line.split()
                        if len(parts) < 6:
                            continue
                        
                        x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                        coordinates.append([x, y, z])
                    
                    # Read box vectors
                    box_line = f.readline().strip()
                    box_vectors = None
                    if box_line:
                        box_parts = box_line.split()
                        if len(box_parts) >= 3:
                            box_vectors = np.array([
                                [float(box_parts[0]), 0, 0],
                                [0, float(box_parts[1]), 0],
                                [0, 0, float(box_parts[2])]
                            ])
                    
                    yield FrameData(
                        frame_number=frame_num,
                        time_ps=frame_num * kwargs.get('dt', 1.0),
                        coordinates=np.array(coordinates),
                        box_vectors=box_vectors
                    )
                    
                    frame_num += 1
                    
                    # For single frame GRO files
                    if kwargs.get('single_frame', True):
                        break
        
        except Exception as e:
            self.logger.error(f"Error reading GRO file {file_path}: {e}")
    
    def _read_pdb_file(self, file_path: Path, **kwargs) -> Iterator[FrameData]:
        """Read PDB format"""
        try:
            with open(file_path, 'r') as f:
                coordinates = []
                frame_num = 0
                
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # Parse PDB atom record
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coordinates.append([x, y, z])
                    
                    elif line.startswith('END') or line.startswith('ENDMDL'):
                        if coordinates:
                            yield FrameData(
                                frame_number=frame_num,
                                time_ps=frame_num * kwargs.get('dt', 1.0),
                                coordinates=np.array(coordinates)
                            )
                            coordinates = []
                            frame_num += 1
                
                # Handle case where file doesn't end with END
                if coordinates:
                    yield FrameData(
                        frame_number=frame_num,
                        time_ps=frame_num * kwargs.get('dt', 1.0),
                        coordinates=np.array(coordinates)
                    )
        
        except Exception as e:
            self.logger.error(f"Error reading PDB file {file_path}: {e}")
    
    def _read_binary_trajectory(self, file_path: Path, **kwargs) -> Iterator[FrameData]:
        """Read binary trajectory files (XTC/TRR)"""
        # This would require specific libraries like MDTraj or MDAnalysis
        # For now, provide a placeholder implementation
        self.logger.warning(f"Binary trajectory reading not implemented for {file_path}")
        return iter([])
    
    def _write_gro_file(self, frames: List[FrameData], file_path: Path, **kwargs) -> bool:
        """Write GROMACS GRO format"""
        try:
            with open(file_path, 'w') as f:
                for frame in frames:
                    # Write title
                    title = kwargs.get('title', f'Frame {frame.frame_number}')
                    f.write(f"{title}\n")
                    
                    # Write number of atoms
                    n_atoms = len(frame.coordinates)
                    f.write(f"{n_atoms}\n")
                    
                    # Write atom coordinates
                    for i, (x, y, z) in enumerate(frame.coordinates):
                        atom_name = kwargs.get('atom_names', ['CA'] * n_atoms)[i]
                        res_name = kwargs.get('residue_names', ['PRO'] * n_atoms)[i]
                        res_id = kwargs.get('residue_ids', list(range(1, n_atoms + 1)))[i]
                        
                        f.write(f"{res_id:5d}{res_name:<5}{atom_name:>5}{i+1:5d}"
                               f"{x:8.3f}{y:8.3f}{z:8.3f}\n")
                    
                    # Write box vectors
                    if frame.box_vectors is not None:
                        box = frame.box_vectors
                        f.write(f"{box[0,0]:10.5f}{box[1,1]:10.5f}{box[2,2]:10.5f}\n")
                    else:
                        # Default box
                        f.write("   0.00000   0.00000   0.00000\n")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing GRO file {file_path}: {e}")
            return False
    
    def _write_pdb_file(self, frames: List[FrameData], file_path: Path, **kwargs) -> bool:
        """Write PDB format"""
        try:
            with open(file_path, 'w') as f:
                for frame_idx, frame in enumerate(frames):
                    if len(frames) > 1:
                        f.write(f"MODEL     {frame_idx + 1}\n")
                    
                    for i, (x, y, z) in enumerate(frame.coordinates):
                        atom_name = kwargs.get('atom_names', ['CA'] * len(frame.coordinates))[i]
                        res_name = kwargs.get('residue_names', ['ALA'] * len(frame.coordinates))[i]
                        res_id = kwargs.get('residue_ids', list(range(1, len(frame.coordinates) + 1)))[i]
                        
                        f.write(f"ATOM  {i+1:5d}  {atom_name:<3} {res_name} A{res_id:4d}    "
                               f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                    
                    if len(frames) > 1:
                        f.write("ENDMDL\n")
                
                f.write("END\n")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing PDB file {file_path}: {e}")
            return False


class EnergyHandler(BaseFormatHandler):
    """Handler for energy data files (XVG, EDR)"""
    
    def __init__(self, format_type: DataFormat = DataFormat.XVG):
        super().__init__(format_type)
    
    def can_handle(self, file_path: Path) -> bool:
        """Check if file is a supported energy format"""
        suffix = file_path.suffix.lower()
        return suffix in ['.xvg', '.edr', '.csv']
    
    def read(self, file_path: Path, **kwargs) -> EnergyData:
        """Read energy data from file"""
        suffix = file_path.suffix.lower()
        
        if suffix == '.xvg':
            return self._read_xvg_file(file_path, **kwargs)
        elif suffix == '.csv':
            return self._read_csv_file(file_path, **kwargs)
        elif suffix == '.edr':
            return self._read_edr_file(file_path, **kwargs)
        else:
            raise ValueError(f"Unsupported energy format: {suffix}")
    
    def write(self, data: EnergyData, file_path: Path, **kwargs) -> bool:
        """Write energy data to file"""
        suffix = file_path.suffix.lower()
        
        try:
            if suffix == '.xvg':
                return self._write_xvg_file(data, file_path, **kwargs)
            elif suffix == '.csv':
                return self._write_csv_file(data, file_path, **kwargs)
            else:
                self.logger.error(f"Writing not supported for format: {suffix}")
                return False
        except Exception as e:
            self.logger.error(f"Failed to write energy data: {e}")
            return False
    
    def _read_xvg_file(self, file_path: Path, **kwargs) -> EnergyData:
        """Read GROMACS XVG format"""
        time_data = []
        energy_columns = {}
        column_names = []
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # Skip comments and empty lines
                    if line.startswith('#') or line.startswith('@') or not line:
                        # Extract column names from legends
                        if 's legend' in line:
                            legend_match = re.search(r'"([^"]+)"', line)
                            if legend_match:
                                column_names.append(legend_match.group(1))
                        continue
                    
                    # Parse data line
                    parts = line.split()
                    if not parts:
                        continue
                    
                    try:
                        values = [float(x) for x in parts]
                        
                        if len(values) >= 1:
                            time_data.append(values[0])
                            
                            # Store additional columns
                            for i, value in enumerate(values[1:]):
                                col_name = f"column_{i+1}"
                                if i < len(column_names):
                                    col_name = column_names[i]
                                
                                if col_name not in energy_columns:
                                    energy_columns[col_name] = []
                                energy_columns[col_name].append(value)
                                
                    except ValueError:
                        continue
            
            # Map to standard energy terms
            potential_energy = energy_columns.get('Potential', energy_columns.get('column_1', [0.0] * len(time_data)))
            kinetic_energy = energy_columns.get('Kinetic-En.', energy_columns.get('column_2', [0.0] * len(time_data)))
            total_energy = energy_columns.get('Total-Energy', energy_columns.get('column_3', [0.0] * len(time_data)))
            temperature = energy_columns.get('Temperature', energy_columns.get('column_4', [300.0] * len(time_data)))
            pressure = energy_columns.get('Pressure', energy_columns.get('column_5', [1.0] * len(time_data)))
            volume = energy_columns.get('Volume', energy_columns.get('column_6', [1000.0] * len(time_data)))
            density = energy_columns.get('Density', energy_columns.get('column_7', [1000.0] * len(time_data)))
            
            return EnergyData(
                time_ps=time_data,
                potential_energy=potential_energy,
                kinetic_energy=kinetic_energy,
                total_energy=total_energy,
                temperature=temperature,
                pressure=pressure,
                volume=volume,
                density=density,
                custom_terms=energy_columns
            )
        
        except Exception as e:
            self.logger.error(f"Error reading XVG file {file_path}: {e}")
            # Return empty data structure
            return EnergyData([], [], [], [], [], [], [], [])
    
    def _read_csv_file(self, file_path: Path, **kwargs) -> EnergyData:
        """Read CSV energy data"""
        try:
            data = {}
            
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f)
                
                for row in reader:
                    for key, value in row.items():
                        if key not in data:
                            data[key] = []
                        try:
                            data[key].append(float(value))
                        except ValueError:
                            data[key].append(0.0)
            
            return EnergyData(
                time_ps=data.get('time', data.get('Time', [])),
                potential_energy=data.get('potential', data.get('Potential', [])),
                kinetic_energy=data.get('kinetic', data.get('Kinetic', [])),
                total_energy=data.get('total', data.get('Total', [])),
                temperature=data.get('temperature', data.get('Temperature', [])),
                pressure=data.get('pressure', data.get('Pressure', [])),
                volume=data.get('volume', data.get('Volume', [])),
                density=data.get('density', data.get('Density', [])),
                custom_terms={k: v for k, v in data.items() 
                             if k not in ['time', 'Time', 'potential', 'Potential', 
                                        'kinetic', 'Kinetic', 'total', 'Total',
                                        'temperature', 'Temperature', 'pressure', 'Pressure',
                                        'volume', 'Volume', 'density', 'Density']}
            )
        
        except Exception as e:
            self.logger.error(f"Error reading CSV file {file_path}: {e}")
            return EnergyData([], [], [], [], [], [], [], [])
    
    def _read_edr_file(self, file_path: Path, **kwargs) -> EnergyData:
        """Read GROMACS EDR binary format"""
        # This would require specific libraries or GROMACS tools
        self.logger.warning(f"EDR reading not implemented for {file_path}")
        return EnergyData([], [], [], [], [], [], [], [])
    
    def _write_xvg_file(self, data: EnergyData, file_path: Path, **kwargs) -> bool:
        """Write GROMACS XVG format"""
        try:
            with open(file_path, 'w') as f:
                # Write header
                f.write("# This file was created by PRISM\n")
                f.write("# Generated on {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
                f.write("@    title \"Energy Data\"\n")
                f.write("@    xaxis  label \"Time (ps)\"\n")
                f.write("@    yaxis  label \"Energy (kJ/mol)\"\n")
                
                # Write column legends
                f.write("@    s0 legend \"Potential\"\n")
                f.write("@    s1 legend \"Kinetic\"\n")
                f.write("@    s2 legend \"Total\"\n")
                f.write("@    s3 legend \"Temperature\"\n")
                
                # Write data
                for i in range(len(data.time_ps)):
                    f.write(f"{data.time_ps[i]:12.3f} {data.potential_energy[i]:12.3f} "
                           f"{data.kinetic_energy[i]:12.3f} {data.total_energy[i]:12.3f} "
                           f"{data.temperature[i]:8.2f}\n")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing XVG file {file_path}: {e}")
            return False
    
    def _write_csv_file(self, data: EnergyData, file_path: Path, **kwargs) -> bool:
        """Write CSV energy data"""
        try:
            with open(file_path, 'w', newline='') as f:
                fieldnames = ['time', 'potential_energy', 'kinetic_energy', 'total_energy',
                             'temperature', 'pressure', 'volume', 'density']
                
                # Add custom terms
                fieldnames.extend(data.custom_terms.keys())
                
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                
                for i in range(len(data.time_ps)):
                    row = {
                        'time': data.time_ps[i],
                        'potential_energy': data.potential_energy[i],
                        'kinetic_energy': data.kinetic_energy[i],
                        'total_energy': data.total_energy[i],
                        'temperature': data.temperature[i],
                        'pressure': data.pressure[i],
                        'volume': data.volume[i],
                        'density': data.density[i]
                    }
                    
                    # Add custom terms
                    for term_name, values in data.custom_terms.items():
                        if i < len(values):
                            row[term_name] = values[i]
                    
                    writer.writerow(row)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing CSV file {file_path}: {e}")
            return False


class PMFDataHandler(BaseFormatHandler):
    """Handler for PMF analysis data"""
    
    def __init__(self, format_type: DataFormat = DataFormat.PMF):
        super().__init__(format_type)
    
    def can_handle(self, file_path: Path) -> bool:
        """Check if file is a supported PMF format"""
        suffix = file_path.suffix.lower()
        return suffix in ['.pmf', '.wham', '.json', '.dat']
    
    def read(self, file_path: Path, **kwargs) -> PMFResult:
        """Read PMF data from file"""
        suffix = file_path.suffix.lower()
        
        if suffix == '.json':
            return self._read_json_pmf(file_path, **kwargs)
        elif suffix in ['.pmf', '.wham', '.dat']:
            return self._read_text_pmf(file_path, **kwargs)
        else:
            raise ValueError(f"Unsupported PMF format: {suffix}")
    
    def write(self, data: PMFResult, file_path: Path, **kwargs) -> bool:
        """Write PMF data to file"""
        suffix = file_path.suffix.lower()
        
        try:
            if suffix == '.json':
                return self._write_json_pmf(data, file_path, **kwargs)
            elif suffix in ['.pmf', '.dat']:
                return self._write_text_pmf(data, file_path, **kwargs)
            else:
                self.logger.error(f"Writing not supported for format: {suffix}")
                return False
        except Exception as e:
            self.logger.error(f"Failed to write PMF data: {e}")
            return False
    
    def _read_json_pmf(self, file_path: Path, **kwargs) -> PMFResult:
        """Read PMF data from JSON format"""
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            return PMFResult(
                reaction_coordinate=data.get('reaction_coordinate', []),
                pmf_values=data.get('pmf_values', []),
                error_estimates=data.get('error_estimates', []),
                binding_energy=data.get('binding_energy', 0.0),
                convergence_status=data.get('convergence_status', False),
                method=data.get('method', 'unknown'),
                temperature=data.get('temperature', 300.0),
                metadata=data.get('metadata', {})
            )
        
        except Exception as e:
            self.logger.error(f"Error reading JSON PMF file {file_path}: {e}")
            return PMFResult([], [], [], 0.0, False, 'unknown', 300.0)
    
    def _read_text_pmf(self, file_path: Path, **kwargs) -> PMFResult:
        """Read PMF data from text format"""
        reaction_coord = []
        pmf_values = []
        error_estimates = []
        metadata = {}
        
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # Skip comments and extract metadata
                    if line.startswith('#') or line.startswith('@'):
                        if 'temperature' in line.lower():
                            temp_match = re.search(r'(\d+\.?\d*)', line)
                            if temp_match:
                                metadata['temperature'] = float(temp_match.group(1))
                        continue
                    
                    if not line:
                        continue
                    
                    # Parse data line
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            coord = float(parts[0])
                            pmf = float(parts[1])
                            error = float(parts[2]) if len(parts) > 2 else 0.5
                            
                            reaction_coord.append(coord)
                            pmf_values.append(pmf)
                            error_estimates.append(error)
                            
                        except ValueError:
                            continue
            
            # Calculate binding energy (minimum PMF)
            binding_energy = min(pmf_values) if pmf_values else 0.0
            
            return PMFResult(
                reaction_coordinate=reaction_coord,
                pmf_values=pmf_values,
                error_estimates=error_estimates,
                binding_energy=binding_energy,
                convergence_status=len(pmf_values) > 10,  # Simple check
                method='WHAM',
                temperature=metadata.get('temperature', 300.0),
                metadata=metadata
            )
        
        except Exception as e:
            self.logger.error(f"Error reading text PMF file {file_path}: {e}")
            return PMFResult([], [], [], 0.0, False, 'unknown', 300.0)
    
    def _write_json_pmf(self, data: PMFResult, file_path: Path, **kwargs) -> bool:
        """Write PMF data to JSON format"""
        try:
            output_data = asdict(data)
            
            with open(file_path, 'w') as f:
                json.dump(output_data, f, indent=2)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing JSON PMF file {file_path}: {e}")
            return False
    
    def _write_text_pmf(self, data: PMFResult, file_path: Path, **kwargs) -> bool:
        """Write PMF data to text format"""
        try:
            with open(file_path, 'w') as f:
                # Write header
                f.write("# PMF Analysis Results\n")
                f.write(f"# Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"# Method: {data.method}\n")
                f.write(f"# Temperature: {data.temperature} K\n")
                f.write(f"# Binding Energy: {data.binding_energy:.3f} kJ/mol\n")
                f.write("# Reaction_Coordinate(nm)  PMF(kJ/mol)  Error(kJ/mol)\n")
                
                # Write data
                for i in range(len(data.reaction_coordinate)):
                    f.write(f"{data.reaction_coordinate[i]:8.3f} {data.pmf_values[i]:12.3f} "
                           f"{data.error_estimates[i]:8.3f}\n")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error writing text PMF file {file_path}: {e}")
            return False


# Format registry for automatic detection
class FormatRegistry:
    """Registry for data format handlers"""
    
    def __init__(self):
        self.handlers: Dict[DataFormat, BaseFormatHandler] = {}
        self.logger = PrismLogger("format_registry")
        
        # Register default handlers
        self.register_handler(DataFormat.XTC, TrajectoryHandler(DataFormat.XTC))
        self.register_handler(DataFormat.GRO, TrajectoryHandler(DataFormat.GRO))
        self.register_handler(DataFormat.PDB, TrajectoryHandler(DataFormat.PDB))
        self.register_handler(DataFormat.XVG, EnergyHandler(DataFormat.XVG))
        self.register_handler(DataFormat.PMF, PMFDataHandler(DataFormat.PMF))
    
    def register_handler(self, format_type: DataFormat, handler: BaseFormatHandler):
        """Register a format handler"""
        self.handlers[format_type] = handler
        self.logger.info(f"Registered handler for {format_type.value} format")
    
    def get_handler(self, file_path: Path) -> Optional[BaseFormatHandler]:
        """Get appropriate handler for file"""
        for handler in self.handlers.values():
            if handler.can_handle(file_path):
                return handler
        return None
    
    def read_data(self, file_path: Path, **kwargs) -> Any:
        """Read data using appropriate handler"""
        handler = self.get_handler(file_path)
        if handler:
            return handler.read(file_path, **kwargs)
        else:
            raise ValueError(f"No handler found for file: {file_path}")
    
    def write_data(self, data: Any, file_path: Path, **kwargs) -> bool:
        """Write data using appropriate handler"""
        handler = self.get_handler(file_path)
        if handler:
            return handler.write(data, file_path, **kwargs)
        else:
            raise ValueError(f"No handler found for file: {file_path}")


# Global format registry
_global_registry = FormatRegistry()

def get_format_registry() -> FormatRegistry:
    """Get global format registry"""
    return _global_registry


# Convenience functions
def read_trajectory(file_path: Union[str, Path]) -> Iterator[FrameData]:
    """Read trajectory file"""
    handler = TrajectoryHandler()
    return handler.read(Path(file_path))


def read_energy_data(file_path: Union[str, Path]) -> EnergyData:
    """Read energy data file"""
    handler = EnergyHandler()
    return handler.read(Path(file_path))


def read_pmf_data(file_path: Union[str, Path]) -> PMFResult:
    """Read PMF data file"""
    handler = PMFDataHandler()
    return handler.read(Path(file_path))


def auto_read_data(file_path: Union[str, Path], **kwargs) -> Any:
    """Automatically detect format and read data"""
    return _global_registry.read_data(Path(file_path), **kwargs)


def auto_write_data(data: Any, file_path: Union[str, Path], **kwargs) -> bool:
    """Automatically detect format and write data"""
    return _global_registry.write_data(data, Path(file_path), **kwargs)