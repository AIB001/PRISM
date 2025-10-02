#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM Data API

Data handling, format conversion, and interoperability with other MD software.
"""

import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
from dataclasses import dataclass, asdict
from enum import Enum
import csv
import xml.etree.ElementTree as ET

from ..utils.logging_system import PrismLogger
from .exceptions import *


class DataFormat(Enum):
    """Supported data formats"""
    GROMACS = "gromacs"
    AMBER = "amber"
    CHARMM = "charmm"
    NAMD = "namd"
    OPENMM = "openmm"
    XVG = "xvg"
    DAT = "dat"
    CSV = "csv"
    JSON = "json"
    XML = "xml"
    HDF5 = "hdf5"


@dataclass
class PMFProfileData:
    """Structure for PMF profile data"""
    reaction_coordinate: List[float]
    pmf_energy: List[float]
    error_estimate: Optional[List[float]] = None
    sampling_count: Optional[List[int]] = None
    metadata: Optional[Dict[str, Any]] = None
    
    def __post_init__(self):
        # Validate data consistency
        if len(self.reaction_coordinate) != len(self.pmf_energy):
            raise ValidationError("Reaction coordinate and PMF energy arrays must have same length")
        
        if self.error_estimate and len(self.error_estimate) != len(self.pmf_energy):
            raise ValidationError("Error estimate array must match PMF energy array length")
        
        if self.sampling_count and len(self.sampling_count) != len(self.pmf_energy):
            raise ValidationError("Sampling count array must match PMF energy array length")


@dataclass
class TrajectoryData:
    """Structure for trajectory data"""
    coordinates: np.ndarray  # Shape: (n_frames, n_atoms, 3)
    topology: Optional[Dict[str, Any]] = None
    box_vectors: Optional[np.ndarray] = None  # Shape: (n_frames, 3, 3)
    time: Optional[List[float]] = None
    metadata: Optional[Dict[str, Any]] = None


@dataclass
class EnergyData:
    """Structure for energy data"""
    time: List[float]
    potential_energy: Optional[List[float]] = None
    kinetic_energy: Optional[List[float]] = None
    total_energy: Optional[List[float]] = None
    temperature: Optional[List[float]] = None
    pressure: Optional[List[float]] = None
    volume: Optional[List[float]] = None
    additional_terms: Optional[Dict[str, List[float]]] = None
    metadata: Optional[Dict[str, Any]] = None


class FormatConverter:
    """Converter for different molecular dynamics data formats"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("format_converter")
        self._registered_converters = {}
        self._register_default_converters()
    
    def _register_default_converters(self):
        """Register default format converters"""
        # PMF profile converters
        self._registered_converters['pmf_to_gromacs'] = self._pmf_to_gromacs_xvg
        self._registered_converters['pmf_to_amber'] = self._pmf_to_amber_dat
        self._registered_converters['pmf_to_json'] = self._pmf_to_json
        self._registered_converters['pmf_to_csv'] = self._pmf_to_csv
        
        # Energy data converters
        self._registered_converters['energy_to_xvg'] = self._energy_to_xvg
        self._registered_converters['energy_to_csv'] = self._energy_to_csv
        self._registered_converters['energy_to_json'] = self._energy_to_json
    
    def convert_pmf_profile(self, data: Union[PMFProfileData, Dict[str, Any]], 
                           target_format: DataFormat, output_file: Optional[str] = None) -> Union[str, Dict]:
        """
        Convert PMF profile data to target format
        
        Args:
            data: PMF profile data
            target_format: Target format for conversion
            output_file: Optional output file path
            
        Returns:
            Converted data as string or dictionary
        """
        # Convert dict to PMFProfileData if needed
        if isinstance(data, dict):
            data = PMFProfileData(**data)
        
        converter_key = f"pmf_to_{target_format.value}"
        
        if converter_key not in self._registered_converters:
            raise DataError(f"No converter available for PMF to {target_format.value}", 
                          data_format=target_format.value)
        
        self.logger.info(f"Converting PMF profile to {target_format.value}")
        
        try:
            converted = self._registered_converters[converter_key](data)
            
            if output_file:
                self._write_output(converted, output_file, target_format)
            
            return converted
            
        except Exception as e:
            self.logger.error(f"PMF conversion failed: {e}")
            raise DataError(f"PMF conversion to {target_format.value} failed: {e}")
    
    def convert_energy_data(self, data: Union[EnergyData, Dict[str, Any]], 
                           target_format: DataFormat, output_file: Optional[str] = None) -> Union[str, Dict]:
        """
        Convert energy data to target format
        
        Args:
            data: Energy data
            target_format: Target format for conversion
            output_file: Optional output file path
            
        Returns:
            Converted data as string or dictionary
        """
        # Convert dict to EnergyData if needed
        if isinstance(data, dict):
            data = EnergyData(**data)
        
        converter_key = f"energy_to_{target_format.value}"
        
        if converter_key not in self._registered_converters:
            raise DataError(f"No converter available for energy to {target_format.value}",
                          data_format=target_format.value)
        
        self.logger.info(f"Converting energy data to {target_format.value}")
        
        try:
            converted = self._registered_converters[converter_key](data)
            
            if output_file:
                self._write_output(converted, output_file, target_format)
            
            return converted
            
        except Exception as e:
            self.logger.error(f"Energy conversion failed: {e}")
            raise DataError(f"Energy conversion to {target_format.value} failed: {e}")
    
    def register_converter(self, name: str, converter_func):
        """Register custom converter function"""
        self._registered_converters[name] = converter_func
        self.logger.info(f"Registered custom converter: {name}")
    
    # PMF Profile converters
    def _pmf_to_gromacs_xvg(self, data: PMFProfileData) -> str:
        """Convert PMF profile to GROMACS XVG format"""
        lines = [
            "# This file was created by PRISM PMF System",
            "# PMF Profile Data",
            "# Data columns:",
            "#  x: Reaction coordinate (nm)",
            "#  y: PMF energy (kJ/mol)",
        ]
        
        if data.error_estimate:
            lines.append("#  dy: Error estimate (kJ/mol)")
        
        lines.extend([
            "@    title \"PMF Profile\"",
            "@    xaxis  label \"Reaction Coordinate (nm)\"",
            "@    yaxis  label \"PMF Energy (kJ/mol)\"",
            "@TYPE xy"
        ])
        
        # Add data points
        for i, (x, energy) in enumerate(zip(data.reaction_coordinate, data.pmf_energy)):
            if data.error_estimate:
                line = f"{x:12.6f} {energy:12.6f} {data.error_estimate[i]:12.6f}"
            else:
                line = f"{x:12.6f} {energy:12.6f}"
            lines.append(line)
        
        return "\n".join(lines)
    
    def _pmf_to_amber_dat(self, data: PMFProfileData) -> str:
        """Convert PMF profile to AMBER DAT format"""
        lines = [
            "# PMF Profile Data (AMBER Format)",
            "# Generated by PRISM PMF System",
            "# Columns: RC(Angstrom) PMF(kcal/mol) Error(kcal/mol)"
        ]
        
        # Convert units: nm to Angstrom, kJ/mol to kcal/mol
        for i, (x, energy) in enumerate(zip(data.reaction_coordinate, data.pmf_energy)):
            x_angstrom = x * 10.0  # nm to Angstrom
            energy_kcal = energy / 4.184  # kJ/mol to kcal/mol
            
            if data.error_estimate:
                error_kcal = data.error_estimate[i] / 4.184
                line = f"{x_angstrom:10.4f} {energy_kcal:12.6f} {error_kcal:12.6f}"
            else:
                line = f"{x_angstrom:10.4f} {energy_kcal:12.6f}"
            
            lines.append(line)
        
        return "\n".join(lines)
    
    def _pmf_to_json(self, data: PMFProfileData) -> Dict[str, Any]:
        """Convert PMF profile to JSON format"""
        result = {
            "type": "pmf_profile",
            "generated_by": "PRISM PMF System",
            "data_points": []
        }
        
        for i, (x, energy) in enumerate(zip(data.reaction_coordinate, data.pmf_energy)):
            point = {
                "reaction_coordinate": x,
                "pmf_energy": energy
            }
            
            if data.error_estimate:
                point["error_estimate"] = data.error_estimate[i]
            
            if data.sampling_count:
                point["sampling_count"] = data.sampling_count[i]
            
            result["data_points"].append(point)
        
        if data.metadata:
            result["metadata"] = data.metadata
        
        return result
    
    def _pmf_to_csv(self, data: PMFProfileData) -> str:
        """Convert PMF profile to CSV format"""
        import io
        
        output = io.StringIO()
        
        # Determine columns
        columns = ["reaction_coordinate", "pmf_energy"]
        if data.error_estimate:
            columns.append("error_estimate")
        if data.sampling_count:
            columns.append("sampling_count")
        
        writer = csv.writer(output)
        writer.writerow(columns)
        
        # Write data rows
        for i in range(len(data.reaction_coordinate)):
            row = [data.reaction_coordinate[i], data.pmf_energy[i]]
            
            if data.error_estimate:
                row.append(data.error_estimate[i])
            if data.sampling_count:
                row.append(data.sampling_count[i])
            
            writer.writerow(row)
        
        return output.getvalue()
    
    # Energy data converters
    def _energy_to_xvg(self, data: EnergyData) -> str:
        """Convert energy data to GROMACS XVG format"""
        lines = [
            "# This file was created by PRISM PMF System",
            "# Energy Data",
            "@    title \"Energy vs Time\"",
            "@    xaxis  label \"Time (ps)\"",
            "@    yaxis  label \"Energy (kJ/mol)\"",
            "@TYPE xy"
        ]
        
        # Determine which energy terms to include
        columns = ["Time"]
        if data.potential_energy:
            columns.append("Potential")
        if data.kinetic_energy:
            columns.append("Kinetic") 
        if data.total_energy:
            columns.append("Total")
        
        lines.append(f"# Columns: {' '.join(columns)}")
        
        # Write data
        for i, time_val in enumerate(data.time):
            row = [str(time_val)]
            
            if data.potential_energy:
                row.append(str(data.potential_energy[i]))
            if data.kinetic_energy:
                row.append(str(data.kinetic_energy[i]))
            if data.total_energy:
                row.append(str(data.total_energy[i]))
            
            lines.append("  ".join(f"{val:12.6f}" for val in row))
        
        return "\n".join(lines)
    
    def _energy_to_csv(self, data: EnergyData) -> str:
        """Convert energy data to CSV format"""
        import io
        
        output = io.StringIO()
        
        # Determine columns
        columns = ["time"]
        if data.potential_energy:
            columns.append("potential_energy")
        if data.kinetic_energy:
            columns.append("kinetic_energy")
        if data.total_energy:
            columns.append("total_energy")
        if data.temperature:
            columns.append("temperature")
        if data.pressure:
            columns.append("pressure")
        
        writer = csv.writer(output)
        writer.writerow(columns)
        
        # Write data rows
        for i, time_val in enumerate(data.time):
            row = [time_val]
            
            if data.potential_energy:
                row.append(data.potential_energy[i])
            if data.kinetic_energy:
                row.append(data.kinetic_energy[i])
            if data.total_energy:
                row.append(data.total_energy[i])
            if data.temperature:
                row.append(data.temperature[i])
            if data.pressure:
                row.append(data.pressure[i])
            
            writer.writerow(row)
        
        return output.getvalue()
    
    def _energy_to_json(self, data: EnergyData) -> Dict[str, Any]:
        """Convert energy data to JSON format"""
        result = {
            "type": "energy_data", 
            "generated_by": "PRISM PMF System",
            "time": data.time
        }
        
        if data.potential_energy:
            result["potential_energy"] = data.potential_energy
        if data.kinetic_energy:
            result["kinetic_energy"] = data.kinetic_energy
        if data.total_energy:
            result["total_energy"] = data.total_energy
        if data.temperature:
            result["temperature"] = data.temperature
        if data.pressure:
            result["pressure"] = data.pressure
        if data.volume:
            result["volume"] = data.volume
        if data.additional_terms:
            result["additional_terms"] = data.additional_terms
        if data.metadata:
            result["metadata"] = data.metadata
        
        return result
    
    def _write_output(self, data: Union[str, Dict], output_file: str, format: DataFormat):
        """Write converted data to output file"""
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if isinstance(data, str):
            with open(output_path, 'w') as f:
                f.write(data)
        elif isinstance(data, dict) and format == DataFormat.JSON:
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2)
        else:
            raise DataError(f"Unsupported data type for writing: {type(data)}")


class DataExporter:
    """Exporter for PRISM calculation results"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("data_exporter")
        self.format_converter = FormatConverter(logger)
    
    def export_pmf_profile(self, pmf_data: Union[PMFProfileData, List[Dict], Dict], 
                          output_file: str, format: DataFormat = DataFormat.XVG,
                          metadata: Optional[Dict[str, Any]] = None) -> bool:
        """
        Export PMF profile to file
        
        Args:
            pmf_data: PMF profile data
            output_file: Output file path
            format: Export format
            metadata: Optional metadata to include
            
        Returns:
            True if export successful
        """
        try:
            # Convert list of dicts to PMFProfileData if needed
            if isinstance(pmf_data, list):
                pmf_data = self._list_to_pmf_profile_data(pmf_data, metadata)
            elif isinstance(pmf_data, dict):
                pmf_data = PMFProfileData(**pmf_data)
            
            # Convert and write
            self.format_converter.convert_pmf_profile(pmf_data, format, output_file)
            
            self.logger.info(f"PMF profile exported to: {output_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"PMF export failed: {e}")
            raise DataError(f"PMF export failed: {e}", file_path=output_file)
    
    def export_energy_data(self, energy_data: Union[EnergyData, Dict], 
                          output_file: str, format: DataFormat = DataFormat.XVG) -> bool:
        """
        Export energy data to file
        
        Args:
            energy_data: Energy data
            output_file: Output file path
            format: Export format
            
        Returns:
            True if export successful
        """
        try:
            # Convert dict to EnergyData if needed
            if isinstance(energy_data, dict):
                energy_data = EnergyData(**energy_data)
            
            # Convert and write
            self.format_converter.convert_energy_data(energy_data, format, output_file)
            
            self.logger.info(f"Energy data exported to: {output_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"Energy export failed: {e}")
            raise DataError(f"Energy export failed: {e}", file_path=output_file)
    
    def export_calculation_results(self, results: Dict[str, Any], output_directory: str,
                                  formats: List[DataFormat] = None) -> Dict[str, List[str]]:
        """
        Export complete calculation results in multiple formats
        
        Args:
            results: Calculation results dictionary
            output_directory: Output directory path
            formats: List of formats to export (default: XVG, JSON, CSV)
            
        Returns:
            Dictionary mapping format names to lists of exported files
        """
        if formats is None:
            formats = [DataFormat.XVG, DataFormat.JSON, DataFormat.CSV]
        
        output_dir = Path(output_directory)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        exported_files = {}
        
        try:
            # Export PMF profile if present
            if 'pmf_profile' in results:
                for format in formats:
                    extension = self._get_file_extension(format)
                    filename = f"pmf_profile.{extension}"
                    filepath = str(output_dir / filename)
                    
                    self.export_pmf_profile(results['pmf_profile'], filepath, format)
                    
                    if format.value not in exported_files:
                        exported_files[format.value] = []
                    exported_files[format.value].append(filepath)
            
            # Export energy data if present
            if 'energy_data' in results:
                for format in formats:
                    extension = self._get_file_extension(format)
                    filename = f"energy_data.{extension}"
                    filepath = str(output_dir / filename)
                    
                    self.export_energy_data(results['energy_data'], filepath, format)
                    
                    if format.value not in exported_files:
                        exported_files[format.value] = []
                    exported_files[format.value].append(filepath)
            
            # Export summary JSON
            summary_file = str(output_dir / "calculation_summary.json")
            with open(summary_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            if 'json' not in exported_files:
                exported_files['json'] = []
            exported_files['json'].append(summary_file)
            
            self.logger.info(f"Calculation results exported to: {output_directory}")
            return exported_files
            
        except Exception as e:
            self.logger.error(f"Results export failed: {e}")
            raise DataError(f"Results export failed: {e}", file_path=output_directory)
    
    def _list_to_pmf_profile_data(self, data_list: List[Dict], metadata: Optional[Dict] = None) -> PMFProfileData:
        """Convert list of dictionaries to PMFProfileData"""
        reaction_coordinate = [point.get('reaction_coordinate', 0.0) for point in data_list]
        pmf_energy = [point.get('pmf_energy', 0.0) for point in data_list]
        
        error_estimate = None
        if all('error_estimate' in point for point in data_list):
            error_estimate = [point['error_estimate'] for point in data_list]
        
        sampling_count = None
        if all('sampling_count' in point for point in data_list):
            sampling_count = [point['sampling_count'] for point in data_list]
        
        return PMFProfileData(
            reaction_coordinate=reaction_coordinate,
            pmf_energy=pmf_energy,
            error_estimate=error_estimate,
            sampling_count=sampling_count,
            metadata=metadata
        )
    
    def _get_file_extension(self, format: DataFormat) -> str:
        """Get file extension for format"""
        extension_map = {
            DataFormat.XVG: "xvg",
            DataFormat.DAT: "dat", 
            DataFormat.CSV: "csv",
            DataFormat.JSON: "json",
            DataFormat.XML: "xml"
        }
        return extension_map.get(format, "txt")


class DataImporter:
    """Importer for external data into PRISM format"""
    
    def __init__(self, logger: Optional[PrismLogger] = None):
        self.logger = logger or PrismLogger("data_importer")
    
    def import_pmf_profile(self, input_file: str, format: DataFormat = None) -> PMFProfileData:
        """
        Import PMF profile from external file
        
        Args:
            input_file: Input file path
            format: File format (auto-detected if None)
            
        Returns:
            PMFProfileData object
        """
        input_path = Path(input_file)
        
        if not input_path.exists():
            raise DataError(f"Input file not found: {input_file}", file_path=input_file)
        
        # Auto-detect format if not specified
        if format is None:
            format = self._detect_format(input_path)
        
        self.logger.info(f"Importing PMF profile from {input_file} ({format.value} format)")
        
        try:
            if format == DataFormat.XVG:
                return self._import_xvg_pmf(input_path)
            elif format == DataFormat.DAT:
                return self._import_dat_pmf(input_path)
            elif format == DataFormat.CSV:
                return self._import_csv_pmf(input_path)
            elif format == DataFormat.JSON:
                return self._import_json_pmf(input_path)
            else:
                raise DataError(f"Unsupported import format: {format.value}", data_format=format.value)
                
        except Exception as e:
            self.logger.error(f"PMF import failed: {e}")
            raise DataError(f"PMF import from {format.value} failed: {e}", file_path=input_file)
    
    def import_energy_data(self, input_file: str, format: DataFormat = None) -> EnergyData:
        """
        Import energy data from external file
        
        Args:
            input_file: Input file path
            format: File format (auto-detected if None)
            
        Returns:
            EnergyData object
        """
        input_path = Path(input_file)
        
        if not input_path.exists():
            raise DataError(f"Input file not found: {input_file}", file_path=input_file)
        
        # Auto-detect format if not specified
        if format is None:
            format = self._detect_format(input_path)
        
        self.logger.info(f"Importing energy data from {input_file} ({format.value} format)")
        
        try:
            if format == DataFormat.XVG:
                return self._import_xvg_energy(input_path)
            elif format == DataFormat.CSV:
                return self._import_csv_energy(input_path)
            elif format == DataFormat.JSON:
                return self._import_json_energy(input_path)
            else:
                raise DataError(f"Unsupported import format: {format.value}", data_format=format.value)
                
        except Exception as e:
            self.logger.error(f"Energy import failed: {e}")
            raise DataError(f"Energy import from {format.value} failed: {e}", file_path=input_file)
    
    def _detect_format(self, file_path: Path) -> DataFormat:
        """Auto-detect file format from extension and content"""
        extension = file_path.suffix.lower()
        
        format_map = {
            '.xvg': DataFormat.XVG,
            '.dat': DataFormat.DAT,
            '.csv': DataFormat.CSV,
            '.json': DataFormat.JSON,
            '.xml': DataFormat.XML
        }
        
        if extension in format_map:
            return format_map[extension]
        
        # Try to detect from content
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                
                if first_line.startswith('@') or first_line.startswith('#'):
                    return DataFormat.XVG
                elif first_line.startswith('{'):
                    return DataFormat.JSON
                elif ',' in first_line:
                    return DataFormat.CSV
                else:
                    return DataFormat.DAT
        except Exception:
            return DataFormat.DAT  # Default fallback
    
    def _import_xvg_pmf(self, file_path: Path) -> PMFProfileData:
        """Import PMF from XVG format"""
        reaction_coordinate = []
        pmf_energy = []
        error_estimate = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and headers
                if line.startswith('#') or line.startswith('@') or not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        reaction_coordinate.append(float(parts[0]))
                        pmf_energy.append(float(parts[1]))
                        
                        if len(parts) >= 3:
                            error_estimate.append(float(parts[2]))
                    except ValueError:
                        continue
        
        return PMFProfileData(
            reaction_coordinate=reaction_coordinate,
            pmf_energy=pmf_energy,
            error_estimate=error_estimate if error_estimate else None
        )
    
    def _import_dat_pmf(self, file_path: Path) -> PMFProfileData:
        """Import PMF from DAT format"""
        # Similar to XVG but may have different units
        return self._import_xvg_pmf(file_path)  # Simplified
    
    def _import_csv_pmf(self, file_path: Path) -> PMFProfileData:
        """Import PMF from CSV format"""
        reaction_coordinate = []
        pmf_energy = []
        error_estimate = []
        sampling_count = []
        
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            
            for row in reader:
                reaction_coordinate.append(float(row.get('reaction_coordinate', 0)))
                pmf_energy.append(float(row.get('pmf_energy', 0)))
                
                if 'error_estimate' in row:
                    error_estimate.append(float(row['error_estimate']))
                
                if 'sampling_count' in row:
                    sampling_count.append(int(row['sampling_count']))
        
        return PMFProfileData(
            reaction_coordinate=reaction_coordinate,
            pmf_energy=pmf_energy,
            error_estimate=error_estimate if error_estimate else None,
            sampling_count=sampling_count if sampling_count else None
        )
    
    def _import_json_pmf(self, file_path: Path) -> PMFProfileData:
        """Import PMF from JSON format"""
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        if 'data_points' in data:
            points = data['data_points']
            reaction_coordinate = [p['reaction_coordinate'] for p in points]
            pmf_energy = [p['pmf_energy'] for p in points]
            
            error_estimate = None
            if all('error_estimate' in p for p in points):
                error_estimate = [p['error_estimate'] for p in points]
            
            sampling_count = None
            if all('sampling_count' in p for p in points):
                sampling_count = [p['sampling_count'] for p in points]
            
            return PMFProfileData(
                reaction_coordinate=reaction_coordinate,
                pmf_energy=pmf_energy,
                error_estimate=error_estimate,
                sampling_count=sampling_count,
                metadata=data.get('metadata')
            )
        else:
            # Direct format
            return PMFProfileData(**data)
    
    def _import_xvg_energy(self, file_path: Path) -> EnergyData:
        """Import energy data from XVG format"""
        time = []
        potential_energy = []
        kinetic_energy = []
        total_energy = []
        
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and headers
                if line.startswith('#') or line.startswith('@') or not line:
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        time.append(float(parts[0]))
                        
                        # Assume columns: Time, Potential, Kinetic, Total
                        if len(parts) > 1:
                            potential_energy.append(float(parts[1]))
                        if len(parts) > 2:
                            kinetic_energy.append(float(parts[2]))
                        if len(parts) > 3:
                            total_energy.append(float(parts[3]))
                    except ValueError:
                        continue
        
        return EnergyData(
            time=time,
            potential_energy=potential_energy if potential_energy else None,
            kinetic_energy=kinetic_energy if kinetic_energy else None,
            total_energy=total_energy if total_energy else None
        )
    
    def _import_csv_energy(self, file_path: Path) -> EnergyData:
        """Import energy data from CSV format"""
        with open(file_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            
            time = []
            potential_energy = []
            kinetic_energy = []
            total_energy = []
            temperature = []
            pressure = []
            
            for row in reader:
                time.append(float(row.get('time', 0)))
                
                if 'potential_energy' in row:
                    potential_energy.append(float(row['potential_energy']))
                if 'kinetic_energy' in row:
                    kinetic_energy.append(float(row['kinetic_energy']))
                if 'total_energy' in row:
                    total_energy.append(float(row['total_energy']))
                if 'temperature' in row:
                    temperature.append(float(row['temperature']))
                if 'pressure' in row:
                    pressure.append(float(row['pressure']))
        
        return EnergyData(
            time=time,
            potential_energy=potential_energy if potential_energy else None,
            kinetic_energy=kinetic_energy if kinetic_energy else None,
            total_energy=total_energy if total_energy else None,
            temperature=temperature if temperature else None,
            pressure=pressure if pressure else None
        )
    
    def _import_json_energy(self, file_path: Path) -> EnergyData:
        """Import energy data from JSON format"""
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        return EnergyData(**data)