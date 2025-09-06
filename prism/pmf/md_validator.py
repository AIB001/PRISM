#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MD Results Validator - Enhanced validation for MD simulation results

This module provides comprehensive validation of MD simulation results,
going beyond simple file existence to check file integrity, completeness,
and compatibility with PMF workflows.
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class ValidationLevel(Enum):
    """Validation levels for MD results"""
    BASIC = "basic"          # File existence only
    STANDARD = "standard"    # File existence + basic integrity
    COMPREHENSIVE = "comprehensive"  # Full validation + recommendations


class FileType(Enum):
    """File types in MD results"""
    STRUCTURE = "structure"    # .gro files
    TOPOLOGY = "topology"      # .top, .itp files  
    TRAJECTORY = "trajectory"  # .xtc, .trr files
    ENERGY = "energy"         # .edr files
    LOG = "log"               # .log files
    RESTART = "restart"       # .cpt files
    PARAMETER = "parameter"   # .mdp files
    INDEX = "index"           # .ndx files


@dataclass
class FileInfo:
    """Information about a file in MD results"""
    path: Path
    file_type: FileType
    size: int
    required: bool = True
    description: str = ""
    
    @property
    def exists(self) -> bool:
        return self.path.exists()
    
    @property
    def size_mb(self) -> float:
        return self.size / (1024 * 1024) if self.size > 0 else 0.0


@dataclass 
class ValidationResult:
    """Result of MD validation"""
    is_valid: bool
    level: ValidationLevel
    score: float  # 0-100
    files_found: List[FileInfo]
    files_missing: List[str]
    warnings: List[str]
    errors: List[str]
    recommendations: List[str]
    summary: str
    
    @property
    def grade(self) -> str:
        """Get validation grade"""
        if self.score >= 90: return "A"
        elif self.score >= 80: return "B" 
        elif self.score >= 70: return "C"
        elif self.score >= 60: return "D"
        else: return "F"


class MDResultsValidator:
    """
    Enhanced MD Results Validator
    
    Provides comprehensive validation of MD simulation results with
    intelligent detection of different PRISM output formats.
    """
    
    def __init__(self, validation_level: ValidationLevel = ValidationLevel.STANDARD):
        self.validation_level = validation_level
        self.reset()
    
    def reset(self):
        """Reset validator state"""
        self.files_found = []
        self.files_missing = []
        self.warnings = []
        self.errors = []
        self.recommendations = []
    
    def validate(self, system_dir: Union[str, Path]) -> ValidationResult:
        """
        Validate MD results directory
        
        Parameters:
        -----------
        system_dir : str or Path
            Directory containing MD results
            
        Returns:
        --------
        ValidationResult : Comprehensive validation result
        """
        self.reset()
        system_path = Path(system_dir).resolve()
        
        logger.info(f"Validating MD results: {system_path}")
        logger.info(f"Validation level: {self.validation_level.value}")
        
        if not system_path.exists():
            return self._create_failed_result("System directory does not exist")
        
        # Locate MD results directory
        md_dir = self._locate_md_directory(system_path)
        if not md_dir:
            return self._create_failed_result("No valid MD results directory found")
        
        # Validate based on level
        if self.validation_level == ValidationLevel.BASIC:
            result = self._validate_basic(md_dir)
        elif self.validation_level == ValidationLevel.STANDARD:
            result = self._validate_standard(md_dir)
        else:  # COMPREHENSIVE
            result = self._validate_comprehensive(md_dir)
        
        # Add final recommendations
        self._add_recommendations(md_dir, result)
        
        return result
    
    def _locate_md_directory(self, system_path: Path) -> Optional[Path]:
        """Intelligently locate MD results directory"""
        
        # Strategy 1: Check if system_path itself is a GMX directory
        if self._is_gmx_directory(system_path):
            logger.debug(f"Found GMX directory: {system_path}")
            return system_path
        
        # Strategy 2: Look for standard GMX_PROLIG_MD subdirectory
        gmx_dir = system_path / "GMX_PROLIG_MD"
        if self._is_gmx_directory(gmx_dir):
            logger.debug(f"Found GMX directory: {gmx_dir}")
            return gmx_dir
        
        # Strategy 3: Search all subdirectories
        for subdir in system_path.iterdir():
            if subdir.is_dir():
                # Check for GMX_PROLIG_MD in subdirectories
                candidate = subdir / "GMX_PROLIG_MD"
                if self._is_gmx_directory(candidate):
                    logger.debug(f"Found GMX directory: {candidate}")
                    return candidate
                
                # Check if subdirectory itself is GMX directory
                if self._is_gmx_directory(subdir):
                    logger.debug(f"Found GMX directory: {subdir}")
                    return subdir
        
        # Strategy 4: Look for any directory with GMX-like structure
        for subdir in system_path.rglob("*"):
            if subdir.is_dir() and self._looks_like_gmx_directory(subdir):
                logger.debug(f"Found potential GMX directory: {subdir}")
                self.warnings.append(f"Non-standard GMX directory structure: {subdir}")
                return subdir
        
        return None
    
    def _is_gmx_directory(self, path: Path) -> bool:
        """Check if directory contains valid GROMACS results"""
        if not path.exists():
            return False
        
        # Check for essential GROMACS files
        essential_files = [
            "topol.top",
            # At least one structure file
        ]
        
        structure_files = [
            "solv_ions.gro",
            "prod/md.gro", 
            "npt/npt.gro",
            "nvt/nvt.gro",
            "em/em.gro"
        ]
        
        # Check essential files
        for file in essential_files:
            if not (path / file).exists():
                return False
        
        # Check at least one structure file exists
        has_structure = any((path / sf).exists() for sf in structure_files)
        if not has_structure:
            return False
        
        return True
    
    def _looks_like_gmx_directory(self, path: Path) -> bool:
        """Less strict check for GROMACS-like directories"""
        if not path.exists():
            return False
        
        gmx_indicators = [
            "topol.top", "*.gro", "*.mdp", "*.edr", "*.log", "*.xtc", "*.trr"
        ]
        
        indicators_found = 0
        for pattern in gmx_indicators:
            if list(path.glob(pattern)) or list(path.glob(f"**/{pattern}")):
                indicators_found += 1
        
        return indicators_found >= 3  # At least 3 indicators
    
    def _validate_basic(self, md_dir: Path) -> ValidationResult:
        """Basic validation - file existence only"""
        logger.debug("Running basic validation")
        
        # Define essential files for basic validation
        essential_files = {
            "topol.top": (FileType.TOPOLOGY, True, "Main topology file"),
        }
        
        # Structure files (at least one required)
        structure_files = {
            "solv_ions.gro": (FileType.STRUCTURE, False, "Solvated system structure"),
            "prod/md.gro": (FileType.STRUCTURE, False, "Production MD final structure"),
            "npt/npt.gro": (FileType.STRUCTURE, False, "NPT equilibration structure"),
        }
        
        # Check essential files
        for filename, (file_type, required, desc) in essential_files.items():
            self._check_file(md_dir / filename, file_type, required, desc)
        
        # Check structure files (at least one required)
        structure_found = False
        for filename, (file_type, required, desc) in structure_files.items():
            if (md_dir / filename).exists():
                self._check_file(md_dir / filename, file_type, True, desc)
                structure_found = True
        
        if not structure_found:
            self.errors.append("No structure files found")
            self.files_missing.extend(list(structure_files.keys()))
        
        # Calculate score
        score = self._calculate_basic_score()
        is_valid = score >= 60  # 60% minimum for basic validation
        
        return ValidationResult(
            is_valid=is_valid,
            level=ValidationLevel.BASIC,
            score=score,
            files_found=self.files_found,
            files_missing=self.files_missing,
            warnings=self.warnings,
            errors=self.errors,
            recommendations=self.recommendations,
            summary=f"Basic validation {'PASSED' if is_valid else 'FAILED'} (Score: {score:.1f}%)"
        )
    
    def _validate_standard(self, md_dir: Path) -> ValidationResult:
        """Standard validation - file existence + basic integrity"""
        logger.debug("Running standard validation")
        
        # Start with basic validation
        basic_result = self._validate_basic(md_dir)
        
        # Additional files for standard validation
        additional_files = {
            "em/em.gro": (FileType.STRUCTURE, False, "Energy minimized structure"),
            "nvt/nvt.gro": (FileType.STRUCTURE, False, "NVT equilibrated structure"),
            "npt/npt.gro": (FileType.STRUCTURE, False, "NPT equilibrated structure"),
            "posre.itp": (FileType.TOPOLOGY, False, "Position restraints"),
            "prod/md.xtc": (FileType.TRAJECTORY, False, "Production trajectory"),
            "prod/md.edr": (FileType.ENERGY, False, "Production energy file"),
            "prod/md.log": (FileType.LOG, False, "Production log file"),
        }
        
        for filename, (file_type, required, desc) in additional_files.items():
            self._check_file(md_dir / filename, file_type, required, desc)
        
        # Check file integrity
        self._check_file_integrity(md_dir)
        
        # Check simulation completeness
        self._check_simulation_completeness(md_dir)
        
        # Calculate enhanced score
        score = self._calculate_standard_score()
        is_valid = score >= 70  # 70% minimum for standard validation
        
        return ValidationResult(
            is_valid=is_valid,
            level=ValidationLevel.STANDARD,
            score=score,
            files_found=self.files_found,
            files_missing=self.files_missing,
            warnings=self.warnings,
            errors=self.errors,
            recommendations=self.recommendations,
            summary=f"Standard validation {'PASSED' if is_valid else 'FAILED'} (Score: {score:.1f}%)"
        )
    
    def _validate_comprehensive(self, md_dir: Path) -> ValidationResult:
        """Comprehensive validation - full analysis"""
        logger.debug("Running comprehensive validation")
        
        # Start with standard validation
        standard_result = self._validate_standard(md_dir)
        
        # Additional comprehensive checks
        self._check_trajectory_quality(md_dir)
        self._check_energy_conservation(md_dir)
        self._check_ligand_topology(md_dir)
        self._check_force_field_consistency(md_dir)
        self._analyze_system_composition(md_dir)
        
        # Calculate comprehensive score
        score = self._calculate_comprehensive_score()
        is_valid = score >= 80  # 80% minimum for comprehensive validation
        
        return ValidationResult(
            is_valid=is_valid,
            level=ValidationLevel.COMPREHENSIVE,
            score=score,
            files_found=self.files_found,
            files_missing=self.files_missing,
            warnings=self.warnings,
            errors=self.errors,
            recommendations=self.recommendations,
            summary=f"Comprehensive validation {'PASSED' if is_valid else 'FAILED'} (Score: {score:.1f}%, Grade: {self._get_grade(score)})"
        )
    
    def _check_file(self, file_path: Path, file_type: FileType, required: bool, description: str):
        """Check individual file"""
        if file_path.exists():
            size = file_path.stat().st_size
            file_info = FileInfo(
                path=file_path,
                file_type=file_type,
                size=size,
                required=required,
                description=description
            )
            self.files_found.append(file_info)
            logger.debug(f"Found {description}: {file_path} ({file_info.size_mb:.1f} MB)")
        else:
            if required:
                self.errors.append(f"Missing required file: {file_path}")
            else:
                self.warnings.append(f"Missing optional file: {file_path}")
            self.files_missing.append(str(file_path))
    
    def _check_file_integrity(self, md_dir: Path):
        """Check basic file integrity"""
        for file_info in self.files_found:
            if file_info.size == 0:
                self.warnings.append(f"Empty file: {file_info.path}")
            elif file_info.file_type == FileType.STRUCTURE and file_info.size < 1000:
                self.warnings.append(f"Unusually small structure file: {file_info.path}")
            elif file_info.file_type == FileType.TRAJECTORY and file_info.size < 10000:
                self.warnings.append(f"Unusually small trajectory file: {file_info.path}")
    
    def _check_simulation_completeness(self, md_dir: Path):
        """Check if simulation appears complete"""
        # Check for complete equilibration chain
        equilibration_files = ["em/em.gro", "nvt/nvt.gro", "npt/npt.gro"]
        completed_steps = sum(1 for f in equilibration_files if (md_dir / f).exists())
        
        if completed_steps == 0:
            self.errors.append("No equilibration steps completed")
        elif completed_steps < 3:
            self.warnings.append(f"Incomplete equilibration chain ({completed_steps}/3 steps)")
        
        # Check production run
        prod_files = ["prod/md.gro", "prod/md.xtc", "prod/md.edr"]
        prod_completed = sum(1 for f in prod_files if (md_dir / f).exists())
        
        if prod_completed == 0:
            self.warnings.append("No production simulation files found")
        elif prod_completed < 3:
            self.warnings.append(f"Incomplete production simulation ({prod_completed}/3 files)")
    
    def _check_trajectory_quality(self, md_dir: Path):
        """Check trajectory file quality (comprehensive only)"""
        trajectory_files = [f for f in self.files_found if f.file_type == FileType.TRAJECTORY]
        
        if not trajectory_files:
            self.warnings.append("No trajectory files found - cannot assess simulation quality")
        else:
            for traj in trajectory_files:
                if traj.size_mb < 10:
                    self.warnings.append(f"Very short trajectory: {traj.path} ({traj.size_mb:.1f} MB)")
                elif traj.size_mb > 10000:  # 10 GB
                    self.warnings.append(f"Very large trajectory: {traj.path} ({traj.size_mb:.1f} MB)")
    
    def _check_energy_conservation(self, md_dir: Path):
        """Check energy file presence (comprehensive only)"""
        energy_files = [f for f in self.files_found if f.file_type == FileType.ENERGY]
        
        if not energy_files:
            self.warnings.append("No energy files found - cannot assess energy conservation")
    
    def _check_ligand_topology(self, md_dir: Path):
        """Check ligand topology files (comprehensive only)"""
        parent_dir = md_dir.parent
        ligand_dirs = ["LIG.amb2gmx", "LIG.openff2gmx"]
        
        ligand_found = False
        for lig_dir in ligand_dirs:
            lig_path = parent_dir / lig_dir
            if lig_path.exists():
                ligand_found = True
                # Check essential ligand files
                lig_files = ["LIG.itp", "LIG.gro", "atomtypes_LIG.itp"]
                missing_lig_files = [f for f in lig_files if not (lig_path / f).exists()]
                
                if missing_lig_files:
                    self.warnings.append(f"Missing ligand files in {lig_dir}: {missing_lig_files}")
                else:
                    logger.debug(f"Complete ligand topology found: {lig_dir}")
        
        if not ligand_found:
            self.warnings.append("No ligand topology directory found")
    
    def _check_force_field_consistency(self, md_dir: Path):
        """Check force field consistency (comprehensive only)"""
        # This could be expanded to actually parse topology files
        self.recommendations.append("Consider validating force field consistency in topology files")
    
    def _analyze_system_composition(self, md_dir: Path):
        """Analyze system composition (comprehensive only)"""
        # This could be expanded to parse structure files and analyze composition
        structure_files = [f for f in self.files_found if f.file_type == FileType.STRUCTURE]
        
        if structure_files:
            largest_structure = max(structure_files, key=lambda x: x.size)
            self.recommendations.append(f"System size analysis available from: {largest_structure.path}")
    
    def _calculate_basic_score(self) -> float:
        """Calculate basic validation score"""
        total_weight = 0
        achieved_weight = 0
        
        # Essential files weight more
        for file_info in self.files_found:
            weight = 50 if file_info.required else 10
            total_weight += weight
            achieved_weight += weight
        
        # Penalty for missing files
        missing_penalty = len(self.files_missing) * 20
        errors_penalty = len(self.errors) * 30
        
        max_possible = 100
        penalty = min(missing_penalty + errors_penalty, max_possible)
        
        return max(0, max_possible - penalty)
    
    def _calculate_standard_score(self) -> float:
        """Calculate standard validation score"""
        base_score = self._calculate_basic_score()
        
        # Bonus for completeness
        structure_bonus = len([f for f in self.files_found if f.file_type == FileType.STRUCTURE]) * 3
        trajectory_bonus = len([f for f in self.files_found if f.file_type == FileType.TRAJECTORY]) * 5
        energy_bonus = len([f for f in self.files_found if f.file_type == FileType.ENERGY]) * 5
        
        # Warning penalty
        warning_penalty = len(self.warnings) * 2
        
        total_score = base_score + structure_bonus + trajectory_bonus + energy_bonus - warning_penalty
        
        return min(100, max(0, total_score))
    
    def _calculate_comprehensive_score(self) -> float:
        """Calculate comprehensive validation score"""
        standard_score = self._calculate_standard_score()
        
        # Additional bonuses for comprehensive analysis
        quality_bonus = 0
        if any(f.file_type == FileType.TRAJECTORY for f in self.files_found):
            quality_bonus += 5
        if any(f.file_type == FileType.ENERGY for f in self.files_found):
            quality_bonus += 5
        
        # Reduce penalty for warnings in comprehensive mode (they're informational)
        warning_adjustment = len(self.warnings) * 1  # Smaller penalty
        
        total_score = standard_score + quality_bonus - warning_adjustment
        
        return min(100, max(0, total_score))
    
    def _get_grade(self, score: float) -> str:
        """Convert score to letter grade"""
        if score >= 90: return "A"
        elif score >= 80: return "B"
        elif score >= 70: return "C" 
        elif score >= 60: return "D"
        else: return "F"
    
    def _add_recommendations(self, md_dir: Path, result: ValidationResult):
        """Add context-specific recommendations"""
        if result.score < 60:
            self.recommendations.insert(0, "System validation failed - consider rebuilding MD simulation")
        elif result.score < 80:
            self.recommendations.insert(0, "System validation marginal - review warnings before proceeding")
        
        # PMF-specific recommendations
        if not any("trajectory" in f.path.name.lower() for f in self.files_found):
            self.recommendations.append("No trajectory files found - PMF analysis may be limited")
        
        if not any("prod" in str(f.path) for f in self.files_found):
            self.recommendations.append("No production simulation files - consider running production MD first")
    
    def _create_failed_result(self, error_message: str) -> ValidationResult:
        """Create a failed validation result"""
        self.errors.append(error_message)
        
        return ValidationResult(
            is_valid=False,
            level=self.validation_level,
            score=0.0,
            files_found=[],
            files_missing=[],
            warnings=self.warnings,
            errors=self.errors,
            recommendations=["Fix critical errors before proceeding"],
            summary=f"Validation FAILED: {error_message}"
        )


# Convenience functions
def validate_md_results(system_dir: Union[str, Path], 
                       level: ValidationLevel = ValidationLevel.STANDARD) -> ValidationResult:
    """
    Convenience function to validate MD results
    
    Parameters:
    -----------
    system_dir : str or Path
        Directory containing MD results
    level : ValidationLevel
        Level of validation to perform
        
    Returns:
    --------
    ValidationResult : Comprehensive validation result
    """
    validator = MDResultsValidator(level)
    return validator.validate(system_dir)


def quick_validate(system_dir: Union[str, Path]) -> bool:
    """
    Quick validation - returns True if system is valid for PMF
    
    Parameters:
    -----------
    system_dir : str or Path
        Directory containing MD results
        
    Returns:
    --------
    bool : True if valid, False otherwise
    """
    result = validate_md_results(system_dir, ValidationLevel.BASIC)
    return result.is_valid


def get_validation_report(system_dir: Union[str, Path], 
                         level: ValidationLevel = ValidationLevel.STANDARD) -> str:
    """
    Get a formatted validation report
    
    Parameters:
    -----------
    system_dir : str or Path
        Directory containing MD results
    level : ValidationLevel
        Level of validation to perform
        
    Returns:
    --------
    str : Formatted validation report
    """
    result = validate_md_results(system_dir, level)
    
    report = []
    report.append(f"MD Results Validation Report")
    report.append(f"{'=' * 40}")
    report.append(f"Directory: {system_dir}")
    report.append(f"Validation Level: {result.level.value}")
    report.append(f"Score: {result.score:.1f}% (Grade: {result.grade})")
    report.append(f"Status: {'VALID' if result.is_valid else 'INVALID'}")
    report.append("")
    
    if result.files_found:
        report.append("Files Found:")
        for file_info in result.files_found:
            report.append(f"  {file_info.description}: {file_info.path.name} ({file_info.size_mb:.1f} MB)")
        report.append("")
    
    if result.files_missing:
        report.append("Missing Files:")
        for file in result.files_missing:
            report.append(f"  ✗ {file}")
        report.append("")
    
    if result.errors:
        report.append("Errors:")
        for error in result.errors:
            report.append(f"  ERROR: {error}")
        report.append("")
    
    if result.warnings:
        report.append("Warnings:")
        for warning in result.warnings:
            report.append(f"  WARNING: {warning}")
        report.append("")
    
    if result.recommendations:
        report.append("Recommendations:")
        for rec in result.recommendations:
            report.append(f"  💡 {rec}")
    
    return "\n".join(report)