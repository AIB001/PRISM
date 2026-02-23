#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAMD FEP analysis for calculating free energy perturbations.

This module provides tools for analyzing NAMD FEP output files (.fepout),
including energy decomposition, convergence analysis, and binding free energy calculations.
"""

import os
import subprocess
import logging
import pickle
import glob
from pathlib import Path
from typing import List, Dict, Optional, Union
import pandas as pd
import numpy as np

from ..core.config import AnalysisConfig
from ..resources import get_namd_fep_script

logger = logging.getLogger(__name__)


class NAMDFEPAnalyzer:
    """
    Analyzer for NAMD Free Energy Perturbation calculations.

    This class provides methods for processing NAMD .fepout files,
    calculating energy decompositions, and analyzing convergence.
    """

    def __init__(self, config: AnalysisConfig):
        """
        Initialize NAMD FEP analyzer.

        Parameters
        ----------
        config : AnalysisConfig
            Analysis configuration object
        """
        self.config = config
        self._cache_dir = Path(config.cache_dir).resolve()  # Make absolute
        self._cache_dir.mkdir(parents=True, exist_ok=True)

        # Get path to shell script
        self._script_path = get_namd_fep_script()
        if not Path(self._script_path).exists():
            raise FileNotFoundError(f"NAMD FEP decomposition script not found: {self._script_path}")
        logger.info(f"NAMD FEP analyzer initialized with script: {self._script_path}")

    def parse_fepout_summary(self, fepout_file: str) -> Dict[str, float]:
        """
        Parse summary information from a NAMD .fepout file using fast grep.

        Parameters
        ----------
        fepout_file : str
            Path to .fepout file

        Returns
        -------
        dict
            Dictionary with 'lambda2', 'dG', and 'net_change' keys
        """
        if not Path(fepout_file).exists():
            raise FileNotFoundError(f"FEP output file not found: {fepout_file}")

        try:
            # Use grep to extract summary lines (much faster than Python parsing)
            cmd = ["grep", "^#Free energy change", fepout_file]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            lines = result.stdout.strip().split("\n")
            if not lines:
                raise ValueError(f"No FEP summary lines found in {fepout_file}")

            # Parse last line (final result)
            last_line = lines[-1].split()
            summary = {"lambda2": float(last_line[8]), "dG": float(last_line[11]), "net_change": float(last_line[18])}

            logger.info(f"Parsed {fepout_file}: ΔG = {summary['net_change']:.2f} kcal/mol")
            return summary

        except subprocess.CalledProcessError as e:
            logger.error(f"Error parsing {fepout_file}: {e}")
            raise
        except Exception as e:
            logger.error(f"Error parsing FEP output: {e}")
            raise

    def run_decomposition(
        self,
        fepout_files: Union[str, List[str]],
        start_frame: int = 10000,
        end_frame: int = 510000,
        stride: int = 5,
        num_points: int = 10,
        output_dir: Optional[str] = None,
        cache_name: Optional[str] = None,
        force_recompute: bool = False,
    ) -> pd.DataFrame:
        """
        Run NAMD FEP decomposition analysis using optimized bash script.

        This method calls the mknamd_fep_decomp_convergence.sh script which uses
        grep/awk for fast parsing of large .fepout files and calculates energy
        decomposition at multiple time points for convergence analysis.

        Parameters
        ----------
        fepout_files : str or list of str
            Path(s) to .fepout file(s). Can be a single file, list of files,
            or glob pattern (e.g., "*.fepout")
        start_frame : int
            Starting frame for analysis (alchEquilSteps)
        end_frame : int
            Ending frame for analysis (numSteps)
        stride : int
            Output frequency (alchOutFreq)
        num_points : int
            Number of convergence points to calculate
        output_dir : str, optional
            Directory for intermediate output. If None, uses 'outdecomp'.
        cache_name : str, optional
            Custom cache name. If None, generates from parameters.
        force_recompute : bool
            If True, ignore cache and recompute

        Returns
        -------
        pd.DataFrame
            Decomposition data with columns:
            - fraction: Simulation fraction (0 to 1)
            - sum_dg: Total free energy change
            - sum_delec: Electrostatic contribution
            - sum_dvdw: van der Waals contribution
            - sum_couple: Coupling term
        """
        # Handle file paths
        if isinstance(fepout_files, str):
            if "*" in fepout_files or "?" in fepout_files:
                # Glob pattern
                fepout_files = glob.glob(fepout_files)
                if not fepout_files:
                    raise ValueError(f"No files found matching pattern")
            else:
                fepout_files = [fepout_files]

        # Validate files exist
        for f in fepout_files:
            if not Path(f).exists():
                raise FileNotFoundError(f"FEP output file not found: {f}")

        # Create cache key
        if cache_name is None:
            file_names = [Path(f).name for f in fepout_files]
            cache_key = f"namd_fep_decomp_{'_'.join(file_names[:3])}_{start_frame}_{end_frame}_{stride}_{num_points}"
            cache_key = cache_key.replace(".", "_").replace("-", "_")[:200]  # Limit length
        else:
            cache_key = cache_name

        cache_file = self._cache_dir / f"{cache_key}.pkl"

        # Check cache
        if not force_recompute and cache_file.exists():
            with open(cache_file, "rb") as f:
                return pickle.load(f)

        # Determine output directory
        if output_dir is None:
            output_dir = "outdecomp"

        # Store original directory
        original_dir = os.getcwd()
        work_dir = Path(fepout_files[0]).parent

        try:
            # Change to working directory
            os.chdir(work_dir)
            logger.info(f"Working directory: {work_dir}")

            # Prepare command - use relative paths in work directory
            fepout_pattern = " ".join([Path(f).name for f in fepout_files])
            cmd = f"bash {self._script_path} {fepout_pattern} {start_frame} {end_frame} {stride} {num_points}"

            # Run bash script
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=600,  # 10 minute timeout
            )

            if result.returncode != 0:
                logger.error(f"Decomposition script failed: {result.stderr}")
                raise RuntimeError(f"NAMD FEP decomposition failed:\n{result.stderr}")

            # Parse results from decompose_summary.dat
            # Note: we're already in work_dir, so use relative path
            summary_file = Path("decompose_summary.dat")
            if not summary_file.exists():
                # Try absolute path as fallback
                summary_file = work_dir / "decompose_summary.dat"
                if not summary_file.exists():
                    raise FileNotFoundError(
                        f"Expected output file not found: {summary_file}\n" "The decomposition script may have failed."
                    )

            # Read decomposition data
            data = pd.read_csv(summary_file, sep=",", skipinitialspace=True)
            data.columns = [c.strip() for c in data.columns]

            logger.info(f"Loaded decomposition data: {len(data)} time points")
            logger.info(f"Final ΔG = {data['sum_dg'].iloc[-1]:.2f} kcal/mol")

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(data, f)

            logger.info(f"Decomposition cached to {cache_file.name}")
            return data

        except subprocess.TimeoutExpired:
            logger.error("Decomposition script timed out after 10 minutes")
            raise TimeoutError("NAMD FEP decomposition timed out")
        except Exception as e:
            logger.error(f"Error running decomposition: {e}")
            raise
        finally:
            # Return to original directory
            os.chdir(original_dir)

    def calculate_system_decomposition(
        self,
        system_dir: str,
        system_type: str = "auto",
        n_repeats: int = 4,
        parallel: bool = True,
        cache_name: Optional[str] = None,
        **decomp_kwargs,
    ) -> Dict[str, np.ndarray]:
        """
        Calculate FEP decomposition for complex and/or ligand systems with multiple repeats.

        This method automatically detects and processes multiple repeat directories
        (e.g., complex1, complex2, ..., ligand1, ligand2, ...) and averages the results.

        Parameters
        ----------
        system_dir : str
            Path to directory containing FEP calculation repeats
        system_type : str
            Type of system to process: 'complex', 'ligand', or 'auto'
            If 'auto', automatically detects both complex and ligand directories
        n_repeats : int
            Expected number of repeat calculations (used for validation)
        parallel : bool
            Whether to process repeats in parallel
        cache_name : str, optional
            Custom cache name for the combined system results
        **decomp_kwargs
            Additional arguments passed to run_decomposition()

        Returns
        -------
        dict
            Dictionary with keys:
            - 'complex': np.array([dG, dElec, dVdW, dCouple]) if complex found
            - 'ligand': np.array([dG, dElec, dVdW, dCouple]) if ligand found
            - 'ddG': difference (complex - ligand) if both found
            Each array contains the final decomposition values averaged across repeats.
        """
        system_dir = Path(system_dir).resolve()
        if not system_dir.exists():
            raise FileNotFoundError(f"System directory not found: {system_dir}")

        logger.info(f"Processing FEP system: {system_dir}")
        if parallel:
            logger.info("Parallel processing requested but not implemented; running serially.")

        # Determine system types to process
        if system_type == "auto":
            system_types = []
            if any(d.name.startswith("complex") for d in system_dir.iterdir() if d.is_dir()):
                system_types.append("complex")
            if any(d.name.startswith("ligand") for d in system_dir.iterdir() if d.is_dir()):
                system_types.append("ligand")
            if not system_types:
                logger.warning("No complex or ligand directories found")
                return {}
        else:
            system_types = [system_type]

        results = {}

        for sys_type in system_types:
            logger.info(f"Processing {sys_type} phase...")

            # Find repeat directories
            repeat_dirs = sorted([d for d in system_dir.iterdir() if d.is_dir() and d.name.startswith(sys_type)])

            if not repeat_dirs:
                logger.warning(f"No {sys_type} directories found in {system_dir}")
                continue

            if n_repeats and len(repeat_dirs) != n_repeats:
                logger.warning(
                    "Expected %s repeat(s) for %s, but found %s.",
                    n_repeats,
                    sys_type,
                    len(repeat_dirs),
                )

            logger.info(f"Found {len(repeat_dirs)} {sys_type} repeat(s): {[d.name for d in repeat_dirs]}")

            # Process each repeat
            decomp_results = []

            for repeat_dir in repeat_dirs:
                # Find .fepout files in this repeat
                fepout_files = list(repeat_dir.glob("*.fepout"))
                if not fepout_files:
                    # Silently skip - it's normal to have incomplete repeats
                    continue

                try:
                    # Run decomposition for this repeat
                    # Create unique cache name for this repeat
                    repeat_cache = (
                        f"{cache_name}_{sys_type}_{repeat_dir.name}" if cache_name else f"{sys_type}_{repeat_dir.name}"
                    )
                    decomp_data = self.run_decomposition(
                        fepout_files=[str(f) for f in fepout_files], cache_name=repeat_cache, **decomp_kwargs
                    )

                    # Save intermediate decomposition data for this repeat
                    intermediate_file = repeat_dir / f"{sys_type}_{repeat_dir.name}_decomposition.csv"
                    decomp_data.to_csv(intermediate_file, index=False)

                    # Extract final values (last row, skip fraction column)
                    final_values = decomp_data.iloc[-1][1:].values
                    decomp_results.append(final_values)

                except Exception as e:
                    print(f"  ✗ {repeat_dir.name}: Error - {str(e)}")
                    continue

            if decomp_results:
                # Average across repeats
                decomp_array = np.array(decomp_results)
                mean_decomp = np.mean(decomp_array, axis=0)
                std_decomp = np.std(decomp_array, axis=0)

                results[sys_type] = mean_decomp
                results[f"{sys_type}_std"] = std_decomp
                results[f"{sys_type}_n"] = len(decomp_results)

                logger.info(
                    f"✓ {sys_type.capitalize()}: ΔG = {mean_decomp[0]:.2f} ± {std_decomp[0]:.2f} kcal/mol "
                    f"(n={len(decomp_results)} repeats)"
                )
            else:
                logger.warning(f"No valid {sys_type} calculations processed")

        # Calculate binding free energy difference if both complex and ligand available
        if "complex" in results and "ligand" in results:
            ddG = results["complex"] - results["ligand"]
            ddG_std = np.sqrt(results["complex_std"] ** 2 + results["ligand_std"] ** 2)

            results["ddG"] = ddG
            results["ddG_std"] = ddG_std

        return results

    def calculate_convergence_data(
        self,
        fepout_files: Union[str, List[str]],
        num_points: int = 10,
        cache_name: Optional[str] = None,
        **decomp_kwargs,
    ) -> pd.DataFrame:
        """
        Calculate FEP convergence over simulation time.

        This method runs the decomposition at multiple time points to show
        how the free energy estimate converges as more simulation data is included.

        Parameters
        ----------
        fepout_files : str or list of str
            Path(s) to .fepout file(s)
        num_points : int
            Number of convergence time points to calculate
        cache_name : str, optional
            Custom cache name
        **decomp_kwargs
            Additional arguments passed to run_decomposition()

        Returns
        -------
        pd.DataFrame
            Convergence data with columns:
            - fraction: Simulation time fraction (0 to 1)
            - sum_dg: Total free energy at this time point
            - sum_delec: Electrostatic contribution
            - sum_dvdw: van der Waals contribution
            - sum_couple: Coupling term
        """
        logger.info(f"Calculating convergence with {num_points} time points")

        # Simply call run_decomposition with the specified num_points
        convergence_data = self.run_decomposition(
            fepout_files=fepout_files, num_points=num_points, cache_name=cache_name, **decomp_kwargs
        )

        logger.info(f"Convergence analysis complete: {len(convergence_data)} data points")
        return convergence_data

    def batch_calculate_systems(
        self, base_dir: str, system_pattern: str, max_systems: Optional[int] = None, parallel: bool = True, **kwargs
    ) -> Dict[str, Dict[str, np.ndarray]]:
        """
        Batch process multiple FEP systems.

        Parameters
        ----------
        base_dir : str
            Base directory containing multiple system directories
        system_pattern : str
            Pattern to match system directories (e.g., "1p-*" for position 1')
        max_systems : int, optional
            Maximum number of systems to process
        parallel : bool
            Whether to process systems in parallel
        **kwargs
            Additional arguments passed to calculate_system_decomposition()

        Returns
        -------
        dict
            Dictionary mapping system names to their decomposition results
        """
        base_dir = Path(base_dir)
        system_dirs = sorted(base_dir.glob(system_pattern))

        if max_systems is not None:
            system_dirs = system_dirs[:max_systems]

        if not system_dirs:
            logger.warning(f"No systems found matching pattern: {system_pattern}")
            return {}

        if parallel:
            logger.info("Parallel processing requested but not implemented; running serially.")

        logger.info(f"Batch processing {len(system_dirs)} systems")

        results = {}

        for system_dir in system_dirs:
            system_name = system_dir.name
            logger.info(f"\n=== Processing system: {system_name} ===")

            try:
                system_results = self.calculate_system_decomposition(str(system_dir), cache_name=system_name, **kwargs)

                if system_results:
                    results[system_name] = system_results
                    logger.info(f"✓ {system_name} completed successfully")
                else:
                    logger.warning(f"✗ {system_name} produced no results")

            except Exception as e:
                logger.error(f"✗ Error processing {system_name}: {e}")
                continue

        logger.info(f"\n✓ Batch processing complete: {len(results)}/{len(system_dirs)} systems succeeded")
        return results

    def parse_fepout_raw(
        self, fepout_file: str, cache_name: Optional[str] = None, force_recompute: bool = False
    ) -> pd.DataFrame:
        """
        Parse raw dG/lambda data from NAMD .fepout file.

        This method extracts the free energy change data for each lambda window
        using fast grep parsing. Useful for creating dG/dλ plots.

        Parameters
        ----------
        fepout_file : str
            Path to .fepout file
        cache_name : str, optional
            Custom cache name
        force_recompute : bool
            If True, ignore cache and recompute

        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
            - start: Lambda window start
            - stop: Lambda window stop
            - dG: Free energy change in this window (kcal/mol)
            - dG_accum: Accumulated free energy (kcal/mol)
        """
        fepout_file = Path(fepout_file)
        if not fepout_file.exists():
            raise FileNotFoundError(f"FEP output file not found: {fepout_file}")

        # Create cache key
        if cache_name is None:
            cache_key = f"namd_fep_raw_{fepout_file.stem}"
        else:
            cache_key = f"namd_fep_raw_{cache_name}"

        cache_file = self._cache_dir / f"{cache_key}.pkl"

        # Check cache
        if not force_recompute and cache_file.exists():
            logger.info(f"✓ Loading cached raw data from {cache_file.name}")
            with open(cache_file, "rb") as f:
                return pickle.load(f)

        logger.info(f"Parsing raw FEP data from {fepout_file.name}")

        try:
            # Use grep for fast extraction (same as reference script)
            cmd = ["grep", "#Free energy", str(fepout_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)

            if not result.stdout:
                raise ValueError(f"No '#Free energy' lines found in {fepout_file}")

            lines = result.stdout.strip().split("\n")

            # Parse data: columns are at indices 7, 8, 11, 18
            data = []
            for line in lines:
                parts = line.split()
                data.append(
                    [
                        float(parts[7]),  # start (lambda1)
                        float(parts[8]),  # stop (lambda2)
                        float(parts[11]),  # dG
                        float(parts[18]),  # dG_accum (net change)
                    ]
                )

            df = pd.DataFrame(data, columns=["start", "stop", "dG", "dG_accum"])

            logger.info(f"Parsed {len(df)} lambda windows, final ΔG = {df['dG_accum'].iloc[-1]:.2f} kcal/mol")

            # Cache results
            with open(cache_file, "wb") as f:
                pickle.dump(df, f)

            return df

        except subprocess.CalledProcessError as e:
            logger.error(f"Grep command failed: {e}")
            raise
        except Exception as e:
            logger.error(f"Error parsing raw FEP data: {e}")
            raise

    def batch_process_fep_system(
        self, system_dir: str, output_dir: Optional[str] = None, save_data: bool = True, **kwargs
    ) -> Dict:
        """
        Comprehensive one-stop processing for a complete FEP system.

        This method processes all complex and ligand repeats, calculates
        decomposition, parses raw data, computes statistics, and optionally
        saves results to files. Designed for production use.

        Parameters
        ----------
        system_dir : str
            Path to FEP system directory containing complex*/ligand* subdirectories
        output_dir : str, optional
            Directory to save output files. If None, saves to system_dir.
        save_data : bool
            Whether to save results to JSON and CSV files
        **kwargs
            Additional arguments passed to run_decomposition()

        Returns
        -------
        dict
            Comprehensive results dictionary containing:
            - 'decomposition': Decomposition results (mean ± std)
            - 'raw_data': Combined raw dG/lambda data for all repeats
            - 'statistics': Statistical summary
            - 'file_paths': Paths to saved output files (if save_data=True)
        """
        system_dir = Path(system_dir).resolve()
        if not system_dir.exists():
            raise FileNotFoundError(f"System directory not found: {system_dir}")

        if output_dir is None:
            output_dir = system_dir
        else:
            output_dir = Path(output_dir).resolve()
            output_dir.mkdir(parents=True, exist_ok=True)

        results = {}

        # 1. Run decomposition analysis
        decomp_results = self.calculate_system_decomposition(str(system_dir), cache_name=system_dir.name, **kwargs)
        results["decomposition"] = decomp_results

        # 2. Parse raw FEP data from all repeats
        raw_data_combined = []

        for phase in ["complex", "ligand"]:
            # Find repeat directories
            repeat_dirs = sorted([d for d in system_dir.iterdir() if d.is_dir() and d.name.startswith(phase)])

            if not repeat_dirs:
                logger.warning(f"No {phase} directories found")
                continue

            phase_data_all = []

            for repeat_dir in repeat_dirs:
                # Find .fepout files
                fepout_files = list(repeat_dir.glob("*-prod-forward.fepout"))
                if not fepout_files:
                    fepout_files = list(repeat_dir.glob("*.fepout"))

                if not fepout_files:
                    # Silently skip - normal to have incomplete repeats
                    continue

                fepout_file = fepout_files[0]  # Use first (should be only one)

                try:
                    # Parse raw data
                    raw_data = self.parse_fepout_raw(str(fepout_file), cache_name=f"{phase}_{repeat_dir.name}")

                    # Add metadata columns
                    raw_data["phase"] = phase
                    raw_data["repeat"] = repeat_dir.name

                    phase_data_all.append(raw_data)
                    raw_data_combined.append(raw_data)

                except Exception as e:
                    logger.error(f"  ✗ Error parsing {repeat_dir.name}: {e}")
                    continue

            # Calculate statistics for this phase
            if phase_data_all:
                # Get mean ± STD across repeats at each lambda point
                n_repeats = len(phase_data_all)

                # Assuming all repeats have same lambda windows
                ref_data = phase_data_all[0]
                n_windows = len(ref_data)

                dG_accum_all = np.array([df["dG_accum"].values for df in phase_data_all])
                dG_accum_mean = np.mean(dG_accum_all, axis=0)
                dG_accum_std = np.std(dG_accum_all, axis=0, ddof=1)  # Use STD, not SEM

                results[f"{phase}_raw_mean"] = dG_accum_mean
                results[f"{phase}_raw_std"] = dG_accum_std  # Changed from _sem to _std
                results[f"{phase}_lambda"] = ref_data[["start", "stop"]].values

        # Combine all raw data
        if raw_data_combined:
            results["raw_data"] = pd.concat(raw_data_combined, ignore_index=True)
        else:
            logger.warning("No raw FEP data could be parsed")
            results["raw_data"] = pd.DataFrame()

        # 3. Compile statistics
        results["statistics"] = {
            "system_name": system_dir.name,
            "n_complex_repeats": decomp_results.get("complex_n", 0),
            "n_ligand_repeats": decomp_results.get("ligand_n", 0),
        }

        if "complex" in decomp_results:
            results["statistics"]["complex_dG"] = float(decomp_results["complex"][0])
            results["statistics"]["complex_dG_std"] = float(decomp_results["complex_std"][0])

        if "ligand" in decomp_results:
            results["statistics"]["ligand_dG"] = float(decomp_results["ligand"][0])
            results["statistics"]["ligand_dG_std"] = float(decomp_results["ligand_std"][0])

        if "ddG" in decomp_results:
            results["statistics"]["binding_ddG"] = float(decomp_results["ddG"][0])
            results["statistics"]["binding_ddG_std"] = float(decomp_results["ddG_std"][0])

        # 4. Save output files
        if save_data:
            results["file_paths"] = {}

            # Save decomposition summary (JSON)
            json_file = output_dir / "namd_fep_summary.json"
            import json

            with open(json_file, "w") as f:
                # Convert numpy arrays to lists for JSON serialization
                json_data = {}
                for key, val in results["statistics"].items():
                    json_data[key] = float(val) if isinstance(val, (np.floating, np.integer)) else val

                # Add decomposition details
                if "decomposition" in results:
                    decomp = results["decomposition"]
                    for phase in ["complex", "ligand", "ddG"]:
                        if phase in decomp:
                            json_data[f"{phase}_decomp"] = {
                                "dG": float(decomp[phase][0]),
                                "dElec": float(decomp[phase][1]),
                                "dVdW": float(decomp[phase][2]),
                                "dCouple": float(decomp[phase][3]),
                            }
                            if f"{phase}_std" in decomp:
                                json_data[f"{phase}_decomp_std"] = {
                                    "dG": float(decomp[f"{phase}_std"][0]),
                                    "dElec": float(decomp[f"{phase}_std"][1]),
                                    "dVdW": float(decomp[f"{phase}_std"][2]),
                                    "dCouple": float(decomp[f"{phase}_std"][3]),
                                }

                json.dump(json_data, f, indent=2)
            results["file_paths"]["summary_json"] = str(json_file)

            # Save detailed summary text file
            text_file = output_dir / "namd_fep_detailed_summary.txt"
            with open(text_file, "w") as f:
                f.write("=" * 60 + "\n")
                f.write(f"NAMD FEP Analysis - {system_dir.name}\n")
                f.write("=" * 60 + "\n\n")

                if "complex" in decomp_results:
                    f.write("Complex Phase:\n")
                    f.write(
                        f"  Final ΔG: {decomp_results['complex'][0]:.2f} ± {decomp_results['complex_std'][0]:.2f} kcal/mol\n"
                    )
                    f.write(f"  Number of repeats: {decomp_results['complex_n']}\n")
                    if "complex_lambda" in results:
                        f.write(f"  Lambda windows: {len(results['complex_lambda'])}\n")
                    f.write(f"  Energy decomposition:\n")
                    f.write(
                        f"    ΔElec:   {decomp_results['complex'][1]:7.2f} ± {decomp_results['complex_std'][1]:.2f} kcal/mol\n"
                    )
                    f.write(
                        f"    ΔVdW:    {decomp_results['complex'][2]:7.2f} ± {decomp_results['complex_std'][2]:.2f} kcal/mol\n"
                    )
                    f.write(
                        f"    ΔCouple: {decomp_results['complex'][3]:7.2f} ± {decomp_results['complex_std'][3]:.2f} kcal/mol\n"
                    )
                    f.write("\n")

                if "ligand" in decomp_results:
                    f.write("Ligand Phase:\n")
                    f.write(
                        f"  Final ΔG: {decomp_results['ligand'][0]:.2f} ± {decomp_results['ligand_std'][0]:.2f} kcal/mol\n"
                    )
                    f.write(f"  Number of repeats: {decomp_results['ligand_n']}\n")
                    if "ligand_lambda" in results:
                        f.write(f"  Lambda windows: {len(results['ligand_lambda'])}\n")
                    f.write(f"  Energy decomposition:\n")
                    f.write(
                        f"    ΔElec:   {decomp_results['ligand'][1]:7.2f} ± {decomp_results['ligand_std'][1]:.2f} kcal/mol\n"
                    )
                    f.write(
                        f"    ΔVdW:    {decomp_results['ligand'][2]:7.2f} ± {decomp_results['ligand_std'][2]:.2f} kcal/mol\n"
                    )
                    f.write(
                        f"    ΔCouple: {decomp_results['ligand'][3]:7.2f} ± {decomp_results['ligand_std'][3]:.2f} kcal/mol\n"
                    )
                    f.write("\n")

                if "ddG" in decomp_results:
                    f.write("Binding Free Energy:\n")
                    f.write(f"  ΔΔG: {decomp_results['ddG'][0]:.2f} ± {decomp_results['ddG_std'][0]:.2f} kcal/mol\n")
                    f.write(f"  Energy decomposition:\n")
                    f.write(
                        f"    ΔΔElec:   {decomp_results['ddG'][1]:7.2f} ± {decomp_results['ddG_std'][1]:.2f} kcal/mol\n"
                    )
                    f.write(
                        f"    ΔΔVdW:    {decomp_results['ddG'][2]:7.2f} ± {decomp_results['ddG_std'][2]:.2f} kcal/mol\n"
                    )
                    f.write(
                        f"    ΔΔCouple: {decomp_results['ddG'][3]:7.2f} ± {decomp_results['ddG_std'][3]:.2f} kcal/mol\n"
                    )

            results["file_paths"]["detailed_summary_txt"] = str(text_file)

            # Save raw dG/lambda data (combined CSV)
            if not results["raw_data"].empty:
                csv_file = output_dir / "namd_fep_raw_data.csv"
                results["raw_data"].to_csv(csv_file, index=False)
                results["file_paths"]["raw_data_csv"] = str(csv_file)

            # Save decomposition data
            if "decomposition" in results and results["decomposition"]:
                decomp_file = output_dir / "namd_fep_decomposition.csv"
                decomp = results["decomposition"]

                decomp_df = pd.DataFrame(
                    {
                        "Phase": [],
                        "dG": [],
                        "dG_std": [],
                        "dElec": [],
                        "dElec_std": [],
                        "dVdW": [],
                        "dVdW_std": [],
                        "dCouple": [],
                        "dCouple_std": [],
                        "n_repeats": [],
                    }
                )

                for phase in ["complex", "ligand", "ddG"]:
                    if phase in decomp:
                        row = {
                            "Phase": phase.capitalize(),
                            "dG": decomp[phase][0],
                            "dElec": decomp[phase][1],
                            "dVdW": decomp[phase][2],
                            "dCouple": decomp[phase][3],
                            "n_repeats": decomp.get(f"{phase}_n", 0),
                        }
                        if f"{phase}_std" in decomp:
                            row["dG_std"] = decomp[f"{phase}_std"][0]
                            row["dElec_std"] = decomp[f"{phase}_std"][1]
                            row["dVdW_std"] = decomp[f"{phase}_std"][2]
                            row["dCouple_std"] = decomp[f"{phase}_std"][3]
                        else:
                            row["dG_std"] = 0
                            row["dElec_std"] = 0
                            row["dVdW_std"] = 0
                            row["dCouple_std"] = 0

                        decomp_df = pd.concat([decomp_df, pd.DataFrame([row])], ignore_index=True)

                decomp_df.to_csv(decomp_file, index=False)
                results["file_paths"]["decomposition_csv"] = str(decomp_file)

        return results
