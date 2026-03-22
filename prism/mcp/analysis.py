"""Trajectory analysis tools (contacts, RMSD, H-bonds, MM/PBSA)."""

import os
import re
import json
import traceback

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def register(mcp):

    @mcp.tool()
    def analyze_trajectory(
        topology: str,
        trajectory: str,
        ligand_resname: str = "LIG",
        output_dir: str = "analysis_results",
        contact_cutoff: float = 0.4,
    ) -> str:
        """Run comprehensive protein-ligand trajectory analysis.

        Performs contact analysis, hydrogen bond analysis, and distance analysis
        on an MD trajectory, generating CSV data files and a summary report.

        Requires MDTraj and MDAnalysis packages.

        Args:
            topology: Absolute path to topology/structure file (.pdb, .gro).
            trajectory: Absolute path to trajectory file (.xtc, .trr, .dcd).
            ligand_resname: Residue name of the ligand in topology. Default: "LIG".
            output_dir: Directory for analysis output files. Default: "analysis_results".
            contact_cutoff: Distance cutoff for contacts in nm. Default: 0.4.

        Returns:
            JSON with analysis summary including top contacts, H-bonds, and output file paths.
        """
        logger.info(f"analyze_trajectory: {topology} + {trajectory}")

        errors = []
        if not os.path.exists(topology):
            errors.append(f"Topology file not found: {topology}")
        if not os.path.exists(trajectory):
            errors.append(f"Trajectory file not found: {trajectory}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.analysis.core.analyzer import analyze_and_visualize
                from prism.analysis.core.config import convert_numpy_types

                results = analyze_and_visualize(
                    topology_file=topology,
                    trajectory_file=trajectory,
                    resname=ligand_resname,
                    distance_cutoff=contact_cutoff,
                    output_dir=output_dir,
                    verbose=True,
                    save_all_data=True,
                    return_detailed=True,
                    generate_report=True,
                )

            # Extract summary for JSON response
            contact_props = results.get("contact_proportions", {})
            hbond_freqs = results.get("hbond_frequencies", [])
            system_info = results.get("system_info", {})

            # Top 10 contacts
            sorted_contacts = sorted(contact_props.items(), key=lambda x: x[1], reverse=True)[:10]
            top_contacts = [{"residue": r, "proportion": float(p)} for r, p in sorted_contacts]

            # Top 10 H-bonds
            top_hbonds = []
            for hb in hbond_freqs[:10]:
                key, freq, avg_dist, avg_angle = hb
                top_hbonds.append({
                    "bond": key,
                    "frequency": float(freq),
                    "avg_distance_A": float(avg_dist),
                    "avg_angle_deg": float(avg_angle),
                })

            # List output files
            output_files = []
            if os.path.isdir(output_dir):
                output_files = sorted(os.listdir(output_dir))

            summary = convert_numpy_types({
                "success": True,
                "output_dir": os.path.abspath(output_dir),
                "system_info": system_info,
                "top_contacts": top_contacts,
                "top_hbonds": top_hbonds,
                "total_contacts": len(contact_props),
                "total_hbonds": len(hbond_freqs),
                "output_files": output_files,
            })
            return json.dumps(summary, indent=2)

        except Exception as e:
            logger.error(f"analyze_trajectory failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}, indent=2)

    @mcp.tool()
    def analyze_contacts(
        topology: str,
        trajectory: str,
        ligand_resname: str = "LIG",
        cutoff: float = 0.4,
        top_n: int = 20,
    ) -> str:
        """Analyze protein-ligand contacts from an MD trajectory.

        Calculates contact proportions (fraction of frames each residue is in
        contact with the ligand) and returns the top contacting residues.

        Args:
            topology: Absolute path to topology/structure file (.pdb, .gro).
            trajectory: Absolute path to trajectory file (.xtc, .trr, .dcd).
            ligand_resname: Residue name of the ligand. Default: "LIG".
            cutoff: Contact distance cutoff in nm. Default: 0.4.
            top_n: Number of top contacts to return. Default: 20.

        Returns:
            JSON with ranked contact residues and their proportions.
        """
        logger.info(f"analyze_contacts: {topology} + {trajectory}")

        errors = []
        if not os.path.exists(topology):
            errors.append(f"Topology file not found: {topology}")
        if not os.path.exists(trajectory):
            errors.append(f"Trajectory file not found: {trajectory}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.analysis.core.analyzer import IntegratedProteinLigandAnalyzer
                from prism.analysis.core.config import AnalysisConfig, convert_numpy_types

                config = AnalysisConfig(
                    ligand_resname=ligand_resname,
                    contact_cutoff_nm=cutoff,
                    contact_enter_threshold_nm=cutoff * 0.75,
                    contact_exit_threshold_nm=cutoff * 1.125,
                )

                analyzer = IntegratedProteinLigandAnalyzer(
                    topology_file=topology,
                    trajectory_file=trajectory,
                    config=config,
                )

                contact_props = analyzer.calculate_contact_proportions()
                top_contacts = analyzer.select_top_contacts(contact_props, top_n=top_n)

            result = convert_numpy_types({
                "success": True,
                "n_frames": analyzer.traj.n_frames if analyzer.traj else 0,
                "cutoff_nm": cutoff,
                "total_contacting_residues": sum(1 for v in contact_props.values() if v > 0),
                "top_contacts": [
                    {"residue": r, "proportion": float(p)}
                    for r, p in top_contacts
                ],
            })
            return json.dumps(result, indent=2)

        except Exception as e:
            logger.error(f"analyze_contacts failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}, indent=2)

    @mcp.tool()
    def analyze_rmsd(
        topology: str,
        trajectory: str,
        selection: str = "protein and name CA",
    ) -> str:
        """Calculate RMSD for an MD trajectory.

        Computes the root-mean-square deviation of selected atoms over
        the trajectory, using the first frame as reference. Returns
        statistics (mean, std, min, max) and per-frame values.

        Args:
            topology: Absolute path to topology/structure file (.pdb, .gro).
            trajectory: Absolute path to trajectory file (.xtc, .trr, .dcd).
            selection: MDAnalysis atom selection string. Default: "protein and name CA".

        Returns:
            JSON with RMSD statistics and per-frame values (in nm).
        """
        logger.info(f"analyze_rmsd: {topology} + {trajectory}, selection='{selection}'")

        errors = []
        if not os.path.exists(topology):
            errors.append(f"Topology file not found: {topology}")
        if not os.path.exists(trajectory):
            errors.append(f"Trajectory file not found: {trajectory}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.analysis.calc.rmsd import RMSDAnalyzer
                from prism.analysis.core.config import AnalysisConfig, convert_numpy_types
                import numpy as np

                config = AnalysisConfig()
                rmsd_analyzer = RMSDAnalyzer(config)

                rmsd_data = rmsd_analyzer.calculate_rmsd(
                    universe=topology,
                    trajectory=trajectory,
                    selection=selection,
                )

                # rmsd_data shape: (n_frames, 3) — [frame, time, rmsd]
                # RMSD values are in Angstroms from MDAnalysis, convert to nm
                rmsd_values = rmsd_data[:, 2] / 10.0  # Å -> nm
                times = rmsd_data[:, 1]  # ps

            result = convert_numpy_types({
                "success": True,
                "selection": selection,
                "n_frames": len(rmsd_values),
                "statistics": {
                    "mean_nm": float(np.mean(rmsd_values)),
                    "std_nm": float(np.std(rmsd_values)),
                    "min_nm": float(np.min(rmsd_values)),
                    "max_nm": float(np.max(rmsd_values)),
                },
                "time_ps": times.tolist(),
                "rmsd_nm": rmsd_values.tolist(),
            })
            return json.dumps(result, indent=2)

        except Exception as e:
            logger.error(f"analyze_rmsd failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}, indent=2)

    @mcp.tool()
    def analyze_hbonds(
        topology: str,
        trajectory: str,
        ligand_resname: str = "LIG",
    ) -> str:
        """Analyze protein-ligand hydrogen bonds from an MD trajectory.

        Identifies hydrogen bonds between protein and ligand, computing
        their frequency, average distance, and average angle.

        Args:
            topology: Absolute path to topology/structure file (.pdb, .gro).
            trajectory: Absolute path to trajectory file (.xtc, .trr, .dcd).
            ligand_resname: Residue name of the ligand. Default: "LIG".

        Returns:
            JSON with hydrogen bond list sorted by frequency.
        """
        logger.info(f"analyze_hbonds: {topology} + {trajectory}")

        errors = []
        if not os.path.exists(topology):
            errors.append(f"Topology file not found: {topology}")
        if not os.path.exists(trajectory):
            errors.append(f"Trajectory file not found: {trajectory}")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.analysis.core.analyzer import IntegratedProteinLigandAnalyzer
                from prism.analysis.core.config import AnalysisConfig, convert_numpy_types

                config = AnalysisConfig(ligand_resname=ligand_resname)
                analyzer = IntegratedProteinLigandAnalyzer(
                    topology_file=topology,
                    trajectory_file=trajectory,
                    config=config,
                )

                hbond_freqs = analyzer.analyze_hydrogen_bonds()

            hbonds = []
            for hb in hbond_freqs:
                key, freq, avg_dist, avg_angle = hb
                hbonds.append({
                    "bond": key,
                    "frequency": float(freq),
                    "avg_distance_A": float(avg_dist),
                    "avg_angle_deg": float(avg_angle),
                })

            result = convert_numpy_types({
                "success": True,
                "n_frames": analyzer.traj.n_frames if analyzer.traj else 0,
                "total_hbonds": len(hbonds),
                "significant_hbonds": sum(1 for h in hbonds if h["frequency"] > 0.05),
                "hbonds": hbonds,
            })
            return json.dumps(result, indent=2)

        except Exception as e:
            logger.error(f"analyze_hbonds failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}, indent=2)

    @mcp.tool()
    def analyze_mmpbsa(
        mmpbsa_dir: str,
    ) -> str:
        """Parse MM/PBSA results from a completed gmx_MMPBSA calculation.

        Reads the FINAL_RESULTS_MMPBSA.dat file produced by gmx_MMPBSA and
        extracts the binding free energy components (van der Waals, electrostatic,
        polar solvation, non-polar solvation, and total).

        Also checks for the CSV output file if available.

        Call this AFTER mmpbsa_run.sh has completed successfully.

        Args:
            mmpbsa_dir: Absolute path to the GMX_PROLIG_MMPBSA directory.

        Returns:
            JSON with binding energy components and total ΔG_bind.
        """
        logger.info(f"analyze_mmpbsa: {mmpbsa_dir}")

        if not os.path.isdir(mmpbsa_dir):
            return json.dumps(
                {"success": False, "error": f"Directory not found: {mmpbsa_dir}"},
                indent=2,
            )

        # Look for the results file
        dat_path = os.path.join(mmpbsa_dir, "FINAL_RESULTS_MMPBSA.dat")
        csv_path = os.path.join(mmpbsa_dir, "FINAL_RESULTS_MMPBSA.csv")

        if not os.path.exists(dat_path):
            return json.dumps(
                {
                    "success": False,
                    "error": f"Results file not found: {dat_path}. Has mmpbsa_run.sh finished?",
                },
                indent=2,
            )

        try:
            with open(dat_path, "r") as f:
                content = f.read()

            result = {
                "success": True,
                "mmpbsa_dir": mmpbsa_dir,
                "results_file": dat_path,
                "energy_components": {},
                "warnings": [],
            }

            # Parse energy terms from gmx_MMPBSA output
            # Typical format:
            #   DELTA TOTAL = -35.2 +/- 2.1
            #   VDWAALS     = -42.1 +/- 1.5
            #   EEL         = -15.3 +/- 3.2
            #   EPB/EGB     =  28.4 +/- 2.8
            #   ESURF/ENPOLAR = -6.2 +/- 0.4
            energy_pattern = re.compile(
                r"^\s*([\w\s/]+?)\s*=\s*([-\d.]+)\s*\+/-\s*([-\d.]+)",
                re.MULTILINE,
            )
            for match in energy_pattern.finditer(content):
                name = match.group(1).strip()
                mean = float(match.group(2))
                std = float(match.group(3))
                result["energy_components"][name] = {
                    "mean_kcal_mol": mean,
                    "std_kcal_mol": std,
                }

            # Extract DELTA TOTAL as the headline binding energy
            delta_total = result["energy_components"].get("DELTA TOTAL")
            if delta_total:
                result["delta_g_bind_kcal_mol"] = delta_total["mean_kcal_mol"]
                result["delta_g_bind_std"] = delta_total["std_kcal_mol"]
            else:
                result["warnings"].append(
                    "Could not find 'DELTA TOTAL' in results. "
                    "File may use a different format — check raw content."
                )
                result["raw_content"] = content[:2000]

            # Check for CSV
            if os.path.exists(csv_path):
                result["csv_file"] = csv_path

            return json.dumps(result, indent=2)

        except Exception as e:
            logger.error(f"analyze_mmpbsa failed: {e}\n{traceback.format_exc()}")
            return json.dumps(
                {"success": False, "error": str(e), "traceback": traceback.format_exc()},
                indent=2,
            )
