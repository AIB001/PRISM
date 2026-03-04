"""PMF (Potential of Mean Force) analysis tools."""

import json
import traceback

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def register(mcp):

    @mcp.tool()
    def analyze_pmf(
        pmf_dir: str,
        generate_plots: bool = True,
    ) -> str:
        """Analyze PMF results from steered MD and umbrella sampling.

        Runs all available analyses on a completed PMF simulation:
        - SMD analysis: rupture force, work
        - Umbrella window validation: completion rate
        - Histogram overlap: sampling quality
        - PMF profile: binding free energy, barriers
        - Convergence: block averaging

        Optionally generates publication-quality plots (force, PMF profile,
        histogram overlap).

        Args:
            pmf_dir: Absolute path to the GMX_PROLIG_PMF directory.
            generate_plots: Whether to generate analysis plots. Default: true.

        Returns:
            JSON with combined analysis results from all available stages.
        """
        logger.info(f"analyze_pmf: {pmf_dir}")

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                from prism.pmf.analysis import PMFAnalyzer
                import numpy as np

                analyzer = PMFAnalyzer(pmf_dir, verbose=True)
                combined = {"success": True, "pmf_dir": pmf_dir}

                # --- SMD analysis ---
                try:
                    smd = analyzer.analyze_smd()
                    combined["smd"] = {
                        "rupture_force": smd["rupture_force"],
                        "rupture_time": smd["rupture_time"],
                        "rupture_position": smd["rupture_position"],
                        "work": smd["work"],
                    }
                    if generate_plots:
                        try:
                            combined["smd"]["plot"] = analyzer.plot_smd_force()
                        except Exception:
                            pass
                except FileNotFoundError:
                    combined["smd"] = {"status": "not_available", "reason": "SMD data files not found"}
                except Exception as e:
                    combined["smd"] = {"status": "error", "reason": str(e)}

                # --- Umbrella window validation ---
                try:
                    validation = analyzer.validate_umbrella_windows()
                    combined["umbrella_validation"] = {
                        "total_windows": validation["total_windows"],
                        "completed_windows": validation["completed_windows"],
                        "completion_rate": validation["completion_rate"],
                        "failed_windows": validation["failed_windows"],
                    }
                except FileNotFoundError:
                    combined["umbrella_validation"] = {"status": "not_available", "reason": "Umbrella directory not found"}
                except Exception as e:
                    combined["umbrella_validation"] = {"status": "error", "reason": str(e)}

                # --- Histogram overlap ---
                try:
                    hist = analyzer.check_histogram_overlap()
                    combined["histogram_overlap"] = {
                        "n_windows": hist["n_windows"],
                        "min_overlap": hist["min_overlap"],
                        "mean_overlap": hist["mean_overlap"],
                        "poorly_sampled_pairs": hist["poorly_sampled"],
                    }
                    if generate_plots:
                        try:
                            combined["histogram_overlap"]["plot"] = analyzer.plot_histogram_overlap()
                        except Exception:
                            pass
                except FileNotFoundError:
                    combined["histogram_overlap"] = {"status": "not_available", "reason": "histo.xvg not found"}
                except Exception as e:
                    combined["histogram_overlap"] = {"status": "error", "reason": str(e)}

                # --- PMF profile ---
                try:
                    pmf = analyzer.analyze_pmf_profile()
                    combined["pmf_profile"] = {
                        "binding_energy": pmf["binding_energy"],
                        "min_pmf": pmf["min_pmf"],
                        "min_distance_nm": pmf["min_distance"],
                        "barrier_height": pmf["barrier_height"],
                        "barrier_distance_nm": pmf["barrier_distance"],
                    }
                    if generate_plots:
                        try:
                            combined["pmf_profile"]["plot"] = analyzer.plot_pmf_profile()
                        except Exception:
                            pass
                except FileNotFoundError:
                    combined["pmf_profile"] = {"status": "not_available", "reason": "pmf.xvg not found"}
                except Exception as e:
                    combined["pmf_profile"] = {"status": "error", "reason": str(e)}

                # --- Convergence ---
                try:
                    conv = analyzer.analyze_convergence()
                    combined["convergence"] = {
                        "converged": conv["converged"],
                        "n_blocks": conv["n_blocks"],
                        "n_windows_analyzed": len(conv["window_means"]),
                    }
                except Exception as e:
                    combined["convergence"] = {"status": "error", "reason": str(e)}

                # --- Generate text report ---
                try:
                    report_path = analyzer.generate_report()
                    combined["report_path"] = report_path
                except Exception:
                    pass

            # Convert any remaining numpy types
            def _convert(obj):
                import numpy as np
                if isinstance(obj, (np.integer,)):
                    return int(obj)
                if isinstance(obj, (np.floating,)):
                    return float(obj)
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                if isinstance(obj, dict):
                    return {k: _convert(v) for k, v in obj.items()}
                if isinstance(obj, list):
                    return [_convert(i) for i in obj]
                return obj

            return json.dumps(_convert(combined), indent=2)

        except Exception as e:
            logger.error(f"analyze_pmf failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}, indent=2)
