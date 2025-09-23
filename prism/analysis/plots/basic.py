#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic plotting utilities for PRISM analysis module.

Moved from visualize.py to separate calculation and visualization concerns.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from typing import List
from ..config import AnalysisConfig

logger = logging.getLogger(__name__)


class BasicPlotter:
    """Basic plotting utilities for trajectory analysis"""

    def plot_distance_time(self, distances: np.ndarray,
                          residue_id: str, config: AnalysisConfig, save_path: str):
        """Plot distance vs time with thresholds"""
        try:
            times = np.arange(len(distances)) * 0.5
            distances_angstrom = distances * 10.0

            if len(distances_angstrom) > 10:
                smoothed_distances = np.convolve(
                    distances_angstrom,
                    np.ones(config.smooth_window)/config.smooth_window,
                    mode='valid'
                )
                smoothed_times = times[config.smooth_window-1:]
                plot_data = smoothed_distances
                plot_times = smoothed_times
            else:
                plot_data = distances_angstrom
                plot_times = times

            fig, ax = plt.subplots(figsize=(12, 6))
            ax.plot(plot_times, plot_data, 'b-', linewidth=2,
               label=f"Distance (window={config.smooth_window})")

            ax.grid(True)

            ax.axhline(y=config.contact_enter_threshold_nm * 10.0,
                      color='g', linestyle='--', label='Contact start')
            ax.axhline(y=config.contact_exit_threshold_nm * 10.0,
                      color='r', linestyle='--', label='Contact end')

            ax.set_xlabel("Time (ns)")
            ax.set_ylabel("Distance (Å)")
            ax.set_title(f"Distance to {residue_id}")
            ax.legend()

            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Distance plot saved: {save_path}")
        except Exception as e:
            logger.error(f"Failed to create distance plot: {e}")
            raise

    def plot_multi_system(self, all_data: List[np.ndarray], all_times: List[np.ndarray],
                         system_names: List[str], weights: List[float], color_mode: str,
                         ylabel: str, title: str, output_file: str, config: AnalysisConfig):
        """Plot multiple systems comparison"""
        try:
            fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
            ax.grid(alpha=0.3)

            numeric_weights = [float(w) for w in weights]
            combined = list(zip(numeric_weights, all_data, all_times, system_names))

            if color_mode == "mw" and len(set(numeric_weights)) > 1:
                combined.sort(key=lambda x: x[0])
                norm = mcolors.Normalize(vmin=min(numeric_weights), vmax=max(numeric_weights))
                cmap = cm.plasma
                smap = cm.ScalarMappable(norm=norm, cmap=cmap)
            else:
                cmap = plt.get_cmap('tab20')
                smap = None
                colors = [cmap(x) for x in np.linspace(0, 1, len(combined))]

            for idx, (weight, data, times, name) in enumerate(combined):
                if len(data) == 0:
                    logger.warning(f"No data for system {name}")
                    continue

                plot_data = data
                plot_times = times

                # Apply smoothing if data is long enough
                if len(data) > config.min_frames_for_smoothing:
                    smoothed_data = np.convolve(
                        data,
                        np.ones(config.smooth_window)/config.smooth_window,
                        mode='valid'
                    )
                    smoothed_times = times[config.smooth_window-1:]
                    if len(smoothed_data) > 0:
                        plot_data = smoothed_data
                        plot_times = smoothed_times

                if weight > 0:
                    label = f"{name} ({int(weight)} Da)"
                else:
                    label = name

                if color_mode == "mw" and smap is not None:
                    color = cmap(norm(weight))
                else:
                    color = colors[idx % len(colors)]

                ax.plot(plot_times, plot_data, linewidth=2.5, alpha=0.9,
                       label=label, color=color)

            # Add 150ns reference line if any trajectory is long enough
            if any(len(times) > 300 for times in all_times):  # 300 * 0.5ns = 150ns
                ax.axvline(x=150, color='green', linestyle='--', linewidth=2, alpha=0.8)
                ax.text(152, ax.get_ylim()[1]*0.95, "150 ns", color='green',
                       fontsize=12, fontweight='bold')

            # Position legend outside plot area
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                     fontsize=10, ncol=1, frameon=True)

            # Add colorbar if using molecular weight coloring
            if color_mode == "mw" and smap is not None:
                cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
                cbar = fig.colorbar(smap, cax=cbar_ax)
                cbar.set_label("Molecular Weight (Da)")

            ax.set_xlabel("Time (ns)", fontweight='bold')
            ax.set_ylabel(ylabel, fontweight='bold')
            ax.set_title(title, fontweight='bold')

            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Multi-system plot saved: {output_file}")
        except Exception as e:
            logger.error(f"Failed to create multi-system plot: {e}")
            raise


# Legacy compatibility
class Visualizer(BasicPlotter):
    """Legacy compatibility class - use BasicPlotter instead"""
    pass