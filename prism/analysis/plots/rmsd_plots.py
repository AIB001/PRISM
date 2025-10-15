#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSD/RMSF plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
from pathlib import Path

# Import publication style
from .publication_utils import apply_publication_style, PUBLICATION_FONTS, get_standard_figsize


def plot_rmsd_simple_timeseries(rmsd_results: Dict,
                               output_path: str,
                               title: str = "") -> bool:
    """
    Create simple RMSD time series plot matching the example format.

    Parameters
    ----------
    rmsd_results : dict
        Dictionary with trajectory names as keys and RMSD data as values
    output_path : str
        Path to save the plot
    title : str
        Plot title (empty by default for publication style)

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("📈 RMSD timeseries analysis showing structural stability over simulation time")
        apply_publication_style()

        # Create single panel figure (matching example)
        fig, ax = plt.subplots(figsize=get_standard_figsize('single'))

        # Use different colors for each trajectory
        colors = ['#4472C4', '#ED7D31', '#70AD47', '#FFC000', '#5B9BD5', '#C55A11']

        # Plot all trajectories with different colors
        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                time = np.arange(len(protein_rmsd)) * 0.5  # 0.5 ns intervals
                color = colors[i % len(colors)]
                ax.plot(time, protein_rmsd, color=color, linewidth=1.5, alpha=0.85, label=traj_name)

        # Format exactly like the example
        ax.set_xlabel('Time (ns)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
        ax.set_ylabel('Ligand RMSD (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')  # Match example label
        ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])

        # Remove grid and titles for clean look
        ax.grid(False)

        # Set clean axis limits
        all_rmsd = []
        for rmsd_data in rmsd_results.values():
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                all_rmsd.extend(rmsd_data['protein'])
            else:
                all_rmsd.extend(rmsd_data)

        if all_rmsd:
            ax.set_ylim(0, max(all_rmsd) * 1.05)

        # Add legend if multiple trajectories
        if len(rmsd_results) > 1:
            ax.legend(loc='best', fontsize=PUBLICATION_FONTS['legend'], framealpha=0.9)

        # Only add title if provided

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in RMSD plotting: {e}")
        return False


def plot_rmsd_analysis(rmsd_results: Dict,
                      output_path: str,
                      title: str = "") -> bool:
    """
    Plot RMSD analysis results for multiple trajectories.

    Parameters
    ----------
    rmsd_results : dict
        Dictionary with trajectory names as keys and RMSD data as values
    output_path : str
        Path to save the plot
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("📈 RMSD 4-panel analysis: timeseries, distribution, averages, and ligand comparison")
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=get_standard_figsize('quad'))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

        # Plot 1: RMSD time series
        ax1 = axes[0, 0]
        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                time = np.arange(len(protein_rmsd)) * 0.5  # Assume 0.5ns interval
                ax1.plot(time, protein_rmsd, color=colors[i % 3],
                        label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel('Time (ns)', fontfamily='Times New Roman')
        ax1.set_ylabel('RMSD (Å)', fontfamily='Times New Roman')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: RMSD distribution
        ax2 = axes[0, 1]
        all_rmsd = []
        for rmsd_data in rmsd_results.values():
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data
            if len(protein_rmsd) > 0:
                all_rmsd.extend(protein_rmsd)

        if all_rmsd:
            ax2.hist(all_rmsd, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
            ax2.set_xlabel('RMSD (Å)', fontfamily='Times New Roman')
            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax2.grid(True, alpha=0.3)

        # Plot 3: Average RMSD per trajectory
        ax3 = axes[1, 0]
        traj_names = []
        avg_rmsd = []
        std_rmsd = []

        for traj_name, rmsd_data in rmsd_results.items():
            if isinstance(rmsd_data, dict) and 'protein' in rmsd_data:
                protein_rmsd = rmsd_data['protein']
            else:
                protein_rmsd = rmsd_data

            if len(protein_rmsd) > 0:
                traj_names.append(traj_name)
                avg_rmsd.append(np.mean(protein_rmsd))
                std_rmsd.append(np.std(protein_rmsd))

        if traj_names:
            bars = ax3.bar(traj_names, avg_rmsd, yerr=std_rmsd,
                          color=colors[:len(traj_names)], alpha=0.7, capsize=5)
            ax3.set_ylabel('Average RMSD (Å)', fontfamily='Times New Roman')
            ax3.grid(True, alpha=0.3, axis='y')

        # Plot 4: Ligand RMSD (if available)
        ax4 = axes[1, 1]
        has_ligand = False

        for i, (traj_name, rmsd_data) in enumerate(rmsd_results.items()):
            if isinstance(rmsd_data, dict) and 'ligand' in rmsd_data:
                ligand_rmsd = rmsd_data['ligand']
                if len(ligand_rmsd) > 0:
                    time = np.arange(len(ligand_rmsd)) * 0.5
                    ax4.plot(time, ligand_rmsd, color=colors[i % 3],
                            label=traj_name, linewidth=2, alpha=0.8)
                    has_ligand = True

        if has_ligand:
            ax4.set_xlabel('Time (ns)', fontfamily='Times New Roman')
            ax4.set_ylabel('Ligand RMSD (Å)', fontfamily='Times New Roman')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(0.5, 0.5, 'No Ligand Data', ha='center', va='center',
                    transform=ax4.transAxes, fontsize=16, fontfamily='Times New Roman')
            ax4.set_xticks([])
            ax4.set_yticks([])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in RMSD plotting: {e}")
        return False


def plot_rmsf_analysis(rmsf_results: Dict,
                      output_path: str,
                      title: str = "") -> bool:
    """
    Plot RMSF analysis results.

    Parameters
    ----------
    rmsf_results : dict
        Dictionary with trajectory names as keys and RMSF data as values
    output_path : str
        Path to save the plot
    title : str
        Plot title

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        print("📈 RMSF 2-panel analysis: per-residue fluctuations and distribution comparison")
        apply_publication_style()

        fig, axes = plt.subplots(1, 2, figsize=get_standard_figsize('horizontal'))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

        # Plot 1: RMSF per residue
        ax1 = axes[0]
        for i, (traj_name, rmsf_data) in enumerate(rmsf_results.items()):
            if isinstance(rmsf_data, tuple) and len(rmsf_data) >= 2:
                residue_indices, rmsf_values = rmsf_data[:2]
            elif isinstance(rmsf_data, dict) and 'rmsf' in rmsf_data:
                rmsf_values = rmsf_data['rmsf']
                residue_indices = np.arange(len(rmsf_values))
            else:
                continue

            if len(rmsf_values) > 0:
                ax1.plot(residue_indices, rmsf_values, color=colors[i % 3],
                        label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel('Residue Index', fontfamily='Times New Roman')
        ax1.set_ylabel('RMSF (Å)', fontfamily='Times New Roman')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: RMSF distribution
        ax2 = axes[1]
        all_rmsf = []
        for rmsf_data in rmsf_results.values():
            if isinstance(rmsf_data, tuple) and len(rmsf_data) >= 2:
                rmsf_values = rmsf_data[1]
            elif isinstance(rmsf_data, dict) and 'rmsf' in rmsf_data:
                rmsf_values = rmsf_data['rmsf']
            else:
                continue

            if len(rmsf_values) > 0:
                all_rmsf.extend(rmsf_values)

        if all_rmsf:
            ax2.hist(all_rmsf, bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
            ax2.set_xlabel('RMSF (Å)', fontfamily='Times New Roman')
            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in RMSF plotting: {e}")
        return False


def plot_rmsf_with_auto_chains(rmsf_results: Dict[str, Tuple[np.ndarray, np.ndarray]],
                                output_path: str,
                                title: str = "",
                                plot_type: str = "multi_panel",
                                figsize: Optional[Tuple[float, float]] = None,
                                save_path: Optional[str] = None,
                                max_chains_per_row: int = 2) -> bool:
    """
    Create RMSF plots with automatic chain detection and flexible layout.

    Parameters
    ----------
    rmsf_results : dict
        Dictionary mapping chain names to (rmsf_values, residue_ids) tuples
        e.g., {"Chain PR1": (rmsf_array, residue_ids), ...}
    output_path : str
        Path to save the plot
    title : str
        Plot title (empty by default for publication style)
    plot_type : str
        Type of plot layout: "multi_panel", "single_panel", "separate_files"
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    save_path : str, optional
        Deprecated: use output_path instead
    max_chains_per_row : int
        Maximum number of chains per row in multi-panel layout

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        apply_publication_style()

        if save_path is not None:
            output_path = save_path

        if not rmsf_results:
            print("No RMSF results provided for plotting")
            return False

        n_chains = len(rmsf_results)
        if n_chains == 0:
            print("No chain data available for RMSF plotting")
            return False

        # Sort chains by name for consistent ordering
        sorted_chains = sorted(rmsf_results.items(), key=lambda x: x[0])
        chain_names = [name for name, _ in sorted_chains]

        print(f"Creating RMSF plot for {n_chains} chains: {chain_names}")

        if plot_type == "multi_panel":
            return _plot_rmsf_multi_panel(sorted_chains, output_path, title, figsize, max_chains_per_row)
        elif plot_type == "single_panel":
            return _plot_rmsf_single_panel(sorted_chains, output_path, title, figsize)
        elif plot_type == "separate_files":
            return _plot_rmsf_separate_files(sorted_chains, output_path, title, figsize)
        else:
            print(f"Unknown plot_type: {plot_type}. Using multi_panel.")
            return _plot_rmsf_multi_panel(sorted_chains, output_path, title, figsize, max_chains_per_row)

    except Exception as e:
        print(f"Error in RMSF plotting with auto chains: {e}")
        return False


def _plot_rmsf_multi_panel(sorted_chains: List[Tuple[str, Tuple[np.ndarray, np.ndarray]]],
                           output_path: str,
                           title: str,
                           figsize: Optional[Tuple[float, float]],
                           max_chains_per_row: int) -> bool:
    """Create multi-panel RMSF plot with automatic chain type detection."""
    try:
        n_chains = len(sorted_chains)
        rows = (n_chains + max_chains_per_row - 1) // max_chains_per_row
        cols = min(n_chains, max_chains_per_row)

        # Auto-calculate figure size for publication quality
        if figsize is None:
            # Make panels similar in size to the original figure
            panel_width = 6  # Each panel width in inches
            panel_height = 4  # Each panel height in inches
            figsize = (panel_width * cols, panel_height * rows)

        fig, axes = plt.subplots(rows, cols, figsize=figsize)
        if rows == 1 and cols == 1:
            axes = [axes]
        elif rows == 1:
            axes = list(axes)
        else:
            axes = axes.flatten()

        # Enhanced chain type detection and colors
        protein_chains = [(name, data) for name, data in sorted_chains if 'Protein' in name or ('Chain' in name and 'Nucleic' not in name)]
        nucleic_chains = [(name, data) for name, data in sorted_chains if 'Nucleic' in name or 'RNA' in name or 'DNA' in name]

        # Standard colors for publication
        protein_color = '#4472C4'  # Professional blue
        nucleic_color = '#E15759'  # Professional red/coral

        for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(sorted_chains):
            ax = axes[i]

            # Determine chain type and appropriate labeling
            is_nucleic = any(keyword in chain_name.lower() for keyword in ['nucleic', 'rna', 'dna'])

            if is_nucleic:
                color = nucleic_color
                x_label = "Residue Number"
                # Clean up chain name for display
                display_name = chain_name.replace("Nucleic Chain", "Nucleic Chain").replace("Chain", "Chain")
                print(f"  📊 Plotting nucleic acid chain: {display_name} ({len(rmsf_values)} P atoms)")
            else:
                color = protein_color
                x_label = "Residue Number"
                # Clean up chain name for display
                display_name = chain_name.replace("Protein Chain", "Chain").replace("Chain", "Chain")
                print(f"  📊 Plotting protein chain: {display_name} ({len(rmsf_values)} CA atoms)")

            # Plot RMSF with filled area
            ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.9)
            ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color=color)

            # Calculate and print statistics (removed from plot to avoid hiding content)
            mean_rmsf = np.mean(rmsf_values)
            std_rmsf = np.std(rmsf_values)
            print(f"  📊 {display_name} RMSF: Mean = {mean_rmsf:.1f} ± {std_rmsf:.2f} Å")

            # Add panel title (matching the original figure style)
            panel_title = display_name.replace("Protein ", "").replace("Nucleic ", "Nucleic ")
            ax.set_title(panel_title, fontsize=PUBLICATION_FONTS['title'], fontweight='bold', pad=10)

            # Formatting with appropriate x-axis label
            ax.set_xlabel(x_label, fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
            ax.set_ylabel('RMSF (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])

            # Set y-axis to start from 0 for better comparison
            ax.set_ylim(bottom=0)

        # Hide unused subplots
        for i in range(n_chains, len(axes)):
            axes[i].set_visible(False)

        # Add main title if provided (but keep minimal per publication style)
        if title:
            fig.suptitle(title, fontsize=PUBLICATION_FONTS['title'] + 2, fontweight='bold', y=0.95)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        print(f"✅ Multi-panel RMSF plot saved: {output_path}")
        return True

    except Exception as e:
        print(f"Error creating multi-panel RMSF plot: {e}")
        return False


def _plot_rmsf_single_panel(sorted_chains: List[Tuple[str, Tuple[np.ndarray, np.ndarray]]],
                             output_path: str,
                             title: str,
                             figsize: Optional[Tuple[float, float]]) -> bool:
    """Create single-panel RMSF plot with all chains."""
    try:
        # Auto-calculate figure size if not provided
        if figsize is None:
            figsize = (12, 8)

        fig, ax = plt.subplots(figsize=figsize)

        # Colors for different chain types
        n_chains = len(sorted_chains)
        colors = plt.cm.tab10(np.linspace(0, 1, n_chains))

        for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(sorted_chains):
            color = colors[i]

            # Plot RMSF
            ax.plot(residue_ids, rmsf_values, color=color, linewidth=2,
                    alpha=0.8, label=chain_name)

        # Formatting
        ax.set_xlabel('Residue Number', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
        ax.set_ylabel('RMSF (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=PUBLICATION_FONTS['legend'])
        ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        print(f"✅ Single-panel RMSF plot saved: {output_path}")
        return True

    except Exception as e:
        print(f"Error creating single-panel RMSF plot: {e}")
        return False


def _plot_rmsf_separate_files(sorted_chains: List[Tuple[str, Tuple[np.ndarray, np.ndarray]]],
                               output_path: str,
                               title: str,
                               figsize: Optional[Tuple[float, float]]) -> bool:
    """Create separate RMSF plots for each chain."""
    try:
        output_base = Path(output_path)
        success_count = 0

        # Auto-calculate figure size if not provided
        if figsize is None:
            figsize = (10, 6)

        for i, (chain_name, (rmsf_values, residue_ids)) in enumerate(sorted_chains):
            try:
                fig, ax = plt.subplots(figsize=figsize)

                # Determine color based on chain type
                if 'Nucleic' in chain_name:
                    color = 'red'
                else:
                    color = 'blue'

                # Plot RMSF
                ax.plot(residue_ids, rmsf_values, color=color, linewidth=2, alpha=0.8)
                ax.fill_between(residue_ids, rmsf_values, alpha=0.3, color=color)

                # Calculate statistics
                mean_rmsf = np.mean(rmsf_values)
                std_rmsf = np.std(rmsf_values)
                # Print statistics (removed from plot to avoid hiding content)
                print(f"  📊 {chain_name} RMSF: Mean = {mean_rmsf:.2f} ± {std_rmsf:.2f} Å")

                # Formatting
                ax.set_xlabel('Residue Number', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
                ax.set_ylabel('RMSF (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
                ax.grid(True, alpha=0.3)
                ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])

                # Save individual file
                chain_output = output_base.parent / f"{output_base.stem}_{chain_name.lower().replace(' ', '_')}{output_base.suffix}"
                plt.tight_layout()
                plt.savefig(chain_output, dpi=300, bbox_inches='tight',
                           facecolor='white', edgecolor='none')
                plt.close()

                print(f"✅ Individual RMSF plot saved: {chain_output}")
                success_count += 1

            except Exception as e:
                print(f"Error creating plot for {chain_name}: {e}")
                continue

        print(f"✅ Created {success_count}/{len(sorted_chains)} individual RMSF plots")
        return success_count > 0

    except Exception as e:
        print(f"Error creating separate RMSF files: {e}")
        return False


def plot_multi_trajectory_rmsf_comprehensive(multi_traj_data: Dict[str, Dict],
                                           output_path: str,
                                           title: str = "",
                                           figsize: Optional[Tuple[float, float]] = None,
                                           plot_individual_trajectories: bool = True,
                                           save_data: bool = True) -> bool:
    """
    Create comprehensive multi-trajectory RMSF plot showing all chains and trajectories.

    Parameters
    ----------
    multi_traj_data : Dict[str, Dict]
        Multi-trajectory RMSF data from calculate_multi_trajectory_rmsf()
    output_path : str
        Path to save the plot
    title : str, optional
        Plot title (empty by default for publication style)
    figsize : Tuple[float, float], optional
        Figure size (width, height)
    plot_individual_trajectories : bool
        Whether to show individual trajectory lines
    save_data : bool
        Whether to save data as CSV files

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        apply_publication_style()

        if not multi_traj_data:
            print("No multi-trajectory RMSF data provided for plotting")
            return False

        n_chains = len(multi_traj_data)
        if n_chains == 0:
            print("No chain data available for multi-trajectory RMSF plotting")
            return False

        # Auto-calculate figure size for publication quality
        if figsize is None:
            if n_chains == 1:
                figsize = (10, 6)
            elif n_chains == 2:
                figsize = (12, 5)
            else:
                figsize = (15, 8)

        # Create subplots
        if n_chains == 1:
            fig, axes = plt.subplots(1, 1, figsize=figsize)
            axes = [axes]
        elif n_chains == 2:
            fig, axes = plt.subplots(1, 2, figsize=figsize)
        else:
            # For more chains, use a grid layout
            rows = (n_chains + 1) // 2
            cols = 2
            fig, axes = plt.subplots(rows, cols, figsize=figsize)
            axes = axes.flatten()

        # Colors for different molecule types
        protein_color = '#4472C4'  # Professional blue
        nucleic_color = '#E15759'  # Professional red/coral
        individual_alpha = 0.3
        mean_alpha = 0.9

        print(f"📊 Creating comprehensive multi-trajectory RMSF plot for {n_chains} chains")

        chain_names = sorted(multi_traj_data.keys())

        for i, chain_name in enumerate(chain_names):
            ax = axes[i]
            chain_data = multi_traj_data[chain_name]

            # Determine chain type and colors
            is_nucleic = any(keyword in chain_name.lower() for keyword in ['nucleic', 'rna', 'dna'])

            if is_nucleic:
                main_color = nucleic_color
                x_label = "Residue Number"
                print(f"  📊 Plotting nucleic acid chain: {chain_name}")
            else:
                main_color = protein_color
                x_label = "Residue Number"
                print(f"  📊 Plotting protein chain: {chain_name}")

            # Get combined data
            combined_data = chain_data["combined"]
            mean_rmsf = combined_data["mean_rmsf"]
            std_rmsf = combined_data["std_rmsf"]
            residue_ids = combined_data["residue_ids"]

            # Plot individual trajectories if requested
            if plot_individual_trajectories and "trajectories" in chain_data:
                traj_names = sorted(chain_data["trajectories"].keys())
                for j, traj_name in enumerate(traj_names):
                    traj_rmsf, traj_residues = chain_data["trajectories"][traj_name]

                    # Ensure residue IDs match (they should if everything is working correctly)
                    if len(traj_rmsf) == len(residue_ids):
                        ax.plot(residue_ids, traj_rmsf,
                               color=main_color, alpha=individual_alpha,
                               linewidth=1, label=f"Trajectory {j+1}" if i == 0 else None)

            # Plot mean with error bars/band
            ax.plot(residue_ids, mean_rmsf, color=main_color, linewidth=3,
                   alpha=mean_alpha, label="Mean" if i == 0 else None)

            # Add standard deviation as filled area
            ax.fill_between(residue_ids, mean_rmsf - std_rmsf, mean_rmsf + std_rmsf,
                           alpha=0.3, color=main_color, label="±1 SD" if i == 0 else None)

            # Add statistics
            stats = chain_data["statistics"]
            overall_mean = stats["overall_mean"]
            overall_std = stats["overall_std"]
            n_trajectories = stats["n_trajectories"]

            # Print statistics (removed from plot to avoid hiding content)
            print(f"  📊 {chain_name} RMSF: Mean = {overall_mean:.1f} ± {overall_std:.2f} Å ({n_trajectories} trajectories)")

            # Panel title
            clean_chain_name = chain_name.replace("Protein ", "").replace("Nucleic ", "Nucleic ")
            ax.set_title(clean_chain_name, fontsize=PUBLICATION_FONTS['title'],
                        fontweight='bold', pad=10)

            # Formatting
            ax.set_xlabel(x_label, fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
            ax.set_ylabel('RMSF (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])
            ax.set_ylim(bottom=0)

            # Add legend only to first subplot
            if i == 0 and plot_individual_trajectories:
                ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

        # Hide unused subplots if any
        for i in range(n_chains, len(axes)):
            axes[i].set_visible(False)

        # Add main title if provided
        if title:
            fig.suptitle(title, fontsize=PUBLICATION_FONTS['title'] + 2,
                        fontweight='bold', y=0.95)

        plt.tight_layout()

        # Save plot
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        print(f"✅ Comprehensive multi-trajectory RMSF plot saved: {output_path}")

        # Save data to CSV if requested
        if save_data:
            output_path_obj = Path(output_path)
            csv_path = output_path_obj.with_suffix('.csv')

            # Create comprehensive data table
            all_data = []
            for chain_name, chain_data in multi_traj_data.items():
                combined_data = chain_data["combined"]
                stats = chain_data["statistics"]

                for i, (rmsf_val, res_id) in enumerate(zip(combined_data["mean_rmsf"],
                                                         combined_data["residue_ids"])):
                    all_data.append({
                        'Chain': chain_name,
                        'Residue_ID': res_id,
                        'Mean_RMSF': rmsf_val,
                        'Std_RMSF': combined_data["std_rmsf"][i],
                        'Overall_Mean': stats["overall_mean"],
                        'Overall_Std': stats["overall_std"],
                        'N_Trajectories': stats["n_trajectories"]
                    })

            import pandas as pd
            df = pd.DataFrame(all_data)
            df.to_csv(csv_path, index=False)
            print(f"✅ RMSF data saved to CSV: {csv_path}")

        return True

    except Exception as e:
        print(f"Error creating comprehensive multi-trajectory RMSF plot: {e}")
        import traceback
        traceback.print_exc()
        return False