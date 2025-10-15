#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact analysis and distribution plots.

Refactored from contact_plots.py for better modularity.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional
from pathlib import Path

# Import publication style
from ..publication_utils import apply_publication_style, PUBLICATION_COLORS, PUBLICATION_FONTS, get_standard_figsize
# Import residue formatting utility
from ....utils.residue import format_residue_list

def plot_contact_analysis(contact_results: Dict,

                         output_path: str,

                         title: str = "",

                         key_residues: Optional[Dict] = None) -> bool:

    """

    Plot contact analysis results for multiple trajectories.



    Parameters

    ----------

    contact_results : dict

        Dictionary with trajectory names as keys and contact data as values

    output_path : str

        Path to save the plot

    title : str

        Plot title

    key_residues : dict, optional

        Dictionary of key residues for analysis



    Returns

    -------

    bool

        True if successful, False otherwise

    """

    try:

        print("📊 4-panel contact analysis: time series, distribution, trajectory averages, and running averages")

        # Apply publication style

        apply_publication_style()



        fig, axes = plt.subplots(2, 2, figsize=get_standard_figsize('quad'))

        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']



        # Plot 1: Contact time series (convert boolean to contact count)

        ax1 = axes[0, 0]

        for i, (traj_name, contacts) in enumerate(contact_results.items()):

            if len(contacts) > 0:

                time = np.arange(len(contacts)) * 0.5  # Assume 0.5ns interval

                # Convert boolean contacts to integers (0 or 1)

                contact_counts = np.array(contacts, dtype=int)

                ax1.plot(time, contact_counts, color=colors[i % 3],

                        label=traj_name, linewidth=2, alpha=0.8)



        ax1.set_xlabel('Time (ns)', fontfamily='Times New Roman')

        ax1.set_ylabel('Contact Count', fontfamily='Times New Roman')

        ax1.legend()

        ax1.grid(True, alpha=0.3)



        # Plot 2: Contact distribution

        ax2 = axes[0, 1]

        all_contacts = []

        for contacts in contact_results.values():

            if len(contacts) > 0:

                # Convert boolean contacts to integers for histogram

                contact_counts = np.array(contacts, dtype=int)

                all_contacts.extend(contact_counts)



        if all_contacts:

            ax2.hist(all_contacts, bins=[0, 0.5, 1], alpha=0.7, color='skyblue', edgecolor='black')

            ax2.set_xlabel('Contact', fontfamily='Times New Roman')

            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')

            ax2.set_xticks([0, 1])

            ax2.grid(True, alpha=0.3)



        # Plot 3: Average contacts per trajectory

        ax3 = axes[1, 0]

        traj_names = []

        avg_contacts = []

        std_contacts = []



        for traj_name, contacts in contact_results.items():

            if len(contacts) > 0:

                # Convert boolean to float for percentage calculation

                contact_percentage = np.mean(np.array(contacts, dtype=float)) * 100

                traj_names.append(traj_name)

                avg_contacts.append(contact_percentage)

                std_contacts.append(np.std(np.array(contacts, dtype=float)) * 100)



        if traj_names:

            bars = ax3.bar(traj_names, avg_contacts, yerr=std_contacts,

                          color=colors[:len(traj_names)], alpha=0.7, capsize=5)

            ax3.set_ylabel('Contact (%)', fontfamily='Times New Roman')

            ax3.grid(True, alpha=0.3, axis='y')



        # Plot 4: Contact stability (running average)

        ax4 = axes[1, 1]

        window_size = min(50, len(list(contact_results.values())[0]) // 10) if contact_results else 10



        for i, (traj_name, contacts) in enumerate(contact_results.items()):

            if len(contacts) > window_size:

                # Convert boolean to float for running average calculation

                contact_floats = np.array(contacts, dtype=float)

                running_avg = np.convolve(contact_floats, np.ones(window_size)/window_size, mode='valid') * 100

                time = np.arange(len(running_avg)) * 0.5

                ax4.plot(time, running_avg, color=colors[i % 3],

                        label=f'{traj_name} (avg)', linewidth=2)



        ax4.set_xlabel('Time (ns)', fontfamily='Times New Roman')

        ax4.set_ylabel('Contact (%)', fontfamily='Times New Roman')

        ax4.legend()

        ax4.grid(True, alpha=0.3)



        plt.tight_layout()

        plt.savefig(output_path, dpi=300, bbox_inches='tight')

        plt.close()



        return True



    except Exception as e:

        print(f"Error in contact plotting: {e}")

        return False







def plot_grouped_contact_bars(contact_data: Dict[str, Dict[str, float]],

                             output_path: str,

                             title: str = "",

                             separate_panels: bool = False,

                             residue_format: str = "1letter") -> bool:

    """

    Plot grouped bar charts for contact analysis across trajectories.



    Parameters

    ----------

    contact_data : dict

        Dictionary with residue names as keys and trajectory data as nested dict

        Format: {'ASP618': {'Repeat1': 98.1, 'Repeat2': 100.0, 'Repeat3': 65.3}, ...}

    output_path : str

        Path to save the main figure

    title : str

        Main plot title

    separate_panels : bool

        If True, also save individual panel figures

    residue_format : str

        Amino acid display format: "1letter" (e.g., D618) or "3letter" (e.g., ASP618)



    Returns

    -------

    bool

        True if successful

    """

    try:

        print("📊 Grouped contact probability bars showing residue-specific interactions across trajectories")

        apply_publication_style()



        # Extract data for plotting

        residues = list(contact_data.keys())



        # Auto-detect number of trajectories from data

        all_traj_keys = set()

        for residue_data in contact_data.values():

            all_traj_keys.update(residue_data.keys())



        # Detect trajectory format (Repeat1/Repeat2/Repeat3 or single trajectory)

        repeat_keys = sorted([k for k in all_traj_keys if k.startswith('Repeat')])

        n_trajectories = len(repeat_keys) if repeat_keys else 1



        if n_trajectories == 1 and not repeat_keys:

            # Single trajectory with custom name (shouldn't happen after test fix)

            trajectories = [list(all_traj_keys)[0]]

            trajectory_nums = [1]

        else:

            # Multiple trajectories or properly formatted single trajectory

            trajectory_nums = [int(k.replace('Repeat', '')) for k in repeat_keys]

            trajectories = [f'Repeat {i}' for i in trajectory_nums]



        colors = PUBLICATION_COLORS['example'][:n_trajectories]



        # Calculate means and standard errors

        means = []

        stderrs = []



        for residue in residues:

            values = [contact_data[residue].get(f'Repeat{i}', 0.0) for i in trajectory_nums]

            means.append(np.mean(values))

            stderrs.append(np.std(values) / np.sqrt(n_trajectories) if n_trajectories > 1 else 0.0)



        # Sort by mean contact probability

        sorted_indices = np.argsort(means)[::-1]

        residues_sorted = [residues[i] for i in sorted_indices]

        means_sorted = [means[i] for i in sorted_indices]

        stderrs_sorted = [stderrs[i] for i in sorted_indices]



        # Format residue names for display

        residues_sorted_display = format_residue_list(residues_sorted, residue_format)



        fig, ax = plt.subplots(figsize=get_standard_figsize('single'))



        # Create grouped bars (adapt width based on number of trajectories)

        x = np.arange(len(residues_sorted))

        width = 0.25 if n_trajectories > 1 else 0.6  # Wider bars for single trajectory



        for i, traj_num in enumerate(trajectory_nums):

            traj_key = f'Repeat{traj_num}'

            values = [contact_data[residues_sorted[j]].get(traj_key, 0.0) for j in range(len(residues_sorted))]



            # Center bars for single trajectory, offset for multiple

            if n_trajectories == 1:

                x_offset = 0

            else:

                x_offset = (i - (n_trajectories-1)/2) * width



            bars = ax.bar(x + x_offset, values, width,

                         label=trajectories[i], color=colors[i],

                         alpha=0.9, edgecolor='#457B9D', linewidth=0.5)



        # Formatting - apply_publication_style() handles fontweight='bold' automatically

        ax.set_xlabel('Amino Acid Residues')

        ax.set_ylabel('Contact Probability (%)')

        ax.set_xticks(x)

        ax.set_xticklabels(residues_sorted_display, rotation=45, ha='center')



        # Only show legend for multiple trajectories

        if n_trajectories > 1:

            ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)



        ax.set_ylim(0, 115)

        ax.set_facecolor('#FAFAFA')



        plt.tight_layout()

        plt.savefig(output_path, dpi=300, bbox_inches='tight',

                   facecolor='white', edgecolor='none')



        # Save separate panels if requested (only if filename doesn't already contain the suffix)

        if separate_panels:

            base_path = Path(output_path)

            # 避免重复添加后缀

            if "_grouped_bars" not in base_path.stem:

                separate_path = base_path.parent / f"{base_path.stem}_grouped_bars{base_path.suffix}"

                plt.savefig(separate_path, dpi=300, bbox_inches='tight',

                           facecolor='white', edgecolor='none')



        plt.close()

        return True



    except Exception as e:

        print(f"Error in grouped contact bars: {e}")

        return False







def plot_key_residue_contacts(key_residue_results: Dict,

                            output_path: str,

                            title: str = "") -> bool:

    """

    Plot contact analysis for key residues.



    Parameters

    ----------

    key_residue_results : dict

        Nested dict: {traj_name: {residue_id: contact_data}}

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

        apply_publication_style()



        # Get all unique residues

        all_residues = set()

        for traj_data in key_residue_results.values():

            all_residues.update(traj_data.keys())

        all_residues = sorted(list(all_residues))



        if not all_residues:

            return False



        n_residues = len(all_residues)

        n_cols = min(4, n_residues)

        n_rows = (n_residues + n_cols - 1) // n_cols



        # 使用标准尺寸计算多panel尺寸

        base_width, base_height = get_standard_figsize('single')

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(base_width*n_cols*0.6, base_height*n_rows*0.6))

        if n_residues == 1:

            axes = [axes]

        elif n_rows == 1:

            axes = [axes]

        else:

            axes = axes.flatten()



        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']



        for i, residue in enumerate(all_residues):

            ax = axes[i]



            for j, (traj_name, traj_data) in enumerate(key_residue_results.items()):

                if residue in traj_data and len(traj_data[residue]) > 0:

                    contacts = traj_data[residue]

                    time = np.arange(len(contacts)) * 0.5

                    ax.plot(time, contacts, color=colors[j % 3],

                           label=traj_name, linewidth=1.5, alpha=0.8)



            ax.set_xlabel('Time (ns)', fontfamily='Times New Roman')

            ax.set_ylabel('Contacts', fontfamily='Times New Roman')

            ax.legend()

            ax.grid(True, alpha=0.3)



        # Hide unused subplots

        for i in range(n_residues, len(axes)):

            axes[i].set_visible(False)



        plt.tight_layout()

        plt.savefig(output_path, dpi=300, bbox_inches='tight')

        plt.close()



        return True



    except Exception as e:

        print(f"Error in key residue contact plotting: {e}")

        return False







def plot_key_residue_contact_distribution(contact_data: Dict[str, Dict[str, float]],

                                        output_path: str,

                                        title: str = "Key Residue Contact Distribution",

                             residue_format: str = "1letter") -> bool:

    """

    Plot grouped bar charts for key residue contacts with RemTP ligand.



    Parameters

    ----------

    contact_data : dict

        Dictionary with residue names as keys and trajectory data as nested dict

    output_path : str

        Path to save the figure

    title : str

        Main plot title



    Returns

    -------

    bool

        True if successful

    """

    try:

        print("📈 Key residue contact distribution bar plot showing probability comparison across trajectories")

        apply_publication_style()



        # Use all residues in the contact_data (dynamic list)

        available_residues = list(contact_data.keys())



        if not available_residues:

            print("No key residue data available for distribution plot")

            return False



        # Auto-detect number of trajectories from data

        all_traj_keys = set()

        for residue_data in contact_data.values():

            all_traj_keys.update(residue_data.keys())



        repeat_keys = sorted([k for k in all_traj_keys if k.startswith('Repeat')])

        n_trajectories = len(repeat_keys) if repeat_keys else 1

        trajectory_nums = [int(k.replace('Repeat', '')) for k in repeat_keys] if repeat_keys else [1]



        # Calculate means and standard errors for key residues

        means = []

        stderrs = []



        for residue in available_residues:

            values = [contact_data[residue].get(f'Repeat{i}', 0.0) for i in trajectory_nums]

            means.append(np.mean(values))

            # Only calculate error bars for multiple trajectories

            if n_trajectories > 1:

                stderrs.append(np.std(values) / np.sqrt(n_trajectories))

            else:

                stderrs.append(0.0)



        # Sort by mean contact probability

        sorted_indices = np.argsort(means)[::-1]

        residues_sorted = [available_residues[i] for i in sorted_indices]

        means_sorted = [means[i] for i in sorted_indices]

        stderrs_sorted = [stderrs[i] for i in sorted_indices]



        # Use distribution figsize for better height proportions

        figsize = get_standard_figsize('distribution')  # (12, 8) - taller for better proportions



        fig, ax = plt.subplots(figsize=figsize)



        # Create bar chart with error bars

        colors = ['#A8DADC' if x > 60 else '#F1FAEE' if x > 30 else '#E9C46A' for x in means_sorted]

        bars = ax.bar(residues_sorted, means_sorted,

                     color=colors,

                     alpha=0.9, edgecolor='#457B9D', linewidth=0.8)



        # Add error bars only for multiple trajectories

        if n_trajectories > 1:

            ax.errorbar(residues_sorted, means_sorted, yerr=stderrs_sorted,

                       fmt='none', capsize=6, capthick=1.5, color='#457B9D', alpha=0.8)



        # Add value labels on bars - simple percentage format only

        for i, (bar, mean, se) in enumerate(zip(bars, means_sorted, stderrs_sorted)):

            if mean > 10:  # Only show labels for significant contacts

                # Simple percentage format - error already shown by error bars

                if mean >= 99.5:  # Essentially 100%

                    label_text = '100%'

                else:  # Show one decimal place for clarity

                    label_text = f'{mean:.1f}%'



                # Adjust y-offset based on whether error bars are present

                if n_trajectories > 1:

                    y_offset = bar.get_height() + se + 5  # Account for error bar

                else:

                    y_offset = bar.get_height() + 3  # Smaller offset without error bar



                ax.text(bar.get_x() + bar.get_width()/2, y_offset,

                       label_text, ha='center', va='bottom',

                       fontsize=PUBLICATION_FONTS['bar_annotation'], color='#457B9D')



        # Formatting with proper tick alignment - apply_publication_style() controls font sizes and bold

        ax.set_ylabel('Contact Probability (%)')

        ax.set_xlabel('Key Residues')



        # Fix tick alignment - center labels with tick positions

        ax.set_xticks(range(len(residues_sorted)))

        # Format residue names for display

        residues_sorted_display = format_residue_list(residues_sorted, residue_format)



        ax.set_xticklabels(residues_sorted_display, rotation=45, ha='center')

        # Remove manual tick_params - let apply_publication_style() handle it



        ax.set_ylim(0, max(means_sorted) * 1.3 if means_sorted else 1)

        ax.set_facecolor('#FAFAFA')

        ax.grid(True, alpha=0.3, axis='y')



        plt.tight_layout()

        plt.savefig(output_path, dpi=300, bbox_inches='tight',

                   facecolor='white', edgecolor='none')

        plt.close()



        return True



    except Exception as e:

        print(f"Error in key residue contact distribution plot: {e}")

        return False









def _create_half_violin(ax, data, position, color, alpha=1.0, width=0.35, side='right'):

    """Create half violin plot using KDE"""

    from scipy import stats



    # Calculate kernel density estimation

    density = stats.gaussian_kde(data)

    y_min, y_max = min(data), max(data)

    y_range = np.linspace(y_min, y_max, 100)

    density_values = density(y_range)



    # Normalize density values

    density_values = density_values / density_values.max() * width



    if side == 'right':

        # Right half: extend from center to right

        x_coords = np.concatenate([

            [position] * len(y_range),  # Center line

            position + density_values[::-1]  # Right contour

        ])

        y_coords = np.concatenate([y_range, y_range[::-1]])

    else:

        # Left half: extend from center to left

        x_coords = np.concatenate([

            position - density_values,  # Left contour

            [position] * len(y_range)  # Center line

        ])

        y_coords = np.concatenate([y_range, y_range[::-1]])



    # Draw half violin

    ax.fill(x_coords, y_coords, color=color, alpha=alpha, zorder=1)

    return ax









def plot_contact_distances_raincloud(distance_data: Dict,

                                    output_path: str,

                                    title: str = "",

                                    distance_cutoff: float = -1.0,

                             residue_format: str = "1letter") -> bool:

    """

    Plot raincloud plots for contact distance distributions.

    Follows reference implementation with half violin, box plot, and density-based scatter.



    Parameters

    ----------

    distance_data : dict

        Dictionary with residue names as keys and distance arrays as values

        Format: {'ASP618': distances_array, 'ASN691': distances_array, ...}

    output_path : str

        Path to save the figure

    title : str

        Main plot title (empty by default for publication)

    distance_cutoff : float, optional

        Maximum distance to include in plot (in Angstroms).

        Default: -1.0 (no cutoff, include all distances)

        Typical values: 6.0 Å for contact analysis



    Returns

    -------

    bool

        True if successful

    """

    try:

        cutoff_msg = f"(cutoff: {distance_cutoff:.1f} Å)" if distance_cutoff > 0 else "(no cutoff)"

        print(f"🌧️ Contact distance raincloud plot {cutoff_msg}")

        apply_publication_style()

        from scipy import stats



        # Apply distance cutoff to all residue data if specified

        filtered_distance_data = {}

        for residue, distances in distance_data.items():

            if distance_cutoff > 0:

                filtered_distances = distances[distances <= distance_cutoff]

            else:

                filtered_distances = distances



            if len(filtered_distances) > 0:

                filtered_distance_data[residue] = filtered_distances



        if not filtered_distance_data:

            print("  ⚠ No data available after applying cutoff")

            return False



        residues = sorted(filtered_distance_data.keys())

        # Format residue names for display

        residues_display = format_residue_list(residues, residue_format)



        fig, ax = plt.subplots(figsize=get_standard_figsize('distribution'))



        # Color scheme matching reference

        fill_colors = [

            '#E8F4F8', '#F0E8F4', '#E8F8F0', '#F8F0E8', '#F4E8E8',

            '#F5F5DC', '#E6E6FA', '#F0FFF0', '#FFF5EE', '#F0F8FF'

        ]

        edge_colors = [

            '#4A90A4', '#8E44AD', '#27AE60', '#E67E22', '#E74C3C',

            '#D2B48C', '#9370DB', '#2E8B57', '#CD853F', '#4682B4'

        ]



        x_positions = np.arange(len(residues))



        # Track max distance for y-axis scaling

        max_distance = 0



        for i, residue in enumerate(residues):

            residue_distances = filtered_distance_data[residue]

            max_distance = max(max_distance, np.max(residue_distances))



            if len(residue_distances) == 0:

                continue



            pos = x_positions[i]

            fill_color = fill_colors[i % len(fill_colors)]

            edge_color = edge_colors[i % len(edge_colors)]



            # 1. Half violin (right side)

            _create_half_violin(ax, residue_distances, pos, fill_color,

                              alpha=1.0, width=0.35, side='right')



            # 2. Box plot (center)

            q1, median, q3 = np.percentile(residue_distances, [25, 50, 75])

            iqr = q3 - q1

            lower_whisker = max(min(residue_distances), q1 - 1.5 * iqr)

            upper_whisker = min(max(residue_distances), q3 + 1.5 * iqr)



            box_width = 0.1

            box = plt.Rectangle((pos - box_width / 2, q1), box_width, q3 - q1,

                              facecolor='white', edgecolor=edge_color,

                              linewidth=2, alpha=0.7, zorder=3)

            ax.add_patch(box)



            # Median line

            ax.plot([pos - box_width / 2, pos + box_width / 2], [median, median],

                   color=edge_color, linewidth=3, zorder=4)



            # Whiskers

            ax.plot([pos, pos], [q3, upper_whisker], color=edge_color, linewidth=2, zorder=3)

            ax.plot([pos, pos], [q1, lower_whisker], color=edge_color, linewidth=2, zorder=3)

            ax.plot([pos - 0.02, pos + 0.02], [upper_whisker, upper_whisker],

                   color=edge_color, linewidth=2, zorder=3)

            ax.plot([pos - 0.02, pos + 0.02], [lower_whisker, lower_whisker],

                   color=edge_color, linewidth=2, zorder=3)



            # 3. Density-based scatter (left side)

            # Sample data if too many points

            max_points = 200

            if len(residue_distances) > max_points:

                sample_indices = np.random.choice(len(residue_distances), max_points, replace=False)

                sampled_distances = residue_distances[sample_indices]

            else:

                sampled_distances = residue_distances



            # Calculate density-based x-offset

            density = stats.gaussian_kde(sampled_distances)

            local_density = density(sampled_distances)

            max_offset = 0.25

            if local_density.max() > local_density.min():

                normalized_density = (local_density - local_density.min()) / (local_density.max() - local_density.min())

            else:

                normalized_density = np.zeros_like(local_density)



            x_offsets = -max_offset * normalized_density - 0.05

            x_offsets += np.random.normal(0, 0.03, len(sampled_distances))  # Add jitter

            x_scatter = pos + x_offsets



            ax.scatter(x_scatter, sampled_distances, s=15, alpha=0.6,

                      color=edge_color, edgecolors='white',

                      linewidth=0.5, zorder=2)



        # Formatting

        ax.set_ylabel('Contact Distance (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')

        ax.set_xlabel('Key Residues', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')

        ax.set_xticks(x_positions)

        ax.set_xticklabels(residues_display, rotation=45, ha='center', fontsize=PUBLICATION_FONTS['tick_label'])

        ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])

        ax.set_xlim(-0.6, len(residues) - 0.4)



        # Set y-axis limits dynamically based on data

        if distance_cutoff > 0:

            ax.set_ylim(0, min(distance_cutoff * 1.1, max_distance * 1.1))

        else:

            ax.set_ylim(0, max_distance * 1.1)



        # Grid and spine styling

        ax.grid(True, alpha=0.3, axis='y')

        ax.set_axisbelow(True)

        ax.spines['top'].set_visible(False)

        ax.spines['right'].set_visible(False)



        plt.tight_layout()

        plt.savefig(output_path, dpi=300, bbox_inches='tight',

                   facecolor='white', edgecolor='none')

        plt.close()



        return True



    except Exception as e:

        print(f"Error in contact distances raincloud plot: {e}")

        import traceback

        traceback.print_exc()

        return False







