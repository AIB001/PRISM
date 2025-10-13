#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact analysis plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional
from pathlib import Path

# Import publication style
from .publication_utils import apply_publication_style, PUBLICATION_COLORS, PUBLICATION_FONTS, get_standard_figsize
# Import residue formatting utility
from ...utils.residue import format_residue_list

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
        trajectories = ['Repeat 1', 'Repeat 2', 'Repeat 3']
        colors = PUBLICATION_COLORS['example'][:3]  # Use example colors

        # Calculate means and standard errors
        means = []
        stderrs = []

        for residue in residues:
            values = [contact_data[residue].get(f'Repeat{i}', 0.0) for i in [1, 2, 3]]
            means.append(np.mean(values))
            stderrs.append(np.std(values) / np.sqrt(3))  # Standard error

        # Sort by mean contact probability
        sorted_indices = np.argsort(means)[::-1]
        residues_sorted = [residues[i] for i in sorted_indices]
        means_sorted = [means[i] for i in sorted_indices]
        stderrs_sorted = [stderrs[i] for i in sorted_indices]

        # Format residue names for display
        residues_sorted_display = format_residue_list(residues_sorted, residue_format)

        fig, ax = plt.subplots(figsize=get_standard_figsize('single'))

        # Create grouped bars
        x = np.arange(len(residues_sorted))
        width = 0.25

        for i, traj in enumerate(['Repeat1', 'Repeat2', 'Repeat3']):
            values = [contact_data[residues_sorted[j]].get(traj, 0.0) for j in range(len(residues_sorted))]
            bars = ax.bar(x + (i-1)*width, values, width,
                         label=trajectories[i], color=colors[i],
                         alpha=0.9, edgecolor='#457B9D', linewidth=0.5)

        # Formatting - apply_publication_style() handles fontweight='bold' automatically
        ax.set_xlabel('Amino Acid Residues')
        ax.set_ylabel('Contact Probability (%)')
        ax.set_xticks(x)
        ax.set_xticklabels(residues_sorted_display, rotation=45, ha='center')
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


def plot_contact_heatmap_annotated(contact_data: Dict[str, Dict[str, float]],
                                  output_path: str,
                                  title: str = "",
                                  separate_panels: bool = False,
                                  residue_format: str = "1letter") -> bool:
    """
    Plot annotated heatmap for contact analysis.

    Parameters
    ----------
    contact_data : dict
        Dictionary with residue names as keys and trajectory data as nested dict
    output_path : str
        Path to save the figure
    title : str
        Main plot title
    separate_panels : bool
        If True, also save individual panel figure
    residue_format : str
        Amino acid display format: "1letter" (e.g., D618) or "3letter" (e.g., ASP618)

    Returns
    -------
    bool
        True if successful
    """
    try:
        print("🔥 Contact probability heatmap with annotations showing residue-trajectory interactions")
        apply_publication_style()
        from .publication_utils import PUBLICATION_FONTS

        # Prepare data for heatmap
        residues = list(contact_data.keys())
        trajectories = ['Repeat 1', 'Repeat 2', 'Repeat 3']

        # Create data matrix
        data_matrix = []
        for traj_num in [1, 2, 3]:
            row = [contact_data[res].get(f'Repeat{traj_num}', 0.0) for res in residues]
            data_matrix.append(row)

        data_matrix = np.array(data_matrix)

        # Use wider figsize for more residues (15+)
        if len(residues) > 12:
            figsize = get_standard_figsize('wide')  # (16, 6)
        elif len(residues) > 6:
            figsize = get_standard_figsize('horizontal')  # (12, 6)
        else:
            figsize = get_standard_figsize('single')  # (8, 6)

        fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap
        im = ax.imshow(data_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=100)

        # Add colorbar with proper font size
        cbar = plt.colorbar(im, label='Contact Probability (%)', shrink=0.8)
        cbar.set_label('Contact Probability (%)', fontweight='bold',
                      fontsize=PUBLICATION_FONTS['colorbar'])
        cbar.ax.tick_params(labelsize=PUBLICATION_FONTS['tick_label'])

        # Format residue names for display
        residues_display = format_residue_list(residues, residue_format)

        # Set ticks and labels
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(trajectories)
        ax.set_xticks(range(len(residues)))
        ax.set_xticklabels(residues_display, rotation=45, ha='center')

        # Add text annotations with proper font size
        for i in range(3):
            for j in range(len(residues)):
                text_color = 'white' if data_matrix[i, j] > 50 else 'black'
                ax.text(j, i, f'{data_matrix[i, j]:.1f}',
                       ha='center', va='center', color=text_color,
                       fontweight='bold', fontsize=PUBLICATION_FONTS['value_text'])

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()
        return True

    except Exception as e:
        print(f"Error in contact heatmap: {e}")
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


def plot_key_residue_contact_violin(contact_data: Dict[str, Dict[str, float]],
                                  output_path: str,
                                  title: str = "Key Residue Contact Analysis",
                             residue_format: str = "1letter") -> bool:
    """
    Plot violin plots for key residue contacts with RemTP ligand using seaborn.

    Creates violin plots showing contact probability distributions for key residues
    interacting with the LIG ligand, with proper tick alignment.

    Parameters
    ----------
    contact_data : dict
        Dictionary with residue names as keys and trajectory data as nested dict
        Format: {'ASP618': {'Repeat1': 98.1, 'Repeat2': 100.0, 'Repeat3': 65.3}, ...}
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
        print("🎻 Key residue contact violin plot showing probability distributions across trajectories")
        apply_publication_style()

        # Use all residues in the contact_data (dynamic list)
        available_residues = list(contact_data.keys())

        if not available_residues:
            print("No key residue data available for violin plot")
            return False

        # Prepare data for seaborn violin plot
        import pandas as pd
        plot_data = []

        for residue in available_residues:
            for i in [1, 2, 3]:
                value = contact_data[residue].get(f'Repeat{i}', 0.0)
                if value > 0:  # Only include non-zero values
                    plot_data.append({
                        'Residue': residue,
                        'Trajectory': f'Repeat {i}',
                        'Contact_Probability': value
                    })
                else:
                    # Add small value for visualization if all are zero
                    plot_data.append({
                        'Residue': residue,
                        'Trajectory': f'Repeat {i}',
                        'Contact_Probability': 0.1
                    })

        if not plot_data:
            return False

        df = pd.DataFrame(plot_data)

        # Use wider figsize for more residues (15+)
        if len(available_residues) > 12:
            figsize = get_standard_figsize('wide')  # (16, 6)
        elif len(available_residues) > 6:
            figsize = get_standard_figsize('horizontal')  # (12, 6)
        else:
            figsize = get_standard_figsize('single')  # (8, 6)

        # Create violin plot with seaborn
        fig, ax = plt.subplots(figsize=figsize)

        # Define colors
        colors = ['#A8DADC', '#F1FAEE', '#E9C46A', '#D4A574', '#E76F51',
                 '#2A9D8F', '#F4A261', '#E76F51', '#457B9D']

        # Create violin plot with proper hue assignment to avoid deprecation warning
        violin_parts = sns.violinplot(
            data=df,
            x='Residue',
            y='Contact_Probability',
            hue='Residue',
            palette=colors[:len(available_residues)],
            inner=None,  # No inner representation to avoid clutter
            legend=False,
            ax=ax
        )

        # Add individual points
        sns.stripplot(
            data=df,
            x='Residue',
            y='Contact_Probability',
            size=8,
            alpha=0.7,
            color='black',
            ax=ax
        )

        # Fix tick alignment - center labels with ticks
        ax.set_xticks(range(len(available_residues)))
        # Format residue names for display
        available_residues_display = format_residue_list(available_residues, residue_format)

        ax.set_xticklabels(available_residues_display, rotation=45, ha='center')

        # Formatting - apply_publication_style() controls font sizes and bold automatically
        ax.set_ylabel('Contact Probability (%)')
        ax.set_xlabel('Key Residues')
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#FAFAFA')

        # Set y-axis limits
        ax.set_ylim(0, max(df['Contact_Probability']) * 1.1)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in key residue contact violin plot: {e}")
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

        # Calculate means and standard errors for key residues
        means = []
        stderrs = []

        for residue in available_residues:
            values = [contact_data[residue].get(f'Repeat{i}', 0.0) for i in [1, 2, 3]]
            means.append(np.mean(values))
            stderrs.append(np.std(values) / np.sqrt(3))  # Standard error

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

        # Add error bars
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

                # Increase y-offset to avoid overlap with error bar caps
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + se + 5,
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


def plot_contact_numbers_timeseries(contact_timeseries_data: Dict,
                                   output_path: str,
                                   title: str = "Contact Numbers Time Series") -> bool:
    """
    Plot contact numbers time series for multiple trajectories.

    Parameters
    ----------
    contact_timeseries_data : dict
        Dictionary with trajectory names as keys and timeseries data as values
        Format: {'Repeat 1': {'contact_numbers': array, 'times': array}, ...}
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
        print("📈 Contact numbers 4-panel timeseries: overlay, distribution, averages, and smoothed trends")
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=get_standard_figsize('quad'))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

        # Plot 1: Time series overlay
        ax1 = axes[0, 0]
        for i, (traj_name, data) in enumerate(contact_timeseries_data.items()):
            contact_numbers = data['contact_numbers']
            times = data['times']
            ax1.plot(times, contact_numbers, color=colors[i % 3],
                    label=traj_name, linewidth=1.5, alpha=0.8)

        ax1.set_xlabel('Time (ns)')
        ax1.set_ylabel('Contact Number')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Contact number distribution
        ax2 = axes[0, 1]
        all_contacts = []
        for data in contact_timeseries_data.values():
            all_contacts.extend(data['contact_numbers'])

        ax2.hist(all_contacts, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        ax2.set_xlabel('Contact Number')
        ax2.set_ylabel('Frequency')
        ax2.grid(True, alpha=0.3)

        # Plot 3: Average contact numbers per trajectory
        ax3 = axes[1, 0]
        traj_names = list(contact_timeseries_data.keys())
        avg_contacts = []
        std_contacts = []

        for data in contact_timeseries_data.values():
            contacts = data['contact_numbers']
            avg_contacts.append(np.mean(contacts))
            std_contacts.append(np.std(contacts))

        bars = ax3.bar(traj_names, avg_contacts, yerr=std_contacts,
                      color=colors[:len(traj_names)], alpha=0.7, capsize=5)
        ax3.set_ylabel('Average Contact Number')
        ax3.grid(True, alpha=0.3, axis='y')

        # Plot 4: Running average
        ax4 = axes[1, 1]
        window_size = 50

        for i, (traj_name, data) in enumerate(contact_timeseries_data.items()):
            contacts = data['contact_numbers']
            times = data['times']

            if len(contacts) > window_size:
                running_avg = np.convolve(contacts, np.ones(window_size)/window_size, mode='valid')
                time_avg = times[:len(running_avg)]
                ax4.plot(time_avg, running_avg, color=colors[i % 3],
                        label=f'{traj_name} (smoothed)', linewidth=2)

        ax4.set_xlabel('Time (ns)')
        ax4.set_ylabel('Running Average Contact Number')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in contact numbers timeseries plot: {e}")
        return False


def plot_residue_contact_numbers_violin(residue_contact_data: Dict,
                                       output_path: str,
                                       title: str = "Residue Contact Numbers",
                                       residue_format: str = "1letter") -> bool:
    """
    Plot violin plots for residue contact numbers (normalized by atom count).

    Parameters
    ----------
    residue_contact_data : dict
        Dictionary with trajectory names as keys and residue data as nested dict
        Format: {'Repeat 1': {'ASP618': {'normalized': array, 'n_atoms': int}, ...}, ...}
    output_path : str
        Path to save the figure
    title : str
        Main plot title
    residue_format : str
        Amino acid display format: "1letter" (e.g., D618) or "3letter" (e.g., ASP618)

    Returns
    -------
    bool
        True if successful
    """
    try:
        print("🎻 Residue contact numbers violin plot showing normalized distributions per residue")
        apply_publication_style()
        import pandas as pd

        # Prepare data for violin plot
        plot_data = []
        all_residues = set()

        for traj_name, residue_data in residue_contact_data.items():
            all_residues.update(residue_data.keys())

        all_residues = sorted(list(all_residues))

        for traj_name, residue_data in residue_contact_data.items():
            for residue in all_residues:
                if residue in residue_data:
                    normalized_contacts = residue_data[residue]['normalized']
                    for contact_val in normalized_contacts:
                        plot_data.append({
                            'Residue': residue,
                            'Trajectory': traj_name,
                            'Normalized_Contact_Number': contact_val
                        })

        if not plot_data:
            return False

        df = pd.DataFrame(plot_data)

        fig, ax = plt.subplots(figsize=get_standard_figsize('single'))

        # Define colors
        colors = ['#A8DADC', '#F1FAEE', '#E9C46A', '#D4A574', '#E76F51',
                 '#2A9D8F', '#F4A261', '#E76F51', '#457B9D']

        # Create violin plot
        violin_parts = sns.violinplot(
            data=df,
            x='Residue',
            y='Normalized_Contact_Number',
            hue='Residue',
            palette=colors[:len(all_residues)],
            inner=None,  # No inner representation to avoid clutter
            legend=False,
            ax=ax
        )

        # Add strip plot for individual points
        sns.stripplot(
            data=df,
            x='Residue',
            y='Normalized_Contact_Number',
            size=4,
            alpha=0.6,
            color='black',
            ax=ax
        )

        # Formatting - apply_publication_style() controls font sizes and bold
        ax.set_ylabel('Normalized Contact Number')
        ax.set_xlabel('Key Residues')
        ax.set_xticks(range(len(all_residues)))
        # Format residue names for display
        all_residues_display = format_residue_list(all_residues, residue_format)

        ax.set_xticklabels(all_residues_display, rotation=45, ha='center')
        # Remove manual tick_params - let apply_publication_style() handle it
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#FAFAFA')

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in residue contact numbers violin plot: {e}")
        return False


def plot_contact_distances_violin(distance_data: Dict,
                                 output_path: str,
                                 title: str = "Contact Distance Distributions",
                                 distance_cutoff: float = -1.0,
                                 residue_format: str = "1letter") -> bool:
    """
    Plot violin plots for contact distance distributions.

    Parameters
    ----------
    distance_data : dict
        Dictionary with trajectory names as keys and residue distance data as nested dict
        Format: {'Repeat 1': {'ASP618': distances_array, ...}, ...}
    output_path : str
        Path to save the figure
    title : str
        Main plot title
    distance_cutoff : float, optional
        Maximum distance to include in plot (in Angstroms).
        Default: -1.0 (no cutoff, include all distances)
        Typical values: 6.0 Å for contact analysis
    residue_format : str
        Amino acid display format: "1letter" (e.g., D618) or "3letter" (e.g., ASP618)

    Returns
    -------
    bool
        True if successful
    """
    try:
        cutoff_msg = f"(cutoff: {distance_cutoff:.1f} Å)" if distance_cutoff > 0 else "(no cutoff)"
        print(f"🎻 Contact distance distributions violin plot {cutoff_msg}")
        apply_publication_style()
        import pandas as pd

        # Prepare data for violin plot
        plot_data = []
        all_residues = set()

        for traj_name, residue_data in distance_data.items():
            all_residues.update(residue_data.keys())

        all_residues = sorted(list(all_residues))

        for traj_name, residue_data in distance_data.items():
            for residue in all_residues:
                if residue in residue_data:
                    distances = residue_data[residue]

                    # Apply distance cutoff if specified
                    if distance_cutoff > 0:
                        distances = distances[distances <= distance_cutoff]

                    if len(distances) == 0:
                        continue

                    # Sample distances to avoid overloading plot
                    if len(distances) > 1000:
                        sample_indices = np.random.choice(len(distances), 1000, replace=False)
                        distances = distances[sample_indices]

                    for distance in distances:
                        plot_data.append({
                            'Residue': residue,
                            'Trajectory': traj_name,
                            'Contact_Distance': distance
                        })

        if not plot_data:
            return False

        df = pd.DataFrame(plot_data)

        fig, ax = plt.subplots(figsize=get_standard_figsize('distribution'))

        # Define colors
        colors = ['#A8DADC', '#F1FAEE', '#E9C46A', '#D4A574', '#E76F51',
                 '#2A9D8F', '#F4A261', '#E76F51', '#457B9D']

        # Create violin plot - no scattering by default
        violin_parts = sns.violinplot(
            data=df,
            x='Residue',
            y='Contact_Distance',
            hue='Residue',
            palette=colors[:len(all_residues)],
            inner=None,  # No inner representation to avoid scattering
            legend=False,
            ax=ax
        )

        # Formatting - apply_publication_style() controls font sizes and bold
        ax.set_ylabel('Contact Distance (Å)')
        ax.set_xlabel('Key Residues')
        ax.set_xticks(range(len(all_residues)))
        # Format residue names for display
        all_residues_display = format_residue_list(all_residues, residue_format)

        ax.set_xticklabels(all_residues_display, rotation=45, ha='center')
        # Remove manual tick_params - let apply_publication_style() handle it
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#FAFAFA')

        # Set y-axis limits dynamically based on data
        if distance_cutoff > 0:
            ax.set_ylim(0, min(distance_cutoff * 1.1, df['Contact_Distance'].max() * 1.1))
        else:
            ax.set_ylim(0, df['Contact_Distance'].max() * 1.1)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in contact distances violin plot: {e}")
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


def plot_residue_distance_timeseries(distance_timeseries: Dict,
                                    residue_name: str,
                                    output_path: str,
                                    contact_start: float = 3.0,
                                    contact_end: float = 4.5,
                                    smooth_window: int = 11) -> bool:
    """
    Plot distance timeseries for a single residue with contact thresholds.

    Parameters
    ----------
    distance_timeseries : dict
        Dictionary containing 'times' and residue distance data
        Format: {'times': time_array, 'ASP618': distance_array, ...}
    residue_name : str
        Residue name to plot (e.g., 'ASP618', 'ARG555')
    output_path : str
        Path to save the figure
    contact_start : float
        Contact start threshold in Angstroms (green line)
    contact_end : float
        Contact end threshold in Angstroms (red line)
    smooth_window : int
        Window size for smoothing (default 11)

    Returns
    -------
    bool
        True if successful
    """
    try:
        print(f"📈 Distance timeseries plot for {residue_name} showing binding/unbinding dynamics")
        apply_publication_style()

        if residue_name not in distance_timeseries:
            print(f"  ✗ Residue {residue_name} not found in timeseries data")
            return False

        times = distance_timeseries['times']
        distances = distance_timeseries[residue_name]

        # Smooth the data
        if len(distances) > smooth_window:
            smoothed_distances = np.convolve(
                distances,
                np.ones(smooth_window) / smooth_window,
                mode='valid'
            )
            smoothed_times = times[smooth_window-1:]
        else:
            smoothed_distances = distances
            smoothed_times = times

        fig, ax = plt.subplots(figsize=get_standard_figsize('single'))

        # Plot distance timeseries
        ax.plot(smoothed_times, smoothed_distances, 'b-', linewidth=2,
               label=f'Distance (window={smooth_window})')

        # Add contact threshold lines
        ax.axhline(y=contact_start, color='g', linestyle='--', linewidth=2,
                  label='Contact start')
        ax.axhline(y=contact_end, color='r', linestyle='--', linewidth=2,
                  label='Contact end')

        # Formatting
        ax.set_xlabel('Time (ns)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
        ax.set_ylabel('Distance (Å)', fontsize=PUBLICATION_FONTS['axis_label'], fontweight='bold')
        ax.tick_params(axis='both', labelsize=PUBLICATION_FONTS['tick_label'])
        ax.legend(loc='best', fontsize=PUBLICATION_FONTS['legend'], framealpha=0.9)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()

        print(f"  ✅ Distance timeseries plot saved: {output_path}")
        return True

    except Exception as e:
        print(f"Error in distance timeseries plot: {e}")
        import traceback
        traceback.print_exc()
        return False