#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact probability and distance violin plots.

Refactored from contact_plots.py for better modularity.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Optional
from pathlib import Path

# Import publication style
from ...core.publication_utils import apply_publication_style, PUBLICATION_COLORS, PUBLICATION_FONTS, get_standard_figsize
# Import residue formatting utility
from .....utils.residue import format_residue_list



def plot_comparison_contact_distances_violin(

    group1_distances: Dict[str, np.ndarray],

    group2_distances: Dict[str, np.ndarray],

    output_path: str,

    group1_name: str = "Group 1",

    group2_name: str = "Group 2",

    key_residues: List[str] = None,

    residue_format: str = "1letter"

) -> bool:

    """

    Plot comparison violin plots for contact distances between two groups.



    Creates split violin plots with Group 1 on the left and Group 2 on the right

    for each residue, enabling direct visual comparison of distance distributions.



    Parameters

    ----------

    group1_distances : dict

        Distance data for group 1 (e.g., RemTP)

        Format: {'ASP618': distances_array, 'ASN691': distances_array, ...}

        Arrays can be from single or multiple trajectories combined

    group2_distances : dict

        Distance data for group 2 (e.g., 1'-CF3)

        Format: {'ASP618': distances_array, 'ASN691': distances_array, ...}

    output_path : str

        Path to save the figure

    group1_name : str

        Display name for group 1 (e.g., "RemTP")

    group2_name : str

        Display name for group 2 (e.g., "1'-CF3")

    key_residues : list, optional

        List of specific residues to plot. If None, uses all common residues

    residue_format : str

        Amino acid display format: "1letter" (e.g., D618) or "3letter" (e.g., ASP618)



    Returns

    -------

    bool

        True if successful



    Examples

    --------

    >>> # Compare RemTP (3 trajs) vs 1'-CF3 (1 traj)

    >>> remtp_combined = {

    ...     'ASP618': np.concatenate([remtp_traj1['ASP618'], remtp_traj2['ASP618'], remtp_traj3['ASP618']]),

    ...     'ASN691': np.concatenate([remtp_traj1['ASN691'], remtp_traj2['ASN691'], remtp_traj3['ASN691']])

    ... }

    >>> cf3_distances = {'ASP618': cf3_array, 'ASN691': cf3_array}

    >>> plot_comparison_contact_distances_violin(

    ...     remtp_combined, cf3_distances,

    ...     'remtp_vs_1pcf3_comparison.png',

    ...     'RemTP', "1'-CF3",

    ...     key_residues=['ASP618', 'ASN691', 'SER759']

    ... )

    """

    try:

        print(f"🎻 Comparison violin plot: {group1_name} vs {group2_name}")

        apply_publication_style()

        import pandas as pd



        # Determine residues to plot

        if key_residues is None:

            # Use all common residues

            common_residues = sorted(set(group1_distances.keys()) & set(group2_distances.keys()))

        else:

            # Use specified residues (filter to those present in both groups)

            common_residues = [r for r in key_residues

                             if r in group1_distances and r in group2_distances]



        if not common_residues:

            print(f"  ⚠ No common residues found between {group1_name} and {group2_name}")

            return False



        print(f"  📊 Plotting {len(common_residues)} residues")



        # Prepare data for split violin plot

        plot_data = []



        for residue in common_residues:

            # Group 1 data

            distances_g1 = group1_distances[residue]

            # Sample if too many points

            if len(distances_g1) > 1000:

                sample_idx = np.random.choice(len(distances_g1), 1000, replace=False)

                distances_g1 = distances_g1[sample_idx]



            for dist in distances_g1:

                plot_data.append({

                    'Residue': residue,

                    'Group': group1_name,

                    'Distance': dist

                })



            # Group 2 data

            distances_g2 = group2_distances[residue]

            if len(distances_g2) > 1000:

                sample_idx = np.random.choice(len(distances_g2), 1000, replace=False)

                distances_g2 = distances_g2[sample_idx]



            for dist in distances_g2:

                plot_data.append({

                    'Residue': residue,

                    'Group': group2_name,

                    'Distance': dist

                })



        df = pd.DataFrame(plot_data)



        # Figure size based on number of residues

        if len(common_residues) > 8:

            figsize = get_standard_figsize('wide')

        else:

            figsize = get_standard_figsize('distribution')



        fig, ax = plt.subplots(figsize=figsize)



        # Color scheme for two groups

        colors = [PUBLICATION_COLORS['example'][0], PUBLICATION_COLORS['example'][1]]



        # Create split violin plot

        sns.violinplot(

            data=df,

            x='Residue',

            y='Distance',

            hue='Group',

            split=True,  # Split violins

            palette=colors,

            inner=None,  # No inner markers

            ax=ax

        )



        # Format residue names for display

        residues_display = format_residue_list(common_residues, residue_format)



        # Formatting

        ax.set_ylabel('Contact Distance (Å)')

        ax.set_xlabel('Amino Acid Residues')

        ax.set_xticks(range(len(common_residues)))

        ax.set_xticklabels(residues_display, rotation=45, ha='center')

        ax.grid(True, alpha=0.3, axis='y')

        ax.set_facecolor('#FAFAFA')



        # Legend

        ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)



        plt.tight_layout()

        plt.savefig(output_path, dpi=300, bbox_inches='tight',

                   facecolor='white', edgecolor='none')

        plt.close()



        print(f"  ✅ Comparison plot saved: {output_path}")

        return True



    except Exception as e:

        print(f"Error in comparison violin plot: {e}")

        import traceback

        traceback.print_exc()

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





