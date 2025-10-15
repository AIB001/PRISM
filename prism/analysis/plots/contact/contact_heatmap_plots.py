#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact heatmap and annotation plots.

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

        from ..publication_utils import PUBLICATION_FONTS



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







