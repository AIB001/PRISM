#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact heatmap and annotation plots.

Refactored from contact_plots.py for better modularity.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, Optional
from pathlib import Path

# Import publication style
from ...core.publication_utils import apply_publication_style, get_standard_figsize

# Import residue formatting utility
from .....utils.residue import format_residue_list


def plot_contact_heatmap_annotated(
    contact_data: Dict[str, Dict[str, float]],
    output_path: str = None,
    title: str = "",
    separate_panels: bool = False,
    residue_format: str = "1letter",
    ax: Optional[plt.Axes] = None,
) -> bool:
    """

    Plot annotated heatmap for contact analysis.



    Parameters

    ----------

    contact_data : dict

        Dictionary with residue names as keys and trajectory data as nested dict

    output_path : str, optional

        Path to save the figure (required in standalone mode, ignored in panel mode)

    title : str

        Main plot title

    separate_panels : bool

        If True, also save individual panel figure

    residue_format : str

        Amino acid display format: "1letter" (e.g., D618) or "3letter" (e.g., ASP618)

    ax : matplotlib.axes.Axes, optional

        Axes object to plot on. If provided, function operates in panel mode.



    Returns

    -------

    bool or matplotlib.axes.Axes

        In standalone mode (ax=None): Returns True if successful

        In panel mode (ax provided): Returns the Axes object

    """

    # Detect mode

    if ax is None:
        # Standalone mode - original behavior

        if output_path is None:
            raise ValueError("output_path is required in standalone mode (when ax is not provided)")

        own_figure = True

    else:
        # Panel mode - use provided axis

        own_figure = False

    try:
        if own_figure:
            print("🔥 Contact probability heatmap with annotations showing residue-trajectory interactions")

            apply_publication_style()

        from ...core.publication_utils import PUBLICATION_FONTS

        # Prepare data for heatmap

        residues = list(contact_data.keys())

        trajectories = ["Repeat 1", "Repeat 2", "Repeat 3"]

        # Create data matrix

        data_matrix = []

        for traj_num in [1, 2, 3]:
            row = [contact_data[res].get(f"Repeat{traj_num}", 0.0) for res in residues]

            data_matrix.append(row)

        data_matrix = np.array(data_matrix)

        # Create figure only in standalone mode

        if own_figure:
            # Use wider figsize for more residues (15+)

            if len(residues) > 12:
                figsize = get_standard_figsize("wide")  # (16, 6)

            elif len(residues) > 6:
                figsize = get_standard_figsize("horizontal")  # (12, 6)

            else:
                figsize = get_standard_figsize("single")  # (8, 6)

            fig, ax = plt.subplots(figsize=figsize)

        # Create heatmap

        im = ax.imshow(data_matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=100)

        # Add colorbar with proper font size

        cbar = plt.colorbar(im, label="Contact Probability (%)", shrink=0.8)

        cbar.set_label("Contact Probability (%)", fontweight="bold", fontsize=PUBLICATION_FONTS["colorbar"])

        cbar.ax.tick_params(labelsize=PUBLICATION_FONTS["tick_label"])

        # Format residue names for display

        residues_display = format_residue_list(residues, residue_format)

        # Set ticks and labels

        ax.set_yticks([0, 1, 2])

        ax.set_yticklabels(trajectories)

        ax.set_xticks(range(len(residues)))

        ax.set_xticklabels(residues_display, rotation=45, ha="center")

        if title:
            ax.set_title(title, fontsize=PUBLICATION_FONTS["title"], fontweight="bold")

        # Add text annotations with proper font size

        for i in range(3):
            for j in range(len(residues)):
                text_color = "white" if data_matrix[i, j] > 50 else "black"

                ax.text(
                    j,
                    i,
                    f"{data_matrix[i, j]:.1f}",
                    ha="center",
                    va="center",
                    color=text_color,
                    fontweight="bold",
                    fontsize=PUBLICATION_FONTS["value_text"],
                )

        # Handle output based on mode

        if own_figure:
            plt.tight_layout()

            plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")

            if separate_panels:
                base_path = Path(output_path)
                panel_path = base_path.parent / f"{base_path.stem}_panel{base_path.suffix}"
                plt.savefig(panel_path, dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")

            plt.close()

            return True

        else:
            # Panel mode - return the axis

            return ax

    except Exception as e:
        print(f"Error in contact heatmap: {e}")

        return False
