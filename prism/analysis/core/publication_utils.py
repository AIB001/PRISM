#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-quality plotting utilities for PRISM analysis.
"""

import matplotlib.pyplot as plt
from typing import Dict, Tuple

# Publication-quality font settings
PUBLICATION_FONTS = {
    "family": "sans-serif",
    "sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "size": 10,
}


def apply_publication_style():
    """Apply publication-quality style to matplotlib."""
    plt.rcParams.update(
        {
            "font.family": PUBLICATION_FONTS["family"],
            "font.sans-serif": PUBLICATION_FONTS["sans-serif"],
            "font.size": PUBLICATION_FONTS["size"],
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.linewidth": 1.0,
            "grid.linewidth": 0.5,
            "lines.linewidth": 1.5,
        }
    )


def get_publication_style() -> Dict:
    """Get publication-quality style dictionary."""
    return {
        "font.family": PUBLICATION_FONTS["family"],
        "font.sans-serif": PUBLICATION_FONTS["sans-serif"],
        "font.size": PUBLICATION_FONTS["size"],
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
    }


def get_standard_figsize(width: str = "single", aspect: float = 0.75) -> Tuple[float, float]:
    """
    Get standard figure size for publications.

    Parameters
    ----------
    width : str
        'single' (3.5 inches) or 'double' (7.0 inches)
    aspect : float
        Height/width ratio

    Returns
    -------
    tuple
        (width, height) in inches
    """
    widths = {"single": 3.5, "double": 7.0}
    w = widths.get(width, 3.5)
    return (w, w * aspect)


def fix_rotated_labels(ax, rotation: float = 45, ha: str = "right"):
    """
    Fix rotated x-axis labels for better readability.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object
    rotation : float
        Rotation angle in degrees
    ha : str
        Horizontal alignment
    """
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation, ha=ha)


def save_publication_figure(
    fig, filepath: str, dpi: int = 300, bbox_inches: str = "tight", transparent: bool = False, **kwargs
):
    """
    Save figure with publication-quality settings.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object to save
    filepath : str
        Output file path
    dpi : int
        Resolution in dots per inch
    bbox_inches : str
        Bounding box setting
    transparent : bool
        Whether to use transparent background
    **kwargs
        Additional arguments passed to fig.savefig()
    """
    fig.savefig(filepath, dpi=dpi, bbox_inches=bbox_inches, transparent=transparent, **kwargs)
