#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Publication-quality plotting utilities for scientific figures.

This module provides styling and utility functions for creating publication-ready
scientific plots that will be used as subfigures or panels in larger publications.
Fonts are sized appropriately for figures that will be scaled down.
"""

import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, Optional, Tuple, Any
import logging

logger = logging.getLogger(__name__)

# Publication-quality font sizes for subfigures/panels
# These are EXTREMELY large because figures will be heavily scaled down in publications
PUBLICATION_FONTS = {
    'title': 36,           # Main title - extremely large for visibility when scaled
    'subtitle': 32,        # Subplot titles
    'axis_label': 32,      # X/Y axis labels - must be readable when tiny
    'tick_label': 26,      # Tick labels - critical for data reading when small
    'legend': 14,          # Legend text - reduced to avoid clutter
    'annotation': 22,      # Statistical annotations
    'colorbar': 24,        # Colorbar labels
    'value_text': 18,      # Heatmap cell values - reduced for better fit in cells
    'bar_annotation': 14    # Bar chart annotations - very small to prevent overlap
}

# 标准panel尺寸定义 - 确保所有单个panel图片字体相对大小一致
STANDARD_PANEL_SIZES = {
    'single': (8, 6),        # 单panel标准尺寸 - 所有单个图片统一使用
    'horizontal': (12, 6),   # 水平排列两panel (1x2)
    'vertical': (8, 12),     # 垂直排列两panel (2x1)
    'quad': (12, 10),        # 2x2四panel
    'wide': (16, 6),         # 宽format单panel (特殊情况)
    'tall': (8, 10),         # 高format单panel (特殊情况)
    'distribution': (12, 8)  # 分布图专用 - 增加高度改善比例
}

# Color palettes for scientific publications
PUBLICATION_COLORS = {
    'default': ['#4A90E2', '#F5A623', '#7ED321', '#D0021B', '#9013FE'],  # Professional
    'nature': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'],   # Nature-style
    'science': ['#0173B2', '#DE8F05', '#029E73', '#CC78BC', '#CA9161'],  # Science-style
    'accessible': ['#0173B2', '#D55E00', '#009E73', '#CC79A7', '#F0E442'], # Colorblind-friendly
    'example': ['#A8DADC', '#D4A574', '#E9C46A', '#457B9D', '#F1FAEE']    # Matching example figure.py
}

def get_publication_style(style_type: str = 'default') -> Dict[str, Any]:
    """
    Get publication-quality matplotlib style parameters optimized for subfigures.

    Parameters
    ----------
    style_type : str
        Style type: 'default', 'nature', 'science', 'accessible'

    Returns
    -------
    dict
        Style parameters for matplotlib rcParams
    """
    base_style = {
        # Font settings - extra large for subfigures, Times New Roman for publication
        'font.family': 'Times New Roman',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
        'font.size': PUBLICATION_FONTS['tick_label'],
        'axes.titlesize': PUBLICATION_FONTS['title'],
        'axes.labelsize': PUBLICATION_FONTS['axis_label'],
        'xtick.labelsize': PUBLICATION_FONTS['tick_label'],
        'ytick.labelsize': PUBLICATION_FONTS['tick_label'],
        'legend.fontsize': PUBLICATION_FONTS['legend'],

        # Additional text elements for better control
        'figure.titlesize': PUBLICATION_FONTS['title'],        # Figure suptitle
        'axes.labelpad': 8,                                   # Space between axis labels and ticks

        # Bold axis labels by default for scientific publications
        'axes.labelweight': 'bold',
        'axes.titleweight': 'bold',

        # High-quality output - reduced DPI to prevent size issues with large fonts
        'figure.dpi': 100,      # Reduced from 300 to prevent extreme pixel sizes
        'savefig.dpi': 200,     # Still high quality for saving but won't cause layout issues
        # Removed 'savefig.bbox': 'tight' as it causes extreme size calculations with large fonts
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'none',

        # Enhanced grid and axes for scientific data - NO GRID GLOBALLY
        'grid.alpha': 0.0,         # No grid opacity by default
        'grid.linewidth': 0.0,     # No grid lines by default
        'axes.linewidth': 2.0,     # Thicker axes for subfigures
        'axes.edgecolor': 'black',
        'axes.grid': False,        # GLOBAL: No grid by default
        'axes.axisbelow': True,

        # Tick styling for professional appearance - larger for subfigures
        'xtick.major.width': 2.0,  # Increased from 1.5
        'ytick.major.width': 2.0,  # Increased from 1.5
        'xtick.major.size': 8,     # Increased from 6
        'ytick.major.size': 8,     # Increased from 6
        'xtick.minor.width': 1.5,  # Increased from 1.0
        'ytick.minor.width': 1.5,  # Increased from 1.0
        'xtick.minor.size': 4,     # Increased from 3
        'ytick.minor.size': 4,     # Increased from 3,

        # Legend styling
        'legend.frameon': True,
        'legend.fancybox': True,
        'legend.shadow': True,
        'legend.framealpha': 0.95,
        'legend.edgecolor': 'black',
        'legend.facecolor': 'white',

        # Error bar styling - enhanced for publication subfigures
        'errorbar.capsize': 12,     # Increased from 10 for better visibility
        'lines.linewidth': 4.0,     # Increased from 3.0 for subfigure scaling
        'lines.markersize': 12,     # Increased from 10 for better visibility

        # Figure background
        'figure.facecolor': 'white',
        'axes.facecolor': 'white'
    }

    return base_style

def get_standard_figsize(panel_type: str = 'single') -> Tuple[float, float]:
    """
    获取标准panel尺寸，确保所有图片字体相对大小一致。

    Parameters
    ----------
    panel_type : str
        Panel类型: 'single', 'horizontal', 'vertical', 'quad', 'wide', 'tall'

    Returns
    -------
    tuple
        标准figure尺寸 (width, height)
    """
    if panel_type in STANDARD_PANEL_SIZES:
        return STANDARD_PANEL_SIZES[panel_type]
    else:
        logger.warning(f"Unknown panel type '{panel_type}', using 'single' as default")
        return STANDARD_PANEL_SIZES['single']

def validate_figure_size(figsize: Tuple[float, float], dpi: int = 150) -> Tuple[float, float]:
    """
    Validate and adjust figure size to prevent matplotlib dimension limits.

    Matplotlib has a limit of 2^16 (65536) pixels in each dimension.
    This function ensures figures don't exceed safe limits while maintaining aspect ratio.

    Parameters
    ----------
    figsize : tuple
        Figure size in inches (width, height)
    dpi : int
        Dots per inch for size calculation

    Returns
    -------
    tuple
        Validated figure size that won't exceed matplotlib limits
    """
    max_pixels = 32000  # Safe limit well below matplotlib's 65536 limit
    width_px = figsize[0] * dpi
    height_px = figsize[1] * dpi

    if width_px > max_pixels or height_px > max_pixels:
        # Calculate scale factor to fit within limits
        scale = min(max_pixels / width_px, max_pixels / height_px)
        new_figsize = (figsize[0] * scale, figsize[1] * scale)
        logger.warning(f"Figure size {figsize} would exceed matplotlib limits. "
                      f"Scaled to {new_figsize} to prevent errors.")
        return new_figsize

    return figsize

def get_color_palette(palette: str = 'default', n_colors: int = 5) -> list:
    """
    Get publication-appropriate color palette.

    Parameters
    ----------
    palette : str
        Palette name: 'default', 'nature', 'science', 'accessible'
    n_colors : int
        Number of colors needed

    Returns
    -------
    list
        List of color hex codes
    """
    if palette not in PUBLICATION_COLORS:
        palette = 'default'

    colors = PUBLICATION_COLORS[palette]

    # Extend palette if more colors needed
    if n_colors > len(colors):
        # Use matplotlib colormap for additional colors
        extra_colors = plt.cm.Set3(np.linspace(0, 1, n_colors - len(colors)))
        colors = colors + [plt.colors.rgb2hex(c) for c in extra_colors]

    return colors[:n_colors]

def setup_publication_figure(figsize: Optional[Tuple[float, float]] = None,
                           panel_type: str = 'single', style_type: str = 'default',
                           title: str = '') -> Tuple[plt.Figure, plt.Axes]:
    """
    Create a publication-ready figure with proper styling and standard sizing.

    Parameters
    ----------
    figsize : tuple, optional
        Figure size in inches. If None, uses standard size for panel_type
    panel_type : str
        Panel type for standard sizing: 'single', 'horizontal', 'vertical', 'quad'
    style_type : str
        Style type for publication
    title : str
        Figure title (empty by default for publication)

    Returns
    -------
    tuple
        (figure, axes) objects
    """
    if figsize is None:
        figsize = get_standard_figsize(panel_type)

    plt.style.use('default')  # Reset first
    with plt.rc_context(get_publication_style(style_type)):
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_title(title)  # Empty by default

    return fig, ax

def apply_publication_style(style_type: str = 'default'):
    """
    Apply publication style globally to all matplotlib figures.

    This function sets the global matplotlib rcParams to use publication-quality
    styling for all subsequent plots. Call this once at the beginning of your
    plotting session.

    Parameters
    ----------
    style_type : str
        Style type: 'default', 'nature', 'science', 'accessible'
    """
    plt.rcParams.update(get_publication_style(style_type))
    logger.info(f"Applied {style_type} publication style globally")

def fix_rotated_labels(ax, labels, positions=None, rotation=45,
                      ha='center', va='top', fontsize=None, manual_alignment=True):
    """
    Fix rotated labels alignment for publication quality with perfect text centering.

    This function ensures that rotated text labels are properly centered on their
    corresponding tick marks, addressing the common issue where rotated labels
    appear misaligned.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object
    labels : list
        List of label strings
    positions : array-like, optional
        Tick positions. If None, uses current tick positions
    rotation : float
        Rotation angle in degrees
    ha : str
        Horizontal alignment (for manual_alignment=True, 'center' works best)
    va : str
        Vertical alignment
    fontsize : int, optional
        Font size. If None, uses publication default
    manual_alignment : bool
        Whether to use manual alignment calculation for perfect centering
    """
    import math

    if positions is None:
        positions = ax.get_xticks()

    if fontsize is None:
        fontsize = PUBLICATION_FONTS['tick_label']

    # Set tick positions
    ax.set_xticks(positions)

    if manual_alignment and rotation != 0:
        # Use standard matplotlib with improved alignment
        # This avoids extreme positioning that causes layout issues
        ax.set_xticklabels(labels, rotation=rotation, ha=ha, va=va,
                          fontsize=fontsize, weight='normal')

        # Ensure proper spacing with tick parameters
        ax.tick_params(axis='x', which='major', pad=8, labelsize=fontsize)
    else:
        # Standard matplotlib alignment (fallback)
        ax.set_xticklabels(labels, rotation=rotation, ha=ha, va=va,
                          fontsize=fontsize, weight='normal')

def add_statistical_annotations(ax, x_positions, statistics,
                              y_offset=5, fontsize=None, color='#2C3E50'):
    """
    Add statistical annotations (mean±std) above bars.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object
    x_positions : array-like
        X positions for annotations
    statistics : dict
        Dictionary with 'mean' and 'std' values for each position
    y_offset : float
        Vertical offset above data
    fontsize : int, optional
        Font size for annotations
    color : str
        Text color
    """
    if fontsize is None:
        fontsize = PUBLICATION_FONTS['annotation']

    for i, (x_pos, stats) in enumerate(zip(x_positions, statistics.values())):
        if 'mean' in stats and 'std' in stats:
            mean_val = stats['mean']
            std_val = stats['std']

            # Calculate y position
            if 'values' in stats:
                max_height = max(stats['values'])
                y_pos = max_height + std_val + y_offset
            else:
                y_pos = mean_val + std_val + y_offset

            # Add annotation with background box
            ax.text(x_pos, y_pos, f"{mean_val:.1f}±{std_val:.1f}",
                   ha='center', va='bottom', fontsize=fontsize,
                   color=color, weight='bold',
                   bbox=dict(boxstyle='round,pad=0.3',
                            facecolor='white', alpha=0.9,
                            edgecolor='gray', linewidth=0.8))

def style_axes_for_publication(ax, spine_width=1.5, tick_width=1.5,
                             tick_length=6, grid_alpha=0.6):
    """
    Apply publication-quality styling to axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object
    spine_width : float
        Width of axis spines
    tick_width : float
        Width of tick marks
    tick_length : float
        Length of tick marks
    grid_alpha : float
        Alpha transparency of grid lines
    """
    # Style spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(spine_width)
    ax.spines['bottom'].set_linewidth(spine_width)

    # Style ticks
    ax.tick_params(axis='both', which='major',
                  width=tick_width, length=tick_length)
    ax.tick_params(axis='both', which='minor',
                  width=tick_width*0.7, length=tick_length*0.5)

    # Style grid
    ax.grid(True, alpha=grid_alpha, linewidth=1.0, linestyle='-')
    ax.set_axisbelow(True)

def save_publication_figure(fig, save_path, dpi=300, format='png',
                          bbox_inches='tight', transparent=False):
    """
    Save figure with publication-quality settings.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save
    save_path : str
        Path to save figure
    dpi : int
        Dots per inch for raster formats (300+ recommended for publication)
    format : str
        File format ('png', 'pdf', 'svg', 'eps')
    bbox_inches : str
        Bounding box setting
    transparent : bool
        Whether to save with transparent background
    """
    # Ensure tight layout before saving
    fig.tight_layout()

    fig.savefig(save_path, dpi=dpi, format=format, bbox_inches=bbox_inches,
               facecolor='white' if not transparent else None,
               edgecolor='none', transparent=transparent)

    logger.info(f"Publication figure saved: {save_path} ({format.upper()}, {dpi} DPI)")

def create_publication_colormap(name='scientific', n_colors=256):
    """
    Create a publication-appropriate colormap.

    Parameters
    ----------
    name : str
        Colormap style: 'scientific', 'heatmap', 'diverging'
    n_colors : int
        Number of colors in map

    Returns
    -------
    matplotlib.colors.ListedColormap
        Custom colormap
    """
    from matplotlib.colors import ListedColormap

    if name == 'scientific':
        # Blue to red through white (good for heatmaps)
        colors = ['#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020']
    elif name == 'heatmap':
        # Yellow to red (classic heatmap)
        colors = ['#ffffcc', '#fed976', '#feb24c', '#fd8d3c', '#f03b20', '#bd0026']
    elif name == 'diverging':
        # Blue to red diverging
        colors = ['#2166ac', '#67a9cf', '#d1e5f0', '#fddbc7', '#ef8a62', '#b2182b']
    else:
        colors = ['#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020']

    return ListedColormap(colors, name=name, N=n_colors)