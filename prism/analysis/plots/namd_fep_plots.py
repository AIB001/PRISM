#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NAMD FEP plotting functions for publication-quality visualizations.

This module provides plotting functions for NAMD Free Energy Perturbation analysis,
including convergence plots, energy decomposition charts, and dG/dλ profiles.
All plots follow PRISM publication standards with Times New Roman fonts and proper sizing.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, MaxNLocator
from typing import Dict, List, Optional, Tuple, Union
import logging

from .publication_utils import (
    get_publication_style,
    apply_publication_style,
    get_standard_figsize,
    save_publication_figure,
    get_color_palette,
    PUBLICATION_FONTS
)

logger = logging.getLogger(__name__)


def plot_namd_fep_convergence(convergence_data: pd.DataFrame,
                               energy_type: str = 'sum_dg',
                               xlabel: str = 'Simulation Time Fraction',
                               ylabel: str = 'ΔG (kcal/mol)',
                               title: str = '',
                               figsize: Optional[Tuple[float, float]] = None,
                               save_path: Optional[str] = None,
                               show_components: bool = False,
                               final_value_line: bool = True,
                               **kwargs) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot NAMD FEP convergence over simulation time.

    This function visualizes how the free energy estimate converges as more
    simulation data is included. Useful for assessing simulation convergence.

    Parameters
    ----------
    convergence_data : pd.DataFrame
        Convergence data from calculate_convergence_data(), with columns:
        'fraction', 'sum_dg', 'sum_delec', 'sum_dvdw', 'sum_couple'
    energy_type : str
        Energy component to plot: 'sum_dg' (default), 'sum_delec', 'sum_dvdw', 'sum_couple'
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    title : str
        Plot title (empty by default per PRISM convention)
    figsize : tuple, optional
        Figure size. If None, uses standard 'single' panel size
    save_path : str, optional
        Path to save the figure
    show_components : bool
        If True, plot all energy components on the same plot
    final_value_line : bool
        If True, show horizontal line at final converged value

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    # Set up figure with publication style
    if figsize is None:
        figsize = get_standard_figsize('single')

    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        if show_components:
            # Plot all energy components
            colors = get_color_palette('default', n_colors=4)
            components = [
                ('sum_dg', 'Total ΔG', colors[0]),
                ('sum_delec', 'ΔElec', colors[1]),
                ('sum_dvdw', 'ΔvdW', colors[2]),
                ('sum_couple', 'ΔCouple', colors[3])
            ]

            for col_name, label, color in components:
                if col_name in convergence_data.columns:
                    ax.plot(convergence_data['fraction'],
                           convergence_data[col_name],
                           marker='o',
                           linewidth=3,
                           markersize=10,
                           label=label,
                           color=color)

                    # Print final values
                    final_val = convergence_data[col_name].iloc[-1]
                    print(f"  Final {label}: {final_val:.2f} kcal/mol")

            ax.legend(frameon=True, loc='best', fontsize=PUBLICATION_FONTS['legend'])
            ylabel = 'Energy (kcal/mol)'

        else:
            # Plot single energy component
            if energy_type not in convergence_data.columns:
                raise ValueError(
                    f"Energy type '{energy_type}' not found in data. "
                    f"Available: {list(convergence_data.columns)}"
                )

            colors = get_color_palette('default', n_colors=1)
            ax.plot(convergence_data['fraction'],
                   convergence_data[energy_type],
                   marker='o',
                   linewidth=4,
                   markersize=12,
                   color=colors[0])

            # Add final value line
            if final_value_line:
                final_value = convergence_data[energy_type].iloc[-1]
                ax.axhline(y=final_value,
                          color='red',
                          linestyle='--',
                          linewidth=2,
                          alpha=0.7,
                          label=f'Final: {final_value:.2f}')
                ax.legend(frameon=True, loc='best', fontsize=PUBLICATION_FONTS['legend'])

                print(f"  Final {energy_type}: {final_value:.2f} kcal/mol")

        # Style axes
        ax.set_xlabel(xlabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel(ylabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        if title:
            ax.set_title(title, fontsize=PUBLICATION_FONTS['title'], weight='bold')

        ax.set_xlim(0, 1)
        ax.xaxis.set_major_locator(MultipleLocator(0.2))

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Set spine width
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)

        fig.tight_layout()

        # Save if requested
        if save_path:
            save_publication_figure(fig, save_path, **kwargs)
            logger.info(f"Convergence plot saved: {save_path}")

        print(f"✓ NAMD FEP convergence plot generated")
        return fig, ax


def plot_namd_fep_decomposition_bar(decomp_data: Union[np.ndarray, Dict[str, np.ndarray]],
                                    labels: Optional[List[str]] = None,
                                    error_data: Optional[Union[np.ndarray, Dict[str, np.ndarray]]] = None,
                                    xlabel: str = 'Energy Component',
                                    ylabel: str = 'ΔG (kcal/mol)',
                                    title: str = '',
                                    figsize: Optional[Tuple[float, float]] = None,
                                    save_path: Optional[str] = None,
                                    show_values: bool = True,
                                    zero_line: bool = True,
                                    **kwargs) -> Tuple[plt.Figure, plt.Axes]:
    """
    Create bar chart for NAMD FEP energy decomposition.

    This function visualizes the decomposition of free energy changes into
    electrostatic, van der Waals, and coupling contributions.

    Parameters
    ----------
    decomp_data : np.ndarray or dict
        Decomposition data. If array, should be [dG, dElec, dVdW, dCouple].
        If dict, keys are system names and values are decomposition arrays.
    labels : list of str, optional
        Labels for multiple datasets
    error_data : np.ndarray or dict, optional
        Error bars (standard deviations). Same structure as decomp_data.
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    title : str
        Plot title (empty by default)
    figsize : tuple, optional
        Figure size. If None, uses standard size
    save_path : str, optional
        Path to save figure
    show_values : bool
        Whether to show numerical values on bars
    zero_line : bool
        Whether to show horizontal line at y=0

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    # Set up figure
    if figsize is None:
        figsize = get_standard_figsize('single')

    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        # Component labels
        component_labels = ['ΔG', 'ΔElec', 'ΔvdW', 'ΔCouple']

        # Handle different input types
        if isinstance(decomp_data, dict):
            datasets = list(decomp_data.values())
            if labels is None:
                labels = list(decomp_data.keys())
        else:
            datasets = [decomp_data]
            if labels is None:
                labels = ['Data']

        # Handle error data
        error_datasets = None
        if error_data is not None:
            if isinstance(error_data, dict):
                error_datasets = list(error_data.values())
            else:
                error_datasets = [error_data]

        n_components = 4
        n_datasets = len(datasets)

        # Get colors
        colors = get_color_palette('default', n_colors=n_datasets)

        # Bar positioning
        width = 0.35 if n_datasets > 1 else 0.6
        x = np.arange(n_components)

        # Plot bars
        for i, (data, label, color) in enumerate(zip(datasets, labels, colors)):
            if len(data) < n_components:
                logger.warning(f"Data for {label} has {len(data)} components, expected {n_components}")
                continue

            # Calculate bar positions
            if n_datasets > 1:
                x_pos = x + width * (i - (n_datasets - 1) / 2)
            else:
                x_pos = x

            # Get error bars if available
            yerr = error_datasets[i] if error_datasets is not None else None

            # Plot bars
            bars = ax.bar(x_pos, data[:n_components], width,
                         label=label, color=color, alpha=0.8,
                         edgecolor='black', linewidth=1.5,
                         yerr=yerr, capsize=8,
                         error_kw={'linewidth': 2, 'capthick': 2})

            # Add value labels
            if show_values:
                for bar, value in zip(bars, data[:n_components]):
                    height = bar.get_height()
                    va = 'bottom' if value >= 0 else 'top'
                    y_offset = 0.5 if value >= 0 else -0.5
                    ax.text(bar.get_x() + bar.get_width() / 2., height + y_offset,
                           f'{value:.1f}',
                           ha='center', va=va,
                           fontsize=PUBLICATION_FONTS['bar_annotation'],
                           weight='bold')

            # Print values
            print(f"  {label}:")
            for comp, val in zip(component_labels, data[:n_components]):
                if yerr is not None and i < len(error_datasets):
                    print(f"    {comp}: {val:7.2f} ± {yerr[component_labels.index(comp)]:.2f} kcal/mol")
                else:
                    print(f"    {comp}: {val:7.2f} kcal/mol")

        # Add zero line
        if zero_line:
            ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.5)

        # Style axes
        ax.set_xticks(x)
        ax.set_xticklabels(component_labels, fontsize=PUBLICATION_FONTS['tick_label'])
        ax.set_xlabel(xlabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel(ylabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        if title:
            ax.set_title(title, fontsize=PUBLICATION_FONTS['title'], weight='bold')

        # Add legend if multiple datasets
        if n_datasets > 1:
            ax.legend(frameon=True, loc='best', fontsize=PUBLICATION_FONTS['legend'])

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)

        fig.tight_layout()

        # Save if requested
        if save_path:
            save_publication_figure(fig, save_path, **kwargs)
            logger.info(f"Decomposition bar plot saved: {save_path}")

        print(f"✓ NAMD FEP decomposition bar plot generated")
        return fig, ax


def plot_namd_fep_dg_lambda(fepout_data: pd.DataFrame,
                            phase: str = 'complex',
                            xlabel: str = 'λ',
                            ylabel_deriv: str = 'dG/dλ (kcal/mol)',
                            ylabel_accum: str = 'ΔG (kcal/mol)',
                            title: str = '',
                            figsize: Optional[Tuple[float, float]] = None,
                            save_path: Optional[str] = None,
                            show_both: bool = True,
                            **kwargs) -> Union[Tuple[plt.Figure, plt.Axes],
                                              Tuple[plt.Figure, Tuple[plt.Axes, plt.Axes]]]:
    """
    Plot dG/dλ derivative and accumulated ΔG vs λ for NAMD FEP.

    This function creates plots showing the free energy derivative (dG/dλ) and
    the accumulated free energy change along the alchemical coordinate λ.

    Parameters
    ----------
    fepout_data : pd.DataFrame
        FEP output data with columns: 'start', 'stop', 'dG', 'dG_accum'
    phase : str
        Phase label (e.g., 'complex', 'ligand')
    xlabel : str
        X-axis label
    ylabel_deriv : str
        Y-axis label for derivative plot
    ylabel_accum : str
        Y-axis label for accumulated plot
    title : str
        Plot title (empty by default)
    figsize : tuple, optional
        Figure size
    save_path : str, optional
        Path to save figure
    show_both : bool
        If True, create 2-panel figure with both derivative and accumulated plots

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes or tuple of axes
        Single axis if show_both=False, otherwise tuple of (ax_deriv, ax_accum)
    """
    # Check required columns
    required_cols = ['start', 'stop', 'dG', 'dG_accum']
    if not all(col in fepout_data.columns for col in required_cols):
        raise ValueError(f"fepout_data must contain columns: {required_cols}")

    # Calculate lambda midpoints and derivatives
    lambda_mid = (fepout_data['start'] + fepout_data['stop']) / 2
    d_lambda = fepout_data['stop'] - fepout_data['start']
    dG_dlambda = fepout_data['dG'] / d_lambda

    # Set up figure
    if figsize is None:
        if show_both:
            figsize = get_standard_figsize('horizontal')  # Wide format for 2 panels
        else:
            figsize = get_standard_figsize('single')

    with plt.rc_context(get_publication_style()):
        if show_both:
            fig, (ax_deriv, ax_accum) = plt.subplots(1, 2, figsize=figsize)
        else:
            fig, ax_deriv = plt.subplots(figsize=figsize)
            ax_accum = None

        colors = get_color_palette('default', n_colors=1)

        # Plot derivative
        ax_deriv.plot(lambda_mid, dG_dlambda,
                     marker='o', linewidth=3, markersize=10,
                     color=colors[0])
        ax_deriv.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

        ax_deriv.set_xlabel(xlabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax_deriv.set_ylabel(ylabel_deriv, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        if title and not show_both:
            ax_deriv.set_title(title, fontsize=PUBLICATION_FONTS['title'], weight='bold')

        ax_deriv.set_xlim(0, 1)
        ax_deriv.spines['top'].set_visible(False)
        ax_deriv.spines['right'].set_visible(False)
        ax_deriv.spines['left'].set_linewidth(2)
        ax_deriv.spines['bottom'].set_linewidth(2)

        # Plot accumulated if requested
        if show_both:
            ax_accum.plot(lambda_mid, fepout_data['dG_accum'],
                         marker='o', linewidth=3, markersize=10,
                         color=colors[0])

            ax_accum.set_xlabel(xlabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
            ax_accum.set_ylabel(ylabel_accum, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

            ax_accum.set_xlim(0, 1)
            ax_accum.spines['top'].set_visible(False)
            ax_accum.spines['right'].set_visible(False)
            ax_accum.spines['left'].set_linewidth(2)
            ax_accum.spines['bottom'].set_linewidth(2)

            # Print final value
            final_dg = fepout_data['dG_accum'].iloc[-1]
            print(f"  {phase.capitalize()} phase final ΔG: {final_dg:.2f} kcal/mol")

        fig.tight_layout()

        # Save if requested
        if save_path:
            save_publication_figure(fig, save_path, **kwargs)
            logger.info(f"dG/dλ plot saved: {save_path}")

        print(f"✓ NAMD FEP dG/dλ plot generated for {phase} phase")

        if show_both:
            return fig, (ax_deriv, ax_accum)
        else:
            return fig, ax_deriv


def plot_namd_fep_multi_comparison(system_results: Dict[str, Dict[str, np.ndarray]],
                                   comparison_type: str = 'ddG',
                                   xlabel: str = 'System',
                                   ylabel: str = 'ΔΔG (kcal/mol)',
                                   title: str = '',
                                   figsize: Optional[Tuple[float, float]] = None,
                                   save_path: Optional[str] = None,
                                   show_components: bool = False,
                                   rotate_labels: bool = True,
                                   **kwargs) -> Tuple[plt.Figure, plt.Axes]:
    """
    Compare NAMD FEP results across multiple systems/compounds.

    This function creates a bar chart comparing binding free energies (ΔΔG)
    or energy components across different molecular systems.

    Parameters
    ----------
    system_results : dict
        Dictionary mapping system names to their FEP results
        Each result should be from calculate_system_decomposition()
    comparison_type : str
        What to compare: 'ddG' (binding free energy), 'complex', 'ligand', or 'decomposition'
    xlabel : str
        X-axis label
    ylabel : str
        Y-axis label
    title : str
        Plot title (empty by default)
    figsize : tuple, optional
        Figure size
    save_path : str, optional
        Path to save figure
    show_components : bool
        If True and comparison_type='decomposition', show all energy components
    rotate_labels : bool
        Whether to rotate x-axis labels for better readability

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : matplotlib.axes.Axes
    """
    if not system_results:
        raise ValueError("No system results provided")

    # Set up figure
    if figsize is None:
        if show_components:
            figsize = get_standard_figsize('wide')
        else:
            figsize = get_standard_figsize('single')

    with plt.rc_context(get_publication_style()):
        fig, ax = plt.subplots(figsize=figsize)

        system_names = list(system_results.keys())
        n_systems = len(system_names)

        if comparison_type == 'ddG':
            # Plot binding free energies
            ddG_values = []
            ddG_errors = []

            for name in system_names:
                result = system_results[name]
                if 'ddG' in result:
                    ddG_values.append(result['ddG'][0])  # First element is total dG
                    if 'ddG_std' in result:
                        ddG_errors.append(result['ddG_std'][0])
                    else:
                        ddG_errors.append(0)
                else:
                    logger.warning(f"No ddG found for {name}")
                    ddG_values.append(0)
                    ddG_errors.append(0)

            colors = get_color_palette('default', n_colors=n_systems)
            x = np.arange(n_systems)

            bars = ax.bar(x, ddG_values, color=colors, alpha=0.8,
                         edgecolor='black', linewidth=1.5,
                         yerr=ddG_errors, capsize=8,
                         error_kw={'linewidth': 2, 'capthick': 2})

            # Add value labels
            for bar, value, error in zip(bars, ddG_values, ddG_errors):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height,
                       f'{value:.1f}',
                       ha='center', va='bottom',
                       fontsize=PUBLICATION_FONTS['bar_annotation'],
                       weight='bold')

            # Print values
            print(f"\nBinding Free Energies (ΔΔG):")
            for name, val, err in zip(system_names, ddG_values, ddG_errors):
                print(f"  {name}: {val:7.2f} ± {err:.2f} kcal/mol")

        elif comparison_type == 'decomposition' and show_components:
            # Plot energy decomposition for all systems
            component_labels = ['ΔG', 'ΔElec', 'ΔvdW', 'ΔCouple']
            n_components = 4
            width = 0.8 / n_components

            for i, name in enumerate(system_names):
                result = system_results[name]
                if 'ddG' in result:
                    data = result['ddG'][:n_components]
                    colors = get_color_palette('default', n_colors=n_components)

                    for j, (val, comp, color) in enumerate(zip(data, component_labels, colors)):
                        x_pos = i + (j - n_components/2 + 0.5) * width
                        ax.bar(x_pos, val, width, color=color,
                              alpha=0.8, edgecolor='black', linewidth=1)

            # Add component legend
            legend_handles = [plt.Rectangle((0,0),1,1, fc=c, ec='black', alpha=0.8)
                            for c in get_color_palette('default', n_colors=n_components)]
            ax.legend(legend_handles, component_labels,
                     frameon=True, loc='best', fontsize=PUBLICATION_FONTS['legend'])

        ax.set_xticks(np.arange(n_systems))
        ax.set_xticklabels(system_names, fontsize=PUBLICATION_FONTS['tick_label'])

        if rotate_labels and n_systems > 5:
            ax.set_xticklabels(system_names, rotation=45, ha='right',
                             fontsize=PUBLICATION_FONTS['tick_label'])

        ax.set_xlabel(xlabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        ax.set_ylabel(ylabel, fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
        if title:
            ax.set_title(title, fontsize=PUBLICATION_FONTS['title'], weight='bold')

        # Add zero line
        ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.5)

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(2)

        fig.tight_layout()

        # Save if requested
        if save_path:
            save_publication_figure(fig, save_path, **kwargs)
            logger.info(f"Multi-system comparison plot saved: {save_path}")

        print(f"✓ NAMD FEP multi-system comparison plot generated ({n_systems} systems)")
        return fig, ax


def plot_namd_fep_full_analysis(results: Dict,
                                 output_path: Optional[str] = None,
                                 endpoint_cutoff_multiplier: float = 3.0,
                                 check_linearity: bool = False,
                                 figsize: Optional[Tuple[float, float]] = None,
                                 **kwargs) -> Tuple[plt.Figure, Tuple]:
    """
    Create comprehensive 4-panel FEP analysis figure matching SciDraft-Studio reference.

    This function creates a publication-quality figure with:
    - dG/dλ derivative plots for complex and ligand phases
    - Accumulated ΔG plots for complex and ligand phases
    - Mean ± STD shaded error bands across repeats
    - End-point catastrophe detection (red stars on problematic windows)
    - Optional linearity checks (R² for VdW/Elec regions)
    - Summary saved to separate text file (not on figure)

    Parameters
    ----------
    results : dict
        Results from batch_process_fep_system() containing:
        - 'raw_data': Combined DataFrame with all repeat data
        - 'complex_raw_mean', 'complex_raw_std': Statistics for complex
        - 'ligand_raw_mean', 'ligand_raw_std': Statistics for ligand
        - 'complex_lambda', 'ligand_lambda': Lambda windows
        - 'statistics': Summary statistics
    output_path : str, optional
        Path to save the figure. Summary will be saved to <basename>_analysis_summary.txt
    endpoint_cutoff_multiplier : float
        Multiplier for dynamic end-point catastrophe threshold (default: 3.0)
    check_linearity : bool
        If True, perform linearity check on dG/dλ for VdW/Elec regions
    figsize : tuple, optional
        Figure size (width, height in inches). If None, uses (18, 12)

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : tuple
        Tuple of (ax_complex_deriv, ax_complex_accum, ax_ligand_deriv, ax_ligand_accum)
    """
    if figsize is None:
        figsize = (18, 12)

    with plt.rc_context(get_publication_style()):
        # Create 2x2 subplot layout
        fig, ((ax_c_deriv, ax_c_accum), (ax_l_deriv, ax_l_accum)) = plt.subplots(
            2, 2, figsize=figsize
        )

        colors = get_color_palette('default', n_colors=2)
        color_complex = colors[0]
        color_ligand = colors[1]

        summary_lines = []
        quality_warnings = []

        # Process each phase
        for phase, color, ax_deriv, ax_accum in [
            ('complex', color_complex, ax_c_deriv, ax_c_accum),
            ('ligand', color_ligand, ax_l_deriv, ax_l_accum)
        ]:
            if f'{phase}_raw_mean' not in results or f'{phase}_lambda' not in results:
                logger.warning(f"No data for {phase} phase")
                ax_deriv.text(0.5, 0.5, f'No {phase} data', ha='center', va='center')
                ax_accum.text(0.5, 0.5, f'No {phase} data', ha='center', va='center')
                continue

            # Get data (now using STD instead of SEM)
            dG_accum_mean = results[f'{phase}_raw_mean']
            dG_accum_std = results[f'{phase}_raw_std']  # Changed from _sem to _std
            lambda_windows = results[f'{phase}_lambda']

            lambda_start = lambda_windows[:, 0]
            lambda_stop = lambda_windows[:, 1]
            lambda_mid = (lambda_start + lambda_stop) / 2
            d_lambda = lambda_stop - lambda_start

            # Calculate dG/dλ derivative
            dG = np.diff(np.concatenate([[0], dG_accum_mean]))  # dG per window
            dG_dlambda_mean = dG / d_lambda

            # Calculate derivative STD (propagation of uncertainty)
            dG_std = np.diff(np.concatenate([[0], dG_accum_std]))
            dG_dlambda_std = dG_std / d_lambda

            n_windows = len(lambda_mid)
            final_dG = dG_accum_mean[-1]
            final_dG_std = dG_accum_std[-1]  # Changed from _sem to _std

            # Plot derivative (dG/dλ) with STD error bands
            ax_deriv.plot(lambda_mid, dG_dlambda_mean,
                         marker='o', linewidth=3, markersize=8,
                         color=color, label=f'{phase.capitalize()} (Step)', linestyle='-')
            ax_deriv.fill_between(lambda_mid,
                                 dG_dlambda_mean - dG_dlambda_std,
                                 dG_dlambda_mean + dG_dlambda_std,
                                 color=color, alpha=0.2)
            ax_deriv.axhline(y=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)

            # Plot accumulated ΔG with STD error bands
            ax_accum.plot(lambda_mid, dG_accum_mean,
                         marker='o', linewidth=3, markersize=8,
                         color=color, label=f'{phase.capitalize()} (Accum)', linestyle='-')
            ax_accum.fill_between(lambda_mid,
                                 dG_accum_mean - dG_accum_std,
                                 dG_accum_mean + dG_accum_std,
                                 color=color, alpha=0.2)

            # End-point catastrophe check
            middle_idx = n_windows // 2
            ref_dG_dlambda = dG_dlambda_mean[max(0, middle_idx - 2)]
            dynamic_cutoff = endpoint_cutoff_multiplier * abs(ref_dG_dlambda)

            epc_windows = []
            for idx in range(n_windows):
                # Check first 5 and last 5 windows
                if (idx < 5 or idx >= n_windows - 5):
                    if abs(dG_dlambda_mean[idx]) > dynamic_cutoff:
                        epc_windows.append(f"λ≈{lambda_mid[idx]:.2f}")
                        # Mark with red star
                        ax_deriv.scatter(lambda_mid[idx], dG_dlambda_mean[idx],
                                       marker='*', s=600, color='red',
                                       zorder=10, edgecolor='black', linewidths=2)

            if epc_windows:
                quality_warnings.append(
                    f"  - {phase.capitalize()}: Potential end-point catastrophe at {', '.join(epc_windows)}"
                )

            # Linearity check (optional)
            if check_linearity:
                for region, name in [(lambda_mid < 0.5, 'VdW'), (lambda_mid >= 0.5, 'Elec')]:
                    if np.sum(region) > 1:
                        x = lambda_mid[region]
                        y = dG_dlambda_mean[region]
                        if len(x) > 1:
                            slope, intercept = np.polyfit(x, y, 1)
                            y_pred = slope * x + intercept
                            ss_res = np.sum((y - y_pred) ** 2)
                            ss_tot = np.sum((y - np.mean(y)) ** 2)
                            r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
                            quality_warnings.append(
                                f"  - {phase.capitalize()} ({name}): slope={slope:.2f}, R²={r2:.2f}"
                            )

            # Style derivative plot
            ax_deriv.set_ylabel('\u0394G/\u0394\u03bb (kcal/mol)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
            ax_deriv.set_xlim(0, 1)
            ax_deriv.grid(True, alpha=0.3)
            # Keep all four spines visible
            ax_deriv.spines['left'].set_linewidth(2)
            ax_deriv.spines['bottom'].set_linewidth(2)
            ax_deriv.spines['top'].set_linewidth(2)
            ax_deriv.spines['right'].set_linewidth(2)

            # Only show x-label on bottom row
            if phase == 'ligand':
                ax_deriv.set_xlabel('\u03bb', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

            # Style accumulated plot
            ax_accum.set_ylabel('\u0394G (kcal/mol)', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')
            ax_accum.set_xlim(0, 1)
            ax_accum.grid(True, alpha=0.3)
            # Keep all four spines visible
            ax_accum.spines['left'].set_linewidth(2)
            ax_accum.spines['bottom'].set_linewidth(2)
            ax_accum.spines['top'].set_linewidth(2)
            ax_accum.spines['right'].set_linewidth(2)

            # Only show x-label on bottom row
            if phase == 'ligand':
                ax_accum.set_xlabel('\u03bb', fontsize=PUBLICATION_FONTS['axis_label'], weight='bold')

            # Add legends instead of text labels
            ax_deriv.legend(loc='upper left', frameon=True, fontsize=PUBLICATION_FONTS['legend'])
            ax_accum.legend(loc='upper left', frameon=True, fontsize=PUBLICATION_FONTS['legend'])

            # Add to summary (no printing)
            n_repeats = results['statistics'].get(f'n_{phase}_repeats', 0)
            summary_lines.append(
                f"{phase.capitalize()}: ΔG = {final_dG:.2f} ± {final_dG_std:.2f} kcal/mol (n={n_repeats})"
            )

        # Add binding free energy if both phases present
        if 'statistics' in results:
            stats = results['statistics']
            if 'binding_ddG' in stats:
                ddG = stats['binding_ddG']
                ddG_std = stats['binding_ddG_std']
                summary_lines.append(f"\nBinding ΔΔG = {ddG:.2f} ± {ddG_std:.2f} kcal/mol")

        # Create summary text
        summary_text = "NAMD FEP Analysis Results\n" + "="*30 + "\n"
        summary_text += "\n".join(summary_lines)

        if quality_warnings:
            summary_text += "\n\nQuality Warnings:\n" + "="*30 + "\n"
            summary_text += "\n".join(quality_warnings)
        else:
            summary_text += "\n\n✓ All quality checks passed"

        # No figure title or text box - save summary to separate file instead
        plt.tight_layout()

        # Save if requested
        if output_path:
            # Save the figure
            save_publication_figure(fig, output_path, **kwargs)
            logger.info(f"Full FEP analysis figure saved: {output_path}")

            # Save summary to separate text file
            from pathlib import Path
            output_base = Path(output_path).with_suffix('')
            summary_file = f"{output_base}_analysis_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(summary_text)
            logger.info(f"Analysis summary saved: {summary_file}")

        return fig, (ax_c_deriv, ax_c_accum, ax_l_deriv, ax_l_accum)
