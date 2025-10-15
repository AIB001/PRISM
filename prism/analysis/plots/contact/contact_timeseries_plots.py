#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact time series and dynamic plots.

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



