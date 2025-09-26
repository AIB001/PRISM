#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contact analysis plotting utilities for PRISM analysis module.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional
from pathlib import Path

# Import publication style
from .publication_utils import apply_publication_style

def plot_contact_analysis(contact_results: Dict,
                         output_path: str,
                         title: str = "Contact Analysis",
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
        # Apply publication style
        apply_publication_style()

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

        # Plot 1: Contact time series
        ax1 = axes[0, 0]
        for i, (traj_name, contacts) in enumerate(contact_results.items()):
            if len(contacts) > 0:
                time = np.arange(len(contacts)) * 0.5  # Assume 0.5ns interval
                ax1.plot(time, contacts, color=colors[i % 3],
                        label=traj_name, linewidth=2, alpha=0.8)

        ax1.set_xlabel('Time (ns)', fontfamily='Times New Roman')
        ax1.set_ylabel('Contact Count', fontfamily='Times New Roman')
        ax1.set_title('RemTP-Protein/RNA Contacts vs Time', fontfamily='Times New Roman', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot 2: Contact distribution
        ax2 = axes[0, 1]
        all_contacts = []
        for contacts in contact_results.values():
            if len(contacts) > 0:
                all_contacts.extend(contacts)

        if all_contacts:
            ax2.hist(all_contacts, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
            ax2.set_xlabel('Contact Count', fontfamily='Times New Roman')
            ax2.set_ylabel('Frequency', fontfamily='Times New Roman')
            ax2.set_title('Contact Distribution', fontfamily='Times New Roman', fontweight='bold')
            ax2.grid(True, alpha=0.3)

        # Plot 3: Average contacts per trajectory
        ax3 = axes[1, 0]
        traj_names = []
        avg_contacts = []
        std_contacts = []

        for traj_name, contacts in contact_results.items():
            if len(contacts) > 0:
                traj_names.append(traj_name)
                avg_contacts.append(np.mean(contacts))
                std_contacts.append(np.std(contacts))

        if traj_names:
            bars = ax3.bar(traj_names, avg_contacts, yerr=std_contacts,
                          color=colors[:len(traj_names)], alpha=0.7, capsize=5)
            ax3.set_ylabel('Average Contacts', fontfamily='Times New Roman')
            ax3.set_title('Average Contacts per Trajectory', fontfamily='Times New Roman', fontweight='bold')
            ax3.grid(True, alpha=0.3, axis='y')

        # Plot 4: Contact stability (running average)
        ax4 = axes[1, 1]
        window_size = min(50, len(list(contact_results.values())[0]) // 10) if contact_results else 10

        for i, (traj_name, contacts) in enumerate(contact_results.items()):
            if len(contacts) > window_size:
                # Calculate running average
                running_avg = np.convolve(contacts, np.ones(window_size)/window_size, mode='valid')
                time = np.arange(len(running_avg)) * 0.5
                ax4.plot(time, running_avg, color=colors[i % 3],
                        label=f'{traj_name} (avg)', linewidth=2)

        ax4.set_xlabel('Time (ns)', fontfamily='Times New Roman')
        ax4.set_ylabel('Running Average Contacts', fontfamily='Times New Roman')
        ax4.set_title('Contact Stability', fontfamily='Times New Roman', fontweight='bold')
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return True

    except Exception as e:
        print(f"Error in contact plotting: {e}")
        return False


def plot_key_residue_contacts(key_residue_results: Dict,
                            output_path: str,
                            title: str = "Key Residue Contacts") -> bool:
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

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
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
            ax.set_title(f'RemTP-{residue}', fontfamily='Times New Roman', fontweight='bold')
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