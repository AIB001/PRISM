#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utility functions for PMF calculations
"""

import os
import numpy as np


def create_mdp_file(template_path, output_path, replacements):
    """Create MDP file from template with replacements"""
    with open(template_path, 'r') as f:
        content = f.read()
    
    for key, value in replacements.items():
        content = content.replace(f"{{{key}}}", str(value))
    
    with open(output_path, 'w') as f:
        f.write(content)


def extract_pull_groups(index_file):
    """Extract pull group indices from index file"""
    groups = {}
    with open(index_file, 'r') as f:
        current_group = None
        for line in f:
            line = line.strip()
            if line.startswith('[') and line.endswith(']'):
                current_group = line[1:-1].strip()
                groups[current_group] = []
            elif current_group and line:
                # Parse atom indices
                indices = [int(x) for x in line.split()]
                groups[current_group].extend(indices)
    return groups


def read_xvg(filename):
    """Read GROMACS XVG file"""
    x_data = []
    y_data = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                x_data.append(float(parts[0]))
                y_data.append(float(parts[1]))
    
    return np.array(x_data), np.array(y_data)