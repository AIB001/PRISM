#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MMPBSA Utilities - Helper functions for AMBER MMPBSA.py workflow.

Standalone usage (fix prmtop files after parmed GROMACS->AMBER conversion):
    python3 utils.py solvated.prmtop complex.prmtop receptor.prmtop ligand.prmtop
"""

import re
import sys


def fix_prmtop_periodicity(prmtop_path):
    """Fix invalid dihedral periodicity values in AMBER prmtop files.

    parmed's GROMACS-to-AMBER conversion can produce non-standard periodicity
    values (e.g., 84) that AMBER's mmpbsa_py_energy rejects. This function
    clamps all periodicity values to the valid AMBER range (1-6).

    Operates directly on the prmtop text format to avoid re-serialization
    issues with parmed round-tripping.

    Parameters
    ----------
    prmtop_path : str
        Path to the AMBER prmtop file to fix in-place.

    Returns
    -------
    int
        Number of periodicity values that were fixed.
    """
    with open(prmtop_path, 'r') as f:
        lines = f.readlines()

    # Find DIHEDRAL_PERIODICITY section
    in_section = False
    fmt_pattern = None
    fixed = 0
    new_lines = []

    i = 0
    while i < len(lines):
        line = lines[i]

        if line.startswith('%FLAG DIHEDRAL_PERIODICITY'):
            in_section = True
            new_lines.append(line)
            i += 1
            # Next line is %FORMAT
            if i < len(lines) and lines[i].startswith('%FORMAT'):
                new_lines.append(lines[i])
                # Parse format: e.g., %FORMAT(5E16.8)
                fmt_match = re.match(r'%FORMAT\((\d+)E(\d+)\.(\d+)\)', lines[i].strip())
                if fmt_match:
                    n_per_line = int(fmt_match.group(1))
                    field_width = int(fmt_match.group(2))
                    decimal = int(fmt_match.group(3))
                    fmt_pattern = (n_per_line, field_width, decimal)
                i += 1
            continue

        if in_section:
            if line.startswith('%FLAG'):
                # End of section
                in_section = False
                new_lines.append(line)
                i += 1
                continue

            if fmt_pattern and line.strip():
                n_per_line, field_width, decimal = fmt_pattern
                # Parse fixed-width fields
                new_fields = []
                for j in range(0, len(line.rstrip('\n')), field_width):
                    field = line[j:j + field_width]
                    if not field.strip():
                        new_fields.append(field)
                        continue
                    try:
                        val = float(field)
                        int_val = abs(int(round(val)))
                        sign = -1.0 if val < 0 else 1.0
                        if int_val > 6 or int_val < 1:
                            # Clamp to valid range 1-6
                            new_val = ((int_val - 1) % 6) + 1
                            new_fields.append(f'{sign * new_val:{field_width}.{decimal}E}')
                            fixed += 1
                        else:
                            new_fields.append(field)
                    except ValueError:
                        new_fields.append(field)
                new_lines.append(''.join(new_fields) + '\n')
                i += 1
                continue

        new_lines.append(line)
        i += 1

    if fixed:
        with open(prmtop_path, 'w') as f:
            f.writelines(new_lines)

    return fixed


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <prmtop_file> [prmtop_file ...]")
        sys.exit(1)

    total_fixed = 0
    for path in sys.argv[1:]:
        n = fix_prmtop_periodicity(path)
        if n:
            print(f"  Fixed {n} invalid periodicity value(s) in {path}")
        else:
            print(f"  {path}: all periodicity values OK")
        total_fixed += n

    if total_fixed:
        print(f"Total: fixed {total_fixed} value(s)")
    else:
        print("No fixes needed.")
