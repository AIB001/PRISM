#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM PMF Alignment Module

This module handles structure alignment for PMF calculations:
1. Parses protein PDB and ligand MOL2/SDF
2. Identifies binding pocket residues (heavy atoms within cutoff of any ligand atom)
3. Optimizes pull direction via Metropolis-Hastings to maximize clearance from pocket
4. Rotates protein+ligand together to align optimized pull vector with Z-axis
5. Generates PyMOL visualization script for inspection
"""

import os
import sys
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List, Dict
import logging

logger = logging.getLogger(__name__)


# ============================================================
# Standalone helper functions for MH-optimized pulling direction
# ============================================================


def _is_heavy(element: str) -> bool:
    """Check if an element is a heavy atom (non-hydrogen)."""
    return element.upper() not in ("H", "D")


def _get_element_from_atom(atom: Dict) -> str:
    """Extract element symbol from a parsed atom dict."""
    # PDB atoms have 'element' field; MOL2 atoms have 'type' field
    elem = atom.get("element", "")
    if not elem and "type" in atom:
        elem = atom["type"].split(".")[0]
    if not elem:
        elem = "".join(c for c in atom.get("name", "C") if c.isalpha())[:1]
    return elem.upper()


def _compute_energy(lig_coords: np.ndarray, pocket_coords: np.ndarray, direction: np.ndarray) -> float:
    """
    Energy = -sum over pocket atoms of (min distance to ligand pulling lines).

    For each pocket atom p, for each ligand atom l, a line passes through l
    in direction v (unit vector). Point-to-line distance:
        d = || (p - l) x v ||

    For each pocket atom, take min over all ligand lines, then sum.
    Negative so lower energy = more clearance = better direction.
    """
    # diff: (N_pocket, N_lig, 3)
    diff = pocket_coords[:, np.newaxis, :] - lig_coords[np.newaxis, :, :]
    # cross product with direction: (N_pocket, N_lig, 3)
    cross = np.cross(diff, direction)
    # distance: (N_pocket, N_lig)
    dist = np.linalg.norm(cross, axis=2)
    # min distance per pocket atom, then sum, negate
    return -np.sum(np.min(dist, axis=1))


def _spherical_to_cartesian(theta: float, phi: float) -> np.ndarray:
    """Convert spherical coordinates to unit vector."""
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def _cartesian_to_spherical(v: np.ndarray) -> Tuple[float, float]:
    """Convert unit vector to spherical coordinates (theta, phi)."""
    theta = np.arccos(np.clip(v[2], -1.0, 1.0))
    phi = np.arctan2(v[1], v[0])
    return theta, phi


def _metropolis_hastings(
    lig_coords: np.ndarray,
    pocket_coords: np.ndarray,
    init_direction: np.ndarray,
    n_iterations: int = 100000,
    initial_temp: float = 2.0,
    cooling_rate: float = 0.99995,
    step_size: float = 0.15,
    conv_window: int = 8000,
    conv_threshold: float = 5e-4,
    verbose: bool = True,
) -> Tuple[np.ndarray, float, list]:
    """
    Minimize energy via Metropolis-Hastings on the unit sphere.

    Energy is negative; lower (more negative) = better clearance.
    Accept lower energy always; accept higher energy with prob exp(-dE/T).

    Returns:
        (best_direction, best_energy, energy_history)
    """
    theta, phi = _cartesian_to_spherical(init_direction)
    current_dir = _spherical_to_cartesian(theta, phi)
    current_energy = _compute_energy(lig_coords, pocket_coords, current_dir)

    best_dir = current_dir.copy()
    best_energy = current_energy
    temp = initial_temp

    energy_history = []
    accept_count = 0
    bar_width = 50

    for i in range(n_iterations):
        # Propose new angles
        new_theta = theta + np.random.normal(0, step_size)
        new_phi = phi + np.random.normal(0, step_size)

        # Wrap theta to [0, pi]
        new_theta = new_theta % (2 * np.pi)
        if new_theta > np.pi:
            new_theta = 2 * np.pi - new_theta
            new_phi += np.pi
        new_phi = new_phi % (2 * np.pi)

        new_dir = _spherical_to_cartesian(new_theta, new_phi)
        new_energy = _compute_energy(lig_coords, pocket_coords, new_dir)

        # Standard Metropolis-Hastings (minimizing energy)
        delta_e = new_energy - current_energy
        if delta_e < 0 or np.random.random() < np.exp(-delta_e / max(temp, 1e-10)):
            theta, phi = new_theta, new_phi
            current_dir = new_dir
            current_energy = new_energy
            accept_count += 1

            if current_energy < best_energy:
                best_energy = current_energy
                best_dir = current_dir.copy()

        energy_history.append(current_energy)
        temp *= cooling_rate

        # Progress bar (update every 500 steps)
        if verbose and ((i + 1) % 500 == 0 or i == n_iterations - 1):
            pct = (i + 1) / n_iterations
            filled = int(bar_width * pct)
            bar = "\u2588" * filled + "\u2591" * (bar_width - filled)
            acc_rate = accept_count / (i + 1)
            sys.stdout.write(
                f"\r  [{bar}] {pct*100:5.1f}% | "
                f"E={current_energy:.2f} | Best={best_energy:.2f} | "
                f"Accept={acc_rate:.3f} | T={temp:.4e}"
            )
            sys.stdout.flush()

        # Convergence check
        if len(energy_history) >= conv_window:
            recent = energy_history[-conv_window:]
            std = np.std(recent)
            mean_abs = abs(np.mean(recent))
            rel_std = std / (mean_abs + 1e-10)
            if rel_std < conv_threshold:
                if verbose:
                    pct = (i + 1) / n_iterations
                    filled = int(bar_width * pct)
                    bar = "\u2588" * filled + "\u2591" * (bar_width - filled)
                    sys.stdout.write(
                        f"\r  [{bar}] {pct*100:5.1f}% | "
                        f"E={current_energy:.2f} | Best={best_energy:.2f} | "
                        f"Accept={accept_count/(i+1):.3f} | T={temp:.4e}"
                    )
                    print(f"\n  Converged at iteration {i+1} " f"(relative std = {rel_std:.2e} < {conv_threshold:.0e})")
                return best_dir, best_energy, energy_history

    if verbose:
        print(f"\n  Reached maximum iterations ({n_iterations})")
    return best_dir, best_energy, energy_history


class PMFAligner:
    """
    Aligns protein-ligand complex for PMF calculations.

    Uses Metropolis-Hastings optimization to find the pulling direction that
    maximizes clearance from binding pocket atoms, then rotates the complex
    so this optimized direction aligns with the Z-axis.
    """

    def __init__(self, pocket_cutoff: float = 4.0, verbose: bool = True):
        """
        Initialize PMF Aligner.

        Parameters:
        -----------
        pocket_cutoff : float
            Distance cutoff (Angstroms) for defining pocket residues around ligand.
            Default: 4.0 A
        verbose : bool
            Whether to print detailed output
        """
        self.pocket_cutoff = pocket_cutoff
        self.verbose = verbose

    def align_for_pmf(
        self, protein_pdb: str, ligand_file: str, output_dir: str, pullvec: Optional[Tuple[int, int]] = None
    ) -> Dict:
        """
        Align protein-ligand complex for PMF calculations.

        When pullvec is None (auto mode), uses Metropolis-Hastings optimization
        to find the pulling direction that maximizes clearance from pocket atoms.

        Parameters:
        -----------
        protein_pdb : str
            Path to protein PDB file
        ligand_file : str
            Path to ligand file (MOL2 or SDF)
        output_dir : str
            Output directory for aligned structures
        pullvec : tuple of (int, int), optional
            User-defined pull vector as (protein_atom_index, ligand_atom_index).
            If None, uses MH-optimized direction.

        Returns:
        --------
        Dict : Alignment results including paths, vectors, and optimization info
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if self.verbose:
            print("\n=== PMF Structure Alignment ===")

        # Step 1: Parse structures
        protein_atoms = self._parse_pdb(protein_pdb)
        ligand_atoms = self._parse_ligand(ligand_file)

        if self.verbose:
            print(f"  Protein atoms: {len(protein_atoms)}")
            print(f"  Ligand atoms: {len(ligand_atoms)}")

        # Step 2: Calculate pull vector
        optimization_info = None

        if pullvec is not None:
            # User-defined pull vector - skip optimization
            prot_atom_idx, lig_atom_idx = pullvec
            pull_start, pull_end = self._get_pullvec_from_atoms(
                protein_atoms, ligand_atoms, prot_atom_idx, lig_atom_idx
            )
            pull_vector = pull_end - pull_start
            pull_vector_normalized = pull_vector / np.linalg.norm(pull_vector)
            pocket_centroid = pull_start

            if self.verbose:
                print(f"  Using user-defined pull vector:")
                print(f"    Protein atom {prot_atom_idx} -> Ligand atom {lig_atom_idx}")
        else:
            # Auto mode: MH-optimized direction
            # Extract coordinates for heavy atoms
            lig_coords = np.array([[a["x"], a["y"], a["z"]] for a in ligand_atoms])

            # Find pocket residues (vectorized: any heavy atom within cutoff of any ligand atom)
            pocket_atom_list, pocket_reskeys = self._find_pocket_residues(protein_atoms, lig_coords)
            pocket_coords = np.array([[a["x"], a["y"], a["z"]] for a in pocket_atom_list])
            pocket_centroid = np.mean(pocket_coords, axis=0) if len(pocket_coords) > 0 else np.zeros(3)

            if self.verbose:
                print(f"  Pocket residues: {len(pocket_reskeys)}")
                print(f"  Pocket heavy atoms: {len(pocket_atom_list)}")
                # List pocket residues
                sorted_keys = sorted(pocket_reskeys, key=lambda x: (x[0], x[1]))
                res_labels = []
                for chain, resid in sorted_keys:
                    resname = next(
                        a["resname"] for a in protein_atoms if a.get("chain", "") == chain and a["resid"] == resid
                    )
                    res_labels.append(f"{resname}{resid}:{chain}")
                print(f"  Residues: {', '.join(res_labels)}")

            # Initial direction: protein heavy-atom COM -> ligand COM
            protein_heavy = [a for a in protein_atoms if _is_heavy(_get_element_from_atom(a))]
            protein_heavy_coords = np.array([[a["x"], a["y"], a["z"]] for a in protein_heavy])
            protein_com = np.mean(protein_heavy_coords, axis=0)
            lig_com = np.mean(lig_coords, axis=0)

            init_vec = lig_com - protein_com
            init_dir = init_vec / np.linalg.norm(init_vec)

            if self.verbose:
                print(f"\n  Initial direction (protein COM -> ligand COM):")
                print(f"    ({init_dir[0]:.6f}, {init_dir[1]:.6f}, {init_dir[2]:.6f})")

            init_energy = _compute_energy(lig_coords, pocket_coords, init_dir)
            if self.verbose:
                print(f"  Initial energy: {init_energy:.4f}")

            # Run MH optimization
            if self.verbose:
                print(f"\n  Running Metropolis-Hastings optimization ...\n")

            best_dir, best_energy, history = _metropolis_hastings(
                lig_coords,
                pocket_coords,
                init_dir,
                n_iterations=100000,
                initial_temp=2.0,
                cooling_rate=0.99995,
                step_size=0.15,
                conv_window=8000,
                conv_threshold=5e-4,
                verbose=self.verbose,
            )

            # Fix direction polarity: energy function is symmetric (||a x v|| == ||a x -v||),
            # so MH may converge to the reverse direction. Ensure it points roughly
            # from protein toward ligand (positive dot product with init_dir).
            if np.dot(best_dir, init_dir) < 0:
                best_dir = -best_dir

            improvement = (best_energy - init_energy) / abs(init_energy) * 100 if abs(init_energy) > 1e-10 else 0.0

            if self.verbose:
                print(f"\n  Optimized direction: ({best_dir[0]:.6f}, {best_dir[1]:.6f}, {best_dir[2]:.6f})")
                print(f"  Energy: {init_energy:.4f} -> {best_energy:.4f} ({improvement:+.2f}%)")

            pull_vector_normalized = best_dir

            optimization_info = {
                "initial_direction": init_dir,
                "initial_energy": float(init_energy),
                "optimized_energy": float(best_energy),
                "improvement_pct": float(improvement),
                "pocket_residues": len(pocket_reskeys),
                "pocket_atoms": len(pocket_atom_list),
            }

            # Generate PyMOL visualization script (before rotation, in original coordinates)
            pml_path = output_dir / f"visualize_{Path(protein_pdb).stem}.pml"
            self._generate_pymol_script(
                str(pml_path),
                lig_com,
                init_dir,
                best_dir,
                os.path.abspath(protein_pdb),
                os.path.abspath(ligand_file),
                pocket_reskeys,
            )
            if self.verbose:
                print(f"  PyMOL script: {pml_path}")

        # Step 3: Rotation matrix to align pull vector with +Z axis
        z_axis = np.array([0.0, 0.0, 1.0])
        rotation_matrix = self._calculate_rotation_matrix(pull_vector_normalized, z_axis)

        if self.verbose:
            print(
                f"\n  Pull vector: ({pull_vector_normalized[0]:.3f}, {pull_vector_normalized[1]:.3f}, {pull_vector_normalized[2]:.3f})"
            )
            print(f"  Aligning to Z-axis...")

        # Step 4: Rotate around the geometric center of the entire complex
        # This minimizes spatial displacement of the structure
        all_atoms = protein_atoms + ligand_atoms
        all_coords_arr = np.array([[a["x"], a["y"], a["z"]] for a in all_atoms])
        system_center = np.mean(all_coords_arr, axis=0)

        if self.verbose:
            print(
                f"  Rotation center (system geometric center): "
                f"({system_center[0]:.2f}, {system_center[1]:.2f}, {system_center[2]:.2f})"
            )

        # Apply rotation: new_coord = R @ (coord - center) + center
        transformed_protein = self._rotate_atoms_around_center(protein_atoms, rotation_matrix, system_center)
        transformed_ligand = self._rotate_atoms_around_center(ligand_atoms, rotation_matrix, system_center)

        # Step 5: Shift all coordinates to be positive (required for GROMACS)
        all_coords = [[a["x"], a["y"], a["z"]] for a in transformed_protein + transformed_ligand]
        min_x = min(c[0] for c in all_coords)
        min_y = min(c[1] for c in all_coords)
        min_z = min(c[2] for c in all_coords)

        margin = 5.0
        shift = np.array([-min_x + margin, -min_y + margin, -min_z + margin])

        for atom in transformed_protein:
            atom["x"] += shift[0]
            atom["y"] += shift[1]
            atom["z"] += shift[2]

        for atom in transformed_ligand:
            atom["x"] += shift[0]
            atom["y"] += shift[1]
            atom["z"] += shift[2]

        if self.verbose:
            protein_z_coords = [a["z"] for a in transformed_protein]
            ligand_z_coords = [a["z"] for a in transformed_ligand]
            print(f"  Final positions (after shift to positive coords):")
            print(f"    Protein Z range: {min(protein_z_coords):.2f} to {max(protein_z_coords):.2f} Å")
            print(f"    Ligand Z range: {min(ligand_z_coords):.2f} to {max(ligand_z_coords):.2f} Å")
            print(
                f"    Ligand is {'above' if np.mean(ligand_z_coords) > np.mean(protein_z_coords) else 'at top of'} protein - ready for +Z pulling"
            )

        # Step 8: Write aligned structures
        aligned_protein_path = output_dir / f"{Path(protein_pdb).stem}_aligned.pdb"
        aligned_ligand_path = output_dir / f"{Path(ligand_file).stem}_aligned{Path(ligand_file).suffix}"

        self._write_pdb(transformed_protein, str(aligned_protein_path))
        self._write_ligand(transformed_ligand, ligand_file, str(aligned_ligand_path))

        if self.verbose:
            print(f"  Aligned protein: {aligned_protein_path}")
            print(f"  Aligned ligand: {aligned_ligand_path}")
            print(f"  System is ready for SMD pulling in +Z direction")

        result = {
            "aligned_protein": str(aligned_protein_path),
            "aligned_ligand": str(aligned_ligand_path),
            "pull_vector": pull_vector_normalized,
            "rotation_matrix": rotation_matrix,
            "rotation_center": system_center,
            "ligand_centroid": lig_com if pullvec is None else pull_end,
            "pocket_centroid": pocket_centroid,
        }

        if optimization_info is not None:
            result["optimization"] = optimization_info
            result["pymol_script"] = str(pml_path)

        return result

    def _parse_pdb(self, pdb_file: str) -> List[Dict]:
        """Parse PDB file and extract atom information."""
        atoms = []
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        atom = {
                            "record": line[:6].strip(),
                            "serial": int(line[6:11].strip()),
                            "name": line[12:16].strip(),
                            "altloc": line[16].strip(),
                            "resname": line[17:20].strip(),
                            "chain": line[21].strip(),
                            "resid": int(line[22:26].strip()),
                            "icode": line[26].strip(),
                            "x": float(line[30:38].strip()),
                            "y": float(line[38:46].strip()),
                            "z": float(line[46:54].strip()),
                            "occupancy": float(line[54:60].strip()) if len(line) > 54 and line[54:60].strip() else 1.0,
                            "tempfactor": float(line[60:66].strip()) if len(line) > 60 and line[60:66].strip() else 0.0,
                            "element": line[76:78].strip() if len(line) > 76 else "",
                            "charge": line[78:80].strip() if len(line) > 78 else "",
                            "original_line": line,
                        }
                        atoms.append(atom)
                    except (ValueError, IndexError) as e:
                        logger.warning(f"Could not parse PDB line: {line.strip()}")
                        continue
        return atoms

    def _parse_ligand(self, ligand_file: str) -> List[Dict]:
        """Parse ligand file (MOL2 or SDF) and extract atom information."""
        suffix = Path(ligand_file).suffix.lower()

        if suffix == ".mol2":
            return self._parse_mol2(ligand_file)
        elif suffix in [".sdf", ".sd"]:
            return self._parse_sdf(ligand_file)
        else:
            raise ValueError(f"Unsupported ligand format: {suffix}")

    def _parse_mol2(self, mol2_file: str) -> List[Dict]:
        """Parse MOL2 file."""
        atoms = []
        in_atom_section = False

        with open(mol2_file, "r") as f:
            for line in f:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atom_section = True
                    continue
                elif line.startswith("@<TRIPOS>"):
                    in_atom_section = False
                    continue

                if in_atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            atom = {
                                "serial": int(parts[0]),
                                "name": parts[1],
                                "x": float(parts[2]),
                                "y": float(parts[3]),
                                "z": float(parts[4]),
                                "type": parts[5],
                                "resid": int(parts[6]) if len(parts) > 6 else 1,
                                "resname": parts[7] if len(parts) > 7 else "LIG",
                                "charge": float(parts[8]) if len(parts) > 8 else 0.0,
                                "original_parts": parts,
                            }
                            atoms.append(atom)
                        except (ValueError, IndexError) as e:
                            logger.warning(f"Could not parse MOL2 line: {line.strip()}")
                            continue
        return atoms

    def _parse_sdf(self, sdf_file: str) -> List[Dict]:
        """Parse SDF file."""
        atoms = []

        with open(sdf_file, "r") as f:
            lines = f.readlines()

        # Find the atom count line (4th line typically)
        if len(lines) < 4:
            raise ValueError("Invalid SDF file: too few lines")

        # Parse counts line
        counts_line = lines[3]
        try:
            num_atoms = int(counts_line[:3].strip())
        except ValueError:
            raise ValueError(f"Cannot parse atom count from SDF: {counts_line}")

        # Parse atom block
        for i in range(4, 4 + num_atoms):
            if i >= len(lines):
                break
            line = lines[i]
            try:
                atom = {
                    "serial": i - 3,
                    "x": float(line[0:10].strip()),
                    "y": float(line[10:20].strip()),
                    "z": float(line[20:30].strip()),
                    "element": line[31:34].strip(),
                    "name": f"{line[31:34].strip()}{i-3}",
                    "resid": 1,
                    "resname": "LIG",
                    "original_line": line,
                }
                atoms.append(atom)
            except (ValueError, IndexError) as e:
                logger.warning(f"Could not parse SDF line: {line.strip()}")
                continue

        return atoms

    def _generate_pymol_script(
        self,
        filename: str,
        origin: np.ndarray,
        init_dir: np.ndarray,
        opt_dir: np.ndarray,
        pdb_path: str,
        ligand_path: str,
        pocket_reskeys: set,
        arrow_length: float = 30.0,
    ):
        """
        Write a PyMOL .pml script that visualizes initial and optimized pull directions.

        Yellow arrow = initial direction (protein COM -> ligand COM)
        Red arrow = MH-optimized pulling direction
        Cyan sticks = binding pocket residues
        """
        shaft_radius = 0.3
        head_radius = 0.6
        head_length = 2.5

        init_shaft_end = origin + init_dir * (arrow_length - head_length)
        init_end = origin + init_dir * arrow_length
        opt_shaft_end = origin + opt_dir * (arrow_length - head_length)
        opt_end = origin + opt_dir * arrow_length

        # Build pocket selection string
        pocket_sele_parts = []
        for chain, resid in sorted(pocket_reskeys):
            pocket_sele_parts.append(f"(chain {chain} and resi {resid})")
        pocket_sele = " or ".join(pocket_sele_parts) if pocket_sele_parts else "none"

        with open(filename, "w") as f:
            f.write(f"""\
# ============================================================
# SMD Pulling Direction Visualization
# Generated by PRISM PMF Aligner
# ============================================================
# Yellow arrow: initial direction (protein COM -> ligand COM)
# Red arrow:    optimized SMD pulling direction
# ============================================================

load {pdb_path}, protein
load {ligand_path}, ligand

hide everything, protein
show cartoon, protein
color gray80, protein

show sticks, ligand
color green, ligand
set stick_radius, 0.15, ligand

# --- Highlight pocket residues ---
select pocket, {pocket_sele}
show sticks, pocket
color cyan, pocket and elem C
set stick_radius, 0.12, pocket
deselect

# --- Ligand COM ---
pseudoatom lig_com, pos=[{origin[0]:.4f}, {origin[1]:.4f}, {origin[2]:.4f}]
show spheres, lig_com
color white, lig_com
set sphere_scale, 0.5, lig_com

python
from pymol.cgo import CYLINDER, CONE
from pymol import cmd

# --- Initial direction arrow (yellow) ---
arrow_init = [
    CYLINDER,
    {origin[0]:.4f}, {origin[1]:.4f}, {origin[2]:.4f},
    {init_shaft_end[0]:.4f}, {init_shaft_end[1]:.4f}, {init_shaft_end[2]:.4f},
    {shaft_radius},
    1.0, 0.9, 0.0,
    1.0, 0.9, 0.0,
    CONE,
    {init_shaft_end[0]:.4f}, {init_shaft_end[1]:.4f}, {init_shaft_end[2]:.4f},
    {init_end[0]:.4f}, {init_end[1]:.4f}, {init_end[2]:.4f},
    {head_radius}, 0.0,
    1.0, 0.9, 0.0,
    1.0, 0.9, 0.0,
    1.0, 0.0,
]
cmd.load_cgo(arrow_init, "vec_initial")

# --- Optimized direction arrow (red) ---
arrow_opt = [
    CYLINDER,
    {origin[0]:.4f}, {origin[1]:.4f}, {origin[2]:.4f},
    {opt_shaft_end[0]:.4f}, {opt_shaft_end[1]:.4f}, {opt_shaft_end[2]:.4f},
    {shaft_radius},
    1.0, 0.2, 0.2,
    1.0, 0.2, 0.2,
    CONE,
    {opt_shaft_end[0]:.4f}, {opt_shaft_end[1]:.4f}, {opt_shaft_end[2]:.4f},
    {opt_end[0]:.4f}, {opt_end[1]:.4f}, {opt_end[2]:.4f},
    {head_radius}, 0.0,
    1.0, 0.2, 0.2,
    1.0, 0.2, 0.2,
    1.0, 0.0,
]
cmd.load_cgo(arrow_opt, "vec_optimized")

print("\\n===== SMD Vector Visualization Loaded =====")
print("  Yellow arrow: initial direction (protein COM -> ligand COM)")
print("  Red arrow:    optimized SMD pulling direction")
print("  Cyan sticks:  binding pocket")
print("============================================")
python end

# --- Arrow endpoint labels ---
pseudoatom init_end, pos=[{init_end[0]:.4f}, {init_end[1]:.4f}, {init_end[2]:.4f}]
pseudoatom opt_end, pos=[{opt_end[0]:.4f}, {opt_end[1]:.4f}, {opt_end[2]:.4f}]
hide everything, init_end
hide everything, opt_end
label init_end, "Initial"
label opt_end, "Optimized"
set label_color, yellow, init_end
set label_color, red, opt_end
set label_size, 20

# --- Camera ---
zoom ligand, 15
bg_color white
set ray_shadow, 0
""")

    def _calculate_centroid(self, atoms: List[Dict]) -> np.ndarray:
        """Calculate the centroid of a list of atoms."""
        if not atoms:
            return np.zeros(3)

        coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
        return np.mean(coords, axis=0)

    def _find_pocket_residues(self, protein_atoms: List[Dict], lig_coords: np.ndarray) -> Tuple[List[Dict], set]:
        """
        Find pocket residues with any heavy atom within cutoff of any ligand atom.

        Uses vectorized distance computation for efficiency.

        Returns:
            (pocket_heavy_atoms, pocket_reskeys) where pocket_reskeys is
            a set of (chain, resid) tuples.
        """
        # Collect heavy protein atoms
        heavy_atoms = []
        for a in protein_atoms:
            if _is_heavy(_get_element_from_atom(a)):
                heavy_atoms.append(a)

        heavy_coords = np.array([[a["x"], a["y"], a["z"]] for a in heavy_atoms])

        # Vectorized: (N_prot, 1, 3) - (1, N_lig, 3) -> (N_prot, N_lig)
        diff = heavy_coords[:, np.newaxis, :] - lig_coords[np.newaxis, :, :]
        dists = np.linalg.norm(diff, axis=2)  # (N_prot, N_lig)
        min_dists = np.min(dists, axis=1)  # (N_prot,)

        # Mark residues within cutoff
        pocket_reskeys = set()
        for i, a in enumerate(heavy_atoms):
            if min_dists[i] <= self.pocket_cutoff:
                pocket_reskeys.add((a.get("chain", ""), a["resid"]))

        # Collect all heavy atoms of pocket residues
        pocket_atoms = [a for a in heavy_atoms if (a.get("chain", ""), a["resid"]) in pocket_reskeys]
        return pocket_atoms, pocket_reskeys

    def _get_pullvec_from_atoms(
        self, protein_atoms: List[Dict], ligand_atoms: List[Dict], prot_atom_idx: int, lig_atom_idx: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Get pull vector from user-specified atom indices."""
        # Find protein atom by serial number
        prot_atom = None
        for atom in protein_atoms:
            if atom["serial"] == prot_atom_idx:
                prot_atom = atom
                break

        if prot_atom is None:
            raise ValueError(f"Protein atom with index {prot_atom_idx} not found")

        # Find ligand atom by serial number
        lig_atom = None
        for atom in ligand_atoms:
            if atom["serial"] == lig_atom_idx:
                lig_atom = atom
                break

        if lig_atom is None:
            raise ValueError(f"Ligand atom with index {lig_atom_idx} not found")

        pull_start = np.array([prot_atom["x"], prot_atom["y"], prot_atom["z"]])
        pull_end = np.array([lig_atom["x"], lig_atom["y"], lig_atom["z"]])

        return pull_start, pull_end

    def _calculate_rotation_matrix(self, v_from: np.ndarray, v_to: np.ndarray) -> np.ndarray:
        """
        Calculate rotation matrix to rotate vector v_from to align with v_to.

        Uses Rodrigues' rotation formula.
        """
        # Normalize vectors
        v_from = v_from / np.linalg.norm(v_from)
        v_to = v_to / np.linalg.norm(v_to)

        # Handle case where vectors are already aligned
        dot = np.dot(v_from, v_to)
        if np.abs(dot - 1.0) < 1e-8:
            return np.eye(3)

        # Handle case where vectors are opposite
        if np.abs(dot + 1.0) < 1e-8:
            # Find perpendicular vector
            if np.abs(v_from[0]) < 0.9:
                perp = np.array([1, 0, 0])
            else:
                perp = np.array([0, 1, 0])
            perp = perp - np.dot(perp, v_from) * v_from
            perp = perp / np.linalg.norm(perp)
            # 180 degree rotation around perpendicular axis
            return 2 * np.outer(perp, perp) - np.eye(3)

        # Cross product gives rotation axis
        cross = np.cross(v_from, v_to)
        cross_norm = np.linalg.norm(cross)

        if cross_norm < 1e-8:
            return np.eye(3)

        # Skew-symmetric cross-product matrix
        k = cross / cross_norm
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])

        # Rodrigues' formula: R = I + sin(theta)*K + (1-cos(theta))*K^2
        sin_theta = cross_norm
        cos_theta = dot

        R = np.eye(3) + sin_theta * K + (1 - cos_theta) * np.dot(K, K)

        return R

    def _rotate_atoms_around_center(
        self, atoms: List[Dict], rotation_matrix: np.ndarray, center: np.ndarray
    ) -> List[Dict]:
        """
        Rotate atoms around a center point.

        new_coord = R @ (coord - center) + center

        This preserves the overall spatial position of the structure
        while aligning the pull vector with the target axis.
        """
        transformed = []
        for atom in atoms:
            coord = np.array([atom["x"], atom["y"], atom["z"]])
            new_coord = np.dot(rotation_matrix, coord - center) + center

            new_atom = atom.copy()
            new_atom["x"] = new_coord[0]
            new_atom["y"] = new_coord[1]
            new_atom["z"] = new_coord[2]
            transformed.append(new_atom)

        return transformed

    def _write_pdb(self, atoms: List[Dict], output_file: str):
        """Write atoms to PDB format."""
        with open(output_file, "w") as f:
            f.write("REMARK   Generated by PRISM PMF Aligner\n")
            f.write("REMARK   Structure aligned with pull vector along Z-axis\n")

            for atom in atoms:
                record = atom.get("record", "ATOM")
                serial = atom["serial"]
                name = atom["name"]
                altloc = atom.get("altloc", "")
                resname = atom["resname"]
                chain = atom.get("chain", "A")
                resid = atom["resid"]
                icode = atom.get("icode", "")
                x, y, z = atom["x"], atom["y"], atom["z"]
                occupancy = atom.get("occupancy", 1.0)
                tempfactor = atom.get("tempfactor", 0.0)
                element = atom.get("element", name[0])

                # Format atom name with proper spacing
                if len(name) < 4:
                    name_fmt = f" {name:<3s}"
                else:
                    name_fmt = f"{name:<4s}"

                line = f"{record:<6s}{serial:5d} {name_fmt}{altloc:1s}{resname:>3s} {chain:1s}{resid:4d}{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{tempfactor:6.2f}          {element:>2s}\n"
                f.write(line)

            f.write("END\n")

    def _write_ligand(self, atoms: List[Dict], original_file: str, output_file: str):
        """Write ligand atoms to original format."""
        suffix = Path(original_file).suffix.lower()

        if suffix == ".mol2":
            self._write_mol2(atoms, original_file, output_file)
        elif suffix in [".sdf", ".sd"]:
            self._write_sdf(atoms, original_file, output_file)

    def _write_mol2(self, atoms: List[Dict], original_file: str, output_file: str):
        """Write atoms to MOL2 format, preserving original structure."""
        # Read original file to preserve non-atom sections
        with open(original_file, "r") as f:
            original_lines = f.readlines()

        with open(output_file, "w") as f:
            in_atom_section = False
            atom_idx = 0

            for line in original_lines:
                if line.startswith("@<TRIPOS>ATOM"):
                    f.write(line)
                    in_atom_section = True
                    continue
                elif line.startswith("@<TRIPOS>"):
                    in_atom_section = False
                    f.write(line)
                    continue

                if in_atom_section and line.strip() and atom_idx < len(atoms):
                    atom = atoms[atom_idx]
                    parts = atom.get("original_parts", line.split())

                    # Update coordinates
                    if len(parts) >= 6:
                        parts[2] = f"{atom['x']:.4f}"
                        parts[3] = f"{atom['y']:.4f}"
                        parts[4] = f"{atom['z']:.4f}"
                        f.write("      ".join(parts[:3]) + "   " + "   ".join(parts[3:]) + "\n")
                    else:
                        f.write(line)
                    atom_idx += 1
                else:
                    f.write(line)

    def _write_sdf(self, atoms: List[Dict], original_file: str, output_file: str):
        """Write atoms to SDF format, preserving original structure."""
        with open(original_file, "r") as f:
            original_lines = f.readlines()

        if len(original_lines) < 4:
            raise ValueError("Invalid SDF file")

        # Parse atom count
        counts_line = original_lines[3]
        num_atoms = int(counts_line[:3].strip())

        with open(output_file, "w") as f:
            # Write header (lines 0-3)
            for i in range(4):
                f.write(original_lines[i])

            # Write atom block with updated coordinates
            for i in range(num_atoms):
                if i < len(atoms):
                    atom = atoms[i]
                    original_line = original_lines[4 + i] if 4 + i < len(original_lines) else ""
                    # Keep original line format but update coordinates
                    rest = (
                        original_line[30:] if len(original_line) > 30 else " C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                    )
                    f.write(f"{atom['x']:10.4f}{atom['y']:10.4f}{atom['z']:10.4f}{rest}")
                else:
                    f.write(original_lines[4 + i])

            # Write rest of file
            for i in range(4 + num_atoms, len(original_lines)):
                f.write(original_lines[i])


def align_complex_for_pmf(
    protein_pdb: str,
    ligand_file: str,
    output_dir: str,
    pullvec: Optional[Tuple[int, int]] = None,
    pocket_cutoff: float = 4.0,
) -> Dict:
    """
    Convenience function to align a protein-ligand complex for PMF calculations.

    Parameters:
    -----------
    protein_pdb : str
        Path to protein PDB file
    ligand_file : str
        Path to ligand file (MOL2 or SDF)
    output_dir : str
        Output directory for aligned structures
    pullvec : tuple of (int, int), optional
        User-defined pull vector as (protein_atom_index, ligand_atom_index)
    pocket_cutoff : float
        Distance cutoff (Angstroms) for defining pocket

    Returns:
    --------
    Dict : Alignment results
    """
    aligner = PMFAligner(pocket_cutoff=pocket_cutoff)
    return aligner.align_for_pmf(protein_pdb, ligand_file, output_dir, pullvec)
