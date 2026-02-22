"""
Topology operations for REST2 setup.

Handles GRO parsing, pocket residue detection, topology merging,
and solute atom marking -- all in pure Python (no MDAnalysis dependency).
Only requires `gmx` (GROMACS) for the grompp -pp step.
"""

import os
import re
import subprocess
import tempfile
from typing import Dict, List, Optional, Set, Tuple

_SECTION_RE = re.compile(r'^\s*\[\s*(\w+)\s*\]')


# ---------------------------------------------------------------------------
# GRO file parsing
# ---------------------------------------------------------------------------

def parse_gro(path: str) -> Tuple[List[dict], List[float]]:
    """Parse a GROMACS .gro file (fixed-width format).

    GRO format (fixed-width columns):
        residue number (5 chars), residue name (5 chars), atom name (5 chars),
        atom number (5 chars), x (8.3f), y (8.3f), z (8.3f)

    Returns:
        atoms: list of dicts with keys: resnr, resname, atomname, atomnr, x, y, z
        box: [bx, by, bz] in nm
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    n_atoms = int(lines[1].strip())
    atoms = []

    for i in range(2, 2 + n_atoms):
        line = lines[i]
        # Fixed-width parsing
        resnr = int(line[0:5].strip())
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atomnr = int(line[15:20].strip())
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())

        atoms.append({
            'resnr': resnr,
            'resname': resname,
            'atomname': atomname,
            'atomnr': atomnr,
            'x': x, 'y': y, 'z': z,
        })

    # Parse box from last line
    box_line = lines[2 + n_atoms].split()
    box = [float(v) for v in box_line[:3]]

    return atoms, box


# ---------------------------------------------------------------------------
# Pocket residue detection
# ---------------------------------------------------------------------------

def _minimum_image_dist(dx: float, box_len: float) -> float:
    """Apply minimum image convention to a single coordinate difference."""
    dx -= box_len * round(dx / box_len)
    return dx


def find_pocket_residues(
    atoms: List[dict],
    box: List[float],
    lig_resname,
    molecules_order: List[Tuple[str, int, int]],
    cutoff_nm: float = 0.5,
) -> Dict[str, Set[int]]:
    """Find protein residues within cutoff of the ligand(s).

    Uses minimum-image convention for periodic boundary conditions.
    Groups results by chain (moleculetype name).

    Args:
        atoms: Parsed GRO atoms
        box: Box dimensions [bx, by, bz] in nm
        lig_resname: Ligand residue name(s) — str or list of str
                     (e.g. 'LIG' or ['LIG_1', 'LIG_2'])
        molecules_order: List of (moltype_name, resname_or_chain, n_atoms) from
                         the topology -- used to assign atoms to moleculetypes
        cutoff_nm: Distance cutoff in nm

    Returns:
        Dict mapping chain moleculetype name -> set of residue numbers
    """
    # Normalize to set for O(1) lookup
    if isinstance(lig_resname, str):
        lig_resnames = {lig_resname}
    else:
        lig_resnames = set(lig_resname)

    cutoff_sq = cutoff_nm * cutoff_nm
    bx, by, bz = box

    # Identify ligand atoms
    lig_atoms = []
    for atom in atoms:
        if atom['resname'] in lig_resnames:
            lig_atoms.append((atom['x'], atom['y'], atom['z']))

    if not lig_atoms:
        raise ValueError(f"No atoms found with resname(s) {sorted(lig_resnames)} in GRO file")

    # Build mapping: atom index -> moleculetype name using molecules_order
    # molecules_order: [(moltype_name, count, n_atoms_per_molecule), ...]
    atom_to_moltype = {}
    atom_idx = 0
    for moltype_name, count, n_atoms_per_mol in molecules_order:
        for _ in range(count):
            for _ in range(n_atoms_per_mol):
                if atom_idx < len(atoms):
                    atom_to_moltype[atom_idx] = moltype_name
                atom_idx += 1

    # Find protein residues near ligand(s)
    pocket: Dict[str, Set[int]] = {}
    for idx, atom in enumerate(atoms):
        if atom['resname'] in lig_resnames:
            continue
        # Only consider protein chains (those with Protein_ prefix)
        moltype = atom_to_moltype.get(idx)
        if moltype is None or not moltype.startswith('Protein_chain_'):
            continue

        ax, ay, az = atom['x'], atom['y'], atom['z']
        for lx, ly, lz in lig_atoms:
            dx = _minimum_image_dist(ax - lx, bx)
            dy = _minimum_image_dist(ay - ly, by)
            dz = _minimum_image_dist(az - lz, bz)
            dist_sq = dx*dx + dy*dy + dz*dz
            if dist_sq <= cutoff_sq:
                if moltype not in pocket:
                    pocket[moltype] = set()
                pocket[moltype].add(atom['resnr'])
                break  # This atom is close enough, move to next

    return pocket


# ---------------------------------------------------------------------------
# Topology merging (gmx grompp -pp)
# ---------------------------------------------------------------------------

def merge_topology(topol: str, gro: str, output: str, gmx: str = 'gmx') -> None:
    """Run gmx grompp -pp to produce a fully-expanded processed topology.

    Creates a temporary MDP file with minimal settings for preprocessing.

    Args:
        topol: Path to input topol.top
        gro: Path to input .gro file
        output: Path for output processed topology
        gmx: Path to gmx executable
    """
    topol_dir = os.path.dirname(os.path.abspath(topol))

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write minimal MDP
        mdp_path = os.path.join(tmpdir, 'minimal.mdp')
        with open(mdp_path, 'w') as f:
            f.write('integrator = steep\nnsteps = 0\n')

        tpr_path = os.path.join(tmpdir, 'topol.tpr')

        cmd = [
            gmx, 'grompp',
            '-f', mdp_path,
            '-c', gro,
            '-p', topol,
            '-pp', output,
            '-o', tpr_path,
            '-maxwarn', '10',
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=topol_dir,
        )

        if not os.path.exists(output):
            raise RuntimeError(
                f"gmx grompp -pp failed.\nstdout: {result.stdout}\nstderr: {result.stderr}"
            )

        print(f"  Processed topology written to: {output}")


# ---------------------------------------------------------------------------
# Parse [ molecules ] section
# ---------------------------------------------------------------------------

def parse_molecules_section(topol_path: str) -> List[Tuple[str, int]]:
    """Parse the [ molecules ] section from a topology file.

    Returns:
        List of (moleculetype_name, count) tuples
    """
    molecules = []
    in_molecules = False

    with open(topol_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            m = _SECTION_RE.match(stripped)
            if m:
                in_molecules = (m.group(1).lower() == 'molecules')
                continue

            if not in_molecules:
                continue

            if not stripped or stripped.startswith(';'):
                continue

            # Remove inline comments
            data = stripped.split(';')[0].split()
            if len(data) >= 2:
                name = data[0]
                count = int(data[1])
                molecules.append((name, count))

    return molecules


# ---------------------------------------------------------------------------
# Count atoms per moleculetype in processed topology
# ---------------------------------------------------------------------------

def count_atoms_per_moltype(processed_top: str) -> Dict[str, int]:
    """Count the number of atoms in each moleculetype from a processed topology.

    Returns:
        Dict mapping moleculetype name -> number of atoms
    """
    counts: Dict[str, int] = {}
    current_moltype: Optional[str] = None
    current_section: Optional[str] = None

    with open(processed_top, 'r') as f:
        for line in f:
            stripped = line.strip()
            m = _SECTION_RE.match(stripped)
            if m:
                current_section = m.group(1).lower()
                if current_section == 'moleculetype':
                    current_moltype = None  # Will be set by next data line
                continue

            if not stripped or stripped.startswith(';') or stripped.startswith('#'):
                continue

            if current_section == 'moleculetype':
                # First data line: name nrexcl
                fields = stripped.split()
                if fields:
                    current_moltype = fields[0]
                    counts[current_moltype] = 0

            elif current_section == 'atoms' and current_moltype is not None:
                counts[current_moltype] = counts.get(current_moltype, 0) + 1

    return counts


# ---------------------------------------------------------------------------
# Mark solute atoms
# ---------------------------------------------------------------------------

def mark_solute_atoms(
    processed_top_path: str,
    output_path: str,
    lig_moltype,
    pocket_residues: Dict[str, Set[int]],
) -> None:
    """Mark solute atoms by adding '_' suffix to their atom types.

    For ligand moleculetype(s): mark ALL atoms.
    For protein chains in pocket_residues: mark only atoms whose resnr is in the set.
    Also duplicates [ atomtypes ] entries for marked types.

    Args:
        processed_top_path: Path to processed topology (from gmx grompp -pp)
        output_path: Path for output marked topology
        lig_moltype: Name(s) of the ligand moleculetype — str or list of str
                     (e.g. 'LIG' or ['LIG_1', 'LIG_2'])
        pocket_residues: Dict of chain_name -> set of residue numbers to mark
    """
    # Normalize to set for O(1) lookup
    if isinstance(lig_moltype, str):
        lig_moltypes = {lig_moltype}
    else:
        lig_moltypes = set(lig_moltype)
    with open(processed_top_path, 'r') as f:
        lines = f.readlines()

    # First pass: identify which types need '_' duplicates
    marked_types: Set[str] = set()
    current_section: Optional[str] = None
    current_moltype: Optional[str] = None

    for line in lines:
        stripped = line.strip()
        m = _SECTION_RE.match(stripped)
        if m:
            current_section = m.group(1).lower()
            if current_section == 'moleculetype':
                current_moltype = None
            continue

        if not stripped or stripped.startswith(';') or stripped.startswith('#'):
            continue

        if current_section == 'moleculetype':
            fields = stripped.split()
            if fields:
                current_moltype = fields[0]

        elif current_section == 'atoms' and current_moltype is not None:
            fields = stripped.split(';')[0].split()
            if len(fields) >= 7:
                atom_type = fields[1]
                resnr = int(fields[2])

                should_mark = False
                if current_moltype in lig_moltypes:
                    should_mark = True
                elif current_moltype in pocket_residues:
                    if resnr in pocket_residues[current_moltype]:
                        should_mark = True

                if should_mark:
                    marked_types.add(atom_type)

    # Second pass: rewrite the topology
    #
    # The processed topology from grompp -pp can have MULTIPLE occurrences of
    # the same type section (e.g. two [ dihedraltypes ] — one from the force
    # field, one from the ligand).  GROMACS 2025.1 requires that all funct=9
    # entries for the same type combination form ONE contiguous block.
    #
    # Strategy: on the FIRST occurrence of each type section, collect entries
    # from ALL occurrences across the whole file, merge them (grouped by
    # type_tuple+funct for contiguity), add REST2 '_' duplicates, and emit
    # the single merged section.  Subsequent occurrences are skipped entirely.

    # Which type sections to merge
    _TYPE_SECTIONS = {'atomtypes', 'bondtypes', 'angletypes', 'dihedraltypes', 'pairtypes'}
    _TYPE_N_FIELDS = {
        'bondtypes': 2, 'angletypes': 3, 'dihedraltypes': 4, 'pairtypes': 2,
    }

    output_lines: List[str] = []
    current_section: Optional[str] = None
    current_moltype: Optional[str] = None
    emitted_type_sections: Set[str] = set()

    for line in lines:
        stripped = line.strip()
        m = _SECTION_RE.match(stripped)
        if m:
            current_section = m.group(1).lower()

            if current_section in _TYPE_SECTIONS:
                if current_section not in emitted_type_sections:
                    # First time seeing this section — emit header + merged content
                    emitted_type_sections.add(current_section)
                    output_lines.append(line)
                    output_lines.extend(
                        _build_merged_type_section(
                            lines, current_section, marked_types,
                            _TYPE_N_FIELDS.get(current_section),
                        )
                    )
                # Skip header (both first — already appended above — and subsequent)
                continue

            if current_section == 'moleculetype':
                current_moltype = None
            output_lines.append(line)
            continue

        # Lines inside a type section are already emitted via the merge — skip
        if current_section in _TYPE_SECTIONS:
            continue

        if not stripped or stripped.startswith(';') or stripped.startswith('#'):
            output_lines.append(line)
            continue

        if current_section == 'moleculetype':
            fields = stripped.split()
            if fields:
                current_moltype = fields[0]
            output_lines.append(line)

        elif current_section == 'atoms' and current_moltype is not None:
            output_lines.append(_mark_atom_line(line, current_moltype, lig_moltypes, pocket_residues))

        else:
            output_lines.append(line)

    with open(output_path, 'w') as f:
        f.writelines(output_lines)

    print(f"  Marked topology written to: {output_path}")
    print(f"  Marked {len(marked_types)} atom types: {sorted(marked_types)}")


# ---------------------------------------------------------------------------
# Helpers: merge type sections and generate REST2 duplicates
# ---------------------------------------------------------------------------

def _build_merged_type_section(
    all_lines: List[str],
    section_name: str,
    marked_types: Set[str],
    n_type_fields: Optional[int],
) -> List[str]:
    """Collect entries from ALL occurrences of *section_name*, add REST2
    duplicates, and return a single merged block of lines.

    For [ atomtypes ]:  deduplicates by type name, adds '_' copies.
    For bonded types:   groups by (type_tuple, funct) for contiguity,
                        then adds all '_' combinations.
    """
    if section_name == 'atomtypes':
        return _merge_atomtypes(all_lines, marked_types)
    else:
        return _merge_bonded_types(all_lines, section_name, marked_types, n_type_fields)


def _merge_atomtypes(all_lines: List[str], marked_types: Set[str]) -> List[str]:
    """Merge all [ atomtypes ] entries and add REST2 '_' duplicates."""
    current_section: Optional[str] = None
    seen: Set[str] = set()
    result: List[str] = []

    # Collect unique entries from all [ atomtypes ] sections
    for line in all_lines:
        stripped = line.strip()
        m = _SECTION_RE.match(stripped)
        if m:
            current_section = m.group(1).lower()
            continue
        if current_section != 'atomtypes':
            continue
        if not stripped or stripped.startswith(';') or stripped.startswith('#'):
            continue

        fields = stripped.split(';')[0].split()
        if not fields:
            continue
        type_name = fields[0]
        if type_name in seen:
            continue
        seen.add(type_name)
        result.append(line)

    # Add '_' duplicates for marked types
    dup_lines: List[str] = []
    for line in result:
        stripped = line.strip()
        data = stripped.split(';')[0]
        comment = ''
        if ';' in stripped:
            comment = stripped[stripped.index(';'):]
        fields = data.split()
        if not fields or fields[0] not in marked_types:
            continue

        new_fields = list(fields)
        new_fields[0] = fields[0] + '_'
        if len(new_fields) > 1:
            try:
                float(new_fields[1])
            except ValueError:
                new_fields[1] = new_fields[1] + '_'

        dup_line = '  '.join(f'{f:>12}' if _looks_numeric(f) else f' {f}' for f in new_fields)
        if comment:
            dup_line += f'  {comment}'
        dup_lines.append(dup_line + '\n')

    if dup_lines:
        result.append('; REST2 solute atom type duplicates\n')
        result.extend(dup_lines)

    return result


def _canonical_bonded_key(type_tuple: tuple) -> tuple:
    """Return canonical form of a bonded type tuple, accounting for symmetry.

    GROMACS treats reversed bonded type combinations as equivalent:
      - bonds/pairs:  (A, B) = (B, A)
      - angles:       (A, B, C) = (C, B, A)
      - dihedrals:    (A, B, C, D) = (D, C, B, A)

    Returns the lexicographically smaller of the forward and reversed tuple.
    """
    rev = tuple(reversed(type_tuple))
    return min(type_tuple, rev)


def _merge_bonded_types(
    all_lines: List[str],
    section_name: str,
    marked_types: Set[str],
    n_type_fields: int,
) -> List[str]:
    """Merge all occurrences of a bonded-type section, group by
    (type_tuple, funct) for contiguity, and add REST2 '_' combinations.
    """
    from collections import OrderedDict

    current_section: Optional[str] = None

    # Collect entries grouped by (type_tuple, funct) across ALL occurrences
    groups: OrderedDict = OrderedDict()

    for line in all_lines:
        stripped = line.strip()
        m = _SECTION_RE.match(stripped)
        if m:
            current_section = m.group(1).lower()
            continue
        if current_section != section_name:
            continue
        if not stripped or stripped.startswith(';') or stripped.startswith('#'):
            continue

        data = stripped.split(';')[0]
        comment = ''
        if ';' in stripped:
            comment = stripped[stripped.index(';'):]
        fields = data.split()
        if len(fields) < n_type_fields + 1:
            continue

        type_tuple = tuple(fields[:n_type_fields])
        funct = fields[n_type_fields]
        key = (type_tuple, funct)
        groups.setdefault(key, []).append((list(fields), comment))

    # Emit original entries — grouped, so multi-term funct=9 blocks are contiguous
    result: List[str] = []
    for (type_tuple, funct), entries in groups.items():
        for fields, comment in entries:
            dup_line = '  '.join(
                f'{f:>12}' if _looks_numeric(f) else f'  {f}'
                for f in fields
            )
            if comment:
                dup_line += f'  {comment}'
            result.append(dup_line + '\n')

    # Generate REST2 '_' duplicate combinations.
    #
    # IMPORTANT: GROMACS treats bonded types as symmetric — reversed type tuples
    # are considered the same type (e.g. A-B-C-D = D-C-B-A for dihedrals).
    # We must deduplicate reverse-equivalent combos so that GROMACS never sees
    # non-contiguous blocks for what it considers the same type.
    #
    # Example: original (C8, C8, C8, C8) generates combos:
    #   combo 1: (C8_, C8, C8, C8) and combo 8: (C8, C8, C8, C8_)
    # These are reverses of each other → same type in GROMACS → keep only one.

    dup_groups: OrderedDict = OrderedDict()  # (canonical_key, funct) -> [(fields, comment)]

    for (type_tuple, funct), entries in groups.items():
        marked_positions = [
            i for i in range(n_type_fields)
            if type_tuple[i] in marked_types and type_tuple[i] != 'X'
        ]
        if not marked_positions:
            continue

        n = len(marked_positions)
        for combo in range(1, 1 << n):
            # Build the new type tuple for this combo
            new_types = list(type_tuple)
            for bit_idx, pos in enumerate(marked_positions):
                if combo & (1 << bit_idx):
                    new_types[pos] = type_tuple[pos] + '_'

            # Canonical key: lexicographically smaller of forward/reverse
            canonical = _canonical_bonded_key(tuple(new_types))
            dup_key = (canonical, funct)

            if dup_key in dup_groups:
                continue  # Reverse-equivalent combo already generated

            # Build fields for all entries (multi-term funct=9) of this combo
            combo_entries = []
            for fields, comment in entries:
                new_fields = list(fields)
                for bit_idx, pos in enumerate(marked_positions):
                    if combo & (1 << bit_idx):
                        new_fields[pos] = fields[pos] + '_'
                combo_entries.append((new_fields, comment))
            dup_groups[dup_key] = combo_entries

    # Emit all duplicate entries — grouped by canonical key for contiguity
    dup_lines: List[str] = []
    for (canonical, funct), combo_entries in dup_groups.items():
        for fields, comment in combo_entries:
            dup_line = '  '.join(
                f'{f:>12}' if _looks_numeric(f) else f'  {f}'
                for f in fields
            )
            if comment:
                dup_line += f'  {comment}'
            dup_lines.append(dup_line + '\n')

    if dup_lines:
        result.append(f'; REST2 solute {section_name} duplicates\n')
        result.extend(dup_lines)

    return result


def _mark_atom_line(
    line: str,
    current_moltype: str,
    lig_moltypes: Set[str],
    pocket_residues: Dict[str, Set[int]],
) -> str:
    """Mark an atom line by suffixing its type with '_' if it should be scaled."""
    stripped = line.strip()
    data_part = stripped.split(';')[0]
    comment_part = ''
    if ';' in stripped:
        comment_part = stripped[stripped.index(';'):]

    fields = data_part.split()
    if len(fields) < 7:
        return line

    atom_type = fields[1]
    resnr = int(fields[2])

    should_mark = False
    if current_moltype in lig_moltypes:
        should_mark = True
    elif current_moltype in pocket_residues:
        if resnr in pocket_residues[current_moltype]:
            should_mark = True

    if should_mark and not atom_type.endswith('_'):
        fields[1] = atom_type + '_'

    # Reconstruct with proper spacing
    # nr(6) type(11) resnr(7) residue(7) atom(7) cgnr(7) charge(11) mass(11) ...
    result = f'{fields[0]:>6} {fields[1]:>11}{fields[2]:>7} {fields[3]:<7} {fields[4]:>6}{fields[5]:>7}{fields[6]:>11}'
    if len(fields) > 7:
        result += f'{fields[7]:>11}'
    for i in range(8, len(fields)):
        result += f'{fields[i]:>11}'

    if comment_part:
        result += f'   {comment_part}'
    return result + '\n'


def _looks_numeric(s: str) -> bool:
    """Check if a string looks like a number."""
    try:
        float(s)
        return True
    except ValueError:
        return False
