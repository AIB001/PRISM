"""
Pure Python replacement for `plumed partial_tempering`.

Scales force-field parameters for atoms whose type ends with '_' (solute atoms),
leaving solvent/environment atoms untouched. This allows REST2 enhanced sampling
by baking the scaling directly into GROMACS topologies.

Usage:
    from prism.rest2.partial_tempering import partial_tempering
    scaled_lines = partial_tempering(topology_lines, lam=0.7)
"""

import math
import re
from typing import Dict, List, Optional, Tuple

# Regex for section headers like [ atoms ], [ bonds ], etc.
_SECTION_RE = re.compile(r'^\s*\[\s*(\w+)\s*\]')


def _strip_comment(line: str) -> Tuple[str, str]:
    """Split a line into data part and comment part (including the ';')."""
    idx = line.find(';')
    if idx == -1:
        return line, ''
    return line[:idx], line[idx:]


def _is_solute_type(atom_type: str) -> bool:
    """Check if an atom type is marked as solute (ends with '_')."""
    return atom_type.endswith('_')


def _find_a_column(header_comment: str, fields: List[str]) -> Tuple[int, int]:
    """Find sigma and epsilon column indices from atomtypes section.

    In GROMACS [ atomtypes ], columns vary. We look for the 'A' ptype column
    to anchor: sigma is the column after 'A', epsilon is the one after that.

    For the processed topology, the format is typically:
    type  bondingtype  at_num  mass  charge  ptype  sigma  epsilon
    or:
    type  at_num  mass  charge  ptype  sigma  epsilon

    We find the field that equals 'A' (ptype), then sigma = next, epsilon = next+1.
    """
    for i, f in enumerate(fields):
        if f == 'A':
            return i + 1, i + 2
    # Fallback: assume last two columns are sigma and epsilon
    return len(fields) - 2, len(fields) - 1


def partial_tempering(lines: List[str], lam: float) -> List[str]:
    """Apply REST2 partial tempering scaling to a GROMACS topology.

    For atoms with types ending in '_' (solute):
    - [ atomtypes ]: epsilon *= lambda
    - [ atoms ]: charge *= sqrt(lambda)
    - [ bonds ] funct=1: force constant scaled (all-solute: *lambda, mixed: *sqrt(lambda))
    - [ angles ] funct=1: force constant scaled similarly
    - [ dihedrals ] funct=9/4/1: kphi scaled similarly
    - [ dihedrals ] funct=3 (RB): C0-C5 scaled similarly
    - [ pairs ] 3-column: no change (derived from scaled atomtypes)

    Args:
        lines: Lines of a processed GROMACS topology (from gmx grompp -pp)
        lam: Scaling factor lambda (T_ref / T_effective), between 0 and 1

    Returns:
        List of scaled topology lines
    """
    sqrt_lam = math.sqrt(lam)

    current_section: Optional[str] = None
    atom_is_solute: Dict[int, bool] = {}  # atom_nr -> is_solute for current moleculetype

    output: List[str] = []

    for line in lines:
        stripped = line.strip()

        # Check for section header
        m = _SECTION_RE.match(stripped)
        if m:
            current_section = m.group(1).lower()
            # Reset atom tracking on new moleculetype
            if current_section == 'moleculetype':
                atom_is_solute = {}
            output.append(line)
            continue

        # Skip blank lines and comments
        if not stripped or stripped.startswith(';'):
            output.append(line)
            continue

        # Handle #ifdef, #endif, #include, #define directives
        if stripped.startswith('#'):
            output.append(line)
            continue

        data_part, comment_part = _strip_comment(stripped)
        fields = data_part.split()

        if not fields:
            output.append(line)
            continue

        # --- Process each section ---

        if current_section == 'atomtypes':
            output.append(_scale_atomtypes(fields, comment_part, lam))

        elif current_section == 'atoms':
            output.append(_scale_atoms(fields, comment_part, sqrt_lam, atom_is_solute))

        elif current_section == 'bonds':
            output.append(_scale_bonds(fields, comment_part, lam, sqrt_lam, atom_is_solute))

        elif current_section == 'angles':
            output.append(_scale_angles(fields, comment_part, lam, sqrt_lam, atom_is_solute))

        elif current_section == 'dihedrals':
            output.append(_scale_dihedrals(fields, comment_part, lam, sqrt_lam, atom_is_solute))

        elif current_section == 'pairs':
            # 3-column pairs (ai, aj, funct): no scaling needed
            # Pairs with explicit parameters would need scaling, but standard
            # AMBER topologies use 3-column pairs derived from atomtypes
            output.append(line)

        else:
            output.append(line)

    return output


def _scale_atomtypes(fields: List[str], comment: str, lam: float) -> str:
    """Scale epsilon in [ atomtypes ] for solute types (ending with '_')."""
    atom_type = fields[0]
    if not _is_solute_type(atom_type):
        return _reconstruct_line(fields, comment)

    # Find sigma/epsilon positions by looking for 'A' column
    sigma_idx, eps_idx = _find_a_column('', fields)

    if eps_idx < len(fields):
        try:
            eps = float(fields[eps_idx])
            fields[eps_idx] = _format_sci(eps * lam)
        except (ValueError, IndexError):
            pass

    return _reconstruct_line(fields, comment)


def _scale_atoms(fields: List[str], comment: str, sqrt_lam: float,
                 atom_is_solute: Dict[int, bool]) -> str:
    """Scale charge in [ atoms ] for solute atoms and record solute status.

    Processed topology format:
    nr  type  resnr  residue  atom  cgnr  charge  mass  [typeB  chargeB  massB]
    """
    try:
        atom_nr = int(fields[0])
        atom_type = fields[1]
        is_solute = _is_solute_type(atom_type)
        atom_is_solute[atom_nr] = is_solute

        if is_solute and len(fields) > 6:
            charge = float(fields[6])
            fields[6] = _format_decimal(charge * sqrt_lam)
            # If there's a chargeB (perturbation), scale that too
            if len(fields) > 9:
                try:
                    chargeB = float(fields[9])
                    fields[9] = _format_decimal(chargeB * sqrt_lam)
                except (ValueError, IndexError):
                    pass
    except (ValueError, IndexError):
        pass

    return _reconstruct_line(fields, comment)


def _classify_interaction(atom_indices: List[int], atom_is_solute: Dict[int, bool]) -> str:
    """Classify a bonded interaction as 'all_solute', 'mixed', or 'all_solvent'.

    Returns:
        'all_solute' if all atoms are solute
        'mixed' if some atoms are solute and some are not
        'all_solvent' if no atoms are solute
    """
    solute_count = sum(1 for idx in atom_indices if atom_is_solute.get(idx, False))
    if solute_count == len(atom_indices):
        return 'all_solute'
    elif solute_count > 0:
        return 'mixed'
    else:
        return 'all_solvent'


def _get_scale_factor(classification: str, lam: float, sqrt_lam: float) -> float:
    """Get the scaling factor based on interaction classification."""
    if classification == 'all_solute':
        return lam
    elif classification == 'mixed':
        return sqrt_lam
    else:
        return 1.0


def _scale_bonds(fields: List[str], comment: str, lam: float, sqrt_lam: float,
                 atom_is_solute: Dict[int, bool]) -> str:
    """Scale force constant in [ bonds ] section.

    Format: ai  aj  funct  [c0  c1  ...]
    For funct=1: c0=r0, c1=kb -> scale kb
    """
    if len(fields) < 3:
        return _reconstruct_line(fields, comment)

    try:
        ai, aj = int(fields[0]), int(fields[1])
        funct = int(fields[2])
    except ValueError:
        return _reconstruct_line(fields, comment)

    classification = _classify_interaction([ai, aj], atom_is_solute)
    if classification == 'all_solvent':
        return _reconstruct_line(fields, comment)

    scale = _get_scale_factor(classification, lam, sqrt_lam)

    if funct == 1 and len(fields) > 4:
        # c0=r0 (col 3), c1=kb (col 4) -> scale kb
        try:
            kb = float(fields[4])
            fields[4] = _format_sci(kb * scale)
        except ValueError:
            pass

    return _reconstruct_line(fields, comment)


def _scale_angles(fields: List[str], comment: str, lam: float, sqrt_lam: float,
                  atom_is_solute: Dict[int, bool]) -> str:
    """Scale force constant in [ angles ] section.

    Format: ai  aj  ak  funct  [c0  c1  ...]
    For funct=1: c0=theta0, c1=ktheta -> scale ktheta (col 5)
    For funct=5 (Urey-Bradley, CHARMM): theta0 ktheta rUB kUB -> scale ktheta (col 5) and kUB (col 7)
    """
    if len(fields) < 4:
        return _reconstruct_line(fields, comment)

    try:
        ai, aj, ak = int(fields[0]), int(fields[1]), int(fields[2])
        funct = int(fields[3])
    except ValueError:
        return _reconstruct_line(fields, comment)

    classification = _classify_interaction([ai, aj, ak], atom_is_solute)
    if classification == 'all_solvent':
        return _reconstruct_line(fields, comment)

    scale = _get_scale_factor(classification, lam, sqrt_lam)

    if funct == 1 and len(fields) > 5:
        # c0=theta0 (col 4), c1=ktheta (col 5) -> scale ktheta
        try:
            ktheta = float(fields[5])
            fields[5] = _format_sci(ktheta * scale)
        except ValueError:
            pass
    elif funct == 5 and len(fields) > 7:
        # Urey-Bradley (CHARMM): theta0(col4) ktheta(col5) rUB(col6) kUB(col7)
        try:
            ktheta = float(fields[5])
            fields[5] = _format_sci(ktheta * scale)
        except ValueError:
            pass
        try:
            kUB = float(fields[7])
            fields[7] = _format_sci(kUB * scale)
        except ValueError:
            pass

    return _reconstruct_line(fields, comment)


def _scale_dihedrals(fields: List[str], comment: str, lam: float, sqrt_lam: float,
                     atom_is_solute: Dict[int, bool]) -> str:
    """Scale force constant in [ dihedrals ] section.

    Format: ai  aj  ak  al  funct  [params...]
    For funct=9/4/1 (proper/improper): phase(col5)  kphi(col6)  pn(col7) -> scale kphi
    For funct=2 (harmonic improper, CHARMM): xi0(col5)  kxi(col6) -> scale kxi
    For funct=3 (RB): C0-C5 (cols 5-10) -> scale all
    """
    if len(fields) < 5:
        return _reconstruct_line(fields, comment)

    try:
        ai, aj, ak, al = int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])
        funct = int(fields[4])
    except ValueError:
        return _reconstruct_line(fields, comment)

    classification = _classify_interaction([ai, aj, ak, al], atom_is_solute)
    if classification == 'all_solvent':
        return _reconstruct_line(fields, comment)

    scale = _get_scale_factor(classification, lam, sqrt_lam)

    if funct in (9, 4, 1):
        # phase (col 5), kphi (col 6), pn (col 7) -> scale kphi
        if len(fields) > 6:
            try:
                kphi = float(fields[6])
                fields[6] = _format_decimal(kphi * scale)
            except ValueError:
                pass
    elif funct == 2:
        # Harmonic improper (CHARMM): xi0(col5) kxi(col6) -> scale kxi
        if len(fields) > 6:
            try:
                kxi = float(fields[6])
                fields[6] = _format_sci(kxi * scale)
            except ValueError:
                pass
    elif funct == 3:
        # Ryckaert-Bellemans: C0-C5 (cols 5-10)
        for i in range(5, min(11, len(fields))):
            try:
                c = float(fields[i])
                fields[i] = _format_sci(c * scale)
            except ValueError:
                pass

    return _reconstruct_line(fields, comment)


def _format_sci(value: float) -> str:
    """Format a float in scientific notation matching GROMACS style."""
    if value == 0.0:
        return '0.00000e+00'
    return f'{value:.5e}'


def _format_decimal(value: float, decimals: int = 5) -> str:
    """Format a float in decimal notation."""
    # Use enough precision to avoid rounding issues
    formatted = f'{value:.6f}'
    # Strip trailing zeros but keep at least one decimal
    if '.' in formatted:
        formatted = formatted.rstrip('0').rstrip('.')
        # Ensure at least some decimal places for readability
        if '.' not in formatted:
            formatted += '.0'
    return formatted


def _reconstruct_line(fields: List[str], comment: str) -> str:
    """Reconstruct a topology line from fields and optional comment."""
    data = '  '.join(f'{f:>12}' if _looks_numeric(f) else f'  {f}' for f in fields)
    if comment:
        return f'{data}  {comment}\n'
    return f'{data}\n'


def _looks_numeric(s: str) -> bool:
    """Check if a string looks like a number."""
    try:
        float(s)
        return True
    except ValueError:
        return s.startswith('-') or s.startswith('+')
