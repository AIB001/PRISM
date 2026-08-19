#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Driver-script generation for membrane (protein-in-bilayer) systems.

The soluble protein-ligand path emits ``GMX_PROLIG_MD/localrun.sh`` from
:meth:`prism.builder.preparation.PRISMBuilder.generate_localrun_script`. That
generator cannot be reused for membranes, for three reasons:

1. it hardcodes the system directory name (``GMX_PROLIG_MD``);
2. it hardcodes the input coordinate name (``solv_ions.gro``);
3. **its ``gmx grompp`` calls pass no ``-n index.ndx``.**

Point 3 is fatal here rather than merely inelegant. The membrane MDPs written by
:func:`prism.membrane.membrane_mdp.write_membrane_mdps` thermostat on
``tc-grps = SOLU MEMB SOLV``. Those are PRISM-defined index groups (written by
the membrane builder), *not* GROMACS built-ins like ``Protein`` or ``Water``, so
``grompp`` aborts with "Group SOLU referenced in the .mdp file was not found in
the index file" unless every single invocation is handed ``-n index.ndx``.

This module therefore emits its own ``localrun.sh`` while keeping the two
properties of the soluble script that matter operationally:

* **runtime hardware detection** -- thread count and GPU offload are resolved on
  the machine that runs the script, never baked in at build time (a hard project
  rule: a script built on a laptop must still be correct on a GPU node);
* **resume/skip logic** -- a stage whose final ``.gro`` exists is skipped, and a
  stage with a ``.tpr`` but no ``.gro`` is restarted from its checkpoint. An
  interrupted overnight membrane equilibration is expensive to repeat. Resuming
  is only safe while the ``.tpr`` still matches its inputs, so a ``.tpr`` older
  than its .mdp/topology/index/coordinates is rebuilt instead of reused, and
  ``PRISM_REDO="<stage>"`` discards a finished stage (and everything after it).

and adds what the membrane protocol needs on top:

* ``-n <index.ndx>`` on *every* ``grompp`` call;
* ``-r <reference.gro>`` on every ``grompp`` call (GROMACS >= 2018 refuses to
  build a .tpr when position restraints are active and no reference coordinate
  file was given -- and the staged protocol has ``-DPOSRES`` live on every
  equilibration stage);
* ``-t <prev>.cpt`` wherever a stage continues the previous stage's velocities
  rather than generating fresh ones.

Every path the script references is relative, and the .mdp files live in the
sibling ``../mdps/`` directory (exactly as in the soluble script), so the
*parent* directory -- system dir plus ``mdps/`` -- is the relocatable unit.

Two layout details are contracts with the rest of PRISM rather than local
choices:

* production writes ``./prod/md.*``. The stage is called ``md`` in the MDP set,
  but ``prism.analysis``, ``prism.sim`` and ``prism.mmpbsa`` all look for
  ``<system dir>/prod/md.xtc``; a membrane trajectory anywhere else is invisible
  to every downstream tool.
* the script is written twice, as ``localrun.sh`` and as ``run_membrane.sh``.
  :class:`prism.sim.gmxsim.GMXSimulator` overwrites ``<system dir>/localrun.sh``
  unconditionally with the soluble script -- which has no ``-n index.ndx`` and so
  cannot grompp a membrane system at all. The second copy, under a name nothing
  else in PRISM writes, turns that clobber into a ``cp`` instead of a rebuild.
  Line 2 of both files carries :data:`RUNSCRIPT_SENTINEL` so a generator can
  recognise a membrane script it did not author.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

__all__ = [
    "INDEX_NAME",
    "RUNSCRIPT_BACKUP_NAME",
    "RUNSCRIPT_NAME",
    "RUNSCRIPT_SENTINEL",
    "SYSTEM_GRO_NAME",
    "TOPOL_NAME",
    "order_membrane_stages",
    "write_membrane_runscript",
]

#: Name of the generated driver script (matches the soluble path's convention).
RUNSCRIPT_NAME = "localrun.sh"

#: Pristine second copy, under a name no other PRISM generator writes. See the
#: module docstring: prism.sim overwrites localrun.sh without asking.
RUNSCRIPT_BACKUP_NAME = "run_membrane.sh"

#: Machine-readable marker on line 2 of both copies, so any other generator can
#: tell "a membrane driver I must not clobber" from "a script I wrote".
RUNSCRIPT_SENTINEL = "# PRISM-RUNSCRIPT: membrane v1"

#: Build outputs the script consumes. Named here rather than spelled as literals
#: at each call site so that prism.sim (which otherwise assumes the soluble
#: ``solv_ions.gro``) can import the membrane contract instead of guessing it.
SYSTEM_GRO_NAME = "system.gro"
TOPOL_NAME = "topol.top"
INDEX_NAME = "index.ndx"

#: Production output location. Not derived from the stage name: every consumer in
#: PRISM (prism/sim/gmxsim.py, prism/mmpbsa/generator.py,
#: prism/analysis/core/analyzer.py) hardcodes ``prod/md.*``.
_PRODUCTION_DIR = "prod"
_PRODUCTION_STEM = "md"


# ---------------------------------------------------------------------------
# Stage ordering
# ---------------------------------------------------------------------------
#
# The caller hands us a mapping of stage-name -> mdp path. We must not hardcode
# the membrane protocol's current stage list (em/nvt/npt/npt2/md): that list is
# being reworked and may gain or lose a stage. What *is* invariant is the
# physics of the pipeline, so we rank stages by the phase their name denotes:
#
#     minimise  ->  heat  ->  fixed volume  ->  fixed pressure  ->  produce
#
# A phase rank is stable across protocol revisions in a way that an explicit
# stage list is not: inserting an "npt3" or a "heat" stage requires no change
# here, and the anchors the rest of the script depends on (a minimisation
# first, because it is the only stage that can consume the freshly built
# coordinates; production last, because it is the unrestrained run everything
# else exists to set up) are enforced structurally rather than by convention.
#
# Gaps between the ranks are deliberate -- new phases can be slotted in without
# renumbering. Unclassifiable names land in _PHASE_UNKNOWN, which sits after all
# recognised equilibration but before production: an unrecognised stage is
# overwhelmingly likely to be extra equilibration, and placing it there cannot
# displace either anchor.

_PHASE_MINIMIZATION = 0
_PHASE_HEAT = 10
_PHASE_NVT = 20
_PHASE_NPT = 30
_PHASE_EQUILIBRATION = 40
_PHASE_UNKNOWN = 50
_PHASE_PRODUCTION = 90

#: Alphabetic stage-name token -> phase rank. Matched exactly first, then by
#: longest prefix, so ``npt2``/``npt_semiiso`` both resolve to the NPT phase.
_PHASE_BY_TOKEN: Dict[str, int] = {
    "em": _PHASE_MINIMIZATION,
    "emin": _PHASE_MINIMIZATION,
    "min": _PHASE_MINIMIZATION,
    "minim": _PHASE_MINIMIZATION,
    "minimize": _PHASE_MINIMIZATION,
    "minimise": _PHASE_MINIMIZATION,
    "minimization": _PHASE_MINIMIZATION,
    "minimisation": _PHASE_MINIMIZATION,
    "steep": _PHASE_MINIMIZATION,
    "cg": _PHASE_MINIMIZATION,
    "heat": _PHASE_HEAT,
    "anneal": _PHASE_HEAT,
    "annealing": _PHASE_HEAT,
    "nvt": _PHASE_NVT,
    "npt": _PHASE_NPT,
    "eq": _PHASE_EQUILIBRATION,
    "equil": _PHASE_EQUILIBRATION,
    "equilibration": _PHASE_EQUILIBRATION,
    "relax": _PHASE_EQUILIBRATION,
    "md": _PHASE_PRODUCTION,
    "prod": _PHASE_PRODUCTION,
    "production": _PHASE_PRODUCTION,
}

#: Stage names become directory names, file stems and shell words, so keep them
#: boring. Anything else is rejected loudly rather than quoted defensively.
_SAFE_STAGE_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")

_LEADING_ALPHA = re.compile(r"^[a-z]+")
_TRAILING_DIGITS = re.compile(r"(\d+)$")


def _phase_rank(stage: str) -> int:
    """Map a stage name onto its physical phase rank (see the table above)."""
    lowered = stage.lower()
    match = _LEADING_ALPHA.match(lowered)
    token = match.group(0) if match else ""
    if token in _PHASE_BY_TOKEN:
        return _PHASE_BY_TOKEN[token]
    # Longest prefix wins so that "minimization" cannot be captured by "min"
    # ahead of a more specific key, and "nptsemiiso" still resolves to NPT.
    for known in sorted(_PHASE_BY_TOKEN, key=len, reverse=True):
        if token.startswith(known):
            return _PHASE_BY_TOKEN[known]
    return _PHASE_UNKNOWN


def _repeat_index(stage: str) -> int:
    """Trailing integer of a stage name (``npt2`` -> 2, ``npt`` -> 0).

    This is what keeps repeated phases in protocol order without naming them:
    ``npt`` < ``npt2`` < ``npt10`` regardless of how the caller's mapping is
    iterated (a lexical sort would put ``npt10`` second).
    """
    match = _TRAILING_DIGITS.search(stage)
    return int(match.group(1)) if match else 0


def order_membrane_stages(stages: Sequence[str]) -> List[str]:
    """Return ``stages`` in execution order.

    Sort key, most significant first:

    1. **phase rank** -- minimisation first, production last, everything else in
       physical order (see ``_PHASE_BY_TOKEN``);
    2. **repeat index** -- the trailing integer, so ``npt`` precedes ``npt2``;
    3. **declaration order** -- the caller's mapping order, which breaks ties
       among stages we cannot classify and makes the result deterministic for
       any input.

    Total and deterministic: key 3 is unique per stage, so no two stages ever
    compare equal and the output never depends on dict iteration luck.
    """
    return [
        name
        for _, _, _, name in sorted(
            (_phase_rank(name), _repeat_index(name), index, name)
            for index, name in enumerate(stages)
        )
    ]


# ---------------------------------------------------------------------------
# MDP introspection
# ---------------------------------------------------------------------------


def _read_mdp_options(path: str) -> Dict[str, str]:
    """Parse ``key = value`` pairs out of an .mdp file.

    Keys are normalised to lowercase with ``-`` separators, because GROMACS
    accepts ``gen_vel`` and ``gen-vel`` interchangeably. Unreadable files yield
    an empty mapping: MDP introspection only *refines* the defaults below, and a
    missing file is reported by the script's own preflight check, not here.
    """
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    except OSError:
        return {}

    options: Dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        options[key.strip().lower().replace("_", "-")] = value.strip()
    return options


def _is_yes(value: Optional[str]) -> bool:
    return (value or "").strip().lower() in {"yes", "true", "on"}


_DEFINE_MACRO = re.compile(r"-D([A-Za-z_][A-Za-z0-9_]*)(=\S*)?")


def _restraint_guards(define: str) -> Tuple[str, ...]:
    """Bare ``-DPOSRES*`` macros in an MDP ``define`` line, in order.

    A macro *with* a value (``-DPOSRES_FC_BB=500.0``) is a force constant
    substituted inside a restraint block; a bare one (``-DPOSRES``,
    ``-DPOSRES_LIPID``) is the ``#ifdef`` guard that switches the block on. Only
    the guards matter downstream: they are what the topology has to consume for
    the restraints to exist at all, and grompp will not say so if it does not.
    """
    guards: List[str] = []
    for name, value in _DEFINE_MACRO.findall(define):
        if not value and name.upper().startswith("POSRES") and name not in guards:
            guards.append(name)
    return tuple(guards)


# ---------------------------------------------------------------------------
# Stage plan
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _Stage:
    """One resolved stage of the run, with every path already relative."""

    name: str
    mdp: str  # relative path to the .mdp
    coord_in: str  # relative -c / -r input (previous stage's output, or the build)
    cpt_in: Optional[str]  # relative -t input, or None
    is_minimization: bool
    is_production: bool
    restrained: bool
    note: str  # short human-readable reason, written into the script
    guards: Tuple[str, ...] = ()  # -DPOSRES* guards the topology must consume

    @property
    def directory(self) -> str:
        # Production is the one stage whose location is not its own business:
        # prism.analysis/sim/mmpbsa discover a run by looking for prod/md.xtc.
        return f"./{_PRODUCTION_DIR}" if self.is_production else f"./{self.name}"

    @property
    def stem(self) -> str:
        return _PRODUCTION_STEM if self.is_production else self.name

    @property
    def prefix(self) -> str:
        return f"{self.directory}/{self.stem}"

    @property
    def tpr(self) -> str:
        return f"{self.prefix}.tpr"

    @property
    def gro(self) -> str:
        return f"{self.prefix}.gro"

    @property
    def cpt(self) -> str:
        return f"{self.prefix}.cpt"


def _relative_to(path: str, start: str) -> str:
    """Path relative to the system directory, for a relocatable script.

    Falls back to the absolute path only when no relative path exists (different
    Windows drives) -- a script that is merely non-portable beats one that is
    wrong.
    """
    try:
        relative = os.path.relpath(os.path.abspath(path), start)
    except ValueError:
        return os.path.abspath(path)
    return relative.replace(os.sep, "/")


def _build_plan(
    ordered: Sequence[str],
    mdp_paths: Mapping[str, str],
    system_dir: str,
    gro: str,
) -> List[_Stage]:
    """Resolve every stage's inputs by chaining the previous stage's outputs."""
    plan: List[_Stage] = []
    previous: Optional[_Stage] = None

    for position, name in enumerate(ordered):
        mdp = _relative_to(os.fspath(mdp_paths[name]), system_dir)
        options = _read_mdp_options(os.path.join(system_dir, mdp))

        rank = _phase_rank(name)
        is_minimization = rank == _PHASE_MINIMIZATION
        # Only the *last* production-ranked stage is redirected to ./prod/md.*:
        # two stages sharing that directory would overwrite each other, and the
        # trajectory downstream tools want is the final one anyway.
        is_production = rank == _PHASE_PRODUCTION and position == len(ordered) - 1

        # Whether the stage makes its own velocities. Prefer the MDP's own
        # answer; fall back to the phase when the file cannot be read, where the
        # convention is that the first dynamics stage after minimisation is the
        # one that heats the system.
        if "gen-vel" in options:
            generates_velocities = _is_yes(options["gen-vel"])
        else:
            generates_velocities = rank in (_PHASE_HEAT, _PHASE_NVT) or (
                previous is not None and previous.is_minimization
            )

        # -DPOSRES* in the MDP's define line means restraints are live. Used to
        # annotate the script and to check the topology actually consumes the
        # guards; -r is passed unconditionally (see below).
        define = options.get("define", "")
        restrained = "POSRES" in define.upper()
        guards = _restraint_guards(define)

        # Chain coordinates: the previous stage's final frame is both the
        # starting structure and the position-restraint reference, so restraints
        # hold atoms where the last stage left them instead of dragging them back
        # to the raw packed geometry.
        coord_in = previous.gro if previous is not None else gro

        # Carry the checkpoint forward only where it means something physically:
        # a checkpoint restores velocities plus thermostat/barostat state, which
        # is exactly what a continuation stage wants, and exactly what a stage
        # that regenerates velocities (or minimises, which has none) does not.
        cpt_in: Optional[str] = None
        if previous is not None and not previous.is_minimization and not is_minimization:
            if not generates_velocities:
                cpt_in = previous.cpt

        if is_minimization:
            note = "relieve packing clashes in the built system"
        elif generates_velocities:
            note = "generates velocities (gen-vel = yes); no checkpoint to inherit"
        elif position == len(ordered) - 1:
            note = "production run; continues the equilibrated state"
        else:
            note = "continues the previous stage's velocities and coupling state"

        stage = _Stage(
            name=name,
            mdp=mdp,
            coord_in=coord_in,
            cpt_in=cpt_in,
            is_minimization=is_minimization,
            is_production=is_production,
            restrained=restrained,
            note=note,
            guards=guards,
        )
        plan.append(stage)
        previous = stage

    return plan


# ---------------------------------------------------------------------------
# Script rendering
# ---------------------------------------------------------------------------

# Hardware detection and failure handling. Deliberately identical in spirit to
# the soluble localrun.sh: nothing about the build host is baked in.
_PREAMBLE = """\
# Abort on the first failing command, on unset variables, and on a failure
# anywhere in a pipeline. A membrane equilibration is an overnight job; limping
# on after a failed grompp would waste the whole night and produce a trajectory
# that looks finished but is not.
set -Eeuo pipefail

CURRENT_STAGE="startup"
trap 'echo "" >&2; echo "ERROR: stage \\"$CURRENT_STAGE\\" failed (line $LINENO). Nothing after this point was run." >&2; echo "       If grompp stopped on a WARNING, read it before overriding: on a membrane" >&2; echo "       system the counted warnings are real (net charge under Ewald, atoms far" >&2; echo "       outside the box). PRISM_GROMPP_MAXWARN=<n> raises the tolerance." >&2' ERR

# GROMACS binary: override for module systems or MPI builds, e.g. GMX=gmx_mpi.
GMX="${GMX:-gmx}"
# grompp warning tolerance. The soluble PRISM script uses 999; that is far too
# permissive here, because grompp is the only automated check a membrane system
# ever gets. A correctly built one passes at 0: the messages a semiisotropic
# bilayer legitimately produces (position restraints, COM-motion removal) are
# NOTEs, which never count against maxwarn. What does count -- a net-charged
# system under Ewald, atoms outside the box after packing -- is a defect worth
# stopping for. So nothing is tolerated by default; raise it deliberately with
# PRISM_GROMPP_MAXWARN=N once you have read what you are suppressing.
MAXWARN="${PRISM_GROMPP_MAXWARN:-0}"

# --- Hardware detection (resolved here, on this host; never hardcoded) -------
# Threads: respect OMP_NUM_THREADS if the caller (or the queue system) set it,
# otherwise use every logical core we can see.
NTOMP="${OMP_NUM_THREADS:-$(nproc 2>/dev/null || echo 1)}"
# Thread-MPI ranks. The flag is built as a whole token, not just a value: a real
# MPI build (GMX=gmx_mpi) rejects -ntmpi outright, so for those the flag has to
# vanish rather than take a number. PRISM_NTMPI=1 (the default) keeps GPU PME on
# a single rank, PRISM_NTMPI=0 lets GROMACS decide, and PRISM_NTMPI= -- set but
# empty -- drops the flag entirely and leaves rank layout to the MPI launcher.
if [ "${PRISM_NTMPI+set}" = "set" ]; then
    NTMPI_ARG="${PRISM_NTMPI:+-ntmpi $PRISM_NTMPI}"
else
    NTMPI_ARG="-ntmpi 1"
fi
# GPU: offload only if an NVIDIA device is actually present and queryable,
# otherwise run CPU-only. On a multi-GPU host, export CUDA_VISIBLE_DEVICES=<id>.
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    GPU_ARGS="-nb gpu -bonded gpu -pme gpu -gpu_id 0"
    # Minimisation takes nonbonded offload only: the bonded and PME GPU kernels
    # are tuned for (and validated against) the MD integrators, and steepest
    # descent is never the bottleneck anyway.
    GPU_ARGS_MIN="-nb gpu -gpu_id 0"
else
    GPU_ARGS=""
    GPU_ARGS_MIN=""
fi
echo "GROMACS driver: $GMX | mdrun: ${NTMPI_ARG:+$NTMPI_ARG }-ntomp $NTOMP | ${GPU_ARGS:+GPU offload}${GPU_ARGS:-CPU-only}"

# --- Stage helpers ----------------------------------------------------------
# A finished stage is skipped, which is what makes an interrupted overnight run
# cheap to restart -- but it also means an edited .mdp never reaches a stage that
# already has its .gro. PRISM_REDO names the stages to throw away and rerun,
# e.g. PRISM_REDO="npt2". Names are matched against the protocol (checked in the
# preflight) and the directory removed is generated below, so nothing the user
# types is ever passed to rm.
#
# The redo cascades: every stage after the named one is discarded too. Each
# stage's inputs are the previous stage's outputs, so keeping a later stage that
# was built from a .gro that no longer exists would silently splice two
# different runs together.
PRISM_REDO_ACTIVE=0
prism_redo() {
    if [ "$PRISM_REDO_ACTIVE" = "1" ]; then
        return 0
    fi
    for requested_stage in ${PRISM_REDO:-}; do
        if [ "$requested_stage" = "$1" ]; then
            PRISM_REDO_ACTIVE=1
            return 0
        fi
    done
    return 1
}

# Discard a .tpr that is older than an input it was built from. Resuming from a
# stale .tpr silently reuses the settings it was built with, and editing an .mdp
# (to lower dt after a LINCS failure, say) is exactly what this script tells the
# user to do -- so the edit has to invalidate the .tpr rather than be ignored.
prism_drop_stale() {
    stage_tpr="$1"
    shift
    [ -f "$stage_tpr" ] || return 0
    for stage_input in "$@"; do
        if [ -e "$stage_input" ] && [ "$stage_input" -nt "$stage_tpr" ]; then
            echo "    NOTE: $stage_input is newer than $stage_tpr - rebuilding it so the current inputs are used."
            rm -f "$stage_tpr" "${stage_tpr%.tpr}.cpt" "${stage_tpr%.tpr}_prev.cpt"
            return 0
        fi
    done
    return 0
}
"""


def _preflight_block(
    required: Sequence[str],
    stage_names: Sequence[str],
    guards: Sequence[str],
    top: str,
) -> str:
    """Check every input exists before burning queue time on stage 1."""
    lines = [
        "# --- Preflight ---------------------------------------------------------",
        "# Fail now, at second zero, rather than after the minimisation has run.",
        'command -v "$GMX" >/dev/null 2>&1 || {',
        '    echo "ERROR: GROMACS binary \'$GMX\' not found in PATH."'
        ' >&2',
        '    echo "       Load your GROMACS module or source GMXRC first." >&2',
        "    exit 1",
        "}",
        "for required_input in \\",
    ]
    lines.extend(f'    "{item}" \\' for item in required)
    lines.extend(
        [
            "    ; do",
            '    [ -f "$required_input" ] || {',
            '        echo "ERROR: required input \'$required_input\' is missing." >&2',
            '        echo "       Run this script from the directory that contains it:" >&2',
            '        echo "         cd <system dir> && bash ' + RUNSCRIPT_NAME + '" >&2',
            "        exit 1",
            "    }",
            "done",
        ]
    )

    if guards:
        # grompp raises nothing at all for a -DPOSRES that no topology consumes:
        # it builds a .tpr with no restraints, and the staged protocol then runs
        # the protein completely free while reporting success. The topology's
        # #ifdef blocks are the only place a failed restraint generation shows.
        lines.extend(
            [
                "for restraint_guard in " + " ".join(guards) + "; do",
                f'    grep -qE "^[[:space:]]*#ifdef[[:space:]]+$restraint_guard[[:space:]]*$" {top} || {{',
                '        echo "ERROR: the MDPs define -D$restraint_guard, but ' + top + ' has no matching" >&2',
                '        echo "       \'#ifdef $restraint_guard\' include. grompp does not reject a define that" >&2',
                '        echo "       nothing consumes, so the equilibration would run entirely unrestrained" >&2',
                '        echo "       and the fold would be gone by production." >&2',
                '        echo "       Rebuild the system: position-restraint generation did not complete." >&2',
                "        exit 1",
                "    }",
                "done",
            ]
        )

    # Reject a mistyped PRISM_REDO here: doing nothing silently looks exactly
    # like the skip the variable exists to override.
    lines.extend(
        [
            "for requested_stage in ${PRISM_REDO:-}; do",
            '    case " ' + " ".join(stage_names) + ' " in',
            '        *" $requested_stage "*) ;;',
            "        *)",
            '            echo "ERROR: PRISM_REDO names unknown stage \'$requested_stage\'." >&2',
            '            echo "       Protocol stages: ' + " ".join(stage_names) + '" >&2',
            "            exit 1",
            "            ;;",
            "    esac",
            "done",
            'echo "Preflight OK: inputs, MDPs and restraint guards present (grompp -maxwarn $MAXWARN)."',
        ]
    )
    return "\n".join(lines) + "\n"


def _stage_block(stage: _Stage, top: str, ndx: str) -> str:
    """Render one stage: skip if done, resume if interrupted, else run it."""
    gpu = "$GPU_ARGS_MIN" if stage.is_minimization else "$GPU_ARGS"
    mdrun_args = f"$NTMPI_ARG -ntomp $NTOMP {gpu}"

    grompp = [
        '    "$GMX" grompp',
        f"-f {stage.mdp}",
        f"-c {stage.coord_in}",
        # -r is unconditional: grompp refuses to build a .tpr when position
        # restraints are active and no reference was given, and it ignores -r
        # when they are not. Since the restraint schedule lives in the MDP's
        # define line (and is being reworked), always supplying it is the only
        # option that cannot break.
        f"-r {stage.coord_in}",
    ]
    if stage.cpt_in is not None:
        grompp.append("$T_ARG")
    grompp.extend(
        [
            f"-p {top}",
            # The reason this module exists: SOLU/MEMB/SOLV are PRISM index
            # groups, so grompp cannot resolve tc-grps without -n.
            f"-n {ndx}",
            f"-o {stage.tpr}",
            "-maxwarn $MAXWARN",
        ]
    )

    header = f"# --- {stage.name} " + "-" * max(4, 66 - len(stage.name))
    restraint_note = (
        "position restraints active (-DPOSRES); reference = "
        f"{stage.coord_in}"
        if stage.restrained
        else "unrestrained"
    )

    lines = [
        header,
        f"# {stage.note}",
        f"# mdp: {stage.mdp} | {restraint_note}",
        f'CURRENT_STAGE="{stage.name}"',
        f"if prism_redo {stage.name}; then",
        f'    echo ">>> [{stage.name}] PRISM_REDO - discarding {stage.directory}'
        ' (a redone stage invalidates every stage after it)."',
        f"    rm -rf {stage.directory}",
        "fi",
        f"mkdir -p {stage.directory}",
        # Inputs, not just the .mdp: a rerun of the previous stage leaves this
        # stage's .tpr built from coordinates that no longer exist on disk. Only
        # for an unfinished stage -- a finished one is skipped either way, and
        # PRISM_REDO is how it gets redone.
        f"[ -f {stage.gro} ] || prism_drop_stale {stage.tpr} {stage.mdp} {top} {ndx} {stage.coord_in}",
        f"if [ -f {stage.gro} ]; then",
        f'    echo ">>> [{stage.name}] already complete ({stage.gro} found) - skipping'
        f' (PRISM_REDO={stage.name} reruns it)."',
        f"elif [ -f {stage.tpr} ]; then",
        # A .tpr without a .gro means the run was interrupted. mdrun -cpi
        # resumes from the checkpoint; if the checkpoint is missing it restarts
        # the stage from the .tpr, which is still the right thing to do.
        f'    echo ">>> [{stage.name}] interrupted run detected - resuming from checkpoint."',
        f'    "$GMX" mdrun -s {stage.tpr} -deffnm {stage.prefix} {mdrun_args} -cpi {stage.cpt} -v',
        "else",
        f'    echo ">>> [{stage.name}] starting from scratch."',
    ]

    if stage.cpt_in is not None:
        lines.extend(
            [
                "    # Inherit velocities and thermostat/barostat state from the previous",
                "    # stage. Guarded so a hand-assembled directory degrades to a warning",
                "    # instead of an abort.",
                "    T_ARG=\"\"",
                f"    if [ -f {stage.cpt_in} ]; then",
                f'        T_ARG="-t {stage.cpt_in}"',
                "    else",
                f'        echo "    WARNING: {stage.cpt_in} not found - continuing from'
                f' {stage.coord_in} coordinates/velocities only."',
                "    fi",
            ]
        )

    lines.append(" ".join(grompp))
    lines.append(
        f'    "$GMX" mdrun -s {stage.tpr} -deffnm {stage.prefix} {mdrun_args} -v'
    )
    lines.append("fi")
    return "\n".join(lines) + "\n"


def _header_block(
    ordered: Sequence[str],
    ndx: str,
    gro: str,
    top: str,
    external_mdps: bool,
) -> str:
    arrow = " -> ".join(name.upper() for name in ordered)
    lines = [
        "#!/bin/bash",
        RUNSCRIPT_SENTINEL,
        "#",
        "# PRISM membrane protocol driver -- semiisotropic protein-in-bilayer MD.",
        f"# Protocol: {arrow}",
        "#",
        "# Generated by prism.membrane.runscript. Edit the .mdp files (or rebuild the",
        "# system) rather than this script -- it is regenerated on every build. An .mdp",
        "# edit takes effect for any stage that has not finished yet; to redo one that",
        '# has, run with PRISM_REDO="<stage> ..." (see the per-stage messages).',
        "#",
        "# Usage:  cd <this directory> && bash " + RUNSCRIPT_NAME,
        "#",
        f"# An identical copy is kept as {RUNSCRIPT_BACKUP_NAME}, because prism.sim",
        f"# overwrites {RUNSCRIPT_NAME} with the soluble-protein script, which passes no",
        f"# '-n {ndx}' and cannot grompp this system. If that happens, restore with",
        f"#   cp {RUNSCRIPT_BACKUP_NAME} {RUNSCRIPT_NAME}",
        "#",
        f"# Every 'gmx grompp' call below passes '-n {ndx}', and that is mandatory, not",
        "# cosmetic: the membrane MDPs thermostat on 'tc-grps = SOLU MEMB SOLV'. Those",
        "# three are PRISM-defined index groups (solute / membrane / solvent+ions), not",
        "# GROMACS built-ins, so without the index file grompp aborts with",
        '#   "Group SOLU referenced in the .mdp file was not found in the index file".',
        "#",
        "# Conventions used throughout:",
        f"#   -c/-r  each stage starts from, and restrains to, the previous stage's .gro",
        "#          (-r is required by GROMACS >= 2018 whenever restraints are active)",
        "#   -t     checkpoint carried forward wherever a stage continues velocities",
        f"#   inputs {gro}, {top}, {ndx}",
        f"#   output equilibration in ./<stage>/, production in ./{_PRODUCTION_DIR}/{_PRODUCTION_STEM}.* --",
        "#          the layout prism.analysis, prism.sim and prism.mmpbsa look for",
        "#",
        "# All paths are relative, so the run can be copied to a cluster and run as-is.",
    ]
    if external_mdps:
        lines.extend(
            [
                "#",
                "# NOTE: the .mdp files live outside this directory (../mdps/...), as they do",
                "#       for a soluble PRISM system. Copy the parent directory, not just this",
                "#       one, when moving the run.",
            ]
        )
    lines.append("#")
    return "\n".join(lines) + "\n"


def _footer_block(final: _Stage, ndx: str) -> str:
    return "\n".join(
        [
            'CURRENT_STAGE="done"',
            "trap - ERR",
            'echo ""',
            'echo "=========================================================="',
            'echo " Membrane protocol complete."',
            f'echo " Trajectory:   {final.prefix}.xtc"',
            f'echo " Final frame:  {final.gro}"',
            'echo ""',
            'echo " Before analysis, fix periodicity so the bilayer is not split"',
            'echo " across the box edge. -center makes trjconv prompt for two groups,"',
            'echo " so they are piped in (centre on SOLU, write the whole System):"',
            f"""echo "   printf 'SOLU\\nSystem\\n' | $GMX trjconv -s {final.tpr} -f {final.prefix}.xtc \\\\" """.rstrip(),
            f'echo "       -n {ndx} -pbc mol -center -ur compact -o {final.stem}_center.xtc"',
            'echo "=========================================================="',
        ]
    ) + "\n"


def _write_executable(path: str, text: str) -> str:
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    os.chmod(path, 0o755)
    return path


def write_membrane_runscript(
    system_dir: str,
    mdp_paths: Mapping[str, str],
    gro: str = SYSTEM_GRO_NAME,
    top: str = TOPOL_NAME,
    ndx: str = INDEX_NAME,
) -> str:
    """Write ``localrun.sh`` for a membrane system and return its path.

    Parameters
    ----------
    system_dir:
        Directory holding the built system (``topol.top``, ``system.gro``,
        ``index.ndx``). The script is written here and is meant to be run here.
    mdp_paths:
        Mapping of stage name -> .mdp path, as returned by
        :func:`prism.membrane.membrane_mdp.write_membrane_mdps`. Paths may be
        absolute or relative to the current working directory; they are rewritten
        relative to ``system_dir``. The stage list is *not* assumed: stages are
        ordered by :func:`order_membrane_stages`.
    gro, top, ndx:
        Input coordinate, topology and index file names, relative to
        ``system_dir``.

    Returns
    -------
    str
        Absolute path to the written script (mode 0755). An identical copy is
        written alongside it as ``run_membrane.sh``; see the module docstring.
    """
    if not mdp_paths:
        raise ValueError(
            "mdp_paths is empty: cannot write a run script with no stages. "
            "Pass the mapping returned by write_membrane_mdps()."
        )

    stage_names = [str(name) for name in mdp_paths]
    for name in stage_names:
        if not _SAFE_STAGE_NAME.match(name):
            raise ValueError(
                f"Unusable stage name {name!r}: stage names become directory names "
                "and shell words, so they must match [A-Za-z0-9][A-Za-z0-9_.-]*."
            )

    system_dir = os.path.abspath(system_dir)
    os.makedirs(system_dir, exist_ok=True)

    ordered = order_membrane_stages(stage_names)
    plan = _build_plan(ordered, {str(k): v for k, v in mdp_paths.items()}, system_dir, gro)

    # Preflight covers the build outputs plus every MDP; a stage's own inputs are
    # produced by the stage before it, so they cannot be checked up front. The
    # restraint .itp files are checked through the topology instead of by name:
    # what matters is not that they exist but that the topology includes them.
    required = [gro, top, ndx] + [stage.mdp for stage in plan]
    guards: List[str] = []
    for stage in plan:
        guards.extend(guard for guard in stage.guards if guard not in guards)
    external_mdps = any(stage.mdp.startswith("../") for stage in plan)

    parts = [
        _header_block(ordered, ndx, gro, top, external_mdps),
        "",
        _PREAMBLE,
        "",
        _preflight_block(required, ordered, guards, top),
        "",
    ]
    for stage in plan:
        parts.append(_stage_block(stage, top, ndx))
        parts.append("")
    parts.append(_footer_block(plan[-1], ndx))

    text = "\n".join(parts)
    script_path = _write_executable(os.path.join(system_dir, RUNSCRIPT_NAME), text)
    # Keep a copy under a name no other PRISM generator writes: prism.sim
    # overwrites localrun.sh with the soluble script, which has no -n index.ndx
    # and cannot grompp this system. Losing that race must cost a cp, not a
    # rebuild of the whole membrane.
    _write_executable(os.path.join(system_dir, RUNSCRIPT_BACKUP_NAME), text)
    return script_path
