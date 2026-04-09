"""Trajectory processing tools (PBC correction, centering)."""

import os
import json
import traceback

from ._common import _StdoutToStderr, _ensure_prism_importable, logger


def register(mcp):
    @mcp.tool()
    def process_trajectory(
        input_trajectory: str,
        output_trajectory: str,
        topology_file: str,
        center_selection: str = "Protein",
        pbc_method: str = "mol",
        overwrite: bool = False,
    ) -> str:
        """Process a trajectory to fix PBC artifacts and center the system.

        Applies GROMACS trjconv processing to fix periodic boundary condition
        artifacts and center the protein-ligand complex. Automatically handles
        DCD format conversion if needed.

        This should be run BEFORE analysis tools when working with raw
        simulation trajectories.

        Args:
            input_trajectory: Absolute path to input trajectory (.xtc, .trr, .dcd).
            output_trajectory: Absolute path for output processed trajectory (.xtc recommended).
            topology_file: Absolute path to topology file (.tpr, .gro, .pdb).
            center_selection: GROMACS selection for centering. Default: "Protein".
            pbc_method: PBC correction method. Default: "mol".
                Options: "mol", "atom", "nojump", "whole".
            overwrite: Overwrite existing output file. Default: false.

        Returns:
            JSON with success status and output file path.
        """
        logger.info(f"process_trajectory: {input_trajectory} -> {output_trajectory}")

        errors = []
        if not os.path.exists(input_trajectory):
            errors.append(f"Input trajectory not found: {input_trajectory}")
        if not os.path.exists(topology_file):
            errors.append(f"Topology file not found: {topology_file}")
        if not overwrite and os.path.exists(output_trajectory):
            errors.append(f"Output file already exists: {output_trajectory}. Set overwrite=true to replace.")
        if errors:
            return json.dumps({"success": False, "errors": errors}, indent=2)

        try:
            with _StdoutToStderr():
                _ensure_prism_importable()
                import prism

                result_path = prism.process_trajectory(
                    input_trajectory=input_trajectory,
                    output_trajectory=output_trajectory,
                    topology_file=topology_file,
                    center_selection=center_selection,
                    pbc_method=pbc_method,
                    overwrite=overwrite,
                )

            return json.dumps(
                {
                    "success": True,
                    "output_trajectory": result_path,
                    "message": f"Trajectory processed successfully: {result_path}",
                },
                indent=2,
            )

        except Exception as e:
            logger.error(f"process_trajectory failed: {e}\n{traceback.format_exc()}")
            return json.dumps({"success": False, "error": str(e), "traceback": traceback.format_exc()}, indent=2)
