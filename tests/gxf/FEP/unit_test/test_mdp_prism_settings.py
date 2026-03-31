#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Test FEP MDP generation with PRISM standard settings
"""

from prism.fep.gromacs.mdp_templates import write_fep_mdps
from prism.fep.config import FEPConfig
from prism.utils.mdp import find_mdp_parameters

from test_utils import FEPTestPaths


def test_fep_mdp_with_prism_settings():
    """Test that FEP MDPs use PRISM standard settings"""

    paths = FEPTestPaths("oMeEtPh-EtPh")
    test_dir = paths.system_dir
    output_dir = paths.get_forcefield_dir("gaff") / "test_mdp_prism_settings"
    output_dir.mkdir(exist_ok=True)

    # Load FEP config
    config = FEPConfig(str(test_dir))
    mdp_config = config.get_all_mdp_params()
    lambda_params = config.get_lambda_params()
    soft_core_params = config.get_soft_core_params()

    print("=" * 80)
    print("Testing FEP MDP Generation with PRISM Standard Settings")
    print("=" * 80)
    print()

    # Generate FEP MDPs
    write_fep_mdps(
        output_dir=str(output_dir),
        lambda_strategy=lambda_params["strategy"],
        lambda_distribution=lambda_params["distribution"],
        lambda_windows=lambda_params["windows"],
        coul_windows=lambda_params["coul_windows"],
        vdw_windows=lambda_params["vdw_windows"],
        soft_core_alpha=soft_core_params["alpha"],
        soft_core_sigma=soft_core_params["sigma"],
        config=mdp_config,
        leg_name="test",
    )

    print(f"✓ Generated MDP files in: {output_dir}")
    print()

    # Check generated files
    mdp_files = {
        "em.mdp": "Energy Minimization",
        "nvt.mdp": "NVT Equilibration",
        "npt.mdp": "NPT Equilibration",
        "prod_00.mdp": "FEP Production (lambda=0)",
        "prod_15.mdp": "FEP Production (mid-window)",
        "prod_31.mdp": "FEP Production (final window)",
    }

    for mdp_file, description in mdp_files.items():
        mdp_path = output_dir / mdp_file
        if mdp_path.exists():
            print(f"✓ {mdp_file:15s} - {description}")

            # Check key parameters using prism.utils.mdp
            checks = []

            # Check for PRISM standard parameters
            params_to_check = ["coulombtype", "DispCorr", "cutoff-scheme"]
            if mdp_file.startswith("prod_"):
                params_to_check.extend(["free-energy", "sc-alpha", "coul-lambdas", "vdw-lambdas"])

            found_params = find_mdp_parameters(str(mdp_path), params_to_check)

            if found_params.get("coulomb_type") == "PME":
                checks.append("PME electrostatics")
            if found_params.get("DispCorr") == "EnerPres":
                checks.append("Dispersion correction")
            if found_params.get("cutoff-scheme") == "Verlet":
                checks.append("Verlet cutoff scheme")

            if mdp_file.startswith("prod_"):
                if found_params.get("free-energy") == "yes":
                    checks.append("FEP enabled")
                if found_params.get("sc-alpha"):
                    checks.append("Soft-core parameters")
                if found_params.get("coul-lambdas") and found_params.get("vdw-lambdas"):
                    checks.append("Lambda vectors")

            if checks:
                print(f"  Parameters: {', '.join(checks)}")
        else:
            print(f"✗ {mdp_file:15s} - NOT FOUND")

    print()
    print("=" * 80)
    print("Detailed Check: prod_00.mdp")
    print("=" * 80)
    print()

    prod_mdp = output_dir / "prod_00.mdp"
    if prod_mdp.exists():
        # Extract key parameters using prism.utils.mdp
        params_to_check = [
            "coulombtype",
            "pme_order",
            "fourierspacing",
            "tcoupl",
            "tau-t",
            "pcoupl",
            "tau-p",
            "compressibility",
            "constraints",
            "constraint-algorithm",
            "sc-alpha",
            "sc-sigma",
            "sc-coul",
        ]

        print("Key parameters found:")
        found_params = find_mdp_parameters(str(prod_mdp), params_to_check)
        for param, value in found_params.items():
            if value is not None:
                print(f"  {param} = {value}")

    print()
    print("=" * 80)
    print("Test completed successfully!")
    print("=" * 80)


if __name__ == "__main__":
    test_fep_mdp_with_prism_settings()
