#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path

from prism.fep.config import FEPConfig
from prism.fep.modeling.core import FEPScaffoldLayout
from prism.fep.modeling.script_writer import write_root_scripts


def _make_layout(tmp_path: Path) -> FEPScaffoldLayout:
    root = tmp_path / "GMX_PROLIG_FEP"
    common_dir = root / "common"
    hybrid_dir = common_dir / "hybrid"
    protein_dir = common_dir / "protein"
    bound_dir = root / "bound"
    unbound_dir = root / "unbound"

    for directory in (hybrid_dir, protein_dir, bound_dir / "input", unbound_dir / "input"):
        directory.mkdir(parents=True, exist_ok=True)

    return FEPScaffoldLayout(
        root=root,
        common_dir=common_dir,
        hybrid_dir=hybrid_dir,
        protein_dir=protein_dir,
        bound_dir=bound_dir,
        unbound_dir=unbound_dir,
    )


def test_write_root_scripts_standard_and_repex_helpers(tmp_path):
    layout = _make_layout(tmp_path)

    write_root_scripts(
        layout,
        {
            "execution": {
                "mode": "standard",
                "num_gpus": 4,
                "parallel_windows": 4,
                "omp_threads": 8,
                "use_gpu_pme": True,
            },
            "fep": {"replicas": 1},
        },
    )

    root_script = (layout.root / "run_fep.sh").read_text()
    standard_script = (layout.bound_dir / "run_prod_standard.sh").read_text()
    repex_script = (layout.bound_dir / "run_prod_repex.sh").read_text()
    dispatch_script = (layout.bound_dir / "run_prod.sh").read_text()

    assert 'EXECUTION_MODE="${PRISM_FEP_MODE:-standard}"' in root_script
    assert "./run_prod_standard.sh" in root_script
    assert "./run_prod_repex.sh" in root_script
    assert "PRISM_FEP_MODE" in root_script

    assert "Concurrent windows: 4" in standard_script
    assert "CUDA_VISIBLE_DEVICES" in standard_script
    assert ".run_${lambda_idx}" in standard_script
    assert "-ntomp 8" in standard_script

    assert "mpirun -oversubscribe -np" in repex_script
    assert "gmx_mpi mdrun" in repex_script
    assert "-multidir" in repex_script
    assert "-replex 1000" in repex_script
    assert "NUM_GPUS=4" in repex_script
    assert "export OMP_NUM_THREADS=8" in repex_script

    assert "run_prod_standard.sh" in dispatch_script


def test_write_root_scripts_repex_dispatch_uses_shared_gpu_pool(tmp_path):
    layout = _make_layout(tmp_path)

    write_root_scripts(
        layout,
        {
            "execution": {
                "mode": "repex",
                "num_gpus": 4,
                "omp_threads": 2,
                "use_gpu_pme": False,
            }
        },
    )

    root_script = (layout.root / "run_fep.sh").read_text()
    repex_script = (layout.unbound_dir / "run_prod_repex.sh").read_text()
    dispatch_script = (layout.unbound_dir / "run_prod.sh").read_text()

    assert 'EXECUTION_MODE="${PRISM_FEP_MODE:-repex}"' in root_script
    assert "Production mode: repex" in repex_script
    assert "NUM_GPUS=4" in repex_script
    assert "export OMP_NUM_THREADS=2" in repex_script
    assert '-gpu_id "${gpu_ids}"' in repex_script
    assert "-replex 1000" in repex_script
    assert "-pme gpu" not in repex_script
    assert "run_prod_repex.sh" in dispatch_script


def test_fep_config_exposes_execution_params(tmp_path):
    work_dir = tmp_path / "case"
    work_dir.mkdir()
    (work_dir / "config.yaml").write_text("forcefield:\n  type: gaff2\n")
    (work_dir / "fep.yaml").write_text("execution:\n  mode: repex\n  num_gpus: 4\n  omp_threads: 2\n")

    config = FEPConfig(str(work_dir))
    execution = config.get_execution_params()

    assert execution["mode"] == "repex"
    assert execution["num_gpus"] == 4
    assert execution["omp_threads"] == 2
    assert "execution" in config.to_dict()
