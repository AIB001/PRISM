from pathlib import Path
from types import SimpleNamespace

from prism.fep.modeling.script_writer import write_fep_master_script, write_root_scripts


def _make_layout(tmp_path: Path):
    root = tmp_path / "GMX_PROLIG_FEP"
    bound = root / "bound" / "repeat1"
    unbound = root / "unbound" / "repeat1"
    bound.mkdir(parents=True)
    unbound.mkdir(parents=True)
    return SimpleNamespace(root=root, bound_dir=bound, unbound_dir=unbound)


def test_write_root_scripts_uses_update_cpu_for_charmm(tmp_path: Path):
    layout = _make_layout(tmp_path)
    config = {
        "forcefield": {"name": "charmm36m-mut"},
        "execution": {"num_gpus": 1, "omp_threads": 4, "use_gpu_pme": True},
        "fep": {"replicas": 1},
    }

    write_root_scripts(layout, config)

    standard = (layout.bound_dir / "run_prod_standard.sh").read_text()
    repex = (layout.bound_dir / "run_prod_repex.sh").read_text()
    master = (layout.root / "run_fep.sh").read_text()

    assert 'UPDATE_MODE="${PRISM_MDRUN_UPDATE_MODE:-cpu}"' in master
    assert "MDRUN_GPU_ARGS=" in master
    assert 'UPDATE_MODE="${PRISM_MDRUN_UPDATE_MODE:-cpu}"' in standard
    assert 'UPDATE_MODE="${PRISM_MDRUN_UPDATE_MODE:-cpu}"' in repex
    assert 'UPDATE_FLAG="-update ${UPDATE_MODE}"' in master
    assert 'UPDATE_FLAG="-update ${UPDATE_MODE}"' in standard
    assert "-update gpu" not in standard


def test_write_fep_master_script_defaults_to_no_gpu_update_for_non_charmm(tmp_path: Path):
    fep_dir = tmp_path / "GMX_PROLIG_FEP"
    fep_dir.mkdir(parents=True)
    config = {
        "forcefield": {"name": "amber99sb"},
        "execution": {"num_gpus": 1, "omp_threads": 4, "use_gpu_pme": True},
        "fep": {"replicas": 1},
    }

    write_fep_master_script(fep_dir, config)

    master = (fep_dir / "run_fep.sh").read_text()
    assert 'UPDATE_MODE="${PRISM_MDRUN_UPDATE_MODE:-none}"' in master
    assert "-update gpu" not in master


def test_write_fep_master_script_uses_update_gpu_for_non_charmm(tmp_path: Path):
    fep_dir = tmp_path / "GMX_PROLIG_FEP"
    fep_dir.mkdir(parents=True)
    config = {
        "forcefield": {"name": "amber99sb"},
        "execution": {"num_gpus": 1, "omp_threads": 4, "use_gpu_pme": True, "use_gpu_update": True},
        "fep": {"replicas": 1},
    }

    write_fep_master_script(fep_dir, config)

    master = (fep_dir / "run_fep.sh").read_text()
    assert 'UPDATE_MODE="${PRISM_MDRUN_UPDATE_MODE:-gpu}"' in master
    assert 'UPDATE_FLAG="-update ${UPDATE_MODE}"' in master
    assert "-update cpu" not in master


def test_write_fep_master_script_all_gpu_branches_share_runtime_args(tmp_path: Path):
    fep_dir = tmp_path / "GMX_PROLIG_FEP"
    fep_dir.mkdir(parents=True)
    config = {
        "forcefield": {"name": "charmm36m-mut"},
        "execution": {"num_gpus": 1, "omp_threads": 4, "use_gpu_pme": True},
        "fep": {"replicas": 1},
    }

    write_fep_master_script(fep_dir, config)

    master = (fep_dir / "run_fep.sh").read_text()
    assert master.count("${MDRUN_GPU_ARGS}") == 6


def test_write_fep_master_script_exports_gpu_and_omp_runtime_controls(tmp_path: Path):
    fep_dir = tmp_path / "GMX_PROLIG_FEP"
    fep_dir.mkdir(parents=True)
    config = {
        "forcefield": {"name": "amber99sb"},
        "execution": {"num_gpus": 4, "total_cpus": 40, "use_gpu_pme": True},
        "fep": {"replicas": 1},
    }

    write_fep_master_script(fep_dir, config)

    master = (fep_dir / "run_fep.sh").read_text()
    assert 'GPU_COUNT="${PRISM_NUM_GPUS:-${DEFAULT_NUM_GPUS}}"' in master
    assert 'PARALLEL_WINDOWS="${PRISM_PARALLEL_WINDOWS:-${DEFAULT_PARALLEL_WINDOWS}}"' in master
    assert 'elif [ -n "${PRISM_TOTAL_CPUS:-}" ]; then' in master
    assert 'export OMP_NUM_THREADS="10"' in master
    assert 'export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"' in master
    assert 'export OMP_PLACES="${OMP_PLACES:-cores}"' in master
    assert "gpu_id=$(resolve_gpu_id)" in master
    assert 'export CUDA_VISIBLE_DEVICES="${gpu_id}"' in master
    assert '-gpu_id 0 -pinoffset "${cpu_offset}" -pinstride 1' in master


def test_generated_scripts_support_runtime_resource_overrides(tmp_path: Path):
    layout = _make_layout(tmp_path)
    config = {
        "forcefield": {"name": "amber99sb"},
        "execution": {"num_gpus": 4, "parallel_windows": 4, "omp_threads": 10, "use_gpu_pme": True},
        "fep": {"replicas": 1},
    }

    write_root_scripts(layout, config)

    standard = (layout.bound_dir / "run_prod_standard.sh").read_text()
    repex = (layout.bound_dir / "run_prod_repex.sh").read_text()
    master = (layout.root / "run_fep.sh").read_text()

    assert 'NUM_GPUS="${PRISM_NUM_GPUS:-${DEFAULT_NUM_GPUS}}"' in standard
    assert 'PARALLEL_WINDOWS="${PRISM_PARALLEL_WINDOWS:-${DEFAULT_PARALLEL_WINDOWS}}"' in standard
    assert 'OMP_THREADS="${PRISM_OMP_THREADS:-${OMP_NUM_THREADS:-${DEFAULT_OMP_THREADS}}}"' in standard
    assert 'GPU_COUNT="${PRISM_NUM_GPUS:-${DEFAULT_NUM_GPUS}}"' in master
    assert 'PARALLEL_WINDOWS="${PRISM_PARALLEL_WINDOWS:-${DEFAULT_PARALLEL_WINDOWS}}"' in master
    assert 'elif [ -n "${PRISM_TOTAL_CPUS:-}" ]; then' in master
    assert "PRISM_FEP_MODE=repex" not in master
    assert 'echo "Configured GPUs: $NUM_GPUS"' in repex
