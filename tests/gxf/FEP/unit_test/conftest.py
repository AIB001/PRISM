def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow (uses real data files)")


from pathlib import Path


def resolve_fep_case_dir(case_name: str) -> Path:
    """Resolve a case directory across the current test-data layout."""
    candidates = [
        Path("tests/gxf/FEP/unit_test") / case_name,
        Path("tests/gxf/FEP/test/hif2a") / case_name,
        Path("tests/gxf/FEP/test") / case_name,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Case directory not found for {case_name}")
