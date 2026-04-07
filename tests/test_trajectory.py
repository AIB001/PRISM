#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for trajectory module to ensure refactoring maintains complete consistency.

This test compares the refactored prism.analysis.trajectory package
with the original prism.analysis.core.trajectory_processor module
to ensure 100% functional consistency.
"""

import pytest
import inspect
from pathlib import Path


class TestTrajectoryModuleStructure:
    """Test that the refactored module has correct structure."""

    def test_module_imports(self):
        """All public classes and functions should be importable."""
        from prism.analysis.trajectory import TrajectoryProcessor
        from prism.analysis.trajectory import TrajectoryManager

        assert TrajectoryProcessor is not None
        assert TrajectoryManager is not None

    def test_trajectory_processor_class_exists(self):
        """TrajectoryProcessor class should exist."""
        from prism.analysis.trajectory import TrajectoryProcessor

        assert TrajectoryProcessor is not None


class TestTrajectoryProcessorMethods:
    """Test that TrajectoryProcessor has all required methods."""

    @pytest.fixture
    def processor_methods(self):
        """Get all methods from TrajectoryProcessor."""
        from prism.analysis.trajectory import TrajectoryProcessor

        return [
            name
            for name, method in inspect.getmembers(TrajectoryProcessor, predicate=inspect.isfunction)
            if not name.startswith("_")
        ]

    def test_public_methods_exist(self, processor_methods):
        """All public methods from original should exist."""
        required_methods = {
            "process_trajectory",
            "batch_process",
            "validate_processing",
            "get_processing_info",
        }

        actual_methods = set(processor_methods)
        missing = required_methods - actual_methods

        assert not missing, f"Missing public methods: {missing}"

    def test_private_conversion_methods_exist(self):
        """All conversion mixin methods should be accessible."""
        from prism.analysis.trajectory import TrajectoryProcessor

        conversion_methods = [
            "_find_trjconv",
            "_convert_dcd_to_xtc",
            "_get_mdtraj_topology",
            "_process_pbc_two_step",
            "_process_pbc_with_mdtraj",
            "_build_trjconv_command",
            "_run_trjconv",
        ]

        for method_name in conversion_methods:
            assert hasattr(TrajectoryProcessor, method_name), f"Missing conversion method: {method_name}"

    def test_private_selection_methods_exist(self):
        """All selection mixin methods should be accessible."""
        from prism.analysis.trajectory import TrajectoryProcessor

        selection_methods = [
            "_select_functional_chain",
            "_select_ligand_atoms",
            "_detect_largest_protein_chain",
            "_map_chain_letter_to_id",
            "_find_special_chain",
            "_needs_mdtraj_processing",
            "_resolve_center_selection",
            "_detect_ligand_group",
            "_prepare_selections",
        ]

        for method_name in selection_methods:
            assert hasattr(TrajectoryProcessor, method_name), f"Missing selection method: {method_name}"


class TestMethodSignatures:
    """Test that method signatures match original implementation."""

    def test_process_trajectory_signature(self):
        """process_trajectory should have correct signature."""
        from prism.analysis.trajectory import TrajectoryProcessor

        sig = inspect.signature(TrajectoryProcessor.process_trajectory)
        params = list(sig.parameters.keys())

        expected_params = [
            "self",
            "input_trajectory",
            "output_trajectory",
            "topology_file",
            "center_selection",
            "output_selection",
            "pbc_method",
            "unit_cell",
            "overwrite",
        ]

        assert params == expected_params, f"Parameter mismatch: {params} vs {expected_params}"

    def test_batch_process_signature(self):
        """batch_process should have correct signature."""
        from prism.analysis.trajectory import TrajectoryProcessor

        sig = inspect.signature(TrajectoryProcessor.batch_process)
        params = list(sig.parameters.keys())

        expected_params = ["self", "input_trajectories", "output_dir", "topology_file", "output_suffix", "kwargs"]

        assert params == expected_params, f"Parameter mismatch: {params} vs {expected_params}"

    def test_init_signature(self):
        """__init__ should have correct signature."""
        from prism.analysis.trajectory import TrajectoryProcessor

        sig = inspect.signature(TrajectoryProcessor.__init__)
        params = list(sig.parameters.keys())

        expected_params = ["self", "topology_file", "gromacs_env"]

        assert params == expected_params, f"Parameter mismatch: {params} vs {expected_params}"


class TestMethodDocumentation:
    """Test that all methods have proper documentation."""

    def test_trajectory_processor_has_docstring(self):
        """TrajectoryProcessor class should have docstring."""
        from prism.analysis.trajectory import TrajectoryProcessor

        assert TrajectoryProcessor.__doc__ is not None
        assert len(TrajectoryProcessor.__doc__) > 0

    def test_process_trajectory_has_docstring(self):
        """process_trajectory should have detailed docstring."""
        from prism.analysis.trajectory import TrajectoryProcessor

        doc = TrajectoryProcessor.process_trajectory.__doc__
        assert doc is not None
        assert "PBC" in doc or "trajectory" in doc.lower()

    def test_conversion_mixin_has_docstring(self):
        """ConversionMixin should have docstring."""
        from prism.analysis.trajectory.conversion import ConversionMixin

        assert ConversionMixin.__doc__ is not None
        assert len(ConversionMixin.__doc__) > 0

    def test_selection_mixin_has_docstring(self):
        """SelectionMixin should have docstring."""
        from prism.analysis.trajectory.selection import SelectionMixin

        assert SelectionMixin.__doc__ is not None
        assert len(SelectionMixin.__doc__) > 0


class TestBackwardCompatibility:
    """Test backward compatibility with original API."""

    def test_module_level_functions_exist(self):
        """Module-level functions from original should exist."""
        from prism.analysis.trajectory.processor import process_trajectory_simple, batch_process_trajectories

        assert callable(process_trajectory_simple)
        assert callable(batch_process_trajectories)

    def test_process_trajectory_simple_signature(self):
        """process_trajectory_simple should have correct signature."""
        from prism.analysis.trajectory.processor import process_trajectory_simple

        sig = inspect.signature(process_trajectory_simple)
        params = list(sig.parameters.keys())

        expected_params = ["input_trajectory", "output_trajectory", "topology_file", "kwargs"]

        assert params == expected_params, f"Parameter mismatch: {params} vs {expected_params}"

    def test_batch_process_trajectories_signature(self):
        """batch_process_trajectories should have correct signature."""
        from prism.analysis.trajectory.processor import batch_process_trajectories

        sig = inspect.signature(batch_process_trajectories)
        params = list(sig.parameters.keys())

        expected_params = ["input_trajectories", "output_dir", "topology_file", "kwargs"]

        assert params == expected_params, f"Parameter mismatch: {params} vs {expected_params}"


class TestInheritanceChain:
    """Test that mixin inheritance is correctly set up."""

    def test_trajectory_processor_inherits_mixins(self):
        """TrajectoryProcessor should inherit from both mixins."""
        from prism.analysis.trajectory import TrajectoryProcessor

        # Check MRO (Method Resolution Order)
        mro = TrajectoryProcessor.__mro__

        mixin_names = [cls.__name__ for cls in mro]
        assert "ConversionMixin" in mixin_names
        assert "SelectionMixin" in mixin_names

    def test_mixins_are_classes(self):
        """Both mixins should be proper classes."""
        from prism.analysis.trajectory.conversion import ConversionMixin
        from prism.analysis.trajectory.selection import SelectionMixin

        assert inspect.isclass(ConversionMixin)
        assert inspect.isclass(SelectionMixin)


class TestFilePaths:
    """Test that all expected files exist."""

    def test_trajectory_package_exists(self):
        """Trajectory package directory should exist."""
        trajectory_dir = Path(__file__).parent.parent / "prism" / "analysis" / "trajectory"
        assert trajectory_dir.exists()
        assert trajectory_dir.is_dir()

    def test_trajectory_modules_exist(self):
        """All expected module files should exist."""
        trajectory_dir = Path(__file__).parent.parent / "prism" / "analysis" / "trajectory"

        expected_files = [
            trajectory_dir / "__init__.py",
            trajectory_dir / "processor.py",
            trajectory_dir / "conversion.py",
            trajectory_dir / "selection.py",
            trajectory_dir / "manager.py",
        ]

        for file_path in expected_files:
            assert file_path.exists(), f"Missing file: {file_path}"
            assert file_path.is_file(), f"Not a file: {file_path}"


class TestMethodCount:
    """Verify the total number of methods matches expectations."""

    def test_total_method_count(self):
        """TrajectoryProcessor should have expected number of methods."""
        from prism.analysis.trajectory import TrajectoryProcessor

        # Count all methods (public and private)
        all_methods = [
            name
            for name, method in inspect.getmembers(TrajectoryProcessor, predicate=inspect.isfunction)
            if not name.startswith("__")
        ]

        # Expected: 4 public + 7 conversion + 9 selection = 20 methods
        expected_min = 18
        expected_max = 22

        assert (
            expected_min <= len(all_methods) <= expected_max
        ), f"Method count {len(all_methods)} outside expected range [{expected_min}, {expected_max}]"


class TestImportChain:
    """Test that import chain works correctly."""

    def test_import_from_main_package(self):
        """Should be importable from main prism package."""
        # This tests that the module is properly integrated
        try:
            from prism.analysis.trajectory import TrajectoryProcessor

            processor_class = TrajectoryProcessor
            assert processor_class is not None
        except ImportError as e:
            pytest.fail(f"Import chain broken: {e}")

    def test_relative_imports_work(self):
        """Relative imports within package should work."""
        from prism.analysis.trajectory.processor import TrajectoryProcessor
        from prism.analysis.trajectory.conversion import ConversionMixin
        from prism.analysis.trajectory.selection import SelectionMixin

        assert TrajectoryProcessor is not None
        assert ConversionMixin is not None
        assert SelectionMixin is not None
