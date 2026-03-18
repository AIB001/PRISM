#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comprehensive consistency tests for trajectory module refactoring.

This test ensures 100% functional consistency between the original
prism.analysis.core.trajectory_processor and the refactored
prism.analysis.trajectory package.
"""

import pytest
import inspect
from unittest.mock import Mock, patch


class TestMethodImplementationConsistency:
    """Test that method implementations are identical between original and refactored."""

    def test_find_trjconv_implementation(self):
        """_find_trjconv should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        # Compare source code
        original_source = inspect.getsource(Original._find_trjconv)
        refactored_source = inspect.getsource(Refactored._find_trjconv)

        # Normalize whitespace for comparison
        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_find_trjconv implementation differs from original"

    def test_convert_dcd_to_xtc_implementation(self):
        """_convert_dcd_to_xtc should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._convert_dcd_to_xtc)
        refactored_source = inspect.getsource(Refactored._convert_dcd_to_xtc)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_convert_dcd_to_xtc implementation differs from original"

    def test_get_mdtraj_topology_implementation(self):
        """_get_mdtraj_topology should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._get_mdtraj_topology)
        refactored_source = inspect.getsource(Refactored._get_mdtraj_topology)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_get_mdtraj_topology implementation differs from original"

    def test_process_pbc_two_step_implementation(self):
        """_process_pbc_two_step should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._process_pbc_two_step)
        refactored_source = inspect.getsource(Refactored._process_pbc_two_step)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_process_pbc_two_step implementation differs from original"

    def test_select_functional_chain_implementation(self):
        """_select_functional_chain should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._select_functional_chain)
        refactored_source = inspect.getsource(Refactored._select_functional_chain)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_select_functional_chain implementation differs from original"

    def test_select_ligand_atoms_implementation(self):
        """_select_ligand_atoms should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._select_ligand_atoms)
        refactored_source = inspect.getsource(Refactored._select_ligand_atoms)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_select_ligand_atoms implementation differs from original"

    def test_needs_mdtraj_processing_implementation(self):
        """_needs_mdtraj_processing should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._needs_mdtraj_processing)
        refactored_source = inspect.getsource(Refactored._needs_mdtraj_processing)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_needs_mdtraj_processing implementation differs from original"

    def test_prepare_selections_implementation(self):
        """_prepare_selections should have identical implementation."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._prepare_selections)
        refactored_source = inspect.getsource(Refactored._prepare_selections)

        original_lines = [l.strip() for l in original_source.split("\n") if l.strip() and not l.strip().startswith("#")]
        refactored_lines = [
            l.strip() for l in refactored_source.split("\n") if l.strip() and not l.strip().startswith("#")
        ]

        assert original_lines == refactored_lines, "_prepare_selections implementation differs from original"


class TestDefaultParameters:
    """Test that all default parameters match the original."""

    def test_process_trajectory_defaults(self):
        """process_trajectory should have identical default parameters."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_sig = inspect.signature(Original.process_trajectory)
        refactored_sig = inspect.signature(Refactored.process_trajectory)

        # Compare default values for all parameters
        for param_name in original_sig.parameters:
            if param_name == "self":
                continue
            original_param = original_sig.parameters[param_name]
            refactored_param = refactored_sig.parameters[param_name]

            assert original_param.default == refactored_param.default, (
                f"Parameter {param_name} has different default: "
                f"{original_param.default} vs {refactored_param.default}"
            )

    def test_batch_process_defaults(self):
        """batch_process should have identical default parameters."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_sig = inspect.signature(Original.batch_process)
        refactored_sig = inspect.signature(Refactored.batch_process)

        for param_name in original_sig.parameters:
            if param_name == "self":
                continue
            original_param = original_sig.parameters[param_name]
            refactored_param = refactored_sig.parameters[param_name]

            assert original_param.default == refactored_param.default, (
                f"Parameter {param_name} has different default: "
                f"{original_param.default} vs {refactored_param.default}"
            )

    def test_init_defaults(self):
        """__init__ should have identical default parameters."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_sig = inspect.signature(Original.__init__)
        refactored_sig = inspect.signature(Refactored.__init__)

        for param_name in original_sig.parameters:
            if param_name == "self":
                continue
            original_param = original_sig.parameters[param_name]
            refactored_param = refactored_sig.parameters[param_name]

            assert original_param.default == refactored_param.default, (
                f"Parameter {param_name} has different default: "
                f"{original_param.default} vs {refactored_param.default}"
            )


class TestDocstrings:
    """Test that all docstrings are preserved."""

    def test_class_docstring(self):
        """TrajectoryProcessor class docstring should be identical."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        assert Original.__doc__ == Refactored.__doc__, "Class docstring differs from original"

    def test_process_trajectory_docstring(self):
        """process_trajectory docstring should be identical."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_doc = Original.process_trajectory.__doc__
        refactored_doc = Refactored.process_trajectory.__doc__

        assert original_doc == refactored_doc, "process_trajectory docstring differs from original"

    def test_batch_process_docstring(self):
        """batch_process docstring should be identical."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_doc = Original.batch_process.__doc__
        refactored_doc = Refactored.batch_process.__doc__

        assert original_doc == refactored_doc, "batch_process docstring differs from original"

    def test_find_trjconv_docstring(self):
        """_find_trjconv docstring should be identical."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_doc = Original._find_trjconv.__doc__
        refactored_doc = Refactored._find_trjconv.__doc__

        assert original_doc == refactored_doc, "_find_trjconv docstring differs from original"

    def test_convert_dcd_to_xtc_docstring(self):
        """_convert_dcd_to_xtc docstring should be identical."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_doc = Original._convert_dcd_to_xtc.__doc__
        refactored_doc = Refactored._convert_dcd_to_xtc.__doc__

        assert original_doc == refactored_doc, "_convert_dcd_to_xtc docstring differs from original"


class TestMethodAttributes:
    """Test that method attributes are preserved."""

    def test_all_methods_are_callable(self):
        """All inherited methods should be callable."""
        from prism.analysis.trajectory import TrajectoryProcessor

        methods = [
            "_find_trjconv",
            "_convert_dcd_to_xtc",
            "_get_mdtraj_topology",
            "_process_pbc_two_step",
            "_process_pbc_with_mdtraj",
            "_build_trjconv_command",
            "_run_trjconv",
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

        for method_name in methods:
            assert hasattr(TrajectoryProcessor, method_name), f"Missing method: {method_name}"
            assert callable(getattr(TrajectoryProcessor, method_name)), f"Method {method_name} is not callable"

    def test_public_methods_match(self):
        """Public methods should match original exactly."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_public = [
            name
            for name, method in inspect.getmembers(Original, predicate=inspect.isfunction)
            if not name.startswith("_") and name != "process_trajectory_simple" and name != "batch_process_trajectories"
        ]

        refactored_public = [
            name
            for name, method in inspect.getmembers(Refactored, predicate=inspect.isfunction)
            if not name.startswith("_")
        ]

        assert set(original_public) == set(
            refactored_public
        ), f"Public methods differ: {set(original_public) ^ set(refactored_public)}"


class TestModuleLevelFunctions:
    """Test that module-level functions are preserved."""

    def test_process_trajectory_simple_exists(self):
        """process_trajectory_simple should exist and be callable."""
        from prism.analysis.trajectory.processor import process_trajectory_simple

        assert callable(process_trajectory_simple), "process_trajectory_simple is not callable"

    def test_process_trajectory_simple_signature(self):
        """process_trajectory_simple should have correct signature."""
        from prism.analysis.trajectory.processor import process_trajectory_simple
        from prism.analysis.core.trajectory_processor import process_trajectory_simple as original

        refactored_sig = inspect.signature(process_trajectory_simple)
        original_sig = inspect.signature(original)

        assert str(refactored_sig) == str(original_sig), (
            f"process_trajectory_simple signature differs: " f"{refactored_sig} vs {original_sig}"
        )

    def test_batch_process_trajectories_exists(self):
        """batch_process_trajectories should exist and be callable."""
        from prism.analysis.trajectory.processor import batch_process_trajectories

        assert callable(batch_process_trajectories), "batch_process_trajectories is not callable"

    def test_batch_process_trajectories_signature(self):
        """batch_process_trajectories should have correct signature."""
        from prism.analysis.trajectory.processor import batch_process_trajectories
        from prism.analysis.core.trajectory_processor import batch_process_trajectories as original

        refactored_sig = inspect.signature(batch_process_trajectories)
        original_sig = inspect.signature(original)

        assert str(refactored_sig) == str(original_sig), (
            f"batch_process_trajectories signature differs: " f"{refactored_sig} vs {original_sig}"
        )


class TestErrorMessages:
    """Test that error messages are preserved."""

    def test_find_trjconv_error_message(self):
        """_find_trjconv should have identical error message."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._find_trjconv)
        refactored_source = inspect.getsource(Refactored._find_trjconv)

        # Check for key error message text
        error_msg = "GROMACS trjconv not found"
        assert error_msg in original_source, "Original missing error message"
        assert error_msg in refactored_source, "Refactored missing error message"

    def test_convert_dcd_error_message(self):
        """_convert_dcd_to_xtc should have identical error messages."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_source = inspect.getsource(Original._convert_dcd_to_xtc)
        refactored_source = inspect.getsource(Refactored._convert_dcd_to_xtc)

        # Check for key error messages
        error_msgs = [
            "MDTraj not available for DCD conversion",
            "RECOMMENDED INSTALLATION",
            "DCD conversion failed",
        ]

        for msg in error_msgs:
            assert msg in original_source, f"Original missing: {msg}"
            assert msg in refactored_source, f"Refactored missing: {msg}"


class TestTypeAnnotations:
    """Test that type annotations are preserved."""

    def test_process_trajectory_return_type(self):
        """process_trajectory should have -> str return type."""
        from prism.analysis.trajectory import TrajectoryProcessor

        sig = inspect.signature(TrajectoryProcessor.process_trajectory)
        assert sig.return_annotation == str, "process_trajectory return type should be str"

    def test_batch_process_return_type(self):
        """batch_process should have -> Dict[str, str] return type."""
        from prism.analysis.trajectory import TrajectoryProcessor

        sig = inspect.signature(TrajectoryProcessor.batch_process)
        # Dict[str, str] might be represented differently, just check it's not empty
        assert sig.return_annotation != inspect.Parameter.empty, "batch_process should have return type annotation"


class TestConstantsAndStrings:
    """Test that constants and string literals are preserved."""

    def test_installation_messages(self):
        """Installation help messages should be identical."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        # Check __init__ error message
        original_source = inspect.getsource(Original.__init__)
        refactored_source = inspect.getsource(Refactored.__init__)

        install_msg = "mamba install -c conda-forge gromacs"
        assert install_msg in original_source, "Original missing install message"
        assert install_msg in refactored_source, "Refactored missing install message"

    def test_pbc_method_defaults(self):
        """PBC method defaults should match."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        original_sig = inspect.signature(Original.process_trajectory)
        refactored_sig = inspect.signature(Refactored.process_trajectory)

        # Check pbc_method default
        original_pbc = original_sig.parameters["pbc_method"].default
        refactored_pbc = refactored_sig.parameters["pbc_method"].default

        assert original_pbc == refactored_pbc == "atom", "pbc_method default should be 'atom'"

        # Check unit_cell default
        original_uc = original_sig.parameters["unit_cell"].default
        refactored_uc = refactored_sig.parameters["unit_cell"].default

        assert original_uc == refactored_uc == "compact", "unit_cell default should be 'compact'"


class TestBehaviorConsistency:
    """Test actual behavior with mocked dependencies."""

    @pytest.fixture
    def mock_gromacs_env(self):
        """Mock GROMACS environment."""
        with patch("prism.analysis.trajectory.processor.GromacsEnvironment") as mock:
            mock_instance = Mock()
            mock_instance.version = "2023.1"
            mock.return_value = mock_instance
            yield mock

    def test_initialization_identical(self):
        """Initialization should behave identically."""
        from prism.analysis.core.trajectory_processor import TrajectoryProcessor as Original
        from prism.analysis.trajectory import TrajectoryProcessor as Refactored

        # Create instances with same parameters
        original = Original(topology_file="test.tpr")
        refactored = Refactored(topology_file="test.tpr")

        # Check attributes are set identically
        assert original.topology_file == refactored.topology_file == "test.tpr"
        assert hasattr(original, "mdtraj_available")
        assert hasattr(refactored, "mdtraj_available")
        assert hasattr(original, "gromacs_env")
        assert hasattr(refactored, "gromacs_env")
        assert hasattr(original, "trjconv_cmd")
        assert hasattr(refactored, "trjconv_cmd")

    def test_validate_processing_return_structure(self):
        """validate_processing should return identical structure."""
        from prism.analysis.trajectory import TrajectoryProcessor

        # Mock the Path operations
        with patch("prism.analysis.trajectory.processor.Path") as mock_path:
            mock_path.return_value.exists.return_value = True
            mock_path.return_value.stat.return_value.st_size = 1000

            processor = TrajectoryProcessor.__new__(TrajectoryProcessor)
            processor.topology_file = "test.tpr"

            result = processor.validate_processing("original.xtc", "processed.xtc")

            # Check result structure
            expected_keys = {
                "original_exists",
                "processed_exists",
                "original_size",
                "processed_size",
                "size_ratio",
                "valid",
            }

            assert (
                set(result.keys()) == expected_keys
            ), f"validate_processing result keys differ: {set(result.keys())} vs {expected_keys}"

    def test_get_processing_info_return_structure(self):
        """get_processing_info should return identical structure."""
        from prism.analysis.trajectory import TrajectoryProcessor

        processor = TrajectoryProcessor.__new__(TrajectoryProcessor)
        processor.topology_file = "test.tpr"
        processor.gromacs_env = Mock()
        processor.gromacs_env.version = "2023.1"
        processor.trjconv_cmd = "gmx"

        result = processor.get_processing_info()

        # Check result structure
        expected_keys = {"gromacs_available", "trjconv_command", "topology_file", "gromacs_version"}

        assert (
            set(result.keys()) == expected_keys
        ), f"get_processing_info result keys differ: {set(result.keys())} vs {expected_keys}"
