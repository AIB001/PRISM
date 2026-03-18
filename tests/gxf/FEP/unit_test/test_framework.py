#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM-FEP Tests

Unit tests for FEP module framework.
"""

import pytest
import numpy as np
from prism.fep.core.mapping import Atom, AtomMapping, DistanceAtomMapper
from prism.fep.core.hybrid_topology import HybridAtom, HybridTopologyBuilder
from prism.fep.gromacs.itp_builder import ITPBuilder
from prism.fep.analysis.xvg_parser import XVGParser
from prism.fep.analysis.estimators import FEstimator


class TestAtomMappingFramework:
    """Test atom mapping data structures"""

    def test_atom_creation(self):
        """Test Atom object creation"""
        atom = Atom(name="C1", element="C", coord=np.array([0.0, 0.0, 0.0]), charge=-0.12, atom_type="ca", index=1)
        assert atom.name == "C1"
        assert atom.element == "C"
        assert atom.charge == -0.12
        assert atom.atom_type == "ca"
        assert atom.index == 1

    def test_atom_mapping_creation(self):
        """Test AtomMapping object creation"""
        mapping = AtomMapping(common=[], transformed_a=[], transformed_b=[], surrounding_a=[], surrounding_b=[])
        assert mapping.common == []
        assert mapping.transformed_a == []
        assert mapping.transformed_b == []
        assert mapping.surrounding_a == []
        assert mapping.surrounding_b == []

    def test_distance_mapper_initialization(self):
        """Test DistanceAtomMapper initialization"""
        mapper = DistanceAtomMapper(dist_cutoff=0.6, charge_cutoff=0.05)
        assert mapper.dist_cutoff == 0.6
        assert mapper.charge_cutoff == 0.05

    def test_distance_mapper_default_values(self):
        """Test DistanceAtomMapper default values"""
        mapper = DistanceAtomMapper()
        assert mapper.dist_cutoff == 0.6
        assert mapper.charge_cutoff == 0.05


class TestDualTopologyFramework:
    """Test dual topology data structures"""

    def test_hybrid_atom_common(self):
        """Test HybridAtom for common atom"""
        atom = HybridAtom(
            name="C1",
            index=1,
            state_a_type="ca",
            state_a_charge=-0.12,
            # Common atom - no state B
            state_b_type=None,
            state_b_charge=None,
        )
        assert atom.name == "C1"
        assert atom.state_a_type == "ca"
        assert atom.state_b_type is None
        assert atom.state_b_charge is None

    def test_hybrid_atom_transformed(self):
        """Test HybridAtom for transformed atom"""
        atom = HybridAtom(
            name="C2A",
            index=2,
            state_a_type="ca",
            state_a_charge=-0.15,
            state_b_type="DUM",  # Dummy in state B
            state_b_charge=0.0,
        )
        assert atom.name == "C2A"
        assert atom.state_a_type == "ca"
        assert atom.state_b_type == "DUM"
        assert atom.state_b_charge == 0.0

    def test_hybrid_topology_builder_initialization(self):
        """Test HybridTopologyBuilder initialization"""
        builder = HybridTopologyBuilder(charge_strategy="mean")
        assert builder.charge_strategy == "mean"

    def test_hybrid_topology_builder_invalid_strategy(self):
        """Test HybridTopologyBuilder with invalid strategy"""
        with pytest.raises(ValueError):
            HybridTopologyBuilder(charge_strategy="invalid")


class TestITPBuilderFramework:
    """Test ITP builder framework"""

    def test_itp_builder_initialization(self):
        """Test ITPBuilder initialization"""
        builder = ITPBuilder(hybrid_atoms=[], hybrid_params={"bonds": [], "angles": []})
        assert builder.hybrid_atoms == []
        assert builder.hybrid_params == {"bonds": [], "angles": []}


class TestAnalysisFramework:
    """Test analysis module framework"""

    def test_xvg_parser_initialization(self):
        """Test XVGParser initialization"""
        parser = XVGParser()
        assert parser.data is None
        assert parser.metadata == {}

    def test_festimator_initialization(self):
        """Test FEstimator initialization"""
        estimator = FEstimator()
        assert estimator.delta_g is None
        assert estimator.error is None


class TestInterfaceSignatures:
    """Test that required methods exist and have correct signatures"""

    def test_distance_mapper_has_map_method(self):
        """Test DistanceAtomMapper has map method"""
        mapper = DistanceAtomMapper()
        assert hasattr(mapper, "map")
        assert callable(mapper.map)

    def test_hybrid_topology_builder_has_build_method(self):
        """Test HybridTopologyBuilder has build method"""
        builder = HybridTopologyBuilder()
        assert hasattr(builder, "build")
        assert callable(builder.build)

    def test_itp_builder_has_write_itp_method(self):
        """Test ITPBuilder has write_itp method"""
        builder = ITPBuilder([], {})
        assert hasattr(builder, "write_itp")
        assert callable(builder.write_itp)

    def test_xvg_parser_has_parse_method(self):
        """Test XVGParser has parse method"""
        parser = XVGParser()
        assert hasattr(parser, "parse")
        assert callable(parser.parse)

    def test_festimator_has_bar_method(self):
        """Test FEstimator has bar method"""
        estimator = FEstimator()
        assert hasattr(estimator, "bar")
        assert callable(estimator.bar)

    def test_festimator_has_mbar_method(self):
        """Test FEstimator has mbar method"""
        estimator = FEstimator()
        assert hasattr(estimator, "mbar")
        assert callable(estimator.mbar)


class TestMethodNotImplemented:
    """Legacy scaffold tests updated to current concrete implementations."""

    def test_distance_mapper_map_handles_empty_inputs(self):
        """DistanceAtomMapper should accept empty inputs and return an empty mapping."""
        mapper = DistanceAtomMapper()
        mapping = mapper.map([], [])
        assert mapping.common == []
        assert mapping.transformed_a == []
        assert mapping.transformed_b == []
        assert mapping.surrounding_a == []
        assert mapping.surrounding_b == []

    def test_hybrid_topology_builder_build_requires_full_arguments(self):
        """HybridTopologyBuilder.build now requires mapping, params, and atom lists."""
        builder = HybridTopologyBuilder()
        with pytest.raises(TypeError):
            builder.build(AtomMapping([], [], [], [], []), {}, {})

    def test_itp_builder_write_itp_creates_file(self, tmp_path):
        """ITPBuilder.write_itp is implemented and should create an output file."""
        builder = ITPBuilder([], {})
        output = tmp_path / "test.itp"
        builder.write_itp(str(output))
        assert output.exists()

    def test_charge_reception_does_not_change_mapping_categories(self):
        """charge_reception should only affect redistribution, not atom classification."""
        common_a = Atom("C1", "C", np.array([0.0, 0.0, 0.0]), 0.1, "ca", 1)
        transformed_a = Atom("C2", "C", np.array([2.0, 0.0, 0.0]), -0.1, "ca", 2)

        common_b = Atom("C1", "C", np.array([0.0, 0.0, 0.0]), 0.1, "ca", 1)

        mapping_surround = DistanceAtomMapper(charge_reception="surround").map(
            [common_a, transformed_a],
            [common_b],
        )
        mapping_unique = DistanceAtomMapper(charge_reception="unique").map(
            [
                Atom("C1", "C", np.array([0.0, 0.0, 0.0]), 0.1, "ca", 1),
                Atom("C2", "C", np.array([2.0, 0.0, 0.0]), -0.1, "ca", 2),
            ],
            [Atom("C1", "C", np.array([0.0, 0.0, 0.0]), 0.1, "ca", 1)],
        )

        assert [atom.name for atom in mapping_surround.transformed_a] == ["C2"]
        assert [atom.name for atom in mapping_unique.transformed_a] == ["C2"]

    def test_xvg_parser_parse_not_implemented(self):
        """Test that XVGParser.parse raises NotImplementedError"""
        parser = XVGParser()
        with pytest.raises(NotImplementedError):
            parser.parse("test.xvg")

    def test_festimator_bar_not_implemented(self):
        """Test that FEstimator.bar raises NotImplementedError"""
        estimator = FEstimator()
        with pytest.raises(NotImplementedError):
            estimator.bar({})

    def test_festimator_mbar_not_implemented(self):
        """Test that FEstimator.mbar raises NotImplementedError"""
        estimator = FEstimator()
        with pytest.raises(NotImplementedError):
            estimator.mbar({})
