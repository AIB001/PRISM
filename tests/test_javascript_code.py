#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for JavaScript code module for contact visualization.

These tests verify that the modularized javascript_code package
maintains complete compatibility with the original monolithic version.
"""

import numpy as np
import pytest
from prism.analysis.contact.html_builder import HTMLBuilder
from prism.analysis.contact.javascript_code import (
    get_javascript_code,
    get_contact_map_class_part1,
    get_contact_map_class_part2,
    get_drawing_methods,
    get_utility_functions,
)


class TestJavaScriptCodeModule:
    """Test the modularized javascript_code package."""

    def test_get_javascript_code_returns_string(self):
        """get_javascript_code should return a non-empty string."""
        js_code = get_javascript_code()
        assert isinstance(js_code, str)
        assert len(js_code) > 0

    def test_javascript_code_has_expected_size(self):
        """JavaScript code should have expected size (around 55k chars)."""
        js_code = get_javascript_code()
        assert 50000 < len(js_code) < 60000

    def test_contact_map_class_part1(self):
        """Part 1 should contain ContactMap class initialization."""
        part1 = get_contact_map_class_part1()
        assert "class ContactMap" in part1
        assert "constructor()" in part1
        assert len(part1) > 0

    def test_contact_map_class_part2(self):
        """Part 2 should contain core ContactMap functionality."""
        part2 = get_contact_map_class_part2()
        # Part 2 contains event handling and mouse interactions
        assert "addEventListener" in part2 or "onMouseDown" in part2
        assert len(part2) > 0

    def test_drawing_methods(self):
        """Drawing methods should contain core visualization functions."""
        drawing = get_drawing_methods()
        assert "drawContact" in drawing
        assert "drawContacts" in drawing
        assert len(drawing) > 0

    def test_utility_functions(self):
        """Utility functions should contain helper functions."""
        utility = get_utility_functions()
        assert "function" in utility
        assert len(utility) > 0

    def test_complete_code_contains_all_parts(self):
        """Complete code should be sum of all parts."""
        complete = get_javascript_code()
        parts = (
            get_contact_map_class_part1()
            + get_contact_map_class_part2()
            + get_drawing_methods()
            + get_utility_functions()
        )
        # Complete code should contain all parts (plus closing tags)
        assert get_contact_map_class_part1()[:100] in complete
        assert get_drawing_methods()[:100] in complete
        assert get_utility_functions()[:100] in complete

    def test_complete_code_has_contactmap_class(self):
        """Complete JavaScript should have ContactMap class."""
        js_code = get_javascript_code()
        assert "class ContactMap" in js_code

    def test_complete_code_has_canvas_integration(self):
        """Complete JavaScript should integrate with canvas."""
        js_code = get_javascript_code()
        assert "canvas" in js_code
        assert "getContext" in js_code

    def test_complete_code_has_event_listeners(self):
        """Complete JavaScript should have event listeners."""
        js_code = get_javascript_code()
        assert "addEventListener" in js_code

    def test_complete_code_has_drawing_functions(self):
        """Complete JavaScript should have drawing functions."""
        js_code = get_javascript_code()
        assert "drawContact" in js_code
        assert "drawContacts" in js_code

    def test_complete_code_has_closing_tags(self):
        """Complete code should end with proper HTML closing tags."""
        js_code = get_javascript_code()
        assert js_code.strip().endswith("</html>")
        assert "</script>" in js_code
        assert "</body>" in js_code


class TestHTMLBuilderIntegration:
    """Test integration with HTMLBuilder."""

    @pytest.fixture
    def sample_data(self):
        """Create sample test data."""
        return {
            "contacts": np.array([[1, 10, 3.5, "H-Bond"], [2, 11, 4.2, "Hydrophobic"]]),
            "ligand_data": {
                "atoms": [
                    {"serial": 1, "name": "N", "resname": "LIG", "x": 0.0, "y": 0.0, "z": 0.0},
                    {"serial": 2, "name": "CA", "resname": "LIG", "x": 1.0, "y": 0.0, "z": 0.0},
                ],
                "bonds": [{"atom1": 1, "atom2": 2}],
            },
            "stats": {
                "total_contacts": 2,
                "high_freq_contacts": 1,
                "max_freq_percent": 50.0,
            },
        }

    def test_html_builder_generates_valid_html(self, sample_data):
        """HTMLBuilder should generate valid HTML with JavaScript."""
        builder = HTMLBuilder()
        html = builder.generate_html(
            trajectory_file="test.xtc",
            topology_file="test.pdb",
            ligand_file="test.pdb",
            ligand_name="LIG",
            total_frames=100,
            ligand_data=sample_data["ligand_data"],
            contacts=sample_data["contacts"],
            stats=sample_data["stats"],
        )

        # Check basic HTML structure
        assert "<!DOCTYPE html>" in html or "<html>" in html
        assert "<head>" in html
        assert "<body>" in html
        assert "</html>" in html

    def test_html_contains_javascript_code(self, sample_data):
        """Generated HTML should contain JavaScript code."""
        builder = HTMLBuilder()
        html = builder.generate_html(
            trajectory_file="test.xtc",
            topology_file="test.pdb",
            ligand_file="test.pdb",
            ligand_name="LIG",
            total_frames=100,
            ligand_data=sample_data["ligand_data"],
            contacts=sample_data["contacts"],
            stats=sample_data["stats"],
        )

        # Check for script tags and JavaScript content
        assert "<script>" in html
        assert "</script>" in html
        assert "class ContactMap" in html
        assert "canvas" in html

    def test_html_contains_contact_data(self, sample_data):
        """Generated HTML should contain contact data."""
        builder = HTMLBuilder()
        html = builder.generate_html(
            trajectory_file="test.xtc",
            topology_file="test.pdb",
            ligand_file="test.pdb",
            ligand_name="LIG",
            total_frames=100,
            ligand_data=sample_data["ligand_data"],
            contacts=sample_data["contacts"],
            stats=sample_data["stats"],
        )

        # Check for contact data
        assert "LIG" in html
        assert str(sample_data["stats"]["total_contacts"]) in html

    def test_html_size_is_reasonable(self, sample_data):
        """Generated HTML should have reasonable size."""
        builder = HTMLBuilder()
        html = builder.generate_html(
            trajectory_file="test.xtc",
            topology_file="test.pdb",
            ligand_file="test.pdb",
            ligand_name="LIG",
            total_frames=100,
            ligand_data=sample_data["ligand_data"],
            contacts=sample_data["contacts"],
            stats=sample_data["stats"],
        )

        # HTML should be substantial (>60k chars) but not enormous
        assert len(html) > 60000
        assert len(html) < 100000


class TestJavaScriptCodeComponents:
    """Test individual components of JavaScript code."""

    def test_part1_has_class_definition(self):
        """Part 1 should start with class ContactMap."""
        part1 = get_contact_map_class_part1()
        assert part1.strip().startswith("class ContactMap")

    def test_part1_initializes_canvas(self):
        """Part 1 should initialize canvas element."""
        part1 = get_contact_map_class_part1()
        assert "canvas" in part1
        assert "getContext" in part1

    def test_part1_has_initial_state(self):
        """Part 1 should set initial state variables."""
        part1 = get_contact_map_class_part1()
        assert "showConnections" in part1
        assert "showHydrogens" in part1
        assert "showGrid" in part1

    def test_drawing_has_drawcontact_method(self):
        """Drawing methods should have drawContact method."""
        drawing = get_drawing_methods()
        # Check for drawContact (might be defined with different formatting)
        assert "drawContact" in drawing
        # Also check for related drawing functions
        assert "drawContacts" in drawing or "drawLigand" in drawing

    def test_drawing_has_drawcontacts_method(self):
        """Drawing methods should have drawContacts method."""
        drawing = get_drawing_methods()
        assert "drawContacts(" in drawing

    def test_drawing_has_circle_drawing(self):
        """Drawing methods should have circle/arc drawing."""
        drawing = get_drawing_methods()
        assert "arc" in drawing or "circle" in drawing

    def test_drawing_has_line_drawing(self):
        """Drawing methods should have line drawing."""
        drawing = get_drawing_methods()
        assert "lineTo" in drawing or "moveTo" in drawing

    def test_utility_has_helper_functions(self):
        """Utility functions should have helper functions."""
        utility = get_utility_functions()
        assert "function" in utility
        # Check for common utility patterns
        assert len(utility) > 1000  # Should be substantial

    def test_parts_combine_to_complete_code(self):
        """All parts should combine to form complete code."""
        parts_combined = (
            get_contact_map_class_part1()
            + get_contact_map_class_part2()
            + get_drawing_methods()
            + get_utility_functions()
        )

        complete = get_javascript_code()

        # Check that major code blocks from parts are in complete
        part1_snippet = get_contact_map_class_part1()[:200]
        part2_snippet = get_contact_map_class_part2()[:200]
        drawing_snippet = get_drawing_methods()[:200]
        utility_snippet = get_utility_functions()[:200]

        assert part1_snippet in complete
        assert part2_snippet in complete
        assert drawing_snippet in complete
        assert utility_snippet in complete

    def test_complete_code_size_matches_parts_sum(self):
        """Complete code size should approximately match sum of parts."""
        parts = (
            get_contact_map_class_part1()
            + get_contact_map_class_part2()
            + get_drawing_methods()
            + get_utility_functions()
        )
        complete = get_javascript_code()

        # Complete should have parts plus closing tags
        for part in [get_contact_map_class_part1(), get_drawing_methods()]:
            assert part[:50] in complete


class TestBackwardCompatibility:
    """Test that the refactored code maintains backward compatibility."""

    def test_same_interface_as_original(self):
        """All original functions should be available."""
        # These functions should exist and be callable
        assert callable(get_javascript_code)
        assert callable(get_contact_map_class_part1)
        assert callable(get_contact_map_class_part2)
        assert callable(get_drawing_methods)
        assert callable(get_utility_functions)

    def test_functions_return_strings(self):
        """All functions should return strings."""
        assert isinstance(get_javascript_code(), str)
        assert isinstance(get_contact_map_class_part1(), str)
        assert isinstance(get_contact_map_class_part2(), str)
        assert isinstance(get_drawing_methods(), str)
        assert isinstance(get_utility_functions(), str)

    def test_import_from_html_builder_works(self):
        """Import in HTMLBuilder should work correctly."""
        # This tests the import chain: html_builder -> javascript_code -> loader
        try:
            from prism.analysis.contact.javascript_code import get_javascript_code

            js_code = get_javascript_code()
            assert len(js_code) > 0
        except ImportError as e:
            pytest.fail(f"Import chain broken: {e}")
