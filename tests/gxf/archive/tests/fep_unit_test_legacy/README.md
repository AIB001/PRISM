# FEP Legacy Test Scripts

This directory stores archived FEP validation scripts that were previously
mixed into `tests/gxf/FEP/unit_test`.

These files were moved here because they are not part of the active pytest
suite:

- some are one-off HTML/ITP generators
- some are manual verification scripts
- some use `test_*.py` names but are not structured as stable pytest tests

Active pytest-backed FEP tests remain in `tests/gxf/FEP/unit_test`.
