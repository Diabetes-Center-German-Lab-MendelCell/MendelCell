"""Minimal smoke test showing the package can be imported.

Real biological tests should be added with small example input tables.
"""

from mendelcell import MendelCellResults, run_mendelcell


def test_imports():
    assert MendelCellResults is not None
    assert run_mendelcell is not None
