"""MendelCell: prioritize candidate genes by tissue-specific single-cell expression."""

from .analysis import MendelCellResults, list_tissues, run_mendelcell, run_mendelcell_from_files

__all__ = [
    "MendelCellResults",
    "list_tissues",
    "run_mendelcell",
    "run_mendelcell_from_files",
]

__version__ = "0.1.0"
