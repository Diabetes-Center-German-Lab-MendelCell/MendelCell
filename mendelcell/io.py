"""Input/output helpers for MendelCell."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


REQUIRED_CLUSTER_COLS = ["Tissue", "Cell type", "Gene name", "nCPM"]
REQUIRED_HPA_COLS = ["Gene name", "Cell type"]
REQUIRED_GENE_COLS = ["Gene Symbol"]


def read_tsv(path_or_buffer) -> pd.DataFrame:
    """Read a TSV file. Zipped TSV files are supported when the path ends in .zip."""
    compression = "zip" if str(path_or_buffer).lower().endswith(".zip") else "infer"
    return pd.read_csv(path_or_buffer, sep="\t", compression=compression)


def validate_columns(df: pd.DataFrame, required_cols: Iterable[str], file_label: str) -> None:
    """Raise a clear error if a dataframe is missing required columns."""
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing column(s) in {file_label}: {', '.join(missing)}")


def choose_expression_column(hpa: pd.DataFrame) -> str:
    """Choose the expression column from the HPA cell-type file."""
    if "nTPM" in hpa.columns:
        return "nTPM"
    if "nCPM" in hpa.columns:
        return "nCPM"
    raise ValueError("Could not find expression column. Expected 'nTPM' or 'nCPM'.")


def clean_input_tables(
    clusters: pd.DataFrame,
    hpa: pd.DataFrame,
    gene_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str], str]:
    """Validate and clean the three MendelCell input tables.

    Returns
    -------
    clusters_clean, hpa_clean, gene_symbols, expression_col
    """
    validate_columns(clusters, REQUIRED_CLUSTER_COLS, "cluster file")
    validate_columns(hpa, REQUIRED_HPA_COLS, "HPA cell-type file")
    validate_columns(gene_table, REQUIRED_GENE_COLS, "gene list file")

    expression_col = choose_expression_column(hpa)

    clusters = clusters.copy()
    hpa = hpa.copy()
    gene_table = gene_table.copy()

    clusters["Tissue"] = clusters["Tissue"].astype(str).str.strip()
    clusters["Cell type"] = clusters["Cell type"].astype(str).str.strip()
    clusters["Gene name"] = clusters["Gene name"].astype(str).str.strip()
    clusters["Gene name clean"] = clusters["Gene name"].str.upper()
    clusters["nCPM"] = pd.to_numeric(clusters["nCPM"], errors="coerce")

    hpa["Cell type"] = hpa["Cell type"].astype(str).str.strip()
    hpa["Gene name"] = hpa["Gene name"].astype(str).str.strip()
    hpa["Gene name clean"] = hpa["Gene name"].str.upper()
    hpa[expression_col] = pd.to_numeric(hpa[expression_col], errors="coerce")

    gene_symbols = (
        gene_table["Gene Symbol"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
    )
    gene_symbols = set(gene_symbols)

    return clusters, hpa, gene_symbols, expression_col


def load_input_files(
    cluster_file: str | Path,
    hpa_file: str | Path,
    gene_file: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str], str]:
    """Load and clean MendelCell input files."""
    clusters = read_tsv(cluster_file)
    hpa = read_tsv(hpa_file)
    gene_table = read_tsv(gene_file)
    return clean_input_tables(clusters, hpa, gene_table)
