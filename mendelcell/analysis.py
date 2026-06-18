"""Core MendelCell analysis functions.

This module contains the reusable scientific core. It does not ask for user input,
write files, or generate figures. That makes it easier to test and reuse.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .io import clean_input_tables, load_input_files


@dataclass
class MendelCellResults:
    """Container for MendelCell analysis results."""

    selected_tissue: str
    threshold: float
    expression_col: str
    gene_symbols: set[str]
    unique_cells: list[str]
    unique_to_tissue: pd.DataFrame
    filtered: pd.DataFrame
    genes_in_cell: dict[str, list[str]]
    cell_dict: dict[str, int]
    ncpm_df: pd.DataFrame

    @property
    def cell_count_df(self) -> pd.DataFrame:
        """Return candidate-gene counts per cell type."""
        return pd.DataFrame(
            {
                "Cell type": list(self.cell_dict.keys()),
                "Candidate gene count": list(self.cell_dict.values()),
            }
        )

    @property
    def genes_in_cell_df(self) -> pd.DataFrame:
        """Return a readable table of genes found in each cell type."""
        return pd.DataFrame(
            [
                {"Cell type": cell_type, "Genes": ", ".join(genes)}
                for cell_type, genes in self.genes_in_cell.items()
            ]
        )

    @property
    def filtered_report(self) -> pd.DataFrame:
        """Return a compact filtered-gene table for reports."""
        cols = [
            col
            for col in ["Gene name", "Cell type", self.expression_col]
            if col in self.filtered.columns
        ]
        if not cols:
            return pd.DataFrame()
        return (
            self.filtered[cols]
            .drop_duplicates()
            .sort_values(["Cell type", self.expression_col], ascending=[True, False])
        )

    @property
    def plot_df(self) -> pd.DataFrame:
        """Return the gene-count table used for plotting."""
        plot_df = pd.DataFrame(
            {
                "Cell type": list(self.cell_dict.keys()),
                "Gene count": list(self.cell_dict.values()),
            }
        )
        return plot_df[plot_df["Gene count"] > 0]


def list_tissues(clusters: pd.DataFrame) -> list[str]:
    """List available tissues from a cleaned or raw cluster dataframe."""
    if "Tissue" not in clusters.columns:
        raise ValueError("Cluster dataframe must contain a 'Tissue' column.")
    return sorted(clusters["Tissue"].dropna().astype(str).str.strip().unique())


def resolve_tissue_name(tissue: str, valid_tissues: list[str]) -> str:
    """Resolve case-insensitive tissue input to the original HPA tissue name."""
    tissue_lookup = {valid_tissue.lower(): valid_tissue for valid_tissue in valid_tissues}
    tissue_input = tissue.strip().lower()
    if tissue_input not in tissue_lookup:
        choices = ", ".join(valid_tissues)
        raise ValueError(f"Tissue '{tissue}' was not found. Available tissues: {choices}")
    return tissue_lookup[tissue_input]


def run_mendelcell(
    clusters: pd.DataFrame,
    hpa: pd.DataFrame,
    gene_table: pd.DataFrame,
    tissue: str,
    threshold: float = 1.0,
) -> MendelCellResults:
    """Run MendelCell on already-loaded dataframes.

    Parameters
    ----------
    clusters:
        HPA cluster-level file, usually ``rna_single_cell_cluster.tsv``.
    hpa:
        HPA cell-type-level file, usually ``rna_single_cell_type.tsv``.
    gene_table:
        Candidate gene table containing a ``Gene Symbol`` column.
    tissue:
        Tissue to analyze. Matching is case-insensitive.
    threshold:
        Expression threshold applied to the HPA cell-type expression column.
    """
    threshold = float(threshold)
    clusters, hpa, gene_symbols, expression_col = clean_input_tables(clusters, hpa, gene_table)

    selected_tissue = resolve_tissue_name(tissue, list_tissues(clusters))

    celltype_tissues = (
        clusters.groupby("Cell type")["Tissue"]
        .apply(lambda x: sorted(set(x.dropna())))
        .reset_index(name="Tissues")
    )

    unique_to_tissue = celltype_tissues[
        celltype_tissues["Tissues"].apply(lambda tissues: tissues == [selected_tissue])
    ].copy()

    unique_cells = unique_to_tissue["Cell type"].tolist()
    if not unique_cells:
        raise ValueError(f"No cell types were unique to tissue '{selected_tissue}'.")

    filtered = hpa[
        (hpa["Cell type"].isin(unique_cells))
        & (hpa["Gene name clean"].isin(gene_symbols))
        & (hpa[expression_col] >= threshold)
    ].copy()

    if filtered.empty:
        raise ValueError("No candidate genes passed the expression threshold.")

    genes_in_cell = (
        filtered.groupby("Cell type")["Gene name"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    cell_dict = {cell: len(genes_in_cell.get(cell, [])) for cell in unique_cells}
    cell_dict = dict(sorted(cell_dict.items(), key=lambda item: item[1], reverse=True))

    expressed_gene_pairs = filtered[["Cell type", "Gene name clean"]].drop_duplicates()

    ncpm_df = clusters[clusters["Tissue"] == selected_tissue].merge(
        expressed_gene_pairs,
        on=["Cell type", "Gene name clean"],
        how="inner",
    )

    ncpm_df = (
        ncpm_df.groupby(["Cell type", "Gene name"], as_index=False)["nCPM"]
        .mean()
        .sort_values(["Cell type", "nCPM"], ascending=[True, False])
    )

    return MendelCellResults(
        selected_tissue=selected_tissue,
        threshold=threshold,
        expression_col=expression_col,
        gene_symbols=gene_symbols,
        unique_cells=unique_cells,
        unique_to_tissue=unique_to_tissue,
        filtered=filtered,
        genes_in_cell=genes_in_cell,
        cell_dict=cell_dict,
        ncpm_df=ncpm_df,
    )


def run_mendelcell_from_files(
    cluster_file: str | Path,
    hpa_file: str | Path,
    gene_file: str | Path,
    tissue: str,
    threshold: float = 1.0,
) -> MendelCellResults:
    """Run MendelCell from file paths."""
    clusters, hpa, gene_symbols, expression_col = load_input_files(
        cluster_file=cluster_file,
        hpa_file=hpa_file,
        gene_file=gene_file,
    )

    # Rebuild a minimal gene table so run_mendelcell can reuse one validation path.
    gene_table = pd.DataFrame({"Gene Symbol": sorted(gene_symbols)})
    return run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue=tissue,
        threshold=threshold,
    )
