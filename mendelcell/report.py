"""Report and output helpers for MendelCell."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from .analysis import MendelCellResults


def safe_filename(text: str) -> str:
    """Create a filesystem-safe string from a tissue or sample name."""
    return (
        str(text)
        .replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
    )


def add_text_page(pdf: PdfPages, title: str, lines: list[str]) -> None:
    """Add a text-only page to a PdfPages report."""
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis("off")

    ax.text(0.05, 0.95, title, fontsize=18, fontweight="bold", va="top")
    y = 0.88

    for line in lines:
        ax.text(0.05, y, line, fontsize=11, va="top")
        y -= 0.04

        if y < 0.05:
            pdf.savefig(fig)
            plt.close(fig)
            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis("off")
            y = 0.95

    pdf.savefig(fig)
    plt.close(fig)


def add_dataframe_pages(
    pdf: PdfPages,
    df: pd.DataFrame,
    title: str,
    rows_per_page: int = 25,
) -> None:
    """Add one or more dataframe table pages to a PdfPages report."""
    if df.empty:
        add_text_page(pdf, title, ["No data available."])
        return

    df_to_show = df.copy()

    for col in df_to_show.columns:
        df_to_show[col] = df_to_show[col].apply(
            lambda x: ", ".join(x) if isinstance(x, list) else x
        )

    df_to_show = df_to_show.astype(str)

    for col in df_to_show.columns:
        df_to_show[col] = df_to_show[col].str.slice(0, 80)

    for start in range(0, len(df_to_show), rows_per_page):
        end = start + rows_per_page
        page_df = df_to_show.iloc[start:end]

        fig, ax = plt.subplots(figsize=(11, 8.5))
        ax.axis("off")
        ax.set_title(
            f"{title} rows {start + 1}-{min(end, len(df_to_show))}",
            fontsize=14,
            fontweight="bold",
            pad=20,
        )

        table = ax.table(
            cellText=page_df.values,
            colLabels=page_df.columns,
            loc="center",
            cellLoc="left",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1, 1.3)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)


def write_tsv_outputs(results: MendelCellResults, output_dir: str | Path) -> dict[str, Path]:
    """Write TSV outputs and return their paths."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "unique_cell_types": output_dir / "unique_cell_types.tsv",
        "filtered_candidate_genes": output_dir / "filtered_candidate_genes.tsv",
        "candidate_gene_ncpm_by_cell_type": output_dir / "candidate_gene_ncpm_by_cell_type.tsv",
    }

    results.unique_to_tissue.to_csv(paths["unique_cell_types"], sep="\t", index=False)
    results.filtered.to_csv(paths["filtered_candidate_genes"], sep="\t", index=False)
    results.ncpm_df.to_csv(paths["candidate_gene_ncpm_by_cell_type"], sep="\t", index=False)
    return paths


def create_pdf_report(
    results: MendelCellResults,
    output_pdf: str | Path | None = None,
) -> Path:
    """Create a MendelCell PDF report and return its path."""
    if output_pdf is None:
        output_pdf = f"MendelCell_report_{safe_filename(results.selected_tissue)}.pdf"

    output_pdf = Path(output_pdf)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(output_pdf) as pdf:
        total_input_genes = len(results.gene_symbols)
        total_unique_cells = len(results.unique_cells)
        total_filtered_genes = results.filtered["Gene name"].nunique()
        total_gene_cell_pairs = len(
            results.filtered[["Cell type", "Gene name"]].drop_duplicates()
        )

        summary_lines = [
            f"Selected tissue: {results.selected_tissue}",
            f"Expression column used for filtering: {results.expression_col}",
            f"Expression threshold: {results.threshold}",
            "",
            f"Number of input candidate genes: {total_input_genes}",
            f"Number of tissue-specific cell types found: {total_unique_cells}",
            f"Number of candidate genes passing threshold: {total_filtered_genes}",
            f"Number of gene-cell type expression pairs: {total_gene_cell_pairs}",
        ]

        add_text_page(pdf, "MendelCell Candidate Gene Expression Report", summary_lines)

        add_dataframe_pages(
            pdf,
            results.unique_to_tissue,
            f"Cell types unique to {results.selected_tissue}",
            rows_per_page=20,
        )
        add_dataframe_pages(
            pdf,
            results.cell_count_df,
            "Number of candidate genes per cell type",
            rows_per_page=30,
        )
        add_dataframe_pages(
            pdf,
            results.genes_in_cell_df,
            "Candidate genes found in each cell type",
            rows_per_page=15,
        )
        add_dataframe_pages(
            pdf,
            results.filtered_report,
            "Filtered candidate genes passing expression threshold",
            rows_per_page=30,
        )
        add_dataframe_pages(
            pdf,
            results.ncpm_df,
            "Mean nCPM values for outputted genes",
            rows_per_page=30,
        )

        if not results.plot_df.empty:
            fig, ax = plt.subplots(figsize=(11, 6))
            ax.bar(results.plot_df["Cell type"], results.plot_df["Gene count"])
            ax.set_xlabel("Cell type")
            ax.set_ylabel("Number of candidate genes")
            ax.set_title(
                f"Candidate genes expressed in {results.selected_tissue}-specific cell types"
            )
            ax.tick_params(axis="x", rotation=45)
            for label in ax.get_xticklabels():
                label.set_ha("right")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        if not results.ncpm_df.empty:
            for cell_type in results.ncpm_df["Cell type"].unique():
                cell_df = results.ncpm_df[results.ncpm_df["Cell type"] == cell_type]
                cell_df = cell_df.sort_values("nCPM", ascending=False)

                fig, ax = plt.subplots(figsize=(11, 6))
                ax.bar(cell_df["Gene name"], cell_df["nCPM"])
                ax.set_xlabel("Gene")
                ax.set_ylabel("Mean nCPM")
                ax.set_title(f"Mean nCPM of candidate genes in {cell_type}")
                ax.tick_params(axis="x", rotation=45)
                for label in ax.get_xticklabels():
                    label.set_ha("right")
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

    return output_pdf
