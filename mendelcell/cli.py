"""Command-line interface for MendelCell."""

from __future__ import annotations

import argparse
from pathlib import Path

from .analysis import run_mendelcell_from_files
from .report import create_pdf_report, safe_filename, write_tsv_outputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run MendelCell candidate gene expression prioritization."
    )
    parser.add_argument(
        "--clusters",
        required=True,
        help="Path to HPA rna_single_cell_cluster.tsv or .zip file.",
    )
    parser.add_argument(
        "--hpa",
        required=True,
        help="Path to HPA rna_single_cell_type.tsv file.",
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Path to candidate gene TSV file containing a 'Gene Symbol' column.",
    )
    parser.add_argument(
        "--tissue",
        required=True,
        help="Tissue to analyze, such as Pancreas.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.0,
        help="Expression threshold applied to nTPM or nCPM. Default: 1.0.",
    )
    parser.add_argument(
        "--outdir",
        default="mendelcell_output",
        help="Directory where TSV and PDF outputs will be written.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    results = run_mendelcell_from_files(
        cluster_file=args.clusters,
        hpa_file=args.hpa,
        gene_file=args.genes,
        tissue=args.tissue,
        threshold=args.threshold,
    )

    tsv_paths = write_tsv_outputs(results, outdir)
    pdf_path = create_pdf_report(
        results,
        outdir / f"MendelCell_report_{safe_filename(results.selected_tissue)}.pdf",
    )

    print("MendelCell analysis complete.")
    print(f"Selected tissue: {results.selected_tissue}")
    print(f"Expression column: {results.expression_col}")
    print(f"Threshold: {results.threshold}")
    print(f"Unique cell types: {len(results.unique_cells)}")
    print(f"Genes passing threshold: {results.filtered['Gene name'].nunique()}")
    print("\nOutput files:")
    for path in tsv_paths.values():
        print(f"- {path}")
    print(f"- {pdf_path}")


if __name__ == "__main__":
    main()
