from pathlib import Path
import argparse
import pandas as pd


def read_tsv_or_zip(path):
    path = Path(path)
    name = path.name.lower()

    if name.endswith(".zip"):
        return pd.read_csv(path, sep="\t", compression="zip", low_memory=False)

    if name.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)

    return pd.read_csv(path, sep="\t", low_memory=False)


def check_columns(df, required_cols, file_label):
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"{file_label} is missing columns: {missing}")


def build_reference(cluster_file, hpa_file, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Reading cluster file...")
    clusters = read_tsv_or_zip(cluster_file)

    print("Reading HPA cell-type file...")
    hpa = read_tsv_or_zip(hpa_file)

    required_cluster_cols = ["Tissue", "Cell type", "Gene name", "nCPM"]
    required_hpa_cols = ["Gene name", "Cell type"]

    check_columns(clusters, required_cluster_cols, "Cluster file")
    check_columns(hpa, required_hpa_cols, "HPA cell-type file")

    if "nTPM" in hpa.columns:
        expression_col = "nTPM"
    elif "nCPM" in hpa.columns:
        expression_col = "nCPM"
    else:
        raise ValueError("HPA file must contain either 'nTPM' or 'nCPM'.")

    print(f"Using expression column: {expression_col}")

    print("Cleaning cluster reference...")
    clusters_ref = clusters[required_cluster_cols].copy()

    clusters_ref["Tissue"] = clusters_ref["Tissue"].astype(str).str.strip()
    clusters_ref["Cell type"] = clusters_ref["Cell type"].astype(str).str.strip()
    clusters_ref["Gene name"] = clusters_ref["Gene name"].astype(str).str.strip()
    clusters_ref["nCPM"] = pd.to_numeric(clusters_ref["nCPM"], errors="coerce")

    clusters_ref = clusters_ref.dropna(subset=["Tissue", "Cell type", "Gene name", "nCPM"])
    clusters_ref = clusters_ref.drop_duplicates()

    print("Cleaning HPA cell-type reference...")
    hpa_ref = hpa[["Gene name", "Cell type", expression_col]].copy()

    hpa_ref["Gene name"] = hpa_ref["Gene name"].astype(str).str.strip()
    hpa_ref["Cell type"] = hpa_ref["Cell type"].astype(str).str.strip()
    hpa_ref[expression_col] = pd.to_numeric(hpa_ref[expression_col], errors="coerce")

    hpa_ref = hpa_ref.dropna(subset=["Gene name", "Cell type", expression_col])
    hpa_ref = hpa_ref.drop_duplicates()

    cluster_out = output_dir / "mendelcell_clusters_reference.parquet"
    hpa_out = output_dir / "mendelcell_celltype_reference.parquet"

    print("Saving Parquet reference files...")
    clusters_ref.to_parquet(cluster_out, index=False)
    hpa_ref.to_parquet(hpa_out, index=False)

    print("\nDone.")
    print(f"Cluster reference: {cluster_out}")
    print(f"HPA cell-type reference: {hpa_out}")
    print(f"Cluster rows: {len(clusters_ref):,}")
    print(f"HPA rows: {len(hpa_ref):,}")
    print(f"Cluster file size MB: {cluster_out.stat().st_size / 1_000_000:.2f}")
    print(f"HPA file size MB: {hpa_out.stat().st_size / 1_000_000:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description="Build MendelCell preprocessed HPA Parquet reference files."
    )

    parser.add_argument(
        "--clusters",
        required=True,
        help="Path to rna_single_cell_cluster.tsv or .zip"
    )

    parser.add_argument(
        "--hpa",
        required=True,
        help="Path to rna_single_cell_type.tsv or .zip"
    )

    parser.add_argument(
        "--outdir",
        default="data",
        help="Output directory for Parquet reference files"
    )

    args = parser.parse_args()

    build_reference(
        cluster_file=args.clusters,
        hpa_file=args.hpa,
        output_dir=args.outdir
    )


if __name__ == "__main__":
    main()