import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from mendelcell import list_tissues, run_mendelcell
from mendelcell.report import create_pdf_report, safe_filename


st.set_page_config(
    page_title="MendelCell",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 MendelCell")
st.subheader("Candidate gene prioritization by tissue-specific single-cell expression")

st.write(
    "MendelCell uses a built-in preprocessed Human Protein Atlas reference. "
    "Upload a candidate gene list, choose a tissue and expression threshold, "
    "then generate tables, plots, and a PDF report."
)


@st.cache_data
def load_reference_data():
    cluster_path = Path("data/mendelcell_clusters_reference.parquet")
    hpa_path = Path("data/mendelcell_celltype_reference.parquet")

    if not cluster_path.exists():
        raise FileNotFoundError(f"Missing reference file: {cluster_path}")

    if not hpa_path.exists():
        raise FileNotFoundError(f"Missing reference file: {hpa_path}")

    clusters = pd.read_parquet(cluster_path)
    hpa = pd.read_parquet(hpa_path)

    return clusters, hpa


@st.cache_data
def read_gene_list(file_name, file_bytes):
    import io

    buffer = io.BytesIO(file_bytes)

    if file_name.lower().endswith(".csv"):
        return pd.read_csv(buffer)

    return pd.read_csv(buffer, sep="\t")


def make_gene_count_plot(results):
    fig, ax = plt.subplots(figsize=(11, 6))

    plot_df = results.plot_df

    ax.bar(plot_df["Cell type"], plot_df["Gene count"])
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Number of candidate genes")
    ax.set_title(
        f"Candidate genes expressed in {results.selected_tissue}-specific cell types"
    )

    ax.tick_params(axis="x", rotation=45)

    for label in ax.get_xticklabels():
        label.set_ha("right")

    fig.tight_layout()
    return fig


def make_ncpm_plot(results, cell_type):
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
    return fig


try:
    clusters, hpa = load_reference_data()

except Exception as e:
    st.error("Could not load MendelCell reference files.")
    st.exception(e)
    st.stop()


valid_tissues = list_tissues(clusters)

st.sidebar.header("1. Upload gene list")

gene_file = st.sidebar.file_uploader(
    "Upload candidate gene list TSV or CSV",
    type=["tsv", "txt", "csv"]
)

st.sidebar.header("2. Choose settings")

selected_tissue = st.sidebar.selectbox(
    "Select tissue",
    valid_tissues
)

threshold = st.sidebar.number_input(
    "Expression threshold",
    min_value=0.0,
    value=1.0,
    step=0.5
)

run_button = st.sidebar.button("Run MendelCell analysis")


if gene_file is None:
    st.info("Upload a candidate gene list to begin.")

    st.markdown(
        """
        Your gene list should contain a column named:

        ```text
        Gene Symbol
        ```

        Example:

        ```text
        Gene Symbol
        INS
        GCG
        PDX1
        CD3D
        PTPRC
        ```
        """
    )

    st.stop()


try:
    gene_table = read_gene_list(gene_file.name, gene_file.getvalue())

except Exception as e:
    st.error(f"Could not read gene list file: {e}")
    st.exception(e)
    st.stop()


if "Gene Symbol" not in gene_table.columns:
    st.error("Gene list must contain a column named 'Gene Symbol'.")
    st.write("Columns found:")
    st.write(list(gene_table.columns))
    st.stop()


if not run_button:
    st.info("Choose a tissue and threshold, then click **Run MendelCell analysis**.")
    st.stop()


try:
    results = run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue=selected_tissue,
        threshold=threshold
    )

except Exception as e:
    st.error(f"MendelCell analysis failed: {e}")
    st.exception(e)
    st.stop()


st.success("Analysis complete.")

col1, col2, col3, col4 = st.columns(4)

col1.metric("Input genes", len(results.gene_symbols))
col2.metric("Tissue-specific cell types", len(results.unique_cells))
col3.metric("Genes passing threshold", results.filtered["Gene name"].nunique())
col4.metric(
    "Gene-cell pairs",
    len(results.filtered[["Cell type", "Gene name"]].drop_duplicates())
)

st.header("Cell types unique to selected tissue")
st.dataframe(results.unique_to_tissue, use_container_width=True)

st.header("Candidate gene count per cell type")
st.dataframe(results.cell_count_df, use_container_width=True)

if not results.plot_df.empty:
    fig = make_gene_count_plot(results)
    st.pyplot(fig)
    plt.close(fig)

st.header("Candidate genes found in each cell type")
st.dataframe(results.genes_in_cell_df, use_container_width=True)

st.header("Filtered candidate genes")
st.dataframe(results.filtered_report, use_container_width=True)

st.header("Mean nCPM values")
st.dataframe(results.ncpm_df, use_container_width=True)

st.header("nCPM plots by cell type")

for cell_type in results.ncpm_df["Cell type"].unique():
    fig = make_ncpm_plot(results, cell_type)
    st.pyplot(fig)
    plt.close(fig)


st.header("Download outputs")

unique_tsv = results.unique_to_tissue.to_csv(sep="\t", index=False)
filtered_tsv = results.filtered.to_csv(sep="\t", index=False)
ncpm_tsv = results.ncpm_df.to_csv(sep="\t", index=False)

safe_tissue = safe_filename(results.selected_tissue)
pdf_filename = f"MendelCell_report_{safe_tissue}.pdf"

with tempfile.TemporaryDirectory() as tmpdir:
    pdf_path = Path(tmpdir) / pdf_filename
    create_pdf_report(results, pdf_path)

    with open(pdf_path, "rb") as f:
        pdf_bytes = f.read()

st.download_button(
    label="Download PDF report",
    data=pdf_bytes,
    file_name=pdf_filename,
    mime="application/pdf"
)

st.download_button(
    label="Download unique cell types TSV",
    data=unique_tsv,
    file_name="unique_cell_types.tsv",
    mime="text/tab-separated-values"
)

st.download_button(
    label="Download filtered candidate genes TSV",
    data=filtered_tsv,
    file_name="filtered_candidate_genes.tsv",
    mime="text/tab-separated-values"
)

st.download_button(
    label="Download nCPM table TSV",
    data=ncpm_tsv,
    file_name="candidate_gene_ncpm_by_cell_type.tsv",
    mime="text/tab-separated-values"
)