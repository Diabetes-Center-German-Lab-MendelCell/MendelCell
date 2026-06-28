from pathlib import Path
import io
import tempfile

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from mendelcell import list_tissues, run_mendelcell
from mendelcell.report import create_pdf_report, safe_filename


# -----------------------------
# File paths
# -----------------------------

ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "data"

CLUSTER_REFERENCE = DATA_DIR / "mendelcell_clusters_reference.parquet"
CELLTYPE_REFERENCE = DATA_DIR / "mendelcell_celltype_reference.parquet"


# -----------------------------
# Page setup
# -----------------------------

st.set_page_config(
    page_title="MendelCell",
    page_icon="🧬",
    layout="wide",
)

# Reduce left and right page margins.
# This makes the app area wider without using the huge default side padding.
st.markdown(
    """
    <style>
    .block-container {
        max-width: 1200px;
        padding-left: 2rem;
        padding-right: 2rem;
        padding-top: 2rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 MendelCell")
st.subheader("Candidate gene prioritization by tissue-specific single-cell expression")

st.write(
    "MendelCell uses a preprocessed Human Protein Atlas single-cell reference. "
    "Upload a candidate gene list, enter a tissue, choose an expression threshold, "
    "and generate ranked tables, plots, and a PDF report."
)


# -----------------------------
# Helper functions
# -----------------------------

def show_dataframe_with_1_index(df, height=400, width=1050):
    """
    Display a dataframe in Streamlit with row numbering starting at 1.

    Uses a fixed pixel width to reduce layout shifting.
    """
    display_df = df.copy()
    display_df.index = range(1, len(display_df) + 1)

    st.dataframe(
        display_df,
        width=width,
        height=height,
    )


def load_reference_data():
    """
    Load preprocessed HPA reference files.

    The reference files are loaded only after the user uploads a gene list
    and clicks Run.
    """
    if not CLUSTER_REFERENCE.exists():
        raise FileNotFoundError(
            f"Missing cluster reference file: {CLUSTER_REFERENCE}"
        )

    if not CELLTYPE_REFERENCE.exists():
        raise FileNotFoundError(
            f"Missing cell-type reference file: {CELLTYPE_REFERENCE}"
        )

    clusters = pd.read_parquet(CLUSTER_REFERENCE)
    hpa = pd.read_parquet(CELLTYPE_REFERENCE)

    return clusters, hpa


def read_gene_list(file_name, file_bytes):
    """
    Read uploaded candidate gene list.

    Supported:
    - .tsv
    - .txt
    - .csv

    Required column:
    - Gene Symbol
    """
    buffer = io.BytesIO(file_bytes)
    file_name_lower = file_name.lower()

    if file_name_lower.endswith(".csv"):
        return pd.read_csv(buffer)

    return pd.read_csv(buffer, sep="\t")


def make_top_ncpm_plot(results, top_n=10):
    """
    Plot the top gene-cell type combinations by average nCPM.

    X-axis:
    - Cell type + gene combination

    Y-axis:
    - Average nCPM

    Bars are sorted so the highest value starts from the left.
    """
    if results.ncpm_df.empty:
        return None, pd.DataFrame()

    required_cols = ["Gene name", "Cell type", "nCPM"]

    missing_cols = [
        col for col in required_cols
        if col not in results.ncpm_df.columns
    ]

    if missing_cols:
        raise ValueError(f"nCPM table is missing columns: {missing_cols}")

    plot_df = results.ncpm_df[required_cols].copy()

    plot_df["nCPM"] = pd.to_numeric(plot_df["nCPM"], errors="coerce")
    plot_df = plot_df.dropna(subset=["Gene name", "Cell type", "nCPM"])

    if plot_df.empty:
        return None, pd.DataFrame()

    top_df = (
        plot_df.groupby(["Cell type", "Gene name"], as_index=False)["nCPM"]
        .mean()
        .rename(columns={"nCPM": "Average nCPM"})
        .sort_values("Average nCPM", ascending=False)
        .head(int(top_n))
        .reset_index(drop=True)
    )

    if top_df.empty:
        return None, top_df

    top_df["Cell type + gene"] = (
        top_df["Cell type"].astype(str)
        + " | "
        + top_df["Gene name"].astype(str)
    )

    # Wider figure for long x-axis labels.
    fig_width = max(10, 0.75 * len(top_df))
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    ax.bar(top_df["Cell type + gene"], top_df["Average nCPM"])

    ax.set_xlabel("Cell type and gene")
    ax.set_ylabel("Average nCPM")
    ax.set_title(f"Top {top_n} cell-gene combinations by average nCPM")

    ax.tick_params(axis="x", rotation=45)

    for label in ax.get_xticklabels():
        label.set_ha("right")

    # Reduce plot left/right margins and keep enough bottom room for labels.
    ax.margins(x=0.01)
    fig.subplots_adjust(
        left=0.06,
        right=0.99,
        bottom=0.36,
        top=0.88,
    )

    return fig, top_df


# -----------------------------
# Sidebar
# -----------------------------

st.sidebar.header("1. Upload gene list")

gene_file = st.sidebar.file_uploader(
    "Upload candidate gene list TSV, TXT, or CSV",
    type=["tsv", "txt", "csv"],
)

st.sidebar.header("2. Choose settings")

selected_tissue = st.sidebar.text_input(
    "Enter tissue name",
    value="Immune cells",
)

threshold = st.sidebar.number_input(
    "Expression threshold",
    min_value=0.0,
    value=1.0,
    step=0.5,
)

top_n = st.sidebar.number_input(
    "Number of gene-cell type combinations to show",
    min_value=1,
    max_value=100,
    value=10,
    step=1,
)

run_button = st.sidebar.button("Run MendelCell analysis")


# -----------------------------
# Reference file status
# -----------------------------

with st.expander("Reference file status"):
    st.write("Expected reference files:")

    st.code(str(CLUSTER_REFERENCE))

    if CLUSTER_REFERENCE.exists():
        st.success(
            f"Cluster reference found "
            f"({CLUSTER_REFERENCE.stat().st_size / 1_000_000:.2f} MB)"
        )
    else:
        st.error("Cluster reference file is missing.")

    st.code(str(CELLTYPE_REFERENCE))

    if CELLTYPE_REFERENCE.exists():
        st.success(
            f"Cell-type reference found "
            f"({CELLTYPE_REFERENCE.stat().st_size / 1_000_000:.2f} MB)"
        )
    else:
        st.error("Cell-type reference file is missing.")


# -----------------------------
# Stop if no gene file uploaded
# -----------------------------

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


# -----------------------------
# Read gene list
# -----------------------------

try:
    gene_table = read_gene_list(gene_file.name, gene_file.getvalue())

except Exception as e:
    st.error(f"Could not read gene list file: {e}")
    st.exception(e)
    st.stop()


# -----------------------------
# Validate gene list
# -----------------------------

if "Gene Symbol" not in gene_table.columns:
    st.error("Gene list must contain a column named 'Gene Symbol'.")
    st.write("Columns found:")
    st.write(list(gene_table.columns))
    st.stop()


st.success(
    f"Gene list uploaded: {gene_table.shape[0]:,} rows and "
    f"{gene_table.shape[1]:,} columns"
)

with st.expander("Preview uploaded gene list"):
    show_dataframe_with_1_index(gene_table.head(20), height=250, width=700)


# -----------------------------
# Wait for Run button
# -----------------------------

if not run_button:
    st.info("Enter a tissue and threshold, then click **Run MendelCell analysis**.")
    st.stop()


# -----------------------------
# Load reference files after Run
# -----------------------------

try:
    with st.spinner("Loading MendelCell reference files..."):
        clusters, hpa = load_reference_data()

    st.success(
        f"Reference loaded: {clusters.shape[0]:,} cluster rows and "
        f"{hpa.shape[0]:,} cell-type rows"
    )

except Exception as e:
    st.error("Could not load MendelCell reference files.")
    st.exception(e)
    st.stop()


# -----------------------------
# Optional available tissues display
# -----------------------------

try:
    available_tissues = list_tissues(clusters)

    with st.expander("Available tissues in reference"):
        st.write(available_tissues)

except Exception as e:
    st.warning(f"Could not list available tissues: {e}")


# -----------------------------
# Run MendelCell analysis
# -----------------------------

try:
    with st.spinner("Running MendelCell analysis..."):
        results = run_mendelcell(
            clusters=clusters,
            hpa=hpa,
            gene_table=gene_table,
            tissue=selected_tissue,
            threshold=threshold,
        )

except Exception as e:
    st.error(f"MendelCell analysis failed: {e}")
    st.exception(e)
    st.stop()


# -----------------------------
# Display results
# -----------------------------

st.success("Analysis complete.")

col1, col2, col3, col4 = st.columns(4)

col1.metric("Input genes", len(results.gene_symbols))
col2.metric("Tissue-specific cell types", len(results.unique_cells))
col3.metric("Genes passing threshold", results.filtered["Gene name"].nunique())
col4.metric(
    "Gene-cell pairs",
    len(results.filtered[["Cell type", "Gene name"]].drop_duplicates()),
)


st.header(f"Top {top_n} cell-gene combinations by average nCPM")

try:
    top_ncpm_fig, top_ncpm_df = make_top_ncpm_plot(results, top_n=top_n)

    if top_ncpm_fig is None or top_ncpm_df.empty:
        st.info(f"No nCPM values available for top-{top_n} plotting.")
    else:
        st.pyplot(top_ncpm_fig, width="content", clear_figure=True)

        show_dataframe_with_1_index(top_ncpm_df, height=400, width=1050)

except Exception as e:
    st.error(f"Could not create top-{top_n} nCPM plot.")
    st.exception(e)


# -----------------------------
# Download outputs
# -----------------------------

st.header("Download outputs")

unique_tsv = results.unique_to_tissue.to_csv(sep="\t", index=False)
filtered_tsv = results.filtered.to_csv(sep="\t", index=False)
ncpm_tsv = results.ncpm_df.to_csv(sep="\t", index=False)

if "top_ncpm_df" in locals() and not top_ncpm_df.empty:
    top_ncpm_tsv = top_ncpm_df.to_csv(sep="\t", index=False)
else:
    top_ncpm_tsv = ""

safe_tissue = safe_filename(results.selected_tissue)
pdf_filename = f"MendelCell_report_{safe_tissue}.pdf"


try:
    with tempfile.TemporaryDirectory() as tmpdir:
        pdf_path = Path(tmpdir) / pdf_filename
        create_pdf_report(results, pdf_path)

        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()

    st.download_button(
        label="Download PDF report",
        data=pdf_bytes,
        file_name=pdf_filename,
        mime="application/pdf",
    )

except Exception as e:
    st.error("Could not create PDF report.")
    st.exception(e)


st.download_button(
    label="Download unique cell types TSV",
    data=unique_tsv,
    file_name="unique_cell_types.tsv",
    mime="text/tab-separated-values",
)

st.download_button(
    label="Download filtered candidate genes TSV",
    data=filtered_tsv,
    file_name="filtered_candidate_genes.tsv",
    mime="text/tab-separated-values",
)

st.download_button(
    label="Download nCPM table TSV",
    data=ncpm_tsv,
    file_name="candidate_gene_ncpm_by_cell_type.tsv",
    mime="text/tab-separated-values",
)

if top_ncpm_tsv:
    st.download_button(
        label=f"Download top {top_n} average nCPM TSV",
        data=top_ncpm_tsv,
        file_name=f"top_{top_n}_cell_gene_average_ncpm.tsv",
        mime="text/tab-separated-values",
    )