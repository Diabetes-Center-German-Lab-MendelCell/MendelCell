from dataclasses import dataclass

import pandas as pd

from .io import clean_input_tables, load_input_files


# -----------------------------
# Special pseudo-tissue: Immune cells
# -----------------------------

SPECIAL_IMMUNE_TISSUE = "Immune cells"

IMMUNE_TISSUE_ALIASES = {
    "immune",
    "immune cell",
    "immune cells",
    "immune system",
    "blood immune cells",
    "blood",
    "hematopoietic",
    "haematopoietic",
}

IMMUNE_CELL_KEYWORDS = [
    "t cell",
    "t-cell",
    "t cells",
    "t-cells",
    "cd4",
    "cd8",
    "helper t",
    "cytotoxic t",
    "regulatory t",
    "treg",
    "b cell",
    "b-cell",
    "b cells",
    "b-cells",
    "plasma cell",
    "plasmablast",
    "nk cell",
    "nk-cell",
    "natural killer",
    "lymphocyte",
    "monocyte",
    "macrophage",
    "dendritic",
    "neutrophil",
    "eosinophil",
    "basophil",
    "mast cell",
    "myeloid",
    "granulocyte",
    "leukocyte",
    "leucocyte",
    "microglia",
    "kupffer",
    "immune cell",
]


# -----------------------------
# Results container
# -----------------------------

@dataclass
class MendelCellResults:
    selected_tissue: str
    threshold: float
    gene_symbols: list[str]
    unique_to_tissue: pd.DataFrame
    unique_cells: list[str]
    filtered: pd.DataFrame
    ncpm_df: pd.DataFrame
    expression_col: str

    @property
    def cell_count_df(self) -> pd.DataFrame:
        """Number of candidate genes found in each cell type."""
        if self.filtered.empty:
            return pd.DataFrame(columns=["Cell type", "Gene count"])

        return (
            self.filtered.groupby("Cell type")["Gene name"]
            .nunique()
            .reset_index(name="Gene count")
            .sort_values("Gene count", ascending=False)
            .reset_index(drop=True)
        )

    @property
    def genes_in_cell_df(self) -> pd.DataFrame:
        """Candidate genes found in each cell type."""
        if self.filtered.empty:
            return pd.DataFrame(columns=["Cell type", "Candidate genes"])

        return (
            self.filtered.groupby("Cell type")["Gene name"]
            .apply(lambda genes: ", ".join(sorted(set(genes))))
            .reset_index(name="Candidate genes")
            .sort_values("Cell type")
            .reset_index(drop=True)
        )

    @property
    def filtered_report(self) -> pd.DataFrame:
        """Filtered candidate-gene results for display/reporting."""
        if self.filtered.empty:
            return pd.DataFrame(
                columns=[
                    "Gene name",
                    "Cell type",
                    self.expression_col,
                ]
            )

        cols = ["Gene name", "Cell type"]

        if self.expression_col in self.filtered.columns:
            cols.append(self.expression_col)

        return (
            self.filtered[cols]
            .drop_duplicates()
            .sort_values(["Cell type", "Gene name"])
            .reset_index(drop=True)
        )

    @property
    def plot_df(self) -> pd.DataFrame:
        """Data used for the candidate-gene count plot."""
        return self.cell_count_df


# -----------------------------
# Tissue helpers
# -----------------------------

def list_tissues(clusters: pd.DataFrame) -> list[str]:
    """List available tissues from the HPA cluster dataframe."""
    if "Tissue" not in clusters.columns:
        raise ValueError("Cluster dataframe must contain a 'Tissue' column.")

    tissues = sorted(
        clusters["Tissue"]
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
    )

    if SPECIAL_IMMUNE_TISSUE not in tissues:
        tissues.append(SPECIAL_IMMUNE_TISSUE)

    return sorted(tissues)


def resolve_tissue_name(tissue: str, valid_tissues: list[str]) -> str:
    """Resolve case-insensitive tissue input to the original tissue name."""
    if not isinstance(tissue, str) or not tissue.strip():
        raise ValueError("Please enter a tissue name.")

    tissue_input = tissue.strip().lower()

    if tissue_input in IMMUNE_TISSUE_ALIASES:
        return SPECIAL_IMMUNE_TISSUE

    tissue_lookup = {
        valid_tissue.lower(): valid_tissue
        for valid_tissue in valid_tissues
    }

    if tissue_input not in tissue_lookup:
        choices = ", ".join(valid_tissues)
        raise ValueError(
            f"Tissue '{tissue}' was not found. Available tissues: {choices}"
        )

    return tissue_lookup[tissue_input]


def find_immune_cell_types(hpa: pd.DataFrame, clusters: pd.DataFrame) -> list[str]:
    """
    Find immune-related cell types by matching common immune-cell keywords.

    This creates a pseudo-tissue called 'Immune cells' by collecting immune
    cell types across the reference data.
    """
    cell_types = set()

    if "Cell type" in hpa.columns:
        cell_types.update(
            hpa["Cell type"]
            .dropna()
            .astype(str)
            .str.strip()
            .unique()
        )

    if "Cell type" in clusters.columns:
        cell_types.update(
            clusters["Cell type"]
            .dropna()
            .astype(str)
            .str.strip()
            .unique()
        )

    immune_cells = sorted(
        cell_type
        for cell_type in cell_types
        if any(keyword in cell_type.lower() for keyword in IMMUNE_CELL_KEYWORDS)
    )

    if not immune_cells:
        raise ValueError("No immune cell types were found in the reference data.")

    return immune_cells


# -----------------------------
# Main analysis
# -----------------------------

def run_mendelcell(
    clusters: pd.DataFrame,
    hpa: pd.DataFrame,
    gene_table: pd.DataFrame,
    tissue: str,
    threshold: float = 1.0,
) -> MendelCellResults:
    """
    Run MendelCell analysis.

    Parameters
    ----------
    clusters:
        HPA single-cell cluster dataframe.

    hpa:
        HPA single-cell cell-type dataframe.

    gene_table:
        Candidate gene table. Must contain a column named 'Gene Symbol'.

    tissue:
        Tissue name, such as 'Pancreas', or the special pseudo-tissue
        'Immune cells'.

    threshold:
        Expression threshold for the HPA cell-type expression value.
    """

    clusters, hpa, gene_symbols, expression_col = clean_input_tables(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
    )

    selected_tissue = resolve_tissue_name(tissue, list_tissues(clusters))

    # -----------------------------
    # Select cell types
    # -----------------------------

    if selected_tissue == SPECIAL_IMMUNE_TISSUE:
        unique_cells = find_immune_cell_types(hpa=hpa, clusters=clusters)

        unique_to_tissue = pd.DataFrame(
            {
                "Cell type": unique_cells,
                "Tissues": [[SPECIAL_IMMUNE_TISSUE] for _ in unique_cells],
            }
        )

    else:
        celltype_tissues = (
            clusters.groupby("Cell type")["Tissue"]
            .apply(lambda x: sorted(set(x.dropna())))
            .reset_index(name="Tissues")
        )

        unique_to_tissue = celltype_tissues[
            celltype_tissues["Tissues"].apply(
                lambda tissues: tissues == [selected_tissue]
            )
        ].copy()

        unique_cells = unique_to_tissue["Cell type"].tolist()

        if not unique_cells:
            raise ValueError(
                f"No cell types were unique to tissue '{selected_tissue}'."
            )

    # -----------------------------
    # Filter HPA cell-type expression
    # -----------------------------

    filtered = hpa[
        (hpa["Cell type"].isin(unique_cells))
        & (hpa["Gene name clean"].isin(gene_symbols))
        & (hpa[expression_col] >= threshold)
    ].copy()

    if filtered.empty:
        filtered = pd.DataFrame(
            columns=[
                "Gene name",
                "Gene name clean",
                "Cell type",
                expression_col,
            ]
        )

    else:
        keep_cols = [
            "Gene name",
            "Gene name clean",
            "Cell type",
            expression_col,
        ]

        filtered = (
            filtered[keep_cols]
            .drop_duplicates()
            .sort_values(["Cell type", "Gene name"])
            .reset_index(drop=True)
        )

    # -----------------------------
    # Build nCPM table for plots/report
    # -----------------------------

    expressed_gene_pairs = filtered[
        ["Cell type", "Gene name clean"]
    ].drop_duplicates()

    if expressed_gene_pairs.empty:
        ncpm_df = pd.DataFrame(
            columns=[
                "Tissue",
                "Cell type",
                "Gene name",
                "Gene name clean",
                "nCPM",
            ]
        )

    else:
        if selected_tissue == SPECIAL_IMMUNE_TISSUE:
            ncpm_source = clusters[
                clusters["Cell type"].isin(unique_cells)
            ].copy()
        else:
            ncpm_source = clusters[
                clusters["Tissue"] == selected_tissue
            ].copy()

        ncpm_df = ncpm_source.merge(
            expressed_gene_pairs,
            on=["Cell type", "Gene name clean"],
            how="inner",
        )

        if "nCPM" in ncpm_df.columns:
            ncpm_df = (
                ncpm_df[
                    [
                        "Tissue",
                        "Cell type",
                        "Gene name",
                        "Gene name clean",
                        "nCPM",
                    ]
                ]
                .drop_duplicates()
                .sort_values(["Cell type", "nCPM"], ascending=[True, False])
                .reset_index(drop=True)
            )

    return MendelCellResults(
        selected_tissue=selected_tissue,
        threshold=threshold,
        gene_symbols=gene_symbols,
        unique_to_tissue=unique_to_tissue,
        unique_cells=unique_cells,
        filtered=filtered,
        ncpm_df=ncpm_df,
        expression_col=expression_col,
    )


def run_mendelcell_from_files(
    cluster_file: str,
    hpa_file: str,
    gene_file: str,
    tissue: str,
    threshold: float = 1.0,
) -> MendelCellResults:
    """
    Run MendelCell analysis from file paths.

    This is useful for command-line or notebook use.
    """

    clusters, hpa, gene_symbols, expression_col = load_input_files(
        cluster_file=cluster_file,
        hpa_file=hpa_file,
        gene_file=gene_file,
    )

    gene_table = pd.DataFrame({"Gene Symbol": gene_symbols})

    return run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue=tissue,
        threshold=threshold,
    )