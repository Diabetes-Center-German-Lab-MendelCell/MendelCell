from dataclasses import dataclass
import re

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


IMMUNE_CELL_PATTERNS = [
    r"\bimmune\s+cells?\b",
    r"\bt\s*-?\s*cells?\b",
    r"\bcd4\b",
    r"\bcd8\b",
    r"\bhelper\s+t\b",
    r"\bcytotoxic\s+t\b",
    r"\bregulatory\s+t\b",
    r"\btreg\b",
    r"\bb\s*-?\s*cells?\b",
    r"\bplasma\s+cells?\b",
    r"\bplasmablasts?\b",
    r"\bnk\s*-?\s*cells?\b",
    r"\bnatural\s+killer\b",
    r"\blymphocytes?\b",
    r"\bmonocytes?\b",
    r"\bmonocyte\s+progenitors?\b",
    r"\bmacrophages?\b",
    r"\bdendritic\b",
    r"\bneutrophils?\b",
    r"\bneutrophil\s+progenitors?\b",
    r"\beosinophils?\b",
    r"\bbasophils?\b",
    r"\bmast\s+cells?\b",
    r"\bmyeloid\b",
    r"\bgranulocytes?\b",
    r"\bleukocytes?\b",
    r"\bleucocytes?\b",
    r"\bmicroglia\b",
    r"\bkupffer\s+cells?\b",
    r"\bkupffer\b",
]


NON_IMMUNE_CELL_PATTERNS = [
    r"\bduct\s+cells?\b",
    r"\bislet\s+cells?\b",
    r"\bgoblet\s+cells?\b",
    r"\bclub\s+cells?\b",
    r"\btuft\s+cells?\b",
    r"\bepithelial\s+cells?\b",
    r"\bendocrine\s+cells?\b",
    r"\bexocrine\s+cells?\b",
    r"\bacinar\s+cells?\b",
    r"\bbeta\s+cells?\b",
    r"\balpha\s+cells?\b",
    r"\bdelta\s+cells?\b",
    r"\bpancreatic\s+duct\s+cells?\b",
    r"\bpancreatic\s+islet\s+cells?\b",
    r"\bsalivary\s+duct\s+cells?\b",
    r"\bprostatic\s+club\s+cells?\b",
    r"\bconjunctival\s+goblet\s+cells?\b",
]


# -----------------------------
# Results container
# -----------------------------

@dataclass
class MendelCellResults:
    selected_tissue: str
    threshold: float
    non_selected_threshold: float
    max_non_selected_cell_types: int
    use_fraction_mean_ncpm_threshold: bool
    threshold_fraction: float
    gene_symbols: list[str]
    unique_to_tissue: pd.DataFrame
    unique_cells: list[str]
    filtered: pd.DataFrame
    ncpm_df: pd.DataFrame
    selective_genes_df: pd.DataFrame
    expression_col: str

    @property
    def cell_count_df(self) -> pd.DataFrame:
        """Number and names of candidate genes found in each cell type."""
        if self.filtered.empty:
            return pd.DataFrame(
                columns=["Cell type", "Gene count", "Candidate genes"]
            )

        return (
            self.filtered.groupby("Cell type")
            .agg(
                **{
                    "Gene count": ("Gene name", "nunique"),
                    "Candidate genes": (
                        "Gene name",
                        lambda genes: ", ".join(sorted(set(genes))),
                    ),
                }
            )
            .reset_index()
            .sort_values(["Gene count", "Cell type"], ascending=[False, True])
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

        for optional_col in [
            "Threshold source cell type",
            "Threshold source tissue count",
            "Threshold source tissues",
            "Threshold source mean nCPM",
            "Selected threshold",
            "Other-cell threshold",
        ]:
            if optional_col in self.filtered.columns:
                cols.append(optional_col)

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
    Find immune-related cell types using strict regex matching.
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

    immune_cells = []

    for cell_type in sorted(cell_types):
        cell_type_lower = cell_type.lower()

        is_immune = any(
            re.search(pattern, cell_type_lower)
            for pattern in IMMUNE_CELL_PATTERNS
        )

        is_non_immune = any(
            re.search(pattern, cell_type_lower)
            for pattern in NON_IMMUNE_CELL_PATTERNS
        )

        if is_immune and not is_non_immune:
            immune_cells.append(cell_type)

    if not immune_cells:
        raise ValueError("No immune cell types were found in the reference data.")

    return immune_cells


# -----------------------------
# Threshold helper
# -----------------------------

def build_fraction_mean_ncpm_thresholds(
    clusters: pd.DataFrame,
    gene_symbols: list[str],
    selected_tissue: str,
    selected_cell_types: list[str],
    threshold_fraction: float = 1 / 3,
) -> pd.DataFrame:
    """
    Build gene-specific thresholds using a fraction of the highest
    selected-cell-type mean nCPM.

    For each gene:
    - Keep selected cell types only.
    - For each selected cell type, calculate mean nCPM across all tissues
      where that gene-cell-type combination is expressed.
    - Choose the selected cell type with the highest mean nCPM.
    - Selected-cell threshold = highest mean nCPM * threshold_fraction.
    - Other-cell threshold = highest mean nCPM * threshold_fraction.
    """
    output_cols = [
        "Gene name clean",
        "Gene name",
        "Threshold source cell type",
        "Threshold source tissue count",
        "Threshold source tissues",
        "Threshold source mean nCPM",
        "Selected threshold",
        "Other-cell threshold",
    ]

    required_cols = ["Tissue", "Cell type", "Gene name", "Gene name clean", "nCPM"]

    missing_cols = [
        col for col in required_cols
        if col not in clusters.columns
    ]

    if missing_cols:
        raise ValueError(f"Cluster dataframe is missing columns: {missing_cols}")

    source_df = clusters[required_cols].copy()

    source_df["nCPM"] = pd.to_numeric(source_df["nCPM"], errors="coerce")

    source_df = source_df.dropna(
        subset=[
            "Tissue",
            "Cell type",
            "Gene name",
            "Gene name clean",
            "nCPM",
        ]
    )

    source_df = source_df[source_df["Gene name clean"].isin(gene_symbols)].copy()
    source_df = source_df[source_df["Cell type"].isin(selected_cell_types)].copy()

    # Match the HPA-style mean expression table:
    # average only across tissues where expression is present.
    source_df = source_df[source_df["nCPM"] > 0].copy()

    if source_df.empty:
        return pd.DataFrame(columns=output_cols)

    cell_mean_df = (
        source_df.groupby(["Gene name clean", "Cell type"])
        .agg(
            **{
                "Gene name": ("Gene name", "first"),
                "Threshold source mean nCPM": ("nCPM", "mean"),
                "Threshold source tissue count": ("Tissue", "nunique"),
                "Threshold source tissues": (
                    "Tissue",
                    lambda tissues: ", ".join(sorted(set(tissues))),
                ),
            }
        )
        .reset_index()
    )

    threshold_df = (
        cell_mean_df.sort_values(
            ["Gene name clean", "Threshold source mean nCPM"],
            ascending=[True, False],
        )
        .drop_duplicates(subset=["Gene name clean"], keep="first")
        .copy()
    )

    threshold_df = threshold_df.rename(
        columns={
            "Cell type": "Threshold source cell type",
        }
    )

    threshold_df["Selected threshold"] = (
        threshold_df["Threshold source mean nCPM"] * threshold_fraction
    )

    threshold_df["Other-cell threshold"] = (
        threshold_df["Threshold source mean nCPM"] * threshold_fraction
    )

    for col in [
        "Threshold source mean nCPM",
        "Selected threshold",
        "Other-cell threshold",
    ]:
        threshold_df[col] = pd.to_numeric(
            threshold_df[col],
            errors="coerce",
        ).round(2)

    threshold_df = (
        threshold_df[output_cols]
        .sort_values("Gene name")
        .reset_index(drop=True)
    )

    return threshold_df


# -----------------------------
# Selectivity helper
# -----------------------------

def build_selective_genes_df(
    hpa: pd.DataFrame,
    gene_symbols: list[str],
    selected_cell_types: list[str],
    expression_col: str,
    selected_threshold: float,
    non_selected_threshold: float,
    max_non_selected_cell_types: int = 3,
    gene_thresholds_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Find genes that are high in selected cell types and limited in other cells.

    If gene_thresholds_df is provided:
    - selected-cell threshold is gene-specific
    - other-cell threshold is gene-specific
    """
    output_cols = [
        "Gene name",
        "Selected cell types passing threshold",
        "Number of other cell types above threshold",
        "Threshold source tissue count",
        "Threshold source tissues",
        "Threshold source mean nCPM",
        "Selected threshold",
        "Other-cell threshold",
        "Number of selected cell types",
        f"Max selected {expression_col}",
        f"Mean selected {expression_col}",
        "Other cell types above threshold",
        f"Max non-selected {expression_col}",
        f"Mean non-selected {expression_col}",
        "Selected/non-selected max ratio",
    ]

    if hpa.empty:
        return pd.DataFrame(columns=output_cols)

    required_cols = ["Gene name", "Gene name clean", "Cell type", expression_col]

    missing_cols = [
        col for col in required_cols
        if col not in hpa.columns
    ]

    if missing_cols:
        raise ValueError(f"HPA dataframe is missing columns: {missing_cols}")

    expr_df = hpa[required_cols].copy()

    expr_df[expression_col] = pd.to_numeric(
        expr_df[expression_col],
        errors="coerce",
    )

    expr_df = expr_df.dropna(
        subset=[
            "Gene name",
            "Gene name clean",
            "Cell type",
            expression_col,
        ]
    )

    expr_df = expr_df[expr_df["Gene name clean"].isin(gene_symbols)].copy()

    if expr_df.empty:
        return pd.DataFrame(columns=output_cols)

    selected_expr = expr_df[
        expr_df["Cell type"].isin(selected_cell_types)
    ].copy()

    non_selected_expr = expr_df[
        ~expr_df["Cell type"].isin(selected_cell_types)
    ].copy()

    if gene_thresholds_df is not None and not gene_thresholds_df.empty:
        threshold_lookup = gene_thresholds_df[
            [
                "Gene name clean",
                "Threshold source cell type",
                "Threshold source tissue count",
                "Threshold source tissues",
                "Threshold source mean nCPM",
                "Selected threshold",
                "Other-cell threshold",
            ]
        ].copy()

        selected_expr = selected_expr.merge(
            threshold_lookup,
            on="Gene name clean",
            how="left",
        )

        non_selected_expr = non_selected_expr.merge(
            threshold_lookup,
            on="Gene name clean",
            how="left",
        )

        selected_expr["Selected threshold"] = (
            selected_expr["Selected threshold"]
            .fillna(selected_threshold)
        )

        selected_expr["Other-cell threshold"] = (
            selected_expr["Other-cell threshold"]
            .fillna(non_selected_threshold)
        )

        non_selected_expr["Selected threshold"] = (
            non_selected_expr["Selected threshold"]
            .fillna(selected_threshold)
        )

        non_selected_expr["Other-cell threshold"] = (
            non_selected_expr["Other-cell threshold"]
            .fillna(non_selected_threshold)
        )

    else:
        selected_expr["Threshold source cell type"] = ""
        selected_expr["Threshold source tissue count"] = pd.NA
        selected_expr["Threshold source tissues"] = ""
        selected_expr["Threshold source mean nCPM"] = pd.NA
        selected_expr["Selected threshold"] = selected_threshold
        selected_expr["Other-cell threshold"] = non_selected_threshold

        non_selected_expr["Threshold source cell type"] = ""
        non_selected_expr["Threshold source tissue count"] = pd.NA
        non_selected_expr["Threshold source tissues"] = ""
        non_selected_expr["Threshold source mean nCPM"] = pd.NA
        non_selected_expr["Selected threshold"] = selected_threshold
        non_selected_expr["Other-cell threshold"] = non_selected_threshold

    selected_pass = selected_expr[
        selected_expr[expression_col] >= selected_expr["Selected threshold"]
    ].copy()

    if selected_pass.empty:
        return pd.DataFrame(columns=output_cols)

    selected_pass_summary = (
        selected_pass.groupby("Gene name clean")
        .agg(
            **{
                "Gene name": ("Gene name", "first"),
                "Threshold source cell type": (
                    "Threshold source cell type",
                    "first",
                ),
                "Threshold source tissue count": (
                    "Threshold source tissue count",
                    "first",
                ),
                "Threshold source tissues": (
                    "Threshold source tissues",
                    "first",
                ),
                "Threshold source mean nCPM": (
                    "Threshold source mean nCPM",
                    "first",
                ),
                "Selected threshold": ("Selected threshold", "first"),
                "Other-cell threshold": ("Other-cell threshold", "first"),
                "Number of selected cell types": ("Cell type", "nunique"),
                "Selected cell types passing threshold": (
                    "Cell type",
                    lambda cells: ", ".join(sorted(set(cells))),
                ),
            }
        )
        .reset_index()
    )

    selected_stats = (
        selected_expr.groupby("Gene name clean")
        .agg(
            **{
                f"Max selected {expression_col}": (expression_col, "max"),
                f"Mean selected {expression_col}": (expression_col, "mean"),
            }
        )
        .reset_index()
    )

    if non_selected_expr.empty:
        non_selected_stats = pd.DataFrame(
            {
                "Gene name clean": selected_pass_summary["Gene name clean"],
                "Number of other cell types above threshold": 0,
                "Other cell types above threshold": "",
                f"Max non-selected {expression_col}": 0.0,
                f"Mean non-selected {expression_col}": 0.0,
            }
        )

    else:
        non_selected_above_threshold = non_selected_expr[
            non_selected_expr[expression_col]
            >= non_selected_expr["Other-cell threshold"]
        ].copy()

        if non_selected_above_threshold.empty:
            non_selected_above_summary = pd.DataFrame(
                {
                    "Gene name clean": selected_pass_summary["Gene name clean"],
                    "Number of other cell types above threshold": 0,
                    "Other cell types above threshold": "",
                }
            )

        else:
            non_selected_above_summary = (
                non_selected_above_threshold.groupby("Gene name clean")
                .agg(
                    **{
                        "Number of other cell types above threshold": (
                            "Cell type",
                            "nunique",
                        ),
                        "Other cell types above threshold": (
                            "Cell type",
                            lambda cells: ", ".join(sorted(set(cells))),
                        ),
                    }
                )
                .reset_index()
            )

        non_selected_stats = (
            non_selected_expr.groupby("Gene name clean")
            .agg(
                **{
                    f"Max non-selected {expression_col}": (expression_col, "max"),
                    f"Mean non-selected {expression_col}": (expression_col, "mean"),
                }
            )
            .reset_index()
        )

        non_selected_stats = non_selected_stats.merge(
            non_selected_above_summary,
            on="Gene name clean",
            how="left",
        )

        non_selected_stats["Number of other cell types above threshold"] = (
            non_selected_stats["Number of other cell types above threshold"]
            .fillna(0)
            .astype(int)
        )

        non_selected_stats["Other cell types above threshold"] = (
            non_selected_stats["Other cell types above threshold"]
            .fillna("")
        )

    summary_df = selected_pass_summary.merge(
        selected_stats,
        on="Gene name clean",
        how="left",
    )

    summary_df = summary_df.merge(
        non_selected_stats,
        on="Gene name clean",
        how="left",
    )

    max_selected_col = f"Max selected {expression_col}"
    mean_selected_col = f"Mean selected {expression_col}"
    max_non_selected_col = f"Max non-selected {expression_col}"
    mean_non_selected_col = f"Mean non-selected {expression_col}"
    ratio_col = "Selected/non-selected max ratio"

    summary_df["Number of other cell types above threshold"] = (
        summary_df["Number of other cell types above threshold"]
        .fillna(0)
        .astype(int)
    )

    summary_df["Other cell types above threshold"] = (
        summary_df["Other cell types above threshold"]
        .fillna("")
    )

    summary_df[max_non_selected_col] = summary_df[max_non_selected_col].fillna(0.0)
    summary_df[mean_non_selected_col] = summary_df[mean_non_selected_col].fillna(0.0)

    denominator = summary_df[max_non_selected_col].where(
        summary_df[max_non_selected_col] != 0
    )

    summary_df[ratio_col] = summary_df[max_selected_col] / denominator
    summary_df[ratio_col] = summary_df[ratio_col].fillna(float("inf"))

    summary_df = summary_df[
        (summary_df[max_selected_col] >= summary_df["Selected threshold"])
        & (
            summary_df["Number of other cell types above threshold"]
            <= max_non_selected_cell_types
        )
    ].copy()

    if summary_df.empty:
        return pd.DataFrame(columns=output_cols)

    for col in [
        "Threshold source tissue count",
        "Threshold source mean nCPM",
        "Selected threshold",
        "Other-cell threshold",
        max_selected_col,
        mean_selected_col,
        max_non_selected_col,
        mean_non_selected_col,
        ratio_col,
    ]:
        summary_df[col] = pd.to_numeric(summary_df[col], errors="coerce").round(2)

    summary_df["Threshold source tissue count"] = (
        summary_df["Threshold source tissue count"]
        .fillna(0)
        .astype(int)
    )

    summary_df = (
        summary_df[output_cols]
        .sort_values(
            [
                "Number of other cell types above threshold",
                ratio_col,
                max_selected_col,
                "Gene name",
            ],
            ascending=[True, False, False, True],
        )
        .reset_index(drop=True)
    )

    return summary_df


# -----------------------------
# Main analysis
# -----------------------------

def run_mendelcell(
    clusters: pd.DataFrame,
    hpa: pd.DataFrame,
    gene_table: pd.DataFrame,
    tissue: str,
    threshold: float = 1.0,
    non_selected_threshold: float | None = None,
    max_non_selected_cell_types: int = 3,
    use_fraction_mean_ncpm_threshold: bool = False,
    threshold_fraction: float = 1 / 3,
) -> MendelCellResults:
    """
    Run MendelCell analysis.
    """
    if non_selected_threshold is None:
        non_selected_threshold = threshold

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
    # Optional 1/3 mean nCPM thresholds
    # -----------------------------

    if use_fraction_mean_ncpm_threshold:
        gene_thresholds_df = build_fraction_mean_ncpm_thresholds(
            clusters=clusters,
            gene_symbols=gene_symbols,
            selected_tissue=selected_tissue,
            selected_cell_types=unique_cells,
            threshold_fraction=threshold_fraction,
        )

    else:
        gene_thresholds_df = pd.DataFrame(
            columns=[
                "Gene name clean",
                "Gene name",
                "Threshold source cell type",
                "Threshold source tissue count",
                "Threshold source tissues",
                "Threshold source mean nCPM",
                "Selected threshold",
                "Other-cell threshold",
            ]
        )

    # -----------------------------
    # Filter HPA cell-type expression
    # -----------------------------

    if use_fraction_mean_ncpm_threshold and not gene_thresholds_df.empty:
        filtered_source = hpa[
            (hpa["Cell type"].isin(unique_cells))
            & (hpa["Gene name clean"].isin(gene_symbols))
        ].copy()

        filtered_source = filtered_source.merge(
            gene_thresholds_df[
                [
                    "Gene name clean",
                    "Threshold source cell type",
                    "Threshold source tissue count",
                    "Threshold source tissues",
                    "Threshold source mean nCPM",
                    "Selected threshold",
                    "Other-cell threshold",
                ]
            ],
            on="Gene name clean",
            how="left",
        )

        filtered = filtered_source[
            filtered_source[expression_col] >= filtered_source["Selected threshold"]
        ].copy()

    else:
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
                "Threshold source cell type",
                "Threshold source tissue count",
                "Threshold source tissues",
                "Threshold source mean nCPM",
                "Selected threshold",
                "Other-cell threshold",
            ]
        )

    else:
        keep_cols = [
            "Gene name",
            "Gene name clean",
            "Cell type",
            expression_col,
        ]

        for optional_col in [
            "Threshold source cell type",
            "Threshold source tissue count",
            "Threshold source tissues",
            "Threshold source mean nCPM",
            "Selected threshold",
            "Other-cell threshold",
        ]:
            if optional_col in filtered.columns:
                keep_cols.append(optional_col)

        filtered = (
            filtered[keep_cols]
            .drop_duplicates()
            .sort_values(["Cell type", "Gene name"])
            .reset_index(drop=True)
        )

    # -----------------------------
    # Build selective genes table
    # -----------------------------

    selective_genes_df = build_selective_genes_df(
        hpa=hpa,
        gene_symbols=gene_symbols,
        selected_cell_types=unique_cells,
        expression_col=expression_col,
        selected_threshold=threshold,
        non_selected_threshold=non_selected_threshold,
        max_non_selected_cell_types=max_non_selected_cell_types,
        gene_thresholds_df=gene_thresholds_df,
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
        non_selected_threshold=non_selected_threshold,
        max_non_selected_cell_types=max_non_selected_cell_types,
        use_fraction_mean_ncpm_threshold=use_fraction_mean_ncpm_threshold,
        threshold_fraction=threshold_fraction,
        gene_symbols=gene_symbols,
        unique_to_tissue=unique_to_tissue,
        unique_cells=unique_cells,
        filtered=filtered,
        ncpm_df=ncpm_df,
        selective_genes_df=selective_genes_df,
        expression_col=expression_col,
    )


def run_mendelcell_from_files(
    cluster_file: str,
    hpa_file: str,
    gene_file: str,
    tissue: str,
    threshold: float = 1.0,
    non_selected_threshold: float | None = None,
    max_non_selected_cell_types: int = 3,
    use_fraction_mean_ncpm_threshold: bool = False,
    threshold_fraction: float = 1 / 3,
) -> MendelCellResults:
    """
    Run MendelCell analysis from file paths.
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
        non_selected_threshold=non_selected_threshold,
        max_non_selected_cell_types=max_non_selected_cell_types,
        use_fraction_mean_ncpm_threshold=use_fraction_mean_ncpm_threshold,
        threshold_fraction=threshold_fraction,
    )