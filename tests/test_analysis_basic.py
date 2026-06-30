from pathlib import Path

import pandas as pd

from mendelcell import run_mendelcell


def make_small_reference_tables():
    """
    Create tiny fake reference tables for testing.

    These are not real HPA data. They are small controlled tables used only
    to check that MendelCell's logic works as expected.
    """

    clusters = pd.DataFrame(
        {
            "Tissue": [
                "Pancreas",
                "Pancreas",
                "Pancreas",
                "Pancreas",
                "Colon",
                "Blood",
                "Blood",
                "Blood",
            ],
            "Cell type": [
                "pancreatic beta cells",
                "pancreatic alpha cells",
                "pancreatic duct cells",
                "pancreatic islet cells",
                "goblet cells",
                "T cells",
                "B cells",
                "macrophages",
            ],
            "Gene name": [
                "INS",
                "GCG",
                "KRT19",
                "INS",
                "MUC2",
                "CD3D",
                "MS4A1",
                "PTPRC",
            ],
            "nCPM": [
                100.0,
                90.0,
                50.0,
                80.0,
                75.0,
                85.0,
                70.0,
                65.0,
            ],
        }
    )

    hpa = pd.DataFrame(
        {
            "Gene name": [
                "INS",
                "GCG",
                "KRT19",
                "INS",
                "MUC2",
                "CD3D",
                "CD8A",
                "MS4A1",
                "PTPRC",
            ],
            "Cell type": [
                "pancreatic beta cells",
                "pancreatic alpha cells",
                "pancreatic duct cells",
                "pancreatic islet cells",
                "goblet cells",
                "T cells",
                "T cells",
                "B cells",
                "macrophages",
            ],
            "nTPM": [
                100.0,
                90.0,
                50.0,
                80.0,
                75.0,
                85.0,
                60.0,
                70.0,
                65.0,
            ],
        }
    )

    return clusters, hpa


def test_example_gene_list_loads_correctly():
    """
    Check that the example gene list exists and has the required column.
    """
    repo_root = Path(__file__).resolve().parents[1]
    example_file = repo_root / "examples" / "example_gene_list.tsv"

    assert example_file.exists(), "Missing examples/example_gene_list.tsv"

    gene_table = pd.read_csv(example_file, sep="\t")

    assert "Gene Symbol" in gene_table.columns
    assert len(gene_table) > 0
    assert "INS" in set(gene_table["Gene Symbol"])


def test_pancreas_analysis_returns_expected_pancreatic_genes_and_cell_types():
    """
    Check that Pancreas analysis returns expected pancreatic genes/cell types.
    """
    clusters, hpa = make_small_reference_tables()

    gene_table = pd.DataFrame(
        {
            "Gene Symbol": [
                "INS",
                "GCG",
                "PDX1",
                "CD3D",
            ]
        }
    )

    results = run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue="Pancreas",
        threshold=1.0,
    )

    filtered_genes = set(results.filtered["Gene name"])
    filtered_cell_types = set(results.filtered["Cell type"])

    assert results.selected_tissue == "Pancreas"

    assert "INS" in filtered_genes
    assert "GCG" in filtered_genes

    assert "pancreatic beta cells" in filtered_cell_types
    assert "pancreatic alpha cells" in filtered_cell_types

    assert "CD3D" not in filtered_genes


def test_immune_cells_analysis_returns_immune_cell_types():
    """
    Check that the Immune cells pseudo-tissue returns immune cell types.
    """
    clusters, hpa = make_small_reference_tables()

    gene_table = pd.DataFrame(
        {
            "Gene Symbol": [
                "CD3D",
                "CD8A",
                "MS4A1",
                "PTPRC",
                "INS",
            ]
        }
    )

    results = run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue="Immune cells",
        threshold=1.0,
    )

    filtered_genes = set(results.filtered["Gene name"])
    immune_cell_types = set(results.unique_cells)

    assert "T cells" in immune_cell_types
    assert "B cells" in immune_cell_types
    assert "macrophages" in immune_cell_types

    assert "CD3D" in filtered_genes
    assert "CD8A" in filtered_genes
    assert "MS4A1" in filtered_genes
    assert "PTPRC" in filtered_genes

    assert "INS" not in filtered_genes


def test_non_immune_cell_types_are_not_counted_as_immune_cells():
    """
    Check that non-immune cell types with names like duct cells, islet cells,
    and goblet cells are not accidentally counted as immune cells.

    This protects against the old loose string-matching issue where words like
    'duct cells', 'islet cells', and 'goblet cells' could accidentally match
    't cells'.
    """
    clusters, hpa = make_small_reference_tables()

    gene_table = pd.DataFrame(
        {
            "Gene Symbol": [
                "KRT19",
                "INS",
                "MUC2",
                "CD3D",
            ]
        }
    )

    results = run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue="Immune cells",
        threshold=1.0,
    )

    immune_cell_types = set(results.unique_cells)
    filtered_cell_types = set(results.filtered["Cell type"])
    filtered_genes = set(results.filtered["Gene name"])

    non_immune_cell_types = {
        "pancreatic duct cells",
        "pancreatic islet cells",
        "goblet cells",
    }

    assert non_immune_cell_types.isdisjoint(immune_cell_types)
    assert non_immune_cell_types.isdisjoint(filtered_cell_types)

    assert "KRT19" not in filtered_genes
    assert "INS" not in filtered_genes
    assert "MUC2" not in filtered_genes

    assert "CD3D" in filtered_genes