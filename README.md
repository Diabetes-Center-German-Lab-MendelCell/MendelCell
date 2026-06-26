---
title: MendelCell
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 8501
pinned: false
---
# MendelCell

MendelCell is a Python package for prioritizing candidate genes by tissue-specific single-cell expression using Human Protein Atlas single-cell expression files.

The package version separates the original single script into:

- `mendelcell/io.py` — read and clean input files
- `mendelcell/analysis.py` — reusable analysis functions
- `mendelcell/report.py` — TSV and PDF report generation
- `mendelcell/cli.py` — command-line interface
- `mendelcell/app.py` — placeholder for an optional Streamlit app

## Installation

From this folder:

```bash
pip install -e .
```

For the optional Streamlit app later:

```bash
pip install -r requirements-app.txt
```

## Required input files

MendelCell expects three tab-separated files:

1. HPA cluster-level file, usually `rna_single_cell_cluster.tsv` or `.zip`
2. HPA cell-type-level file, usually `rna_single_cell_type.tsv`
3. Candidate gene file containing a column named `Gene Symbol`

Required columns:

### Cluster file

- `Tissue`
- `Cell type`
- `Gene name`
- `nCPM`

### HPA cell-type file

- `Gene name`
- `Cell type`
- `nTPM` or `nCPM`

### Gene list file

- `Gene Symbol`

## Command-line use

```bash
mendelcell \
  --clusters rna_single_cell_cluster.tsv \
  --hpa rna_single_cell_type.tsv \
  --genes Family16-Set1.tsv \
  --tissue Pancreas \
  --threshold 1.0 \
  --outdir mendelcell_output
```

This creates:

- `unique_cell_types.tsv`
- `filtered_candidate_genes.tsv`
- `candidate_gene_ncpm_by_cell_type.tsv`
- `MendelCell_report_<TISSUE>.pdf`

## Python use

```python
from mendelcell import run_mendelcell_from_files
from mendelcell.report import create_pdf_report, write_tsv_outputs

results = run_mendelcell_from_files(
    cluster_file="rna_single_cell_cluster.tsv",
    hpa_file="rna_single_cell_type.tsv",
    gene_file="Family16-Set1.tsv",
    tissue="Pancreas",
    threshold=1.0,
)

write_tsv_outputs(results, "mendelcell_output")
create_pdf_report(results, "mendelcell_output/MendelCell_report_Pancreas.pdf")
```

## Current analysis logic

The current version identifies cell types that are unique to the selected tissue. A cell type is considered tissue-specific if it appears only in the selected tissue in the HPA cluster table.

Then MendelCell filters candidate genes that:

1. Occur in the user-provided candidate gene list
2. Are expressed in those tissue-specific cell types
3. Have expression greater than or equal to the selected threshold

## Next steps before journal submission

- Add tests with small known example datasets
- Add biological validation examples
- Add threshold sensitivity analysis
- Add documentation for HPA data version and provenance
- Add a real Streamlit app using `mendelcell.analysis` and `mendelcell.report`
- Add a `CITATION.cff` file and Zenodo DOI
