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

MendelCell is a Streamlit app that helps prioritize candidate genes by looking at their expression across tissue-specific and cell-type-specific single-cell reference data.

The app was built to support exploratory analysis of candidate genes from rare disease, monogenic disease, and immune-related genetics projects.

**Disclaimer:** MendelCell is intended for research and exploratory analysis only. It is not a diagnostic tool and should not be used to make clinical decisions.

## Live App

MendelCell is available on Hugging Face Spaces:

```text
https://landonchamberlain-mendelcell.hf.space
```


## What MendelCell Does

MendelCell allows users to upload a candidate gene list and identify which genes are expressed in relevant cell types.

The app can:

* Upload a candidate gene list as a TSV, TXT, or CSV file
* Analyze a selected tissue, such as `Pancreas`
* Analyze the special pseudo-tissue `Immune cells`
* Apply an expression threshold
* Show candidate genes found in each cell type
* Generate summary tables
* Generate nCPM expression plots
* Show the top 10 gene-cell type combinations by average nCPM
* Download TSV results and a PDF report

## Input File Format

The uploaded file must contain a column named:

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
IL2RA
CTLA4
```

## Example Input File

A small example gene list is included in the repository:

```text
examples/example_gene_list.tsv
```

Users can download or open this file and upload it directly into the MendelCell app to test that the app is working.

Example contents:

```text
Gene Symbol
INS
GCG
PDX1
CD3D
PTPRC
IL2RA
CTLA4
```


## Tissue Options

You can enter tissue names such as:

```text
Pancreas
Liver
Kidney
Lung
```

You can also enter:

```text
Immune cells
```

or:

```text
immune
```

This runs the analysis on immune-related cell types across the reference data.

## Validation Examples

MendelCell includes small example gene lists that can be used to check whether the app returns biologically expected results.

### Pancreas Example

File:

```text
examples/validation_pancreas_genes.tsv
```

Example genes:

```text
INS
GCG
PDX1
SST
```

Expected result: these genes should prioritize pancreatic endocrine or islet-related cell types when the selected tissue is `Pancreas`.

### Immune Cell Example

File:

```text
examples/validation_immune_genes.tsv
```

Example genes:

```text
CD3D
CD4
CD8A
PTPRC
IL2RA
CTLA4
```

Expected result: these genes should prioritize immune-related cell types, especially T-cell-associated cell types, when the selected tissue is `Immune cells`.

### Cell-Line-Oriented Sanity Check

File:

```text
examples/validation_cell_line_genes.tsv
```

This gene list includes markers associated with T-cell, monocyte/macrophage, and pancreatic beta-cell contexts.

Expected result: MendelCell should return immune-cell-associated results for immune markers and pancreatic cell-type-associated results for pancreatic endocrine markers when appropriate tissues are selected.

These examples are intended as simple biological sanity checks. They do not establish diagnostic performance and should not be interpreted as clinical validation.


## Output

MendelCell generates:

* Cell types associated with the selected tissue
* Candidate gene counts per cell type
* Candidate genes detected in each cell type
* Filtered gene-cell type expression results
* Mean nCPM expression values
* A top 10 average nCPM plot
* Downloadable TSV files
* A downloadable PDF report

## Reference Data

MendelCell uses preprocessed Human Protein Atlas single-cell expression data.

Human Protein Atlas version: 25.1  
Download date: June 2026  
Download page: https://www.proteinatlas.org/about/download  

Raw files used:
- rna_single_cell_cluster.tsv
- rna_single_cell_type.tsv

