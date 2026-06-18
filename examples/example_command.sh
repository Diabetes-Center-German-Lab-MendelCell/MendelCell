#!/usr/bin/env bash

mendelcell \
  --clusters rna_single_cell_cluster.tsv \
  --hpa rna_single_cell_type.tsv \
  --genes Family16-Set1.tsv \
  --tissue Pancreas \
  --threshold 1.0 \
  --outdir mendelcell_output
