#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
}))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript join_annotations_dplyr.R annotation.txt expression.txt output.txt")
}

anno_file <- args[1]
expr_file <- args[2]
out_file  <- args[3]

# Read annotation file (no header)
anno <- read_tsv(
  anno_file,
  col_names = c(
    "Chromosome", "Source", "Type", "Start", "Stop", "Score",
    "Strand", "Phase", "Name", "GeneID", "BioType", "Function"
  ),
  quote = "",
  show_col_types = FALSE
)

# Read expression file (has header)
expr <- read_tsv(expr_file, show_col_types = FALSE)

# Left join on Name
merged <- left_join(anno, expr, by = "Name")

# Write output
write_tsv(merged, out_file)

