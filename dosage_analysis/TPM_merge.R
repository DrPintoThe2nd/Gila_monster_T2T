#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript script.R file1.tsv file2.tsv output.tsv")
}

file1 <- args[1]
file2 <- args[2]
outfile <- args[3]

# Extract sample names from filenames (before first underscore)
sample1 <- sub("_.*", "", basename(file1))
sample2 <- sub("_.*", "", basename(file2))

# Read data
df1 <- read.delim(file1, sep = "\t", check.names = FALSE)
df2 <- read.delim(file2, sep = "\t", check.names = FALSE)

# Subset columns
df1_subset <- df1[, c(1:13, 15)]
df2_subset <- df2[, 15, drop = FALSE]

# Rename TPM columns
colnames(df1_subset)[ncol(df1_subset)] <- sample1
colnames(df2_subset)[1] <- sample2

# Combine like `paste`
out <- cbind(df1_subset, df2_subset)

# Write output
write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
