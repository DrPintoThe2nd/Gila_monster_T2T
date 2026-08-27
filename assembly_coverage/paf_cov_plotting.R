#!/usr/bin/env Rscript

# =========================================================================
# paf_cov_plotting.R
#
# Produce coverage plots from a PAF file, WITHOUT the pafr package.
# For each of two plot styles (continuous depth + binary covered/not),
# writes PNG, SVG, and PDF.
#
# Code authored by Claude (Anthropic) with user direction. Review and validate before use.
#
# Usage:
#   Rscript paf_cov_plotting.R <file.paf> [output_prefix] [bin_size]
#
# Examples:
#   Rscript paf_cov_plotting.R hap1_20x.paf
#   Rscript paf_cov_plotting.R hap1_20x.paf hap1_cov 5000
#
# Outputs (default prefix = input filename without .paf):
#   <prefix>_continuous.png / .svg / .pdf
#   <prefix>_binary.png     / .svg / .pdf
# =========================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ---- 0. Parse command-line arguments ------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript paf_cov_plotting.R <file.paf> [output_prefix] [bin_size]",
       call. = FALSE)
}

paf_file <- args[1]
if (!file.exists(paf_file)) {
  stop(sprintf("PAF file not found: %s", paf_file), call. = FALSE)
}

# default prefix: input name minus a trailing .paf (case-insensitive)
prefix <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  sub("\\.paf$", "", basename(paf_file), ignore.case = TRUE)
}

bin_size <- if (length(args) >= 3 && nzchar(args[3])) {
  as.numeric(args[3])
} else {
  1e4   # 10 kb bins
}
if (is.na(bin_size) || bin_size <= 0) {
  stop("bin_size must be a positive number", call. = FALSE)
}

message(sprintf("Input:    %s", paf_file))
message(sprintf("Prefix:   %s", prefix))
message(sprintf("Bin size: %s bp", format(bin_size, scientific = FALSE)))

# ---- 1. Read the PAF file manually --------------------------------------
paf_cols <- c("qname", "qlen", "qstart", "qend", "strand",
              "tname", "tlen", "tstart", "tend",
              "nmatch", "alen", "mapq")

ali1 <- read.table(paf_file, sep = "\t", header = FALSE,
                   fill = TRUE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE,
                   col.names = paste0("V", 1:50))   # over-allocate for tags

ali1 <- ali1[, 1:12]
colnames(ali1) <- paf_cols

num_cols <- c("qlen", "qstart", "qend", "tlen", "tstart", "tend",
              "nmatch", "alen", "mapq")
ali1[num_cols] <- lapply(ali1[num_cols], as.numeric)

# ---- 2. Filter: long, confidently-mapped alignments ----------------------
long_ali1 <- subset(ali1, alen > 1e4 & mapq > 30)

if (nrow(long_ali1) == 0) {
  stop("No alignments pass the filter (alen > 1e4 & mapq > 30).", call. = FALSE)
}

# ---- 3. Compute per-bin coverage depth over each target sequence ---------
cov_list <- lapply(split(long_ali1, long_ali1$tname), function(df) {
  tlen  <- df$tlen[1]
  breaks <- seq(0, tlen, by = bin_size)
  if (tail(breaks, 1) < tlen) breaks <- c(breaks, tlen)
  starts <- head(breaks, -1)
  ends   <- tail(breaks, -1)
  
  depth <- vapply(seq_along(starts), function(i) {
    sum(df$tstart < ends[i] & df$tend > starts[i])
  }, numeric(1))
  
  data.frame(tname = df$tname[1], xmin = starts, xmax = ends, depth = depth)
})

cov <- do.call(rbind, cov_list)

# binary covered / not-covered flag
cov$covered <- factor(ifelse(cov$depth > 0, "Covered", "Not covered"),
                      levels = c("Not covered", "Covered"))

# ---- 4. Robust chromosome ordering --------------------------------------
# Rules:
#   * Autosomes sort numerically:  chr1 < chr2 < ... < chr17
#   * "_alt" contigs form a second group, also numeric, placed after them.
#   * Sex chromosomes (Z then W) always go to the very bottom.
#   * Anything unrecognized is parked above the sex chromosomes.
chrom_rank <- function(x) {
  x   <- as.character(x)
  alt <- grepl("_alt", x, ignore.case = TRUE)
  num <- suppressWarnings(as.integer(gsub("\\D", "", x)))
  
  is_z <- grepl("(^|[^a-z])z$", x, ignore.case = TRUE) & is.na(num)
  is_w <- grepl("(^|[^a-z])w$", x, ignore.case = TRUE) & is.na(num)
  
  tier <- ifelse(!is.na(num) & !alt, 1L,
                 ifelse(!is.na(num) &  alt, 2L,
                        ifelse(is_z, 4L,
                               ifelse(is_w, 5L, 3L))))
  
  within <- ifelse(!is.na(num), num, NA_integer_)
  
  data.frame(name = x, tier = tier, num = within, stringsAsFactors = FALSE)
}

uniq   <- unique(as.character(cov$tname))
rk     <- chrom_rank(uniq)
ord_df <- rk[order(rk$tier, rk$num, rk$name, na.last = TRUE), ]
ord    <- ord_df$name

# facet_grid stacks the FIRST level at the TOP, so ranked order puts
# chr1 at top and chrW at the bottom.
cov$tname <- factor(cov$tname, levels = ord)

# ---- 5. Plot builders ----------------------------------------------------
base_theme <- theme_minimal() +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid   = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1))

plot_continuous <- function() {
  ggplot(cov, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
                  fill = depth)) +
    geom_rect() +
    facet_grid(tname ~ ., switch = "y") +
    scale_fill_gradient(low = "#e5f5e0", high = "forestgreen",
                        name = "Coverage") +
    labs(x = "Position (bp)", y = NULL) +
    base_theme
}

plot_binary <- function() {
  ggplot(cov, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
                  fill = covered)) +
    geom_rect() +
    facet_grid(tname ~ ., switch = "y") +
    scale_fill_manual(values = c("Not covered" = "#f0f0f0",
                                 "Covered"     = "forestgreen"),
                      name = NULL) +
    labs(x = "Position (bp)", y = NULL) +
    base_theme
}

# ---- 6. Save each plot in PNG, SVG, and PDF ------------------------------
# height scales with contig count so rows stay legible.
n_contigs <- nlevels(cov$tname)
plot_w <- 10
plot_h <- max(3, n_contigs * 0.4)

#save as svg also, but doesn't work
#save_all <- function(p, tag) {
#  for (ext in c("png", "svg", "pdf")) {
#    fn <- sprintf("%s_%s.%s", prefix, tag, ext)
#    ggsave(fn, plot = p, width = plot_w, height = plot_h,
#           units = "in", dpi = 300, limitsize = FALSE)
#    message(sprintf("Wrote %s", fn))
#  }
#}

save_all <- function(p, tag) {
  for (ext in c("png", "pdf")) {
    fn <- sprintf("%s_%s.%s", prefix, tag, ext)
    ggsave(fn, plot = p, width = plot_w, height = plot_h,
           units = "in", dpi = 300, limitsize = FALSE)
    message(sprintf("Wrote %s", fn))
  }
}

save_all(plot_continuous(), "continuous")
save_all(plot_binary(),     "binary")

message("Done.")