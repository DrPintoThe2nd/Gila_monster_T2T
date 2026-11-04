#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript plot_bed_chromosomes.R <input.bed[.gz]> [output_prefix]")
}

input_file <- args[1]

# Derive output prefix from input filename if not provided
if (length(args) == 2) {
  prefix <- args[2]
} else {
  prefix <- sub("\\.bed(\\.gz)?$", "", basename(input_file),
                ignore.case = TRUE)
}

# Read BED file, automatically handling .gz files
dat <- read.table(
  input_file,
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

if (ncol(dat) < 4) {
  stop("Input BED file must contain at least 4 columns.")
}

# Preserve chromosome order as it appears in the input
chromosomes <- unique(dat[[1]])

# Output PDF and PNG
pdf(
  paste0(prefix, ".chromosome_plots.pdf"),
  width = 14,
  height = 9
)

png(
  paste0(prefix, ".chromosome_plots.png"),
  width = 14,
  height = 9,
  units = "in",
  res = 300
)

# Plot 18 chromosomes per page
for (start in seq(1, length(chromosomes), by = 18)) {
  
  chr_page <- chromosomes[start:min(start + 17, length(chromosomes))]
  
  for (device in c("pdf", "png")) {
    
    if (device == "png") {
      # Open a separate PNG for each page
      dev.off()
      png(
        paste0(prefix, ".chromosome_plots_page",
               ceiling(start / 18), ".png"),
        width = 14,
        height = 9,
        units = "in",
        res = 300
      )
    }
    
    par(
      mfrow = c(3, 6),
      mar = c(2.5, 2.5, 2, 0.5),
      oma = c(3, 3, 1, 1)
    )
    
    for (chr in chr_page) {
      
      chr_dat <- dat[dat[[1]] == chr, ]
      
      values <- suppressWarnings(as.numeric(chr_dat[[4]]))
      values[!is.finite(values) | values <= 0] <- NA
      
      plot(
        log(values),
        type = "l",
        main = chr,
        xlab = "",
        ylab = "",
        xaxt = "n"
      )
    }
    
    mtext("500-kb windows", side = 1, outer = TRUE, line = 1)
    mtext("log(value)", side = 2, outer = TRUE, line = 1)
    
    if (device == "pdf") {
      # Close PDF after all pages have been written
      if (start == max(seq(1, length(chromosomes), by = 18))) {
        dev.off()
      }
    }
  }
}

# Close final PNG device
dev.off()

cat("Plots saved with prefix:", prefix, "\n")
