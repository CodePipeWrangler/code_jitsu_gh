#!/usr/bin/ Rscript
#' @title Ridgeline Plot of Chromosome-Wise Frequencies
#'
#' @description
#' Generates a ridgeline plot showing the frequency distribution of genomic positions
#' across chromosomes. The plot uses a continuous viridis color gradient to represent
#' frequency intensity and produces both a visual output and an image file.
#'
#' @usage
#' Rscript ridgeline_plot.R <input_file>
#'
#' @param input_file Character. Path to the input file containing three columns:
#' \describe{
#'   \item{frequency}{Numeric frequency or density value for each position.}
#'   \item{position}{Numeric genomic coordinate (e.g., base pair position).}
#'   \item{chromosome}{Chromosome identifier (treated as a factor).}
#' }
#'
#' @details
#' The script reads a tab- or space-delimited file and creates a ridgeline plot
#' where each ridge corresponds to one chromosome. The height and color of each ridge
#' reflect the frequency of observations across genomic positions.
#'
#' The resulting plot is saved automatically to a PNG file in the same directory as
#' the input, using the naming convention:
#' \preformatted{
#' <input_file>_ridgeline_colored.png
#' }
#'
#' The viridis color scale (option = "C") ensures perceptual uniformity and accessibility
#' for color-blind viewers.
#'
#' @output
#' PNG image file saved to the same directory as the input file, along with a plot
#' object printed to the console.
#'
#' @author Brandon Jordan
#' @date 2025-10-01
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[ggridges]{geom_density_ridges_gradient}}
#'
#' @example
#' # Example command-line usage:
#' # Rscript ridgeline_plot.R results/chromosome_frequency_data.tsv
#'
#' @keywords visualization genomics ggplot ridgeline
#' 
#' @note
#' Documentation follows Roxygen2 conventions for standalone scripts.
#' See `tools/doc_guidelines.R` for an explanation of each section.
#' 
#' @license MIT
# ------------------------------------------------------------------------------

# Load required libraries
library(ggplot2)
library(ggridges)
library(viridis)  # For better color mapping

# Check if a file argument is provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript ridgeline_plot.R <input_file>", call. = FALSE)
}

# Read the file path from command-line arguments
input_file <- args[1]

# Read the data (assuming tab-separated file)
chrom_data <- read.table(input_file, header = FALSE, sep = " ", col.names = c("frequency", "position", "chromosome"))

# Ensure chromosome column is treated as a factor
chrom_data$chromosome <- factor(chrom_data$chromosome)

# Generate the ridgeline plot with a color gradient based on frequency
plot <- ggplot(chrom_data, aes(x = position, y = chromosome, height = frequency, fill = after_stat(x), group = chromosome)) +
  geom_density_ridges_gradient(stat = "identity", scale = 2, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "Frequency", option = "C") +
  theme_minimal() +
  labs(title = "Chromosome-Wise Frequency Ridgeline Plot",
       x = "Genomic Position",
       y = "Chromosome")

# Define output file name
output_file <- paste0(tools::file_path_sans_ext(input_file), "_ridgeline_colored.png")

# Save the plot
ggsave(output_file, plot = plot, width = 10, height = 6, dpi = 300)

# Print the plot
print(plot)

# Print completion message
cat("Plot saved as:", output_file, "\n")
