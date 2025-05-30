#!/usr/bin/ Rscript

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
