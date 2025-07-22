# usr/bin/env Rscript

#setwd("C:/Users/bdjor/Desktop/temp/apios/blastout") # Set your working directory
setwd(apios_path)

# Load variables from .env file
readRenviron(".env_r")
apios_path <- Sys.getenv("PWD")
apiam_bln_path <- Sys.getenv("APIAM_BLAST")
apipr_bln_path <- Sys.getenv("APIPR_BLAST")
apiam_len <- Sys.getenv("APIAM_LEN")
apipr_len <- Sys.getenv("APIPR_LEN")

# Side-by-side vertical comparison of Apiam and Apipr BLAST results
# Read BLAST output
blast1 <- read.table(apiam_bln_path, header=FALSE, sep="\t",
  col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore","sseq"))
blast2 <- read.table(apipr_bln_path, header=FALSE, sep="\t",
  col.names=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore","sseq"))

# Read chromosome lengths   
chr_lengths_1 <- read.table(apiam_len, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths_1) <- c("sseqid", "chr_length")
chr_lengths_2 <- read.table(apipr_len, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chr_lengths_2) <- c("sseqid", "chr_length")

# Change chromosome naming 
blast1$sseqid <- sub(".*\\.", "", blast1$sseqid)
chr_lengths_1$sseqid <- sub(".*\\.", "", chr_lengths_1$sseqid)
blast2$sseqid <- sub(".*\\.", "", blast2$sseqid)
chr_lengths_2$sseqid <- sub(".*\\.", "", chr_lengths_2$sseqid)


# merge
blast1 <- merge(blast1, chr_lengths_1, by = "sseqid")
blast2 <- merge(blast2, chr_lengths_2, by = "sseqid")

# Get maximum chromosome length for setting x-axis limits
max_chr_length <- max(c(blast1$chr_length, blast2$chr_length), na.rm = TRUE)

# Order chromosomes
chrom_order <- paste0("Chr", sprintf("%02d", 1:11))
blast1$sseqid <- factor(blast1$sseqid, levels = chrom_order)
blast2$sseqid <- factor(blast2$sseqid, levels = chrom_order)

# Define alternating ISU colors (Cardinal and Gold)
#custom_colors <- c("#C8102E", "#F1BE48")
# other color options:
#custom_colors <- c("#1B9E77", "#4C5C68", "#D95F02")  # Clean Genomics set
# or
#custom_colors <- c("#B22222", "#4B0082", "#E69F00")  # Bold Academic
# or
custom_colors <- c("#556B2F", "#8B4000", "#6C7A89")  # Earthy Modern
# a method to add the bar color to the plot should be implemented...

# Load necessary libraries
library(ggplot2) # for plotting
library(patchwork) # for combining plots
library(grid)  # for unit()

# Define a helper plotting function
plot_blast <- function(df, fill_color = "#C8102E", bar_color = "#6C7A89", max_chr_length,
                       individual_title = NULL, show_title = FALSE) {
  ggplot(df, aes(x = sstart)) +
    
    # Chromosome bar under the jitter points
    # geom_segment(
    #     aes(x = 0, xend = chr_length, y = -0.4, yend = -0.4),
    #     inherit.aes = FALSE,
    #     color = "navy",
    #     size = 0.8,
    #     alpha = 0.9) + # adjust this value for transparency
    
    # Or if you want a filled rectangle
    geom_rect(
        aes(xmin = 0, xmax = chr_length, ymin = -0.45, ymax = -0.35),
        inherit.aes = FALSE,
        fill = bar_color,
        alpha = 0.2) +

    # Plot jittered hits
    geom_jitter(aes(y = 0), height = 0.3, size = 0.8, alpha = 0.7, color = fill_color) +
    
    # Force facet scaling to fit chromosome length
    geom_blank(aes(x = chr_length)) +
    
    facet_wrap(~ sseqid, scales = "fixed", ncol = 1, drop = FALSE) +
    coord_cartesian(xlim = c(0, max_chr_length), ylim = c(-0.5, 0.5)) +
    
    labs(
      x = NULL,
      y = "193 bp",
      title = if (show_title) individual_title else NULL
    ) +
    
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(angle = 90, size = 10),
      strip.text = element_text(size = 9, hjust = 0),
      panel.grid = element_blank(),
      panel.spacing = unit(0.4, "lines"),  # tighter vertical space
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
}




# Use color 1 for blast1, color 2 for blast2
p1 <- plot_blast(blast1, fill_color = custom_colors[1], bar_color = custom_colors[3], max_chr_length, individual_title = "Apiam", show_title = TRUE)
p2 <- plot_blast(blast2, fill_color = custom_colors[2], bar_color = custom_colors[3], max_chr_length, individual_title = "Apios", show_title = TRUE)

# Combine them
final_plot <- p1 | p2

# OR add global title
final_plot <- (p1 | p2) +
  patchwork::plot_annotation(
    title = "Locations of Putative 193 bp Centromeric Repeat",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# Save
ggsave("patchwork_dual_blasts.pdf", plot = final_plot, width = 14, height = 14)
