# usr/bin/env Rscript

#' plot_blast
#'
#' Creates a faceted jitter plot showing the distribution of BLAST hits along chromosomes.
#' Each facet represents a chromosome, and each hit is shown as a jittered point at y = 0,
#' with an optional horizontal bar indicating the chromosome span.
#'
#' @param df A data frame of BLAST hits. Must contain the following columns:
#'   - `sstart`: Start position of the hit on the chromosome
#'   - `chr_length`: Total length of the chromosome (used for scaling)
#'   - `sseqid`: Chromosome name (used for faceting)
#'
#' @param fill_color Color for the jittered points (e.g., "#C8102E")
#' @param bar_color Color for the horizontal chromosome bar (e.g., "#4C5C68")
#' @param max_chr_length Numeric value indicating the maximum chromosome length across all datasets
#' @param individual_title Optional title for the full plot (used when faceted or patchworked)
#' @param show_title Logical; if TRUE, displays the individual_title on the plot
#'
#' @return A ggplot object ready for display or composition with patchwork
#'
#' @examples #1
#' # Use color 1 for blast1, color 2 for blast2
#' p1 <- plot_blast(blast1, fill_color = custom_colors[1], bar_color = custom_colors[3], max_chr_length, individual_title = "Apiam", show_title = TRUE)
#' p2 <- plot_blast(blast2, fill_color = custom_colors[2], bar_color = custom_colors[3], max_chr_length, individual_title = "Apios", show_title = TRUE)
#' # Combine them
#' final_plot <- p1 | p2
#'
#' # OR add global title
#' final_plot <- (p1 | p2) +
#'   patchwork::plot_annotation(
#'     title = "Locations of Putative 193 bp Centromeric Repeat",
#'     theme = theme(
#'       plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
#'     )
#'   )
#'
#' # Save
#' ggsave("patchwork_dual_blasts.pdf", plot = final_plot, width = 14, height = 14)


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

