#!/usr/bin/env Rscript
#' @title BLAST Hit Plotting Helpers and Grid Compositor
#'
#' @description
#' Plot helpers for visualizing BLAST hit positions per chromosome and a compositor
#' to combine multiple datasets into a patchwork grid. Intended for **sourcing**
#' in an interactive R session or from other scripts.
#'
#' @usage
#' # Source the file, then call the functions:
#' source("genomics/centromere_characterization/scripts/plot_blasts.R")
#' prep <- readRDS("results/demo/prep.rds")  # named list of data.frames
#' p <- plot_blast_grid(prep, ids = c("Apiam","Apios"), max_chr_length = 1.9e8)
#' ggsave("results/demo/blast_grid.pdf", plot = p, width = 14, height = 14)
#'
#' @details
#' Each data frame needs columns: `sstart` (numeric), `chr_length` (numeric),
#' and `sseqid` (chromosome name). `max_chr_length` enforces a common x-scale.
#'
#' @output
#' Returns `ggplot`/`patchwork` objects; save with `ggsave()`.
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}},
#' \code{\link[ggplot2]{geom_jitter}},
#' \code{\link[patchwork]{plot_annotation}},
#' \code{\link[grDevices]{hcl}},
#' \code{\link[grid]{unit}}
#'
#' @note
#' This is a **standalone script** (not an R package). Roxygen here is for humans
#' and for the repo’s script indexer. Dependencies: ggplot2, patchwork, grid.
#' To view just this header: in R
#'   `cat(paste(readLines(".../plot_blasts.R")[grepl("^#'", readLines(".../plot_blasts.R"))],
#'              collapse = "\n"))`
#' or from shell: `awk '/^#'\''/{print}' .../plot_blasts.R`
#'
#' @author Brandon Jordan
#' @date 2025-10-14
#' @license MIT

# --- Imports (script) ---------------------------------
library(ggplot2) # for plotting
library(patchwork) # for combining plots
library(grid)  # for unit()

#' Plot BLAST hits per chromosome
#'
#' @param df Data frame with columns: `sstart`, `chr_length`, `sseqid`.
#' @param fill_color Character hex color for points.
#' @param bar_color  Character hex color for chromosome bar.
#' @param max_chr_length Numeric max chromosome length for x-scale.
#' @param individual_title Character; optional facet title.
#' @param show_title Logical; show `individual_title` if TRUE.
#' @return A ggplot object.
#' @examples
#' # after source("plot_blasts.R"):
#' # plot_blast(df_apiam, "#C8102E", "#6C7A89", 1.9e8, "Apiam", TRUE)
#' @seealso \code{\link{plot_blast_grid}}
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
      y = NULL,
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

#' Compose multiple BLAST plots into a grid
#'
#' @param blast_list Named list of data frames containing columns:
#'   `sstart`, `chr_length`, and `sseqid`.
#' @param ids Character vector specifying which datasets (and in what order)
#'   to include. Defaults to `names(blast_list)`.
#' @param max_chr_length Numeric maximum chromosome length across datasets.
#'   Ensures a consistent x-scale across panels.
#' @param colors Character vector of hex color codes used when
#'   `fill_method = "custom"`. Defaults to
#'   `c("#556B2F", "#8B4000", "#6C7A89")`.
#' @param title Character. Optional overall title for the combined plot.
#' @param fill_method Character; either `"wheel"` for automatically spaced hues
#'   (via `grDevices::hcl()`) or `"custom"` to use supplied colors.
#' @param bar_color Character or `NULL`. Color for chromosome bars
#'   (defaults to `"#6C7A89"` if not specified).
#' @param alt_titles Optional character vector (same length as `ids`) for
#'   alternate subplot titles.
#'
#' @return A patchwork plot object combining multiple ggplot panels.
#' @examples
#' # after source("plot_blasts.R"):
#' # p <- plot_blast_grid(prep, ids = c("Apiam","Apios"), max_chr_length = 1.9e8)
#' # ggsave("results/demo/blast_grid.pdf", plot = p, width = 14, height = 14)
#' @seealso \code{\link{plot_blast}},
#'   \code{\link[patchwork]{plot_annotation}},
#'   \code{\link[grDevices]{hcl}}
plot_blast_grid <- function(blast_list, ids = names(blast_list),
                            max_chr_length,
                            colors = c("#556B2F", "#8B4000", "#6C7A89"),
                            title = NULL,
                            fill_method = c("wheel", "custom"),
                            bar_color = NULL,
                            alt_titles = NULL) {
  fill_method <- match.arg(fill_method)

  # subset/reorder exactly as requested
  bl <- blast_list[ids]
  n  <- length(ids)

  # optional replacement titles (must align 1:1 with ids if supplied)
  if (!is.null(alt_titles)) {
    if (length(alt_titles) != n) {
      stop("alt_titles must have the same length as ids.")
    }
    titles <- alt_titles
  } else {
    titles <- ids
  }

  # fills
  if (fill_method == "wheel") {
    fills <- grDevices::hcl(h = seq(15, 375, length.out = n + 1)[1:n], c = 70, l = 45)
  } else {
    if (length(colors) == 0) stop("When fill_method='custom', 'colors' must be non-empty.")
    k <- length(colors)
    order_idx <- c(seq(1, k, by = 2), seq(2, k, by = 2))
    fills <- rep(colors[order_idx], length.out = n)
  }

  # bar color default
  barc <- if (!is.null(bar_color)) bar_color else "#6C7A89"

  # build one plot per dataset (iterate by index so we can pair ids/fills/titles)
  plots <- Map(
    f = function(i) {
      df <- bl[[ ids[i] ]]
      plot_blast(
        df               = df,
        fill_color       = fills[i],
        bar_color        = barc,
        max_chr_length   = max_chr_length,
        individual_title = titles[i],
        show_title       = TRUE
      )
    },
    i = seq_len(n)
  )

  p <- Reduce(`|`, plots)
  if (!is.null(title)) {
    p <- p + patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  }
  p
}



