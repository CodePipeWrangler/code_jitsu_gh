# usr/bin/env Rscript

# Load necessary libraries
library(ggplot2) # for plotting
library(patchwork) # for combining plots
library(grid)  # for unit()

#' Plot BLAST hits per chromosome
#'
#' Creates a faceted jitter plot showing the distribution of BLAST hits along
#' chromosomes, one facet per chromosome. Each hit is drawn as a jittered point
#' at y = 0, with an optional horizontal bar indicating the chromosome span.
#'
#' @description
#' **Workflow role:** This function renders a single dataset (one species/sample)
#' returned by `load_blast_data()` after filtering/formatting. Use it on its own
#' to inspect one dataset, or let `plot_blast_grid()` call it repeatedly to
#' assemble multi-panel comparisons.
#'
#' @param df A data frame of BLAST hits with at least:
#'   \itemize{
#'     \item \code{sstart} Numeric start position of the hit on the chromosome.
#'     \item \code{chr_length} Numeric total chromosome length (for scaling).
#'     \item \code{sseqid} Factor/character chromosome name (for faceting).
#'   }
#' @param fill_color Character. Color for the jittered points (e.g., \code{"#C8102E"}).
#' @param bar_color Character. Color for the chromosome bar/rectangle (e.g., \code{"#6C7A89"}).
#' @param max_chr_length Numeric. Maximum chromosome length across all datasets; used
#'   to fix a comparable x-scale across facets/plots.
#' @param individual_title Character. Optional title to display on this plot.
#' @param show_title Logical. If \code{TRUE}, shows \code{individual_title}.
#'
#' @return A \code{ggplot} object (suitable for display or patchwork composition).
#'
#' @examples
#' \dontrun{
#' # Single dataset view
#' p <- plot_blast(df = blast_apiam,
#'                 fill_color = "#C8102E",
#'                 bar_color = "#6C7A89",
#'                 max_chr_length = 1.9e8,
#'                 individual_title = "Apiam",
#'                 show_title = TRUE)
#' }
#'
#' @seealso \code{\link{plot_blast_grid}}, \code{\link[patchwork]{plot_annotation}}
#' @import ggplot2
#' @importFrom grid unit
#' @export
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
#' Iterates over a list of pre-filtered BLAST data frames (e.g., produced by
#' \code{load_blast_data()}) and builds one plot per element via \code{plot_blast()},
#' then combines them with \pkg{patchwork}. Supports automatic color wheels or
#' custom color sequences, a consistent bar color, optional global title, and
#' optional alternate titles for each subplot.
#'
#' @description
#' **Workflow role:** This is the top-level compositor. Typical flow is:
#' \enumerate{
#'   \item Use \code{load_blast_data()} to create a named list of filtered data frames
#'         plus metadata (outside this file).
#'   \item Use \code{plot_blast()} to render each dataset.
#'   \item Use \code{plot_blast_grid()} (this function) to assemble all panels and
#'         apply a global title.
#' }
#' All functions are robust and can be used independently as needed.
#'
#' @param blast_list Named list of data frames. Each element must contain
#'   \code{sstart}, \code{chr_length}, and \code{sseqid}.
#' @param ids Character vector of names (subset/order) to select from \code{blast_list}.
#'   Defaults to \code{names(blast_list)}.
#' @param max_chr_length Numeric. Maximum chromosome length across all datasets; ensures
#'   a consistent x-scale across panels.
#' @param colors Character vector of colors used when \code{fill_method = "custom"}.
#'   Defaults to \code{c("#556B2F", "#8B4000", "#6C7A89")}.
#' @param title Character. Optional global title for the combined plot.
#' @param fill_method Character. Either \code{"wheel"} for evenly spaced hues
#'   via \code{grDevices::hcl()} or \code{"custom"} to use \code{colors} (recycled as needed).
#' @param bar_color Character or \code{NULL}. If \code{NULL}, defaults to \code{"#6C7A89"}.
#'   Passed to \code{plot_blast()} for the chromosome bar in each panel.
#' @param alt_titles Optional character vector of the same length as \code{ids} to override
#'   subplot titles (e.g., display friendly names while retaining list element names for lookup).
#'
#' @return A \code{patchwork} plot object (the result of combining per-dataset \code{ggplot}s).
#'
#' @examples
#' \dontrun{
#' # Assume 'prep' was created by load_blast_data(), giving a named list of dfs.
#' ids <- c("Apiam", "Apios")
#' p <- plot_blast_grid(
#'   blast_list     = prep,
#'   ids            = ids,
#'   max_chr_length = 1.9e8,
#'   colors         = c("#556B2F", "#8B4000", "#6C7A89"),
#'   title          = "Locations of Putative 193 bp Centromeric Repeat",
#'   fill_method    = "custom",
#'   bar_color      = NULL,              # defaults to "#6C7A89"
#'   alt_titles     = c("Apiam LA2127", "Apios")
#' )
#' # ggsave("patchwork_dual_blasts.pdf", plot = p, width = 14, height = 14)
#' }
#'
#' @seealso \code{\link{plot_blast}}, \code{\link[patchwork]{plot_annotation}}
#' @importFrom grDevices hcl
#' @import patchwork
#' @export
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



