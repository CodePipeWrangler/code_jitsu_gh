#!/usr/bin/env Rscript
#' @title [SCRIPT_TITLE]
#'
#' @description
#' [Briefly describe the purpose of the script. One to two sentences max.
#' Explain *what* the script does and *why* it is useful.]
#'
#' @usage
#' source("[relative_path]/[script].R")
#' result <- my_function(arg1 = ..., arg2 = ...)
#'
#' @param [param_name] [Type]. [Description of the parameter and its expected format.]
#' @param [param_name2] [Type]. [Additional parameters if applicable.]
#'
#' @details
#' [Provide longer-form details about how the script works internally.
#' Mention assumptions about data format, algorithmic approach, or dependencies.
#' This section can span multiple lines, and should be technical but concise.]
#'
#' Input files are expected to follow the structure:
#' \preformatted{
#' column1    column2    column3
#' ...        ...        ...
#' }
#'
#' Output files will be generated in:
#' \preformatted{
#' results/[subdirectory]/[date]/[output_files]
#' }
#'
#' @output
#' [Summarize outputs produced — files, plots, or printed summaries.
#' Include naming conventions or file extensions, e.g. *_summary.tsv, *_plot.png.]
#'
#' @author Brandon Jordan
#' @date [YYYY-MM-DD]
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}},
#' \code{\link[dplyr]{mutate}},
#' \code{\link[optparse]{parse_args}}
#'
#' @example
#' # Example command-line usage:
#' # Rscript [script_name].R --input data.tsv --outdir results/
#'
#' @keywords [genomics] [data_viz] [automation]
#'
#' @note
#' Documentation follows Roxygen2 conventions for standalone scripts.
#' See `tools/doc_guidelines.R` for an explanation of each section.
#' This script is intended to be **sourced**, not run as a CLI. Roxygen docs are
#' for readers and for the script index. Package-only tags like @export/@import
#' are intentionally omitted in standalone scripts.
#'
#' @license MIT
# ------------------------------------------------------------------------------

# Load required libraries
# library(ggplot2)
# library(dplyr)
# library(optparse)

# Parse arguments ---------------------------------------------------------------
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 0) {
#   stop("Usage: Rscript [script_name].R <input_file>", call. = FALSE)
# }

# input_file <- args[1]
# outdir <- ifelse(length(args) > 1, args[2], "results/default_output")

# Main script logic -------------------------------------------------------------
# data <- read.table(input_file, header = TRUE, sep = "\t")
# [Add analysis code here...]

# Output results ----------------------------------------------------------------
# output_file <- file.path(outdir, paste0(tools::file_path_sans_ext(basename(input_file)), "_output.png"))
# ggsave(output_file, plot = plot, width = 10, height = 6, dpi = 300)
# cat("Output saved to:", output_file, "\n")
