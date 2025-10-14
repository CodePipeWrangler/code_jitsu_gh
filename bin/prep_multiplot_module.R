#!/usr/bin/env Rscript
#' @title Load and Prepare BLAST Datasets from Environment Manifest
#'
#' @description
#' Reads a user-supplied `.env_r` file to locate and load paired BLAST result
#' files (`*.bln`) and chromosome length files (`*_chr_len.txt`). The data are
#' merged, chromosome IDs normalized, and returned as a ready-to-plot list for
#' use with downstream visualization functions such as `plot_blast()` and
#' `plot_blast_grid()`.
#'
#' @usage
#' # Source the file, then call the function:
#' source("genomics/centromere_characterization/scripts/prep_multiplot_module.R")
#' res <- load_blast_data(".env_r")
#' names(res$blast_list)
#' res$max_chr_length
#'
#' @details
#' The `.env_r` file is parsed with `readRenviron()` and must contain environment
#' variables whose values end with:
#' \itemize{
#'   \item `.bln` — paths to BLAST output files
#'   \item `_chr_len.txt` — paths to chromosome length files
#' }
#' The function detects these automatically, merges them, and returns normalized
#' chromosome names (`Chr##`, `ChrMT`, `ChrPT`). It also computes a global
#' `max_chr_length` for consistent plotting scales.
#'
#' @output
#' A list containing:
#'   \itemize{
#'     \item `blast_list` — named list of merged BLAST tables
#'     \item `max_chr_length` — maximum chromosome length
#'     \item `ids` — dataset identifiers
#'     \item `path_map` — input path metadata
#'   }
#'
#' @seealso
#' \code{\link{plot_blast}}, \code{\link{plot_blast_grid}}
#'
#' @note
#' This is a **standalone script** (not an R package). Roxygen headers are for
#' readability and for inclusion in the repository’s script index. To inspect
#' the documentation directly:
#'   In R:
#'     `cat(paste(readLines(".../prep_multiplot_module.R")[grepl("^#'", readLines(".../prep_multiplot_module.R"))], collapse="\n"))`
#'   From shell:
#'     `awk '/^#'\''/{print}' .../prep_multiplot_module.R`
#'
#' @author Brandon Jordan
#' @date 2025-10-14
#' @license MIT
# --- Imports (script) ---------------------------------
library(utils)   # read.table, readRenviron
library(stats)   # setNames

#' Load and prepare BLAST datasets from a `.env_r` manifest
#'
#' @param env_file Character path to a `.env_r` file read via `readRenviron()`.
#'   The environment must include values ending in `.bln` (BLAST hits) and
#'   `_chr_len.txt` (chromosome lengths). Dataset IDs are inferred from the
#'   `_chr_len.txt` basenames (with `_chr_len.*` removed).
#'
#' @return
#' A list with:
#' \describe{
#'   \item{blast_list}{Named list of merged data frames with columns
#'     `qseqid`, `sseqid` (normalized to `Chr##`/`ChrMT`/`ChrPT`), `pident`,
#'     `length`, `sstart`, `chr_length`.}
#'   \item{max_chr_length}{Numeric maximum chromosome length across datasets.}
#'   \item{ids}{Character vector of dataset IDs.}
#'   \item{path_map}{List with named vectors `$bln` and `$len` (input paths).}
#' }
#'
#' @examples
#' # after source("prep_multiplot_module.R"):
#' # res <- load_blast_data(".env_r")
#' # names(res$blast_list); res$max_chr_length
#'
#' @seealso \code{\link{plot_blast}}, \code{\link{plot_blast_grid}}
#' @note Rows with non-standard chromosome IDs are dropped after normalization
#'   (warnings emitted if any are removed).
load_blast_data <- function(env_file) {
  readRenviron(env_file)

  # Gather all environment variables starting with "bln_path_"
  env_vars <- Sys.getenv()

  # Collect values that end with the desired filenames (case-insensitive)
  vals <- unname(env_vars)

  bln_paths <- vals[grepl("\\.bln$", vals, ignore.case = TRUE)]
  len_paths <- vals[grepl("_chr_len\\.txt$", vals, ignore.case = TRUE)]

  # 
  ids <- sub("_chr_len.*$", "", basename(len_paths))
  names(bln_paths) <- ids
  names(len_paths) <- ids


  # (optional but sensible) make order deterministic
  # bln_paths <- bln_paths[order(basename(bln_paths))]
  # len_paths <- len_paths[order(basename(len_paths))]

  if (length(bln_paths) != length(len_paths)) {
    stop("Mismatch: found ", length(bln_paths), " .bln files but ",
        length(len_paths), " _chr_len.txt files in the environment.")
  }

  blast_list <- setNames(vector("list", length(ids)), ids)
  all_chr_lengths <- c()

  # Helper to normalize chromosome IDs to "Chr%02d"
  normalize_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(gsub("\r", "", x))
  low <- tolower(x)

  # first run of digits anywhere (Chr01, chr1, Gm07, etc.)
  num <- suppressWarnings(as.integer(sub(".*?(\\d+).*", "\\1", low, perl = TRUE)))

  # optional special labels (kept silent)
  special <- ifelse(grepl("\\b(mt|chrm|mtdna)\\b", low), "ChrMT",
             ifelse(grepl("\\b(pt|cp|chrc|plastid|chloroplast)\\b", low), "ChrPT",
                    NA_character_))

  ifelse(!is.na(num), sprintf("Chr%02d", num), special)
  }
  # # More verbose version for debugging
  # normalize_chr <- function(x) {
  # cat("RAW:\n"); print(head(x, 3))

  # x <- as.character(x)
  # cat("AS CHARACTER:\n"); print(head(x, 3))

  # x <- gsub("\r", "", x)          # remove carriage returns
  # cat("REMOVE CR:\n"); print(head(x, 3))

  # x <- trimws(x)                  # remove surrounding whitespace
  # cat("TRIMMED:\n"); print(head(x, 3))

  # low <- tolower(x)               # lowercase for easier matching
  # cat("LOWER:\n"); print(head(low, 3))

  # # extract first number run
  # num <- suppressWarnings(as.integer(sub(".*?(\\d+).*", "\\1", low, perl = TRUE)))
  # cat("NUMERIC:\n"); print(head(num, 3))

  # # special cases for mt/pt
  # special <- ifelse(grepl("\\b(mt|chrm|mtdna)\\b", low), "ChrMT",
  #            ifelse(grepl("\\b(pt|cp|chrc|plastid|chloroplast)\\b", low), "ChrPT",
  #                   NA_character_))
  # cat("SPECIAL:\n"); print(head(special, 3))

  # out <- ifelse(!is.na(num), sprintf("Chr%02d", num), special)
  # cat("FINAL:\n"); print(head(out, 3))

  # out
  # }

  overall_max_chr_length <- -Inf
  keep <- c(1, 2, 3, 4, 9)
  for (i in seq_along(bln_paths)) {
    # Read blast file; Keep only needed BLAST columns: 1,2,3,4,9
    blast_all <- read.table(
      bln_paths[i],
      header = FALSE, sep = "\t",
      quote = "", comment.char = "", fill = TRUE,
      stringsAsFactors = FALSE
    )

    if (ncol(blast_all) < max(keep)) {
      stop("BLAST file has only ", ncol(blast_all), " columns; need at least ", max(keep),
          ". File: ", basename(bln_paths[i]))
    }

    blast <- blast_all[keep]
    names(blast) <- c("qseqid","sseqid","pident","length","sstart")

    # Read chr lengths
    chr_len <- read.table(len_paths[i], header = FALSE, sep = "", stringsAsFactors = FALSE,
                  comment.char = "", quote = "", fill = TRUE)

    if (ncol(chr_len) < 2) stop("chr_len file must have ≥2 columns: ", len_paths[i])
    chr_len <- chr_len[, 1:2, drop = FALSE]
    colnames(chr_len) <- c("sseqid", "chr_length")

    # Simplify chromosome names
    blast$sseqid <- sub(".*\\.", "", blast$sseqid)
    chr_len$sseqid <- sub(".*\\.", "", chr_len$sseqid)

    # Normalize BEFORE merge (quiet)
    blast$sseqid   <- normalize_chr(blast$sseqid)
    chr_len$sseqid <- normalize_chr(chr_len$sseqid)

    # If everything failed in a table, hard-stop with a concise message
    if (all(is.na(chr_len$sseqid))) {
      stop("All chromosome IDs in chr_len failed to normalize for file: ", basename(len_paths[i]))
    }
    if (all(is.na(blast$sseqid))) {
      stop("All chromosome IDs in BLAST failed to normalize for file: ", basename(bln_paths[i]))
    }

    # Drop rows that couldn't normalize; warn concisely only if we actually drop any
    dropped_len   <- sum(is.na(chr_len$sseqid))
    dropped_blast <- sum(is.na(blast$sseqid))
    if (dropped_len > 0)   warning("Dropped ", dropped_len, " chr_len rows (non-standard IDs) in ", basename(len_paths[i]))
    if (dropped_blast > 0) warning("Dropped ", dropped_blast, " BLAST rows (non-standard IDs) in ", basename(bln_paths[i]))

    # Merge, then compute factor levels without keeping extra objects
    blast <- merge(blast, chr_len, by = "sseqid", all.x = TRUE)

    per_max <- max(blast$chr_length, na.rm = TRUE)
    overall_max_chr_length <- max(overall_max_chr_length, per_max, na.rm = TRUE)

    levels_chr <- {
      nums <- as.integer(sub("^Chr", "", chr_len$sseqid))
      num_levels <- sprintf("Chr%02d", sort(unique(nums[!is.na(nums)])))
      c(num_levels, setdiff(unique(as.character(chr_len$sseqid)), num_levels))
    }

    blast$sseqid <- factor(blast$sseqid, levels = levels_chr)

    blast_list[[ ids[i] ]] <- blast # store df in list and name by id
  }

  max_chr_length <- overall_max_chr_length

  return(list(
    blast_list     = blast_list,     # named by ids
    max_chr_length = max_chr_length,
    ids            = ids,            # optional, convenience
    path_map       = list(           # optional, tiny metadata
      bln = setNames(bln_paths, ids),
      len = setNames(len_paths, ids)
    )
  ))
}

