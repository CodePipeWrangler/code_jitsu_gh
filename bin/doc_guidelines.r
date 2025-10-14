# Documentation Guidelines for Standalone Scripts
# ------------------------------------------------
# This repository uses Roxygen2-style headers to document standalone R scripts.
# They serve as structured, human-readable documentation and metadata.
#
# Sections and their purposes:
#   @title       — Short name for the script.
#   @description — One- to two-line summary of what the script does.
#   @usage       — Command-line example showing how to run it.
#   @param       — Description of arguments (if script takes inputs).
#   @details     — Technical details, assumptions, and background.
#   @output      — Files or objects created by the script.
#   @author      — Author name.
#   @seealso     — Related R functions or libraries.
#   @example     — Example usage snippet.
#   @note        — Additional context or repository-level notes.
#   @license     — Licensing declaration (MIT).
#
# Accessing Documentation
# -----------------------
# View only the documentation header without loading the script:
#   In R:
#     view_doc <- function(f) cat(paste(readLines(f)[grepl("^#'", readLines(f))], collapse="\n"))
#     view_doc("path/to/script.R")
#
#   Or from shell:
#     awk '/^#'"'"'/{print}' path/to/script.R
#
# For quick inspection (approximate), the first ~40 lines usually cover
# the complete header since Roxygen comments conventionally appear at
# the top of the file.

