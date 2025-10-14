#!/usr/bin/env bash
#===============================================================================
#: TITLE:        <script_name>.sh
#: DESCRIPTION:  <1–2 lines on what the tool does.>
#: USAGE:        <script_name>.sh [options] <input(s)>
#: EXAMPLES:     <script_name>.sh input1.ext
#:               <script_name>.sh -o out/ *.ext
#: OPTIONS:
#:   -h, --help      Show this help and exit
#:   -o, --outdir D  Write outputs to directory D (default: current dir)
#: INPUTS:       Describe required inputs / formats.
#: OUTPUTS:      Describe files produced and naming pattern.
#: DEPENDENCIES: bash, awk, sed, perl, jq (optional), coreutils
#: EXIT CODES:   0 success; 1 usage error; 2 runtime error
#: AUTHOR:       Brandon Jordan
#: LICENSE:      MIT
#: DATE:         2025-10-14
#===============================================================================

print_help() { grep -E '^#:' "$0" | sed 's/^#:[ ]\?//'; }

# ...rest of the script...
