#!/usr/bin/env bash

# ULTRA_plot_hist
# Whole-genome histograms (no per-chromosome loop)
# Modes:
#   -m periods (default)
#       -x N or L-U (optional). If omitted, includes all repeat sizes ($4).
#       Plots monomer/period size distribution (histogram of $4).
#   -m arrays
#       -x N or L-U (required). Filters by $4, plots array/polymer sizes (int($3/1000) kb).
#       Auto bar scale to 10 unless -b is provided.
#
# Common:
#   -f FILE (required)
#   -b BINS (default 100 for periods; arrays uses 10 unless overridden)

ULTRA_plot_hist() {
    local file=""
    local x_spec=""
    local bin_divisor=100
    local bin_divisor_set=0
    local mode="periods"

    OPTIND=1
    while getopts "f:x:b:m:" opt "$@"; do
        case $opt in
            f) file="$OPTARG" ;;
            x) x_spec="$OPTARG" ;;
            b) bin_divisor="$OPTARG"; bin_divisor_set=1 ;;
            m) mode="$OPTARG" ;;
            *)
                echo "Usage: ULTRA_plot_hist -f file [-m periods|arrays] [-x N|L-U] [-b bins]"
                return 1 ;;
        esac
    done

    if [[ -z "$file" ]]; then
        echo "Error: -f file is required."
        echo "Usage: ULTRA_plot_hist -f file [-m periods|arrays] [-x N|L-U] [-b bins]"
        return 1
    fi
    [[ -r "$file" ]] || { echo "Error: cannot read file: $file"; return 1; }

    # Parse -x into lb/ub (if provided). Allow N or L-U; strip whitespace; normalize order.
    local rep_single=""
    local lb=""
    local ub=""
    if [[ -n "$x_spec" ]]; then
        x_spec="${x_spec//[[:space:]]/}"
        if [[ "$x_spec" =~ ^[0-9]+-[0-9]+$ ]]; then
            lb="${x_spec%-*}"
            ub="${x_spec#*-}"
        elif [[ "$x_spec" =~ ^[0-9]+$ ]]; then
            rep_single="$x_spec"
            lb="$x_spec"
            ub="$x_spec"
        else
            echo "Error: -x must be N or L-U (e.g., 165 or 155-156). Got: $x_spec"
            return 1
        fi
        if (( lb > ub )); then
            local tmp="$lb"; lb="$ub"; ub="$tmp"
        fi
    fi

    case "$mode" in
        periods)
            # If -x omitted → include all sizes; else filter by [lb,ub]
            if [[ -n "$x_spec" ]]; then
                awk -v lb="$lb" -v ub="$ub" '$4 >= lb && $4 <= ub {print $4}' "$file" \
                | sort -n | uniq -c \
                | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' \
                | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
            else
                awk '{print $4}' "$file" \
                | sort -n | uniq -c \
                | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' \
                | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
            fi
            ;;

        arrays)
            # Require -x (single or range). Default bar scale to 10 if not set.
            if [[ -z "$x_spec" ]]; then
                echo "Error: -x is required for -m arrays (use N or L-U)."
                return 1
            fi
            if [[ $bin_divisor_set -eq 0 ]]; then bin_divisor=10; fi

            awk -v lb="$lb" -v ub="$ub" '$4 >= lb && $4 <= ub {print int($3/1000)}' "$file" \
            | sort -n | uniq -c \
            | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' \
            | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
            ;;

        *)
            echo "Error: unknown mode '$mode' (use: periods|arrays)."
            return 1 ;;
    esac
}
