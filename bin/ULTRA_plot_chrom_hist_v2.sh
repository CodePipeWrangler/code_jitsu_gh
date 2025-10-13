#!/usr/bin/env bash

# ULTRA_plot_chrom_hist
# -m periods|arrays|specific   (default: periods)
# -x N or L-U   (required for arrays/specific; optional for periods)
# -p chr prefix (default: Chr) used in periods/arrays pattern
# -b bar scale  (default: 100; arrays defaults to 10 if -b not set)
# -f file, -s start chr, -e end chr   (required)

ULTRA_plot_chrom_hist() {
    local file=""
    local x_spec=""           # raw -x (single or range)
    local chr_start=""
    local chr_end=""
    local bin_divisor=100
    local bin_divisor_set=0
    local chr_prefix="Chr"
    local mode="periods"

    OPTIND=1
    while getopts "f:x:s:e:b:p:m:" opt "$@"; do
        case $opt in
            f) file="$OPTARG" ;;
            x) x_spec="$OPTARG" ;;
            s) chr_start="$OPTARG" ;;
            e) chr_end="$OPTARG" ;;
            b) bin_divisor="$OPTARG"; bin_divisor_set=1 ;;
            p) chr_prefix="$OPTARG" ;;
            m) mode="$OPTARG" ;;
            *)
                echo "Usage: ULTRA_plot_chrom_hist -f file -s start -e end [-m periods|arrays|specific] [-x N|L-U] [-p prefix] [-b bins]"
                return 1 ;;
        esac
    done

    # Basic required args
    if [[ -z "$file" || -z "$chr_start" || -z "$chr_end" ]]; then
        echo "Error: missing -f, -s, or -e."
        echo "Usage: ULTRA_plot_chrom_hist -f file -s start -e end [-m periods|arrays|specific] [-x N|L-U] [-p prefix] [-b bins]"
        return 1
    fi

    # Parse -x into lower/upper (if provided)
    local rep_single=""
    local lb=""
    local ub=""
    if [[ -n "$x_spec" ]]; then
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
    fi

    # Mode-specific requirements/defaults
    case "$mode" in
        periods)
            # -x optional; if omitted we don't filter by $4
            :
            ;;
        arrays)
            # Need range or single value; arrays default bar scale to 10 if -b not set
            if [[ -z "$x_spec" ]]; then
                echo "Error: -x is required for -m arrays (use N or L-U)."
                return 1
            fi
            if [[ $bin_divisor_set -eq 0 ]]; then bin_divisor=10; fi
            ;;
        specific)
            # Must be single value, not a range
            if [[ -z "$x_spec" ]]; then
                echo "Error: -x is required for -m specific (single value)."
                return 1
            fi
            if [[ -z "$rep_single" ]]; then
                echo "Error: -m specific expects a single -x value (e.g., -x 104), not a range."
                return 1
            fi
            ;;
        *)
            echo "Error: unknown mode '$mode' (use: periods|arrays|specific)."
            return 1 ;;
    esac

    for i in $(seq -w "$chr_start" "$chr_end"); do
        chr_id="${chr_prefix}${i}"
        echo "$chr_id"

        if [[ "$mode" == "periods" ]]; then
            # Pattern uses chr_id; optionally filter by rep_single if provided
            if [[ -n "$rep_single" ]]; then
                # positions (Mb) for all repeat sizes
                awk '/'".$chr_id"'/ {print $4}' "$file" \
                | sort -n | uniq -c \
                | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' \
                | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
            else # duplicated logic FIXME
                # positions (Mb) for all repeat sizes
                awk '/'".$chr_id"'/ {print $4}' "$file" \
                | sort -n | uniq -c \
                | awk '$1>=10 && $2>=60 {print $1 "\t" $2}' \
                | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
            fi

        elif [[ "$mode" == "arrays" ]]; then
            # array/polymer sizes (kb) within [lb,ub]
            awk '/'".$chr_id"'/ && $4>='"$lb"' && $4<='"$ub"' {print int($3/1000)}' "$file" \
            | sort -n | uniq -c \
            | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' \
            | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"

        else # specific
            # literal .[Cc]hr$i pattern + single repeat size; bins by position (Mb)
            awk ''/'.[Cc]hr'$i'/ && $4=='"$rep_single"' {print int($2/1000000)}' "$file" \
            | sort -n | uniq -c \
            | awk '$1>=1 && $2>=1 {print $1 "\t" $2}' \
            | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
        fi

        echo
    done
}
