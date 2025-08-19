#!/usr/bin/env bash

# Composite of functions useful for analyzing ULTRA data
#
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
                | awk '$1>=10 && $2>=1 {print $1 "\t" $2}' \
                | perl -lane "print \$F[1], \"\t\", \".\" x int(\$F[0]/$bin_divisor)"
            else
                awk '{print $4}' "$file" \
                | sort -n | uniq -c \
                | awk '$1>=10 && $2>=1 {print $1 "\t" $2}' \
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

# ----------------------------
# largest_bin: Print the largest 1 Mb bin (by count) per chromosome for a repeat
# ----------------------------
largest_bin() {
    local file="" x="" chr_start=1 chr_end=20 chr_prefix="Chr"
    OPTIND=1
    while getopts "f:x:s:e:p:" opt; do
        case $opt in
            f) file="$OPTARG" ;;
            x) x="$OPTARG" ;;
            s) chr_start="$OPTARG" ;;
            e) chr_end="$OPTARG" ;;
            p) chr_prefix="$OPTARG" ;;
        esac
    done
    [[ -z "$file" || -z "$x" ]] && { echo "Usage: largest_bin -f file -x repeat [-s start -e end -p prefix]"; return 1; }

    for i in $(seq -w "$chr_start" "$chr_end"); do
        local chr_id="${chr_prefix}${i}"
        echo "$chr_id"
        awk '/'".$chr_id"'/ && $4=='"$x"' {print int($2/1000000)}' "$file" \
        | sort -n | uniq -c \
        | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{if(max>0)print want; else print "NA"}'
        echo
    done
}

# ----------------------------
# arrays_topN: Take top N arrays per chromosome from a pre-sorted tmp file (chr>arraySize)
# Expects input already sorted by chr then array size (like your tmp.*tsv)
# ----------------------------
arrays_topN() {
    local in="" out="" chr_start=1 chr_end=20 chr_prefix="Chr" N=200
    OPTIND=1
    while getopts "i:o:s:e:p:n:" opt; do
        case $opt in
            i) in="$OPTARG" ;;
            o) out="$OPTARG" ;;
            s) chr_start="$OPTARG" ;;
            e) chr_end="$OPTARG" ;;
            p) chr_prefix="$OPTARG" ;;
            n) N="$OPTARG" ;;
        esac
    done
    [[ -z "$in" || -z "$out" ]] && { echo "Usage: arrays_topN -i tmp.tsv -o out.tsv [-s start -e end -p prefix -n N]"; return 1; }

    : > "$out"
    for i in $(seq -w "$chr_start" "$chr_end"); do
        local chr_id="${chr_prefix}${i}"
        echo "$chr_id"
        awk '/'".$chr_id"'/ {print}' "$in" | tail -"$N" >> "$out"
    done
}

# ----------------------------
# arrays_pick: Select top N (or all) arrays per chromosome, optional period filter, TSV/FASTA output.
# Input is expected pre-sorted by chr then array size (desc) so "top N" = last N rows per chr.
# ----------------------------
arrays_pick() {
    local in="" out="" chr_start=1 chr_end=20 chr_prefix="Chr"
    local N=10 all=0 x="" fmt="tsv"   # fmt: tsv|fasta

    OPTIND=1
    while getopts "i:o:s:e:p:n:ax:f:" opt; do
        case $opt in
            i) in="$OPTARG" ;;
            o) out="$OPTARG" ;;
            s) chr_start="$OPTARG" ;;
            e) chr_end="$OPTARG" ;;
            p) chr_prefix="$OPTARG" ;;
            n) N="$OPTARG" ;;
            a) all=1 ;;                       # take all rows (ignore -n)
            x) x="$OPTARG" ;;                 # repeat size filter (column 4)
            f) fmt="$OPTARG" ;;               # tsv (default) or fasta
        esac
    done
    [[ -z "$in" || -z "$out" ]] && {
        echo "Usage: arrays_pick -i tmp.tsv -o out.{tsv|fasta} [-s start -e end -p prefix] [-n N | -a] [-x period] [-f tsv|fasta]"
        return 1
    }
    [[ "$fmt" != "tsv" && "$fmt" != "fasta" ]] && { echo "fmt must be 'tsv' or 'fasta'"; return 1; }

    : > "$out"
    shopt -s nullglob

    for i in $(seq -w "$chr_start" "$chr_end"); do
        local chr_id="${chr_prefix}${i}"
        echo "$chr_id"

        # Filter by chromosome (string match) and optional period (-x) on column 4
        # If you need stricter chr matching, adapt '$0 ~ chr' to target a specific column.
        if (( all )); then
            awk -v chr="$chr_id" -v x="$x" '
                $0 ~ chr && (x=="" || $4==x) { print }
            ' "$in"
        else
            awk -v chr="$chr_id" -v x="$x" '
                $0 ~ chr && (x=="" || $4==x) { print }
            ' "$in" | tail -"$N"
        fi | {
            if [[ "$fmt" == "fasta" ]]; then
                # Build FASTA from columns (adjust if your schema differs):
                # header: >col1_col2_col3_col4; sequence: col9
                awk '{print ">"$1"_"$2"_"$3"_"$4"\n"$9}'
            else
                cat
            fi
        } >> "$out"
    done
}


# ----------------------------
# fasta_window: append selected repeat sequences into FASTA
# -c accepts a single chromosome (e.g., 3 or 03) or a range (e.g., 1-7 or 01-07)
# -l/-r positional window is OPTIONAL and only applied for a single chromosome
# -x (period) is OPTIONAL; omit to take all periods
# Output is ALWAYS appended; file is created if missing.
# Columns assumed: chr(anywhere on line), pos=$2, period=$4, seq=$9
# ----------------------------
fasta_window() {
    local file="" out="" x="" chr_spec="" chr_prefix="Chr" pos1="" pos2=""
    OPTIND=1
    while getopts "f:o:x:c:p:l:r:" opt; do
        case $opt in
            f) file="$OPTARG" ;;
            o) out="$OPTARG"  ;;
            x) x="$OPTARG"    ;;
            c) chr_spec="$OPTARG" ;;
            p) chr_prefix="$OPTARG" ;;
            l) pos1="$OPTARG" ;;
            r) pos2="$OPTARG" ;;
        esac
    done

    # Required args
    [[ -z "$file" || -z "$out" || -z "$chr_spec" ]] && {
        echo "Usage: fasta_window -f file -o out.fa -c <N | N-M> [-x period] [-p prefix] [-l left -r right]"
        return 1
    }

    # If only one of l/r is provided, that's an error
    if { [[ -n "$pos1" ]] && [[ -z "$pos2" ]]; } || { [[ -z "$pos1" ]] && [[ -n "$pos2" ]]; }; then
        echo "Error: -l and -r must be given together (or both omitted)."
        return 1
    fi

    # Ensure output file exists; never truncate
    [[ -e "$out" ]] || : > "$out"

    # Build list of chromosome IDs from single or range
    local ids=()
    if [[ "$chr_spec" == *-* ]]; then
        local cstart="${chr_spec%-*}"
        local cend="${chr_spec#*-}"
        for i in $(seq -w "$cstart" "$cend"); do ids+=("${chr_prefix}${i}"); done
    else
        local cnum; printf -v cnum "%02d" "$chr_spec"
        ids+=("${chr_prefix}${cnum}")
    fi

    # Alternation pattern: Chr01|Chr02|...
    local pat=""
    for id in "${ids[@]}"; do
        [[ -n "$pat" ]] && pat+="|"
        pat+="$id"
    done

    # Positional filter applies ONLY when selecting a single chromosome and both -l/-r given
    local use_pos=0
    if (( ${#ids[@]} == 1 )) && [[ -n "$pos1" && -n "$pos2" ]]; then use_pos=1; fi

    awk -v PAT="$pat" -v X="$x" -v USE_POS="$use_pos" -v START="$pos1" -v STOP="$pos2" '
        $0 ~ PAT && (X=="" || $4==X) && (USE_POS==0 || ($2>=START && $2<=STOP)) {
            print ">"$1"_"$2"_"$3"_"$4"\n"$9
        }' "$file" >> "$out"
}


# ----------------------------
# fasta_centromere_bins_multi: For many genomes (ULTRA*.tsv), pick max 1 Mb bin per chr and dump FASTA
# Output file pattern: ref.x.chrNN.cent_<bin>M.fn
# ----------------------------
fasta_centromere_bins_multi() {
    local glob="ULTRA*tsv" x="" chr_start=1 chr_end=20 chr_prefix="Chr"
    OPTIND=1
    while getopts "g:x:s:e:p:" opt; do
        case $opt in
            g) glob="$OPTARG" ;;
            x) x="$OPTARG" ;;
            s) chr_start="$OPTARG" ;;
            e) chr_end="$OPTARG" ;;
            p) chr_prefix="$OPTARG" ;;
        esac
    done
    [[ -z "$x" ]] && { echo "Usage: fasta_centromere_bins_multi -x period [-g 'ULTRA*tsv' -s start -e end -p prefix]"; return 1; }

    for filename in $glob; do
        local ref
        ref=$(echo "$filename" | perl -pe 's/ultra\.(.+)\.p3000\.tsv/$1/')
        echo -e "\n$ref"
        for i in $(seq -w "$chr_start" "$chr_end"); do
            local chr_id="${chr_prefix}${i}"
            local maxBin
            maxBin=$(awk -v CUT="$x" '/'".$chr_id"'/ && $4==CUT {print int($2/1000000)}' "$filename" \
                     | sort -n | uniq -c \
                     | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}')
            [[ -z "$maxBin" ]] && maxBin=0
            local pos1=$((maxBin*1000000))
            local pos2=$((pos1+1000000))
            echo -e "$chr_id\t$x\t$pos1\t$pos2"
            awk -v CUT="$x" -v START="$pos1" -v STOP="$pos2" '/'".$chr_id"'/ && $4==CUT && $2>=START && $2<=STOP {print ">"$1"_"$2"_"$3"_"$4"\n"$9}' "$filename" >> "$ref.$x.$chr_id.cent_${maxBin}M.fn"
        done
    done
}

# ----------------------------
# fasta_centromere_top5: Same as above but only top 5 arrays per chr (by $3 desc) for a single file
# ----------------------------
fasta_centromere_top5() {
    local file="" x="" chr_start=1 chr_end=20 chr_prefix="Gm" ref_override=""
    OPTIND=1
    while getopts "f:x:s:e:p:r:" opt; do
        case $opt in
            f) file="$OPTARG" ;;
            x) x="$OPTARG" ;;
            s) chr_start="$OPTARG" ;;
            e) chr_end="$OPTARG" ;;
            p) chr_prefix="$OPTARG" ;;
            r) ref_override="$OPTARG" ;;
        esac
    done
    [[ -z "$file" || -z "$x" ]] && { echo "Usage: fasta_centromere_top5 -f file -x period [-s start -e end -p prefix -r ref]"; return 1; }

    local ref
    if [[ -n "$ref_override" ]]; then
        ref="$ref_override"
    else
        ref=$(echo "$file" | perl -pe 's/ultra\.(.+)\.p1000\.json\.tsv/$1/')
    fi

    echo -e "\n$ref"
    for i in $(seq -w "$chr_start" "$chr_end"); do
        local chr_id="${chr_prefix}${i}"
        local maxBin
        maxBin=$(awk -v CUT="$x" '/'".$chr_id"'/ && $4==CUT {print int($2/1000000)}' "$file" \
                 | sort -n | uniq -c \
                 | awk -v max=0 '{if($1>max){want=$2; max=$1}}END{print want}')
        [[ -z "$maxBin" ]] && maxBin=0
        local pos1=$((maxBin*1000000))
        local pos2=$((pos1+1000000))
        echo -e "$chr_id\t$x\t$pos1\t$pos2"
        awk -v CUT="$x" -v START="$pos1" -v STOP="$pos2" '/'".$chr_id"'/ && $4==CUT && $2>=START && $2<=STOP {print $0}' "$file" \
        | sort -k3,3nr | head -n 5 \
        | awk '{print ">"$1"_"$2"_"$3"_"$4"\n"$9}' >> "$ref.$x.$chr_id.cent_${maxBin}M.fn"
    done
}

# ----------------------------
# fasta_merge: Combine single-sequence files (no headers) into one FASTA with >filename headers
# ----------------------------
fasta_merge() {
    local pattern="*ext" out="combined_sequences.fasta"
    OPTIND=1
    while getopts "i:o:" opt; do
        case $opt in
            i) pattern="$OPTARG" ;;
            o) out="$OPTARG" ;;
        esac
    done

    : > "$out"
    for f in $pattern; do
        local id
        id=$(basename "$f"); id="${id%.*}"
        echo ">$id" >> "$out"
        cat "$f" >> "$out"
    done
}