#!/usr/bin/env bash
# ULTRA related batch alignments + consensus utilities for FASTA (*.fna)
# deps: clustalo, famsa, emboss 'cons', grep, awk, perl (optional)

set -euo pipefail

# ---------------------- helpers ----------------------
die(){ echo "ERROR: $*" >&2; exit 1; }
has(){ command -v "$1" >/dev/null 2>&1; }

need_tools(){
  local missing=()
  for t in "$@"; do has "$t" || missing+=("$t"); done
  ((${#missing[@]}==0)) || die "Missing tools: ${missing[*]}"
}

# ceil(x * frac); frac given as 0.30 style
ceil_frac(){
  # usage: ceil_frac <count> <fraction>
  awk -v n="$1" -v f="$2" '
    BEGIN{
      val = n*f;
      c = (val == int(val)) ? val : int(val)+1;
      if (c < 1) c=1;
      print c
    }'
}

# ---------------------- 1) clustalo batch ----------------------
# ### Generate clustalo alignments for batch of FASTA
# name: aln_clustalo
aln_clustalo(){
  local glob="${1:-*.fna}" outfmt="${2:-clu}"    # clu|fa|msf|phylip
  need_tools clustalo
  shopt -s nullglob
  for f in $glob; do
    local ref="${f%.fna}"
    echo -e "\n[clustalo] $ref"
    clustalo -i "$f" -o "${ref}.clust.${outfmt}" --outfmt="$outfmt" --force
  done
}

# ---------------------- 2) famsa batch ----------------------
# ### Generate famsa alignments for batch of FASTA
# name: aln_famsa
aln_famsa(){
  local glob="${1:-*.fna}" threads="${2:-2}"
  need_tools famsa
  shopt -s nullglob
  for f in $glob; do
    local ref="${f%.fna}"
    echo -e "\n[famsa] $ref"
    famsa -t "$threads" "$f" "${ref}.famsa.aln"
  done
}

# ---------------------- 3) consensus from clustalo ----------------------
# ### All in One: Get consensus sequences from batch of clustal omega alignments
# name: cons_clust
# Notes:
# - Produces a clustal-omega alignment in FASTA, then runs EMBOSS cons.
# - plurality can be integer (e.g., 5) OR a fraction (e.g., 0.30) of seq count.
# - setcase uses the same threshold number (upper vs lower case).
# - You can cut the alignment with: extra CONS args via --cons-args "-sbegin1 1 -send1 200"
cons_clust(){
  local glob="${1:-*.fna}"
  local plurality="${2:-0.30}"       # 0.30 means 30% of sequences; or give 5 for absolute
  local cons_args="${3:-}"           # optional: pass emboss cons extras: "-sbegin1 1 -send1 200"
  need_tools clustalo cons grep

  shopt -s nullglob
  for f in $glob; do
    local ref="${f%.fna}"
    echo -e "\n[cons_clust] $ref"

    # Count sequences in input FASTA
    local nseq
    nseq=$(grep -c '^>' "$f" || true)
    (( nseq > 0 )) || { echo "  skip: no sequences in $f"; continue; }

    # Derive threshold
    local thr
    if [[ "$plurality" == *.* ]]; then
      thr=$(ceil_frac "$nseq" "$plurality")
    else
      thr="$plurality"
    fi
    (( thr < 1 )) && thr=1

    # 1) Align with clustalo to FASTA
    local aln="${ref}.aln.fa"
    clustalo -i "$f" -o "$aln" --outfmt=fa --force

    # 2) Consensus with EMBOSS cons
    local out="${ref}.cons.fa"
    echo "  nseq=$nseq plurality(thr)=$thr -> $out"
    # Using same threshold for -plurality and -setcase to match your behavior
    cons -sequence "$aln" -outseq "$out" -plurality "$thr" -setcase "$thr" ${cons_args:+$cons_args}
  done
}

# ---------------------- 4) clustalo pairwise distances ----------------------
# ### Get pairwaise comparions with clustalo
# name: dist_clustalo
# If given a multi-FASTA, writes a pairwise distance matrix (percent ID).
dist_clustalo(){
  local infile="${1:?Give input FASTA}"
  local outfile="${2:-distances.txt}"
  need_tools clustalo
  echo "[dist] $infile -> $outfile"
  clustalo -i "$infile" --distmat-out="$outfile" --full --percent-id --force
}

# ---------------------- CLI wrapper ----------------------
usage(){
cat <<'EOF'
ULTRA alignment/consensus module

Functions (call after 'source ultra_align_cons.sh'):

  aln_clustalo [glob=*.fna] [outfmt=clu]
      # clustalo batch, output formats: clu|fa|msf|phylip

  aln_famsa   [glob=*.fna] [threads=2]
      # famsa batch, writes *.famsa.aln

  cons_clust  [glob=*.fna] [plurality=0.30] [cons_args=""]
      # clustalo -> EMBOSS cons to *.cons.fa
      # plurality: integer (e.g., 5) or fraction (0.30 = 30% of sequences)
      # cons_args: EMBOSS extra args, e.g. "-sbegin1 1 -send1 200"

  dist_clustalo <input.fasta> [distances.txt]
      # write pairwise percent-ID matrix

Examples:
  source ultra_align_cons.sh
  aln_clustalo "*.fna" fa
  aln_famsa "*.fna" 4
  cons_clust "*.fna" 0.30 "-sbegin1 1 -send1 300"
  dist_clustalo merged.fna matrix.tsv
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
fi
