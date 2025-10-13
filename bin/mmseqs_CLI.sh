#!/usr/bin/env bash

# Obj: mmseqs2 clustering CLI wrapper
# Author: Brandon D. Jordan
# Date: 2025-09-06

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  # Full clustering workflow -> creates: BASE/{db/,clustMode_X/}
  mmseqs_CLI.sh -f INPUT.fasta [-o BASE] [-m MODE] [-t THREADS] [--keep-tmp] [--dry-run]

  # Extract-only (NO clustering): find clustered_pairs.tsv and extract cluster
  mmseqs_CLI.sh --extract-only --rep-id SEQID -f INPUT.fasta
                [--pairs PATH | --find-pairs] [-o BASE] [--dry-run]

Required:
  -f, --fasta PATH     Input FASTA (original source sequences).

Options (clustering):
  -o, --outdir DIR     Base directory (default: .)
  -m, --mode INT       --cluster-mode 0..3 (default: 1)
  -t, --threads INT    Threads (default: all)
      --keep-tmp       Keep tmp dir
      --dry-run        Print commands; don't execute

Options (extract-only):
      --extract-only   Skip clustering; just extract the cluster for --rep-id.
      --rep-id ID      Representative ID to extract cluster for (required).
      --pairs PATH     Explicit clustered_pairs.tsv path (else use --find-pairs).
      --find-pairs     Recursively find newest clustered_pairs.tsv under BASE (or .).

Directory layout (clustering):
  BASE/
    db/
      modeX/
        inputDB/ clusteredDB/ repsDB/ tmp*/
    clustMode_X/
      reps.fasta
      clustered_pairs.tsv
      summary.num_clusters.txt
      summary.rep_counts.tsv

USAGE
}

# -------- defaults --------
FASTA=""
BASEDIR="."
MODE=1
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
KEEP_TMP=0
DRY=0
EXTRACT_ONLY=0

REP_ID=""
PAIRS_PATH=""
FIND_PAIRS=0

# -------- arg parsing --------
LONG_OPTS=fasta:,outdir:,mode:,threads:,keep-tmp,dry-run,extract-only,rep-id:,pairs:,find-pairs,help
SHORT_OPTS=f:o:m:t:h
PARSED=$(getopt -o "$SHORT_OPTS" --long "$LONG_OPTS" -n "$0" -- "$@") || { usage; exit 2; }
eval set -- "$PARSED"
while true; do
  case "$1" in
    -f|--fasta)   FASTA="$2"; shift 2 ;;
    -o|--outdir)  BASEDIR="$2"; shift 2 ;;
    -m|--mode)    MODE="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    --keep-tmp)   KEEP_TMP=1; shift ;;
    --dry-run)    DRY=1; shift ;;
    --extract-only) EXTRACT_ONLY=1; shift ;;
    --rep-id)     REP_ID="$2"; shift 2 ;;
    --pairs)      PAIRS_PATH="$2"; shift 2 ;;
    --find-pairs) FIND_PAIRS=1; shift ;;
    -h|--help)    usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Internal parsing error"; exit 2 ;;
  esac
done

# -------- helpers --------
run() { echo "+ $*"; [[ "$DRY" -eq 0 ]] && "$@"; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH."; exit 1; }; }
find_latest_pairs() {
  local root="$1"
  find "$root" -type f -name 'clustered_pairs.tsv' -printf '%T@ %p\n' 2>/dev/null \
    | sort -nr | awk 'NR==1{ $1=""; sub(/^ /,""); print }'
}

# -------- checks --------
[[ -z "$FASTA" ]] && { echo "ERROR: --fasta is required."; usage; exit 1; }
require_cmd awk
require_cmd seqkit
require_cmd mmseqs

if [[ "$EXTRACT_ONLY" -eq 1 ]]; then
  [[ -z "$REP_ID" ]] && { echo "ERROR: --extract-only requires --rep-id."; exit 1; }
  if [[ -z "$PAIRS_PATH" && "$FIND_PAIRS" -eq 0 ]]; then
    echo "ERROR: --extract-only requires --pairs PATH or --find-pairs."; exit 1
  fi
else
  # Full clustering path must NOT accept --rep-id
  [[ -n "$REP_ID" ]] && { echo "ERROR: --rep-id is only valid with --extract-only."; exit 1; }
  [[ "$MODE" =~ ^[0-3]$ ]] || { echo "ERROR: --mode must be 0..3 (got $MODE)"; exit 1; }
  [[ "$THREADS" =~ ^[0-9]+$ ]] || { echo "ERROR: --threads must be integer"; exit 1; }
fi

# ===========================
# EXTRACT-ONLY MODE
# ===========================
if [[ "$EXTRACT_ONLY" -eq 1 ]]; then
  SEARCH_ROOT="$BASEDIR"
  if [[ -n "$PAIRS_PATH" ]]; then
    [[ -f "$PAIRS_PATH" ]] || { echo "ERROR: --pairs not found: $PAIRS_PATH"; exit 1; }
  else
    echo "Searching for clustered_pairs.tsv under: $SEARCH_ROOT"
    PAIRS_PATH="$(find_latest_pairs "$SEARCH_ROOT")" || true
    [[ -n "$PAIRS_PATH" ]] || { echo "ERROR: No clustered_pairs.tsv found."; exit 1; }
    echo "Using newest: $PAIRS_PATH"
  fi

  IDS_FILE="$(mktemp -u ./tmp_ids.XXXXXX.txt)"
  OUT_DIR="$(dirname "$PAIRS_PATH")"
  OUT_FASTA="$OUT_DIR/clust_${REP_ID}_seqs.fasta"

  run bash -c "awk -v rep='$REP_ID' '\$1==rep {print \$2}' '$PAIRS_PATH' > '$IDS_FILE'"
  run bash -c "echo '$REP_ID' >> '$IDS_FILE'"
  run seqkit grep -f "$IDS_FILE" "$FASTA" > "$OUT_FASTA"
  [[ "$DRY" -eq 0 ]] && rm -f "$IDS_FILE"

  echo "==> Extracted cluster for $REP_ID"
  echo "    Pairs:   $PAIRS_PATH"
  echo "    FASTA:   $FASTA"
  echo "    Output:  $OUT_FASTA"
  exit 0
fi

# ===========================
# CLUSTERING MODE
# ===========================
DB_ROOT="$BASEDIR/db/mode${MODE}"
RES_ROOT="$BASEDIR/clustMode_${MODE}"
mkdir -p "$DB_ROOT" "$RES_ROOT"

INPUT_DB="$DB_ROOT/inputDB"
CLUST_DB="$DB_ROOT/clusteredDB"
REPS_DB="$DB_ROOT/repsDB"
TMPDIR="$(mktemp -d "$DB_ROOT/tmp.XXXXXX")"

PAIRS_TSV="$RES_ROOT/clustered_pairs.tsv"
REPS_FASTA="$RES_ROOT/reps.fasta"
NUM_FILE="$RES_ROOT/summary.num_clusters.txt"
COUNTS_TSV="$RES_ROOT/summary.rep_counts.tsv"

echo "==> Base:        $BASEDIR"
echo "==> DB dir:      $DB_ROOT"
echo "==> Results dir: $RES_ROOT"
echo "==> Mode:        $MODE"
echo "==> Threads:     $THREADS"

if [[ ! -d "$INPUT_DB" ]]; then
  run mmseqs createdb "$FASTA" "$INPUT_DB"
else
  echo "Skipping createdb: $INPUT_DB exists."
fi

run mmseqs cluster "$INPUT_DB" "$CLUST_DB" "$TMPDIR" --cluster-mode "$MODE" --threads "$THREADS"
run mmseqs createtsv "$INPUT_DB" "$INPUT_DB" "$CLUST_DB" "$PAIRS_TSV"
run mmseqs createsubdb "$CLUST_DB" "$INPUT_DB" "$REPS_DB"
run mmseqs convert2fasta "$REPS_DB" "$REPS_FASTA"

run bash -c "grep -c '^>' '$REPS_FASTA' > '$NUM_FILE'"
run bash -c "awk '{counts[\$1]++} END{for (r in counts) print r\"\t\"counts[r]}' '$PAIRS_TSV' \
  | sort -k1,1 > '$COUNTS_TSV'"

if [[ "$KEEP_TMP" -eq 0 ]]; then
  run rm -rf "$TMPDIR"
else
  echo "Keeping tmp: $TMPDIR"
fi

echo "==> Done."
printf "  • %s\n" \
  "$DB_ROOT" "$RES_ROOT" \
  "$INPUT_DB" "$CLUST_DB" "$REPS_DB" \
  "$REPS_FASTA" "$PAIRS_TSV" "$NUM_FILE" "$COUNTS_TSV"

