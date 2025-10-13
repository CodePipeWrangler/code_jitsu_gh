#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Modular centromeric repeat filtering and visualization.

Features:
- Prefilter ULTRA TSVs by chromosome/position/period (replaces AWK if desired)
- Representative sequence selection (auto or manual)
- Feature extraction (GC, entropy, indel variability, distance metrics)
- Alignment scoring vs representative
- Outlier detection (IsolationForest, LOF)
- Optional t-SNE embeddings & plots
- Filtering + reporting
- FASTA export from TSVs (Consensus -> sequence)

Input ULTRA TSV is expected as tab-delimited with the following columns:
0: SequenceName
1: Start
2: Length
3: Period
4: Score
5: Substitutions
6: Insertions
7: Deletions
8: Consensus
9: Sequence

Author: you
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

# Optional plotting imports are guarded inside functions so the script works headless.


# ----------------------------- I/O & schema ----------------------------- #

ULTRA_COLS = [
    "SequenceName",
    "Start",
    "Length",
    "Period",
    "Score",
    "Substitutions",
    "Insertions",
    "Deletions",
    "Consensus",
    "Sequence",
]


def read_ultra_tsv(path: Path) -> pd.DataFrame:
    """Read ULTRA-format TSV (no header) into a DataFrame with named columns."""
    df = pd.read_csv(path, sep="\t", header=None, names=ULTRA_COLS, usecols=range(10))
    # Basic types
    for c in ["Start", "Length", "Period", "Score", "Substitutions", "Insertions", "Deletions"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


# ----------------------------- Prefilter ----------------------------- #

def prefilter_df(
    df: pd.DataFrame,
    chrom_substr: Optional[str],
    start: Optional[int],
    end: Optional[int],
    period: Optional[int],
) -> pd.DataFrame:
    """
    Filter by:
      - chrom_substr: substring to match in SequenceName (e.g., 'Gm15' or '.Gm15')
      - start..end   : inclusive genomic start window (compares to 'Start')
      - period       : exact match on Period
    """
    m = pd.Series(True, index=df.index)
    if chrom_substr:
        m &= df["SequenceName"].astype(str).str.contains(chrom_substr, regex=False)
    if start is not None:
        m &= df["Start"] >= start
    if end is not None:
        m &= df["Start"] <= end
    if period is not None:
        m &= df["Period"] == period
    return df.loc[m].copy()


# ----------------------------- Representative ----------------------------- #

def pick_representative_by_similarity(consensus_list: List[str]) -> Tuple[int, str]:
    """
    Choose the sequence with highest mean global alignment similarity to all others.
    Returns (index, representative_seq).
    """
    if not consensus_list:
        raise ValueError("No sequences provided to choose representative.")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    n = len(consensus_list)
    score_matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        si = consensus_list[i]
        for j in range(i + 1, n):
            sj = consensus_list[j]
            s = aligner.score(si, sj) / max(len(si), len(sj))
            score_matrix[i, j] = s
            score_matrix[j, i] = s

    mean_scores = score_matrix.mean(axis=1)
    idx = int(np.argmax(mean_scores))
    return idx, consensus_list[idx]


# ----------------------------- Features ----------------------------- #

def entropy(seq: str) -> float:
    """Shannon entropy (bits) using symbol frequencies in seq."""
    if not seq:
        return 0.0
    # Faster than set(...) in tight loops
    counts = {}
    for ch in seq:
        counts[ch] = counts.get(ch, 0) + 1
    ent = 0.0
    n = len(seq)
    for c in counts.values():
        p = c / n
        ent -= p * math.log2(p)
    return ent


def add_features(df: pd.DataFrame, centromere_midpoint: Optional[int]) -> pd.DataFrame:
    out = df.copy()
    out["GC_Content"] = out["Consensus"].str.count("G") + out["Consensus"].str.count("C")
    out["GC_Content"] = out["GC_Content"] / out["Consensus"].str.len().clip(lower=1)

    out["Entropy"] = out["Consensus"].apply(entropy)
    out["Indel_Variability"] = out["Substitutions"].fillna(0) + out["Insertions"].fillna(0) + out["Deletions"].fillna(0)

    if centromere_midpoint is not None:
        out["DistanceFromCentromere"] = (out["Start"] - centromere_midpoint).abs()
        out["NormalizedDistance"] = out["DistanceFromCentromere"] / out["Length"].clip(lower=1)
    else:
        out["DistanceFromCentromere"] = np.nan
        out["NormalizedDistance"] = np.nan

    return out


# ----------------------------- Alignment scoring ----------------------------- #

def add_alignment_scores(
    df: pd.DataFrame,
    rep_seq: str,
    repeat_extend: int = 1,
) -> pd.DataFrame:
    """
    Score alignment of each Consensus to representative.
    Optionally "extend" both by repeating N times (helps w/ arrays).
    """
    if repeat_extend < 1:
        repeat_extend = 1

    rep = rep_seq * repeat_extend

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0

    scores = []
    for s in df["Consensus"].astype(str):
        q = s * repeat_extend
        sc = aligner.score(rep, q) / max(len(rep), len(q))
        scores.append(sc)

    out = df.copy()
    out["AlignmentScore"] = scores
    return out


# ----------------------------- Outliers ----------------------------- #

@dataclass
class OutlierParams:
    contamination: float = 0.2
    lof_neighbors: int = 10
    random_state: int = 42


def detect_outliers(
    df: pd.DataFrame,
    feature_cols_basic: List[str],
    feature_cols_enhanced: List[str],
    params: OutlierParams,
) -> pd.DataFrame:
    out = df.copy()

    X_basic = out[feature_cols_basic].to_numpy()
    X_enh = out[feature_cols_enhanced].to_numpy()

    if_model_basic = IsolationForest(
        contamination=params.contamination, random_state=params.random_state
    ).fit(X_basic)
    out["IF_Basic"] = if_model_basic.predict(X_basic)  # 1=inlier, -1=outlier

    lof_basic = LocalOutlierFactor(
        n_neighbors=params.lof_neighbors, contamination=params.contamination
    )
    out["LOF_Basic"] = lof_basic.fit_predict(X_basic)

    if_model_enh = IsolationForest(
        contamination=params.contamination, random_state=params.random_state
    ).fit(X_enh)
    out["IF_Enhanced"] = if_model_enh.predict(X_enh)

    lof_enh = LocalOutlierFactor(
        n_neighbors=params.lof_neighbors, contamination=params.contamination
    )
    out["LOF_Enhanced"] = lof_enh.fit_predict(X_enh)

    return out


# ----------------------------- Filtering & reporting ----------------------------- #

def filter_rows(
    df: pd.DataFrame,
    min_align: float = 0.90,
    use_enhanced: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Keep rows that are NOT outliers and pass alignment threshold.
    """
    if use_enhanced:
        if_mask = df["IF_Enhanced"] == 1
        lof_mask = df["LOF_Enhanced"] == 1
    else:
        if_mask = df["IF_Basic"] == 1
        lof_mask = df["LOF_Basic"] == 1

    align_mask = df["AlignmentScore"] >= min_align
    keep = df[if_mask & lof_mask & align_mask].copy()
    drop = df.loc[~keep.index.isin(df.index[if_mask & lof_mask & align_mask])].copy()
    return keep, drop


def print_breakdown(df: pd.DataFrame, kept: pd.DataFrame, min_align: float) -> None:
    print(f"Total rows: {len(df)}")
    print(f"Filtered rows retained: {len(kept)}")
    print("Removed as IsolationForest (enhanced) outlier:", (df["IF_Enhanced"] == -1).sum())
    print("Removed as LOF (enhanced) outlier:", (df["LOF_Enhanced"] == -1).sum())
    print(f"Removed for low alignment (< {min_align:.2f}):", (df["AlignmentScore"] < min_align).sum())


# ----------------------------- Embeddings & plots ----------------------------- #

def embed_and_plot(
    features_basic: np.ndarray,
    features_enh: np.ndarray,
    labels_if_basic: np.ndarray,
    labels_if_enh: np.ndarray,
    labels_lof_enh: np.ndarray,
    align_scores: np.ndarray,
    out_prefix: Path,
) -> None:
    """
    Create side-by-side plots similar to your draft. Saved to disk.
    """
    # Lazy import to avoid hard dependency when plotting is disabled.
    import matplotlib.pyplot as plt
    from sklearn.manifold import TSNE

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    tsne = TSNE(n_components=2, random_state=42)
    emb_basic = tsne.fit_transform(features_basic)
    emb_enh = tsne.fit_transform(features_enh)

    # IF basic
    plt.figure()
    plt.scatter(emb_basic[:, 0], emb_basic[:, 1], c=(labels_if_basic == -1), cmap="coolwarm", s=10)
    plt.title("Isolation Forest (Basic Features)")
    plt.tight_layout()
    plt.savefig(out_prefix.with_name(out_prefix.stem + "_IF_basic.png"))
    plt.close()

    # IF enhanced
    plt.figure()
    plt.scatter(emb_enh[:, 0], emb_enh[:, 1], c=(labels_if_enh == -1), cmap="coolwarm", s=10)
    plt.title("Isolation Forest (Enhanced Features)")
    plt.tight_layout()
    plt.savefig(out_prefix.with_name(out_prefix.stem + "_IF_enhanced.png"))
    plt.close()

    # LOF enhanced
    plt.figure()
    plt.scatter(emb_enh[:, 0], emb_enh[:, 1], c=(labels_lof_enh == -1), cmap="coolwarm", s=10)
    plt.title("LOF (Enhanced Features)")
    plt.tight_layout()
    plt.savefig(out_prefix.with_name(out_prefix.stem + "_LOF_enhanced.png"))
    plt.close()

    # Alignment gradients
    for name, emb in [("basic", emb_basic), ("enhanced", emb_enh)]:
        plt.figure()
        plt.scatter(emb[:, 0], emb[:, 1], c=align_scores, cmap="viridis", s=10)
        plt.title(f"Alignment Score Gradient ({name})")
        plt.colorbar(label="AlignmentScore")
        plt.tight_layout()
        plt.savefig(out_prefix.with_name(out_prefix.stem + f"_align_{name}.png"))
        plt.close()


# ----------------------------- FASTA export ----------------------------- #

def tsv_to_fasta(
    tsv: Path,
    fasta_out: Path,
    name_cols: Tuple[str, ...] = ("SequenceName", "Start", "Length", "Period"),
    seq_col: str = "Consensus",
) -> None:
    """
    Convert a TSV (ULTRA-like) to FASTA using selected columns in the header.
    """
    df = read_ultra_tsv(tsv)
    with fasta_out.open("w") as fh:
        for _, row in df.iterrows():
            hdr = "_".join(str(row[c]) for c in name_cols)
            seq = str(row[seq_col])
            fh.write(f">{hdr}\n{seq}\n")


# ----------------------------- CLI ----------------------------- #

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Modular centromeric repeat processing pipeline (supervised/stepwise)."
    )
    p.add_argument("-i", "--input", type=Path, required=True, help="Input ULTRA TSV")
    p.add_argument("-o", "--outdir", type=Path, default=Path("results"), help="Output directory")

    # Prefilter (Python alternative to your AWK lines)
    p.add_argument("--prefilter-chr", type=str, default=None, help="Substring to match in SequenceName (e.g., '.Gm15')")
    p.add_argument("--prefilter-start", type=int, default=None, help="Start ≥ this position")
    p.add_argument("--prefilter-end", type=int, default=None, help="Start ≤ this position")
    p.add_argument("--prefilter-period", type=int, default=None, help="Exact Period (e.g., 91)")

    # Representative
    g = p.add_mutually_exclusive_group()
    g.add_argument("--rep-seq", type=str, default=None, help="Manual representative sequence")
    g.add_argument("--rep-auto", action="store_true", help="Pick representative by mean similarity")

    p.add_argument("--repeat-extend", type=int, default=3, help="Repeat factor for array-like alignment (default: 3)")

    # Features
    p.add_argument("--centromere-mid", type=int, default=None, help="Centromere midpoint for distance features")

    # Outliers
    p.add_argument("--contamination", type=float, default=0.2, help="Outlier contamination (0..0.5)")
    p.add_argument("--lof-nn", type=int, default=10, help="LOF n_neighbors")
    p.add_argument("--min-align", type=float, default=0.90, help="Min alignment score to keep")
    p.add_argument("--basic-only", action="store_true", help="Use basic features instead of enhanced")

    # Plots
    p.add_argument("--plots", action="store_true", help="Generate t-SNE plots")

    # FASTA export (standalone)
    p.add_argument("--fasta-from", type=Path, default=None, help="TSV to convert to FASTA (skips rest)")
    p.add_argument("--fasta-out", type=Path, default=None, help="Output FASTA path (used with --fasta-from)")

    return p


def main():
    args = build_arg_parser().parse_args()

    # Standalone FASTA mode
    if args.fasta_from is not None:
        out = args.fasta_out or args.fasta_from.with_suffix(".fna")
        tsv_to_fasta(args.fasta_from, out)
        print(f"FASTA written: {out}")
        return

    outdir: Path = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Load
    df = read_ultra_tsv(args.input)

    # Optional prefilter (Python replacement for AWK lines)
    df_pref = prefilter_df(
        df,
        chrom_substr=args.prefilter_chr,
        start=args.prefilter_start,
        end=args.prefilter_end,
        period=args.prefilter_period,
    )
    pre_name = "prefiltered.tsv" if any(
        x is not None for x in (args.prefilter_chr, args.prefilter_start, args.prefilter_end, args.prefilter_period)
    ) else "input_clean.tsv"
    write_tsv(df_pref, outdir / pre_name)

    # Representative
    if args.rep_seq:
        rep = args.rep_seq
    elif args.rep_auto:
        _, rep = pick_representative_by_similarity(df_pref["Consensus"].astype(str).tolist())
    else:
        raise SystemExit(
            "You must specify either --rep-seq (manual) or --rep-auto (automatic representative)."
        )

    # Features
    df_feat = add_features(df_pref, args.centromere_mid)

    # Alignment
    df_align = add_alignment_scores(df_feat, rep_seq=rep, repeat_extend=args.repeat_extend)

    # Outliers (basic vs enhanced)
    basic_cols = ["GC_Content", "Entropy", "Indel_Variability"]
    enh_cols = basic_cols + ["NormalizedDistance", "AlignmentScore"]

    params = OutlierParams(
        contamination=args.contamination,
        lof_neighbors=args.lof_nn,
        random_state=42,
    )
    df_out = detect_outliers(df_align, basic_cols, enh_cols, params)

    # Optional plots
    if args.plots:
        try:
            X_basic = df_out[basic_cols].to_numpy()
            X_enh = df_out[enh_cols].to_numpy()
            embed_and_plot(
                X_basic,
                X_enh,
                df_out["IF_Basic"].to_numpy(),
                df_out["IF_Enhanced"].to_numpy(),
                df_out["LOF_Enhanced"].to_numpy(),
                df_out["AlignmentScore"].to_numpy(),
                out_prefix=outdir / "tsne",
            )
            print(f"t-SNE plots written to: {outdir}")
        except Exception as e:
            print(f"[warn] Plotting failed (skipping): {e}")

    # Filter + save
    kept, dropped = filter_rows(df_out, min_align=args.min_align, use_enhanced=not args.basic_only)
    write_tsv(df_out, outdir / "scored_full.tsv")
    write_tsv(kept, outdir / "kept.tsv")
    write_tsv(dropped, outdir / "outliers.tsv")

    # Report
    print_breakdown(df_out, kept, args.min_align)
    print(f"Representative length (post-repeat x{args.repeat_extend}): {len(rep) * args.repeat_extend if args.repeat_extend>1 else len(rep)}")
    print(f"All outputs in: {outdir.resolve()}")

if __name__ == "__main__":
    main()
