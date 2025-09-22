#!/usr/bin/env python3
"""
ultra_plot_hist (CLI)

- periods mode: histogram of period sizes (column 'period').
- arrays mode: violin plots of array sizes (len/1000 kb) per chromosome,
  filtered by period -x, with medians, IQR bars, and Tukey whiskers.

Examples
--------
# Period-size histogram (all periods)
python ultra_plot_hist.py -f repeats.tsv -m periods -b 100 --out periods_hist.png

# Period-size histogram for 155–156 bp monomers (log x)
python ultra_plot_hist.py -f repeats.tsv -m periods -x 155-156 -b 60 --logx --out periods_155_156.png

# Arrays: violins per chromosome (165 bp), colored by cmap, y in bp at 1000-bp steps
python ultra_plot_hist.py -f repeats.tsv -m arrays -x 165 --cmap tab20 --yunits bp --ytick-step 1000

# Same, but single color and full labels
python ultra_plot_hist.py -f repeats.tsv -m arrays -x 165 --cmap none --facecolor "#1976d2" --full-labels

# Custom spaced positions (must match #chromosomes)
python ultra_plot_hist.py -f repeats.tsv -m arrays -x 165 --positions 1,5,7,9,12
"""

import argparse
import re
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FuncFormatter, MaxNLocator

# ---------- Global styling ----------
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['font.family']  = 'Arial'


# ---------- Parsing & I/O ----------
def parse_x_spec(x_spec: Optional[str]) -> Optional[Tuple[int, int]]:
    """Parse -x as N or L-U. Returns (lb, ub) or None."""
    if not x_spec:
        return None
    s = re.sub(r"\s+", "", str(x_spec))
    if re.fullmatch(r"\d+-\d+", s):
        l, u = s.split("-")
        lb, ub = int(l), int(u)
        if lb > ub:
            lb, ub = ub, lb
        return lb, ub
    if re.fullmatch(r"\d+", s):
        n = int(s)
        return n, n
    raise ValueError(f"-x must be N or L-U (e.g., 165 or 155-156). Got: {x_spec!r}")


def parse_positions(s: Optional[str]) -> Optional[List[float]]:
    """Parse --positions '1,5,7' -> [1.0, 5.0, 7.0]."""
    if not s:
        return None
    try:
        return [float(tok) for tok in s.split(",") if tok.strip() != ""]
    except Exception as e:
        raise ValueError(f"Could not parse --positions: {s!r}") from e


def read_ultra_file(file: str) -> pd.DataFrame:
    # Use exactly your signature:
    df = pd.read_csv(
        file,
        sep="\t",
        names=['chr', 'pos', 'len', 'period', 'score', 'sub', 'ins', 'del', 'cons', 'seq'],
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    )
    return df


# ---------- Helpers ----------
def apply_period_filter(df: pd.DataFrame, lbub: Optional[Tuple[int, int]]) -> pd.DataFrame:
    if lbub is None:
        return df
    lb, ub = lbub
    return df[(df['period'] >= lb) & (df['period'] <= ub)].copy()


def natural_chr_sort(chr_values: List[str]) -> List[str]:
    """Sort chr labels 'chr1', 'chr2', ..., 'chrX' in a natural-ish way."""
    def keyfun(s: str):
        s = str(s)
        m = re.search(r"(\d+)$", s)
        if m:
            return (0, int(m.group(1)))
        special = {"X": 1, "Y": 2, "M": 3, "MT": 3}
        tag = s.replace("chr", "").upper()
        return (1, special.get(tag, 999), tag)
    return sorted(chr_values, key=keyfun)


def tukey_whiskers(q1: float, q3: float, data: np.ndarray) -> tuple[float, float]:
    """Return (lower, upper) Tukey bounds clipped to data range."""
    iqr = q3 - q1
    upper = np.clip(q3 + 1.5 * iqr, q3, np.max(data))
    lower = np.clip(q1 - 1.5 * iqr, np.min(data), q1)
    return float(lower), float(upper)


def set_length_ticks(ax, yunits: str = "kb", step: float | None = None):
    """
    Configure y-axis ticks for length in kb or bp on a LINEAR axis.
    The plotted data are in kb; for bp labels we reformat tick labels.
    """
    # Default step if not provided
    if step is None:
        step = 1.0 if yunits == "kb" else 1000.0

    ax.set_yscale("linear")

    if yunits == "kb":
        # major ticks every 'step' kb; labels as integers
        ax.yaxis.set_major_locator(MultipleLocator(step))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(round(v))}"))
        ax.set_ylabel("Array length (kb)")
        if abs(step - 1.0) < 1e-9:
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    else:  # "bp"
        step_kb = step / 1000.0  # convert bp step to kb spacing
        ax.yaxis.set_major_locator(MultipleLocator(step_kb))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(round(v * 1000))}"))
        ax.set_ylabel("Array length (bp)")


def shorten_after_last_dot(labels: List[str]) -> List[str]:
    """Keep only the part after the final '.' in each label."""
    out = []
    for s in labels:
        s = str(s)
        if "." in s:
            out.append(s.rsplit(".", 1)[-1])
        else:
            out.append(s)
    return out


def colors_from_cmap(n: int, cmap_name: str = "tab20") -> List:
    """Return n RGBA colors sampled evenly from the given colormap."""
    if cmap_name.lower() == "none":
        return []
    cmap = plt.get_cmap(cmap_name)
    if n == 1:
        return [cmap(0.5)]
    xs = np.linspace(0.0, 1.0, n, endpoint=True)
    return [cmap(x) for x in xs]


# ---------- Plotting ----------
def plot_periods_histogram(df: pd.DataFrame, bins: int, lbub: Optional[Tuple[int, int]],
                           out: Optional[str], logx: bool) -> None:
    vals = df['period'].dropna().astype(float).values
    if vals.size == 0:
        raise SystemExit("No data to plot after filtering periods.")

    plt.figure(figsize=(10, 4))
    plt.hist(vals, bins=bins, density=False)
    if logx:
        plt.xscale('log')
        plt.xlabel("Period size (bp) [log scale]")
    else:
        plt.xlabel("Period size (bp)")
    sel = ""
    if lbub:
        lb, ub = lbub
        sel = f"  [{lb}-{ub}]" if lb != ub else f"  [{lb}]"
    plt.ylabel("Count")
    plt.title(f"Monomer/Period Size Distribution{sel}")
    plt.tight_layout()

    if out:
        plt.savefig(out, dpi=300)
    else:
        plt.show()


def plot_arrays_violin_by_chr(
    df: pd.DataFrame,
    lbub: Tuple[int, int],
    out: Optional[str],
    logy: bool,  # kept for backward compatibility; ignored when yunits specified
    positions: Optional[List[float]] = None,
    facecolor: str = "red",
    edgecolor: str = "black",
    yunits: str = "kb",
    ytick_step: float | None = None,
    cmap_name: str = "tab20",
    full_labels: bool = False,
) -> None:
    if lbub is None:
        raise SystemExit("Error: -x is required for -m arrays (use N or L-U).")

    # Filter on period
    dff = apply_period_filter(df, lbub)
    if dff.empty:
        raise SystemExit("No rows matched the requested period filter for arrays.")

    # Array size in kb
    dff = dff.assign(array_kb=dff['len'].astype(float) / 1000.0)

    # Prepare groups by chromosome with natural sorting
    chrs = natural_chr_sort(dff['chr'].astype(str).unique().tolist())
    dataset = [dff.loc[dff['chr'].astype(str) == c, 'array_kb'].values for c in chrs]

    # Positions: user-specified or sequential 1..N
    if positions is None:
        positions = list(np.arange(1, len(chrs) + 1, dtype=float))
    else:
        if len(positions) != len(chrs):
            raise SystemExit(
                f"--positions has {len(positions)} values but there are {len(chrs)} chromosomes. "
                "They must match so each violin has a location."
            )

    # X labels: short or full
    xtick_labels = chrs if full_labels else shorten_after_last_dot(chrs)

    fig, ax = plt.subplots(figsize=(max(10, min(2 + 0.4 * len(chrs), 20)), 6))

    # Violin plot (keep extrema visible for styling)
    vp = ax.violinplot(dataset=dataset, positions=positions, showmeans=False, showextrema=True, widths=0.9)

    # Colors: per-violin from cmap (or single facecolor if cmap=none)
    colors = colors_from_cmap(len(dataset), cmap_name=cmap_name)
    for i, body in enumerate(vp['bodies']):
        if colors:
            body.set_facecolor(colors[i])
        else:
            body.set_facecolor(facecolor)
        body.set_edgecolor(edgecolor)
        body.set_alpha(1.0)

    # Style extrema bars
    vp['cmaxes'].set_color('black')
    vp['cmins'].set_color('black')
    vp['cbars'].set_color('black')

    # Quantiles and Tukey whiskers
    q = [np.percentile(arr, [25, 50, 75]) for arr in dataset]
    whiskers = [tukey_whiskers(qq[0], qq[2], arr) for qq, arr in zip(q, dataset)]

    # Median markers (white dots)
    ax.scatter(
        positions, [qq[1] for qq in q],
        marker='o', color='white', s=30, zorder=3, edgecolor='black', linewidth=0.6
    )

    # IQR bars
    ax.vlines(
        positions, [qq[0] for qq in q], [qq[2] for qq in q],
        color='black', linestyle='-', lw=5, zorder=2
    )

    # Whiskers
    ax.vlines(
        positions, [w[0] for w in whiskers], [w[1] for w in whiskers],
        color='black', linestyle='-', lw=2, zorder=2
    )

    # Axes/labels
    ax.set_xticks(positions)
    ax.set_xticklabels(xtick_labels, rotation=45, ha='right')
    ax.set_xlabel("Chromosome")

    # Linear, unit-aware ticks (overrides logy for clarity)
    set_length_ticks(ax, yunits=yunits, step=ytick_step)

    lb, ub = lbub
    sel = f"period {lb}-{ub} bp" if lb != ub else f"period {lb} bp"
    ax.set_title(f"Array Size by Chromosome (violins with median/IQR/whiskers), {sel}")

    fig.tight_layout()
    if out:
        fig.savefig(out, dpi=300)
    else:
        plt.show()


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="ULTRA_plot_hist (Python CLI)")
    ap.add_argument("-f", "--file", required=True, help="Input ULTRA TSV file")
    ap.add_argument("-m", "--mode", choices=["periods", "arrays"], default="periods",
                    help="Plot mode: periods (histogram) or arrays (violins per chromosome)")
    ap.add_argument("-x", "--x-spec",
                    help="N or L-U period filter (e.g. 165 or 155-156). Required for arrays.")
    ap.add_argument("-b", "--bins", type=int, default=100,
                    help="Bins for histogram (periods mode). Default 100.")
    ap.add_argument("--out", help="Output image path. If omitted, opens an interactive window.")
    ap.add_argument("--logx", action="store_true", help="Use log-scale on x-axis for periods histogram.")
    ap.add_argument("--logy", action="store_true",
                    help="(Deprecated for arrays when using unit ticks) Log-scale on y-axis.")
    # Existing customization flags
    ap.add_argument("--positions", help="Comma-separated x-positions for violins (must match #chromosomes).")
    ap.add_argument("--facecolor", default="red", help="Violin face color if --cmap none (e.g., 'red', '#1976d2').")
    ap.add_argument("--yunits", choices=["kb", "bp"], default="kb",
                    help="Y-axis units for arrays plot (default: kb).")
    ap.add_argument("--ytick-step", type=float, default=None,
                    help="Tick step (default: 1 for kb, 1000 for bp).")
    # NEW: colormap & label-shortening controls
    ap.add_argument("--cmap", default="tab20",
                    help="Matplotlib colormap name for per-violin colors (e.g., tab20, rainbow, viridis). "
                         "Use 'none' for a single facecolor.")
    ap.add_argument("--full-labels", action="store_true",
                    help="Use the full chromosome labels (disable shortening after last '.').")
    args = ap.parse_args()

    lbub = parse_x_spec(args.x_spec) if args.x_spec else None
    df = read_ultra_file(args.file)

    if args.mode == "periods":
        plot_periods_histogram(
            df if lbub is None else apply_period_filter(df, lbub),
            bins=args.bins, lbub=lbub, out=args.out, logx=args.logx
        )
    elif args.mode == "arrays":
        if lbub is None:
            raise SystemExit("Error: -x is required for -m arrays (use N or L-U).")
        pos = parse_positions(args.positions)
        plot_arrays_violin_by_chr(
            df, lbub, out=args.out,
            logy=False,  # ignore logy to ensure clean unit ticks
            positions=pos,
            facecolor=args.facecolor,
            yunits=args.yunits,
            ytick_step=args.ytick_step,
            cmap_name=args.cmap,
            full_labels=args.full_labels,
        )
    else:
        raise SystemExit(f"Unknown mode: {args.mode}")


if __name__ == "__main__":
    main()

