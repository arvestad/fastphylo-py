"""
Grid-size sensitivity: bias vs true distance for different n_points values.

Usage: plot_npts_sensitivity.py <npts_N.tsv> ... <output.pdf>

Each positional argument except the last is a TSV produced by estimate_distances
with method=expected and a specific n_points value.  The n_points value is read
from the filename (expected_npts{N}_all.tsv → N).
"""
import sys
import os
import re
import math
import numpy as np


def _load(path, col="estimated_distance"):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    df = df.dropna(subset=["true_distance", col])
    df = df[df[col].apply(lambda v: math.isfinite(v))]
    return df


def _rolling_bias(true, est, n_bins=30):
    edges = np.percentile(true, np.linspace(0, 100, n_bins + 1))
    cx, cy = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (true >= lo) & (true < hi)
        if mask.sum() >= 3:
            cx.append(float(np.median(true[mask])))
            cy.append(float(np.mean(est[mask] - true[mask])))
    return np.array(cx), np.array(cy)


def _npts_from_path(path):
    m = re.search(r"npts(\d+)", os.path.basename(path))
    return int(m.group(1)) if m else None


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <npts_N.tsv> ... <output.pdf>")
        sys.exit(1)

    *tsv_paths, out_pdf = sys.argv[1:]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig, ax = plt.subplots(figsize=(7, 5))

    colors = cm.viridis(np.linspace(0.1, 0.9, len(tsv_paths)))

    for path, color in zip(sorted(tsv_paths, key=_npts_from_path), colors):
        n = _npts_from_path(path)
        df = _load(path)
        true = df["true_distance"].values
        est  = df["estimated_distance"].values
        cx, cy = _rolling_bias(true, est)
        ax.plot(cx, cy, lw=1.5, color=color, label=f"n={n}")

    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("True distance")
    ax.set_ylabel("Mean bias (estimated − true)")
    ax.set_title("WAG — expected distance: grid-size sensitivity")
    ax.legend(title="n_points", fontsize=8)

    os.makedirs(os.path.dirname(out_pdf) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
