"""
Produce scatter and error plots from an aggregated TSV.

Usage:
  plot.py scatter <all.tsv> <output.pdf> <model_name>
  plot.py error   <all.tsv> <output.pdf> <model_name>
"""
import sys
import math
import os


def _load(tsv_path):
    import pandas as pd
    df = pd.read_csv(tsv_path, sep="\t")
    df = df.dropna(subset=["true_distance", "estimated_distance"])
    df = df[df["estimated_distance"].apply(lambda v: math.isfinite(v))]
    return df


def _pearson_r(x, y):
    import numpy as np
    if len(x) < 2:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _rolling_mean(x, y, n_bins=30):
    """Bin x into equal-width bins and return (bin_center, mean_y) arrays."""
    import numpy as np
    x = np.asarray(x)
    y = np.asarray(y)
    xmin, xmax = x.min(), x.max()
    bins = np.linspace(xmin, xmax, n_bins + 1)
    centers, means = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (x >= lo) & (x < hi)
        if mask.sum() > 0:
            centers.append(0.5 * (lo + hi))
            means.append(float(y[mask].mean()))
    return centers, means


def plot_scatter(tsv_path, out_pdf, model_name):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = _load(tsv_path)
    x  = df["true_distance"].values
    y  = df["estimated_distance"].values
    r  = _pearson_r(x, y)

    fig, ax = plt.subplots(figsize=(5, 5))
    if len(x) == 0:
        ax.set_title(f"{model_name} — no finite data")
        os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        return
    ax.scatter(x, y, s=2, alpha=0.3, rasterized=True, label=f"r = {r:.4f}")
    lim = max(float(x.max()), float(y.max())) * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=0.8, label="y = x")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("True distance")
    ax.set_ylabel("Estimated distance")
    ax.set_title(f"{model_name}")
    ax.legend(fontsize=8)
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_error(tsv_path, out_pdf, model_name):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = _load(tsv_path)
    x   = df["true_distance"].values
    err = df["estimated_distance"].values - x

    centers, means = _rolling_mean(x, err)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, err, s=2, alpha=0.25, rasterized=True, color="steelblue")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    if centers:
        ax.plot(centers, means, color="red", lw=1.5, label="bin mean")
        ax.legend(fontsize=8)
    ax.set_xlabel("True distance")
    ax.set_ylabel("Estimated − True")
    ax.set_title(f"{model_name} — bias / variance")
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} scatter|error <tsv> <output.pdf> <model>")
        sys.exit(1)
    mode, tsv_path, out_pdf, model_name = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    if mode == "scatter":
        plot_scatter(tsv_path, out_pdf, model_name)
    elif mode == "error":
        plot_error(tsv_path, out_pdf, model_name)
    else:
        print(f"Unknown mode: {mode}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
