"""
Summarise and plot the NJ / FNJ / BioNJ reconstruction-accuracy benchmark.

Usage: plot_benchmark.py <reconstruction.tsv> <plot_dir> <summary.tsv>

Produces:
  <plot_dir>/rf_distance.pdf  — normalised RF distance vs n, one panel per
                                 sigma, one line per method.
  <plot_dir>/path_rmse.pdf    — leaf-pair path-distance RMSE vs n (NJ and
                                 BioNJ only; FNJ has no branch lengths).
  <summary.tsv>                — mean/std of both metrics per
                                 (n, sigma, method).
"""
import os
import sys

METHOD_COLORS = {"nj": "#1f77b4", "fnj": "#2ca02c", "bionj": "#d62728"}
METHOD_LABELS = {"nj": "NJ", "fnj": "FNJ", "bionj": "BioNJ"}


def _grid_plot(df, value_col, methods, ylabel, out_pdf):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    sigmas = sorted(df["sigma"].unique())
    fig, axes = plt.subplots(1, len(sigmas), figsize=(4.2 * len(sigmas), 4), sharey=True)
    if len(sigmas) == 1:
        axes = [axes]

    n_vals_all = sorted(df["n"].unique())
    for ax, sigma in zip(axes, sigmas):
        sub = df[df["sigma"] == sigma]
        for method in methods:
            msub = sub[sub["method"] == method]
            grp = msub.groupby("n")[value_col].agg(["mean", "std", "count"]).dropna(subset=["mean"])
            if grp.empty:
                continue
            n_vals = grp.index.values
            means  = grp["mean"].values
            sem    = (grp["std"] / np.sqrt(grp["count"])).values
            ax.errorbar(n_vals, means, yerr=sem, marker="o", ms=4, lw=1.2,
                        capsize=2, label=METHOD_LABELS[method],
                        color=METHOD_COLORS[method])
        ax.set_xscale("log", base=2)
        ax.set_xticks(n_vals_all)
        ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
        ax.set_xlabel("n (leaves)")
        ax.set_title(f"sigma = {sigma}")

    axes[0].set_ylabel(ylabel)
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_pdf) or ".", exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def _summary(df, out_tsv):
    summary = (
        df.groupby(["n", "sigma", "method"])[["rf", "path_rmse"]]
        .agg(["mean", "std"])
    )
    summary.columns = ["_".join(c) for c in summary.columns]
    summary = summary.reset_index()
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    summary.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")
    return summary


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <reconstruction.tsv> <plot_dir> <summary.tsv>")
        sys.exit(1)
    tsv_path, plot_dir, summary_tsv = sys.argv[1], sys.argv[2], sys.argv[3]

    import pandas as pd
    df = pd.read_csv(tsv_path, sep="\t")

    _grid_plot(df, "rf", ["nj", "fnj", "bionj"],
               "Normalised Robinson-Foulds distance",
               os.path.join(plot_dir, "rf_distance.pdf"))

    bl_df = df[df["method"].isin(["nj", "bionj"])]
    _grid_plot(bl_df, "path_rmse", ["nj", "bionj"],
               "Leaf-pair path-distance RMSE",
               os.path.join(plot_dir, "path_rmse.pdf"))

    summary = _summary(df, summary_tsv)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
