"""
Side-by-side scatter: ML vs true (left) and expected vs true (right).

Usage: plot_comparison.py <ml_all.tsv> <exp_all.tsv> <output.pdf> <model_name>
"""
import sys
import os
import math
import numpy as np


def _load(path, col="estimated_distance"):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    df = df.dropna(subset=["true_distance", col])
    df = df[df[col].apply(lambda v: math.isfinite(v))]
    return df


def _stats(x, y):
    if len(x) < 2:
        return float("nan"), float("nan")
    r    = float(np.corrcoef(x, y)[0, 1])
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    return r, rmse


def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <ml.tsv> <exp.tsv> <output.pdf> <model>")
        sys.exit(1)

    ml_path, exp_path, out_pdf, model_name = sys.argv[1:]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df_ml  = _load(ml_path)
    df_exp = _load(exp_path)

    true_ml  = df_ml["true_distance"].values
    est_ml   = df_ml["estimated_distance"].values
    true_exp = df_exp["true_distance"].values
    est_exp  = df_exp["estimated_distance"].values

    r_ml,  rmse_ml  = _stats(true_ml,  est_ml)
    r_exp, rmse_exp = _stats(true_exp, est_exp)

    lim = max(float(true_ml.max()),  float(est_ml.max()),
              float(true_exp.max()), float(est_exp.max())) * 1.05

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    fig.suptitle(f"{model_name} — ML vs expected distance", fontsize=12)

    for ax, x, y, label, r, rmse in [
        (axes[0], true_ml,  est_ml,  "ML",       r_ml,  rmse_ml),
        (axes[1], true_exp, est_exp, "Expected",  r_exp, rmse_exp),
    ]:
        ax.scatter(x, y, s=2, alpha=0.3, rasterized=True,
                   label=f"r={r:.4f}  RMSE={rmse:.4f}")
        ax.plot([0, lim], [0, lim], "k--", lw=0.8)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.set_xlabel("True distance")
        ax.set_ylabel("Estimated distance")
        ax.set_title(label)
        ax.legend(fontsize=8)

    os.makedirs(os.path.dirname(out_pdf) or ".", exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
