"""
Plot the Gascuel (1997) recreation benchmark: mean RF/2 (topological error)
for NJ vs BioNJ, one panel per rate condition, grouped by model tree and
sequence length.

Usage: plot_gascuel.py <reconstruction.tsv> <output_dir>
"""
import os
import sys

METHOD_COLORS = {"nj": "#1f77b4", "bionj": "#d62728"}
CONDITION_ORDER = ["low", "middle", "high_low", "high"]
TREE_ORDER = ["A", "B", "E", "F", "C", "D"]  # constant -> intermediate -> highly-varying


def _rf_plot(df, out_pdf):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_taxa = 12
    rf_scale = (2 * n_taxa - 6) / 2
    conditions = [c for c in CONDITION_ORDER if c in df["condition"].unique()]
    lengths = sorted(df["length"].unique())

    fig, axes = plt.subplots(len(lengths), len(conditions),
                              figsize=(3.2 * len(conditions), 3.2 * len(lengths)),
                              sharey=True, squeeze=False)

    trees = [t for t in TREE_ORDER if t in df["tree"].unique()]
    x = np.arange(len(trees))
    width = 0.35

    for row, length in enumerate(lengths):
        for col, condition in enumerate(conditions):
            ax = axes[row][col]
            sub = df[(df["length"] == length) & (df["condition"] == condition)]
            for i, method in enumerate(("nj", "bionj")):
                means = []
                for t in trees:
                    g = sub[sub["tree"] == t][f"rf_{method}"]
                    means.append(g.mean() * rf_scale if len(g) else float("nan"))
                ax.bar(x + (i - 0.5) * width, means, width,
                       label=method.upper(), color=METHOD_COLORS[method])
            ax.set_xticks(x)
            ax.set_xticklabels(trees)
            ax.set_title(f"{condition}, L={length}", fontsize=9)
            if col == 0:
                ax.set_ylabel("mean RF / 2")
            if row == 0 and col == len(conditions) - 1:
                ax.legend(fontsize=8)

    fig.suptitle("NJ vs BioNJ topological error by model tree, rate condition, sequence length")
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_pdf) or ".", exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <reconstruction.tsv> <output_dir>")
        sys.exit(1)
    tsv_path, out_dir = sys.argv[1], sys.argv[2]

    import pandas as pd
    df = pd.read_csv(tsv_path, sep="\t")
    _rf_plot(df, os.path.join(out_dir, "rf_comparison.pdf"))


if __name__ == "__main__":
    main()
