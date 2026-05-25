"""
Build a summary table and combined comparison plots across all models.

Usage:
  summarize.py [--output-tsv FILE] [--output-md FILE]
               [--comp-dna FILE] [--comp-prot FILE]
               dtype:model:path.tsv ...

Each positional argument is a colon-separated triple:
  dtype  — "dna" or "protein"
  model  — fastphylo model name
  path   — path to the aggregated all.tsv for that model
"""
import argparse
import math
import os
import sys


def _load(path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t")
    return df


def _stats(df):
    import numpy as np
    x   = df["true_distance"].values
    y   = df["estimated_distance"].values
    nan_frac = float(np.isnan(y).mean())
    mask = ~np.isnan(y)
    if mask.sum() < 2:
        return dict(rmse=float("nan"), mae=float("nan"),
                    pearson_r=float("nan"), mean_bias=float("nan"),
                    nan_frac=nan_frac, n_pairs=int(mask.sum()))
    xm, ym = x[mask], y[mask]
    err    = ym - xm
    rmse   = float(np.sqrt(np.mean(err**2)))
    mae    = float(np.mean(np.abs(err)))
    bias   = float(np.mean(err))
    r      = float(np.corrcoef(xm, ym)[0, 1])
    return dict(rmse=rmse, mae=mae, pearson_r=r, mean_bias=bias,
                nan_frac=nan_frac, n_pairs=int(mask.sum()))


def _comparison_plot(entries, out_pdf, dtype_label):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 7))
    all_vals = []
    for model, df in entries:
        import math as _math
        df2 = df.dropna(subset=["true_distance", "estimated_distance"])
        df2 = df2[df2["estimated_distance"].apply(lambda v: _math.isfinite(v))]
        x   = df2["true_distance"].values
        y   = df2["estimated_distance"].values
        r   = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else float("nan")
        ax.scatter(x, y, s=1, alpha=0.15, rasterized=True, label=f"{model} r={r:.3f}")
        all_vals.extend(x.tolist())
        all_vals.extend(y.tolist())

    finite_vals = [v for v in all_vals if _math.isfinite(v)]
    lim = max(finite_vals) * 1.05 if finite_vals else 1.0
    ax.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("True distance")
    ax.set_ylabel("Estimated distance")
    ax.set_title(f"Model comparison — {dtype_label}")
    ax.legend(fontsize=6, markerscale=3, loc="upper left")
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def _md_table(rows, headers):
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-tsv",  default="summary/summary.tsv")
    ap.add_argument("--output-md",   default="summary/summary.md")
    ap.add_argument("--comp-dna",    default="summary/comparison_dna.pdf")
    ap.add_argument("--comp-prot",   default="summary/comparison_protein.pdf")
    ap.add_argument("inputs", nargs="+", help="dtype:model:path triples")
    args = ap.parse_args()

    import pandas as pd

    records   = []
    dna_data  = []
    prot_data = []

    for triple in args.inputs:
        dtype, model, path = triple.split(":", 2)
        df = _load(path)
        s  = _stats(df)
        records.append(dict(dtype=dtype, model=model, **s))
        if dtype == "dna":
            dna_data.append((model, df))
        else:
            prot_data.append((model, df))

    summary = pd.DataFrame(records)
    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    summary.to_csv(args.output_tsv, sep="\t", index=False, float_format="%.6f")

    headers = ["dtype", "model", "n_pairs", "rmse", "mae", "pearson_r", "mean_bias", "nan_frac"]
    rows = []
    for _, row in summary.iterrows():
        rows.append([
            row["dtype"],
            row["model"],
            int(row["n_pairs"]),
            f"{row['rmse']:.5f}",
            f"{row['mae']:.5f}",
            f"{row['pearson_r']:.5f}",
            f"{row['mean_bias']:.5f}",
            f"{row['nan_frac']:.4f}",
        ])
    md = "# Distance estimation summary\n\n" + _md_table(rows, headers)
    with open(args.output_md, "w") as fh:
        fh.write(md + "\n")

    if dna_data:
        _comparison_plot(dna_data, args.comp_dna, "DNA")
    if prot_data:
        _comparison_plot(prot_data, args.comp_prot, "protein")


if __name__ == "__main__":
    main()
