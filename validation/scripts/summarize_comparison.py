"""
Quantitative comparison of ML vs expected distance estimators.

Usage:
  summarize_comparison.py --output-tsv FILE
                          model:ml_path:exp_path ...
                          [--wag-noprior path]

Each positional argument is a colon-separated triple:
  model     — protein model name
  ml_path   — path to results/protein/{model}/all.tsv
  exp_path  — path to results/protein/{model}/expected_all.tsv

--wag-noprior adds a third set of columns for the flat-prior variant (WAG only).
"""
import argparse
import os
import sys


def _stats(true, est):
    import numpy as np
    import math
    mask = ~(np.isnan(est) | ~np.isfinite(est))
    if mask.sum() < 2:
        return dict(r=float("nan"), rmse=float("nan"), bias=float("nan"), n=0)
    t, e = true[mask], est[mask]
    err  = e - t
    return dict(
        r    = float(np.corrcoef(t, e)[0, 1]),
        rmse = float(np.sqrt(np.mean(err**2))),
        bias = float(np.mean(err)),
        n    = int(mask.sum()),
    )


def _load_pair(ml_path, exp_path):
    import pandas as pd
    import numpy as np
    ml  = pd.read_csv(ml_path,  sep="\t")
    exp = pd.read_csv(exp_path, sep="\t")
    merged = ml.merge(exp, on=["taxon_i", "taxon_j", "true_distance"],
                      suffixes=("_ml", "_exp"))
    true   = merged["true_distance"].values
    est_ml = merged["estimated_distance_ml"].values
    est_exp= merged["estimated_distance_exp"].values
    return true, est_ml, est_exp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-tsv",    required=True)
    ap.add_argument("--wag-noprior",   default=None,
                    help="path to WAG expected_noprior_all.tsv")
    ap.add_argument("inputs", nargs="+")
    args = ap.parse_args()

    import pandas as pd
    import numpy as np

    rows = []
    for triple in args.inputs:
        model, ml_path, exp_path = triple.split(":", 2)
        true, est_ml, est_exp = _load_pair(ml_path, exp_path)
        s_ml  = _stats(true, est_ml)
        s_exp = _stats(true, est_exp)
        row = dict(
            model     = model,
            n_pairs   = s_ml["n"],
            r_ml      = s_ml["r"],
            rmse_ml   = s_ml["rmse"],
            bias_ml   = s_ml["bias"],
            r_exp     = s_exp["r"],
            rmse_exp  = s_exp["rmse"],
            bias_exp  = s_exp["bias"],
        )
        rows.append(row)

    df = pd.DataFrame(rows)

    if args.wag_noprior:
        import pandas as _pd
        np_df  = _pd.read_csv(args.wag_noprior, sep="\t")
        wag_ml = _pd.read_csv(
            next(t.split(":")[1] for t in args.inputs if t.startswith("WAG:")),
            sep="\t",
        )
        merged = wag_ml.merge(np_df, on=["taxon_i", "taxon_j", "true_distance"],
                              suffixes=("_ml", "_np"))
        true   = merged["true_distance"].values
        est_np = merged["estimated_distance_np"].values
        s_np   = _stats(true, est_np)
        df.loc[df["model"] == "WAG", "r_noprior"]    = s_np["r"]
        df.loc[df["model"] == "WAG", "rmse_noprior"] = s_np["rmse"]
        df.loc[df["model"] == "WAG", "bias_noprior"] = s_np["bias"]

    os.makedirs(os.path.dirname(args.output_tsv) or ".", exist_ok=True)
    df.to_csv(args.output_tsv, sep="\t", index=False, float_format="%.6f")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
