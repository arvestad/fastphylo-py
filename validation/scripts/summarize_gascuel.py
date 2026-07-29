"""
Summarize the Gascuel (1997) recreation benchmark into a table matching
the style of the paper's Tables 1-6: for each (tree, condition, length),
NJ and BioNJ's %ME< / %ME> (fraction of replicates where the reconstructed
tree's OLS/ME length is shorter/longer than the true tree's), mean RF/2
(on the paper's raw 0-to-(2n-6) scale, n=12 here so 0-18), %RF=0 (exact
topology recovery), and the BioNJ-vs-NJ head-to-head comparison (%RF<,
%RF>, and the RF error reduction BioNJ achieves over NJ).

Usage: summarize_gascuel.py <reconstruction.tsv> <summary.tsv>
"""
import sys


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <reconstruction.tsv> <summary.tsv>")
        sys.exit(1)
    tsv_path, out_path = sys.argv[1], sys.argv[2]

    import pandas as pd

    df = pd.read_csv(tsv_path, sep="\t")
    n_taxa = 12
    rf_scale = (2 * n_taxa - 6) / 2  # paper reports RF/2; my rf_distance is RF/(2n-6)

    rows = []
    for (tree, condition, length), g in df.groupby(["tree", "condition", "length"]):
        row = dict(tree=tree, condition=condition, length=length, n=len(g))

        for method in ("nj", "bionj"):
            me_est = g[f"me_{method}"]
            me_true = g["me_true"]
            rf = g[f"rf_{method}"]
            row[f"{method}_pct_me_lt"] = 100 * (me_est < me_true).mean()
            row[f"{method}_pct_me_gt"] = 100 * (me_est > me_true).mean()
            row[f"{method}_rf_half"] = (rf * rf_scale).mean()
            row[f"{method}_pct_rf0"] = 100 * (rf == 0).mean()

        rf_nj, rf_bionj = g["rf_nj"], g["rf_bionj"]
        row["pct_bionj_better"] = 100 * (rf_bionj < rf_nj).mean()
        row["pct_nj_better"] = 100 * (rf_bionj > rf_nj).mean()
        mean_rf_nj = rf_nj.mean()
        mean_rf_bionj = rf_bionj.mean()
        row["rf_error_reduction_pct"] = (
            100 * (mean_rf_nj - mean_rf_bionj) / mean_rf_nj if mean_rf_nj > 0 else float("nan")
        )
        rows.append(row)

    result = pd.DataFrame(rows).sort_values(["tree", "condition", "length"])
    result.to_csv(out_path, sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(result.to_string(index=False))


if __name__ == "__main__":
    main()
