"""
Render summary.tsv as a formatted table PDF using matplotlib.

Usage: render_summary.py <summary.tsv> <output.pdf>
"""
import sys
import os


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <summary.tsv> <output.pdf>")
        sys.exit(1)

    tsv_path, out_pdf = sys.argv[1], sys.argv[2]

    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(tsv_path, sep="\t")

    # Format float columns
    float_cols = ["rmse", "mae", "pearson_r", "mean_bias", "nan_frac"]
    for col in float_cols:
        df[col] = df[col].apply(lambda v: f"{v:.4f}" if pd.notna(v) else "—")
    df["n_pairs"] = df["n_pairs"].astype(int)

    col_labels = ["dtype", "model", "n_pairs", "RMSE", "MAE", "Pearson r", "mean bias", "NaN frac"]
    cell_data  = df.values.tolist()

    n_rows = len(cell_data)
    fig_h  = max(2.5, 0.35 * (n_rows + 2))
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.axis("off")

    table = ax.table(
        cellText=cell_data,
        colLabels=col_labels,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.auto_set_column_width(list(range(len(col_labels))))

    # Style header row
    for col_idx in range(len(col_labels)):
        cell = table[0, col_idx]
        cell.set_facecolor("#2c3e50")
        cell.set_text_props(color="white", fontweight="bold")

    # Alternate row shading, highlight best RMSE per dtype
    rmse_vals = df["rmse"].tolist()
    dna_rows    = [i for i, r in enumerate(df["dtype"]) if r == "dna"]
    prot_rows   = [i for i, r in enumerate(df["dtype"]) if r == "protein"]

    def _best(rows):
        vals = [(float(rmse_vals[i]), i) for i in rows]
        return min(vals)[1] if vals else -1

    best_dna  = _best(dna_rows)
    best_prot = _best(prot_rows)

    for row_idx in range(n_rows):
        bg = "#eaf4fb" if row_idx % 2 == 0 else "white"
        if row_idx in (best_dna, best_prot):
            bg = "#d5f5e3"   # green highlight for best RMSE per dtype
        for col_idx in range(len(col_labels)):
            table[row_idx + 1, col_idx].set_facecolor(bg)

    ax.set_title("Distance estimation summary", fontsize=12, fontweight="bold", pad=12)

    os.makedirs(os.path.dirname(out_pdf) or ".", exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
