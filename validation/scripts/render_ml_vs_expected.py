"""
Render ml_vs_expected.tsv as a formatted table PDF.

Usage: render_ml_vs_expected.py <ml_vs_expected.tsv> <output.pdf>
"""
import sys
import os


# Column groups: (header_label, bg_color, columns_in_tsv)
_GROUPS = [
    ("ML",                  "#2c3e50", ["r_ml",       "rmse_ml",       "bias_ml"]),
    ("Expected + prior",    "#1a5276", ["r_exp",      "rmse_exp",      "bias_exp"]),
    ("Expected flat prior", "#4a235a", ["r_noprior",  "rmse_noprior",  "bias_noprior"]),
]
_INFO_COLS = ["model", "n_pairs"]
_COL_LABELS = ["model", "n"] + ["r", "RMSE", "bias"] * 3


def _fmt(v):
    try:
        f = float(v)
        import math
        if math.isnan(f):
            return "—"
        return f"{f:.4f}"
    except (TypeError, ValueError):
        return str(v)


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <ml_vs_expected.tsv> <output.pdf>")
        sys.exit(1)

    tsv_path, out_pdf = sys.argv[1], sys.argv[2]

    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    df = pd.read_csv(tsv_path, sep="\t")

    # Build cell data in display column order
    display_cols = _INFO_COLS + [c for _, _, cols in _GROUPS for c in cols]
    rows = []
    for _, row in df.iterrows():
        cells = [str(int(row["n_pairs"])) if c == "n_pairs" else _fmt(row.get(c))
                 for c in display_cols]
        rows[len(rows):] = [cells]

    n_rows = len(rows)
    n_cols = len(_COL_LABELS)

    # Figure sizing
    fig_h = max(3.5, 0.42 * (n_rows + 3))
    fig, ax = plt.subplots(figsize=(13, fig_h))
    ax.axis("off")

    # --- draw the data table ---
    table = ax.table(
        cellText=rows,
        colLabels=_COL_LABELS,
        loc="center",
        cellLoc="center",
        bbox=[0, 0, 1, 0.88],   # leave 12 % at top for the group header band
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.auto_set_column_width(list(range(n_cols)))

    # Column colour bands (light versions for data cells, dark for header)
    _LIGHT = {"ML": "#d6eaf8", "Expected + prior": "#d0e8f5",
              "Expected flat prior": "#e8daef"}
    col_group = {}   # col_idx → group label
    idx = len(_INFO_COLS)
    for label, _, cols in _GROUPS:
        for _ in cols:
            col_group[idx] = label
            idx += 1

    # Style header row (row 0)
    info_header_bg = "#34495e"
    for col_idx in range(n_cols):
        cell = table[0, col_idx]
        if col_idx < len(_INFO_COLS):
            cell.set_facecolor(info_header_bg)
        else:
            g_label = col_group[col_idx]
            cell.set_facecolor(next(c for l, c, _ in _GROUPS if l == g_label))
        cell.set_text_props(color="white", fontweight="bold")

    # Style data rows: alternate shading + column tinting
    # Highlight cells where expected beats ML (lower RMSE / less absolute bias)
    rmse_ml_col  = display_cols.index("rmse_ml")
    rmse_exp_col = display_cols.index("rmse_exp")
    bias_ml_col  = display_cols.index("bias_ml")
    bias_exp_col = display_cols.index("bias_exp")

    for row_idx, row_data in enumerate(rows):
        alt = row_idx % 2 == 0
        try:
            rmse_ml  = float(row_data[rmse_ml_col])
            rmse_exp = float(row_data[rmse_exp_col])
            bias_ml  = abs(float(row_data[bias_ml_col]))
            bias_exp = abs(float(row_data[bias_exp_col]))
            exp_wins_rmse = rmse_exp < rmse_ml
            exp_wins_bias = bias_exp < bias_ml
        except ValueError:
            exp_wins_rmse = exp_wins_bias = False

        for col_idx in range(n_cols):
            cell = table[row_idx + 1, col_idx]
            if col_idx < len(_INFO_COLS):
                cell.set_facecolor("#eaecee" if alt else "white")
            else:
                base = _LIGHT.get(col_group[col_idx], "white")
                # Slightly de-saturate the alternating row
                cell.set_facecolor(base if alt else "white")
                # Green tint where expected estimator wins
                c_name = display_cols[col_idx]
                if exp_wins_rmse and c_name in ("rmse_exp",):
                    cell.set_facecolor("#a9dfbf")
                if exp_wins_bias and c_name in ("bias_exp",):
                    cell.set_facecolor("#a9dfbf")

    # --- group header band above the table ---
    # Use axes-fraction coordinates matching the bbox above
    band_y  = 0.90
    band_h  = 0.08
    x_left  = 0.0

    # info columns width fraction (approx 2/n_cols)
    info_frac  = len(_INFO_COLS) / n_cols
    group_frac = 3 / n_cols   # each group has 3 columns

    ax.add_patch(mpatches.FancyBboxPatch(
        (x_left, band_y), info_frac, band_h,
        boxstyle="square,pad=0", transform=ax.transAxes,
        facecolor="#34495e", edgecolor="none", clip_on=False,
    ))
    x = x_left + info_frac
    for label, color, _ in _GROUPS:
        has_data = any(
            any(r[display_cols.index(c)] != "—" for r in rows)
            for c in next(cols for l, _, cols in _GROUPS if l == label)
        )
        fc = color if has_data else "#aaaaaa"
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, band_y), group_frac, band_h,
            boxstyle="square,pad=0", transform=ax.transAxes,
            facecolor=fc, edgecolor="white", linewidth=0.5, clip_on=False,
        ))
        ax.text(x + group_frac / 2, band_y + band_h / 2, label,
                transform=ax.transAxes, ha="center", va="center",
                fontsize=8.5, fontweight="bold", color="white", clip_on=False)
        x += group_frac

    ax.set_title("ML vs expected distance — quantitative comparison",
                 fontsize=11, fontweight="bold", pad=4)

    os.makedirs(os.path.dirname(out_pdf) or ".", exist_ok=True)
    fig.savefig(out_pdf, bbox_inches="tight", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
