"""
Bland-Altman figure for the report, from the *_bland_altman.csv files written by the module or by run_headless.py.

  python scripts/plot_bland_altman.py                 # SIMONE expected results -> docs/bland_altman.png and .pdf
  python scripts/plot_bland_altman.py --dir RESULTS --out figure.png

Left panel: the elastic modules (Torso, Femoral). Right panel: the rigid module (Cranial).
x = mean of printed and design diameter, y = printed - design. Each module gets its own bias line (solid) and
95% limits of agreement (band, bias +- 1.96 SD), because the modules do not share one bias.
Needs numpy and matplotlib (not part of Slicer's Python).
"""
import argparse, csv, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
PANELS = [("Elastic modules", ["Torso", "Femoral"]), ("Rigid module", ["Cranial"])]
COLOUR = {"Torso": "#2a78d6", "Femoral": "#eb6834", "Cranial": "#1baf7a"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def load(folder, name):
    rows = list(csv.DictReader(open(os.path.join(folder, f"{name}_bland_altman.csv"))))
    return np.array([float(r["mean_D_mm"]) for r in rows]), np.array([float(r["diff_mm"]) for r in rows])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(ROOT, "SIMONE_validation_files", "expected_results"))
    ap.add_argument("--out", default=os.path.join(ROOT, "docs", "bland_altman.png"))
    a = ap.parse_args()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True, gridspec_kw=dict(width_ratios=[1.45, 1], wspace=0.08))
    for ax, (title, names) in zip(axes, PANELS):
        ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
        ax.grid(axis="y", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
        data = {n: load(a.dir, n) for n in names}
        x_hi = max(m.max() for m, _ in data.values()); x_lo = min(m.min() for m, _ in data.values()); pad = 0.04 * (x_hi - x_lo)
        for n, (m, d) in data.items():
            c = COLOUR[n]; bias, sd = d.mean(), d.std(ddof=1); lo, hi = bias - 1.96 * sd, bias + 1.96 * sd
            ax.fill_between([m.min(), m.max()], lo, hi, color=c, alpha=0.10, lw=0, zorder=2)
            for y in (lo, hi): ax.plot([m.min(), m.max()], [y, y], color=c, lw=0.8, ls=(0, (4, 3)), zorder=3)
            lab = f"{n}, n = {len(d)}\nbias {bias:+.2f} mm ({lo:+.2f} to {hi:+.2f})".replace("-", "−")
            ax.scatter(m, d, s=7, color=c, alpha=0.55, lw=0, zorder=4, label=lab)
            ax.plot([m.min(), m.max()], [bias, bias], color=c, lw=1.6, zorder=5)
            print(f"{n:8s} n={len(d):4d}  bias {bias:+.2f} mm  SD {sd:.2f}  LoA {lo:+.2f} to {hi:+.2f} mm")
        ax.set_xlim(x_lo - pad, x_hi + pad)
        ax.set_title(title, loc="left", fontsize=10, color=INK, pad=64)
        leg = ax.legend(loc="lower left", bbox_to_anchor=(-0.02, 1.0), frameon=False, fontsize=7.5, handletextpad=0.2, markerscale=2.2, labelcolor=INK,
                        labelspacing=0.8, borderaxespad=0.0, ncol=1)
        for h in leg.legend_handles: h.set_alpha(1)
    axes[0].set_ylabel("Printed − design diameter (mm)")
    axes[0].set_ylim(-1.6, 0.8)
    fig.supxlabel("Mean of printed and design diameter (mm)", fontsize=9, color=INK, y=-0.05)
    fig.savefig(a.out, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.splitext(a.out)[0] + ".pdf", bbox_inches="tight", facecolor="white")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
