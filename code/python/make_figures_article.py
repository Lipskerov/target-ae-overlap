"""Figures for the DDO target article. Values copied from out/all_0124.txt and out/03_examples.txt
(re-measured 2026-11-12: out_20261112/) and, for the as-of row, asof_20261112/RESULT.txt (2026-11-12).
Run: python3 make_figures.py  -> fig1.png/.svg, fig2.png/.svg"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#838179", "#e6e5e0"
ACCENT = "#2a78d6"
plt.rcParams.update({"font.family": "Arial", "font.size": 10, "axes.edgecolor": GRID,
                     "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK})

def frame(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=0, colors=INK2)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)

# ---- Figure 1: dumbbell, same vs different target, under each control
rows = [  # label, same, different, pairs(same)
    ("All pairs", 4.13, 1.11),
    ("Molecules with a single target", 4.38, 1.15),
    ("First approved vs later molecules", 4.11, 1.14),
    ("First molecule's record as of the successor's approval", 2.38, 0.44),
    ("Different FDA pharmacologic class", 4.48, 1.24),
    ("Disease and reporting terms removed", 3.38, 0.91),
    ("No shared indication", 2.27, 0.71),
    ("Three filters combined", 2.54, 0.63),
]
fig, ax = plt.subplots(figsize=(7.2, 4.3), dpi=300)
fig.patch.set_facecolor(SURFACE); frame(ax)
y = list(range(len(rows)))[::-1]
for yi, (lab, s, d) in zip(y, rows):
    ax.plot([d, s], [yi, yi], color=GRID, linewidth=3, solid_capstyle="round", zorder=1)
    ax.scatter([d], [yi], s=70, color=MUTED, edgecolor=SURFACE, linewidth=2, zorder=3)
    ax.scatter([s], [yi], s=70, color=ACCENT, edgecolor=SURFACE, linewidth=2, zorder=3)
    ax.text(5.75, yi, f"{s/d:.1f}×", va="center", ha="right", fontsize=10, color=INK, fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
ax.set_xlim(0, 5.8); ax.set_xticks([0, 1, 2, 3, 4, 5])
ax.set_xlabel("Mean top-20 FAERS signals shared by a pair of molecules", fontsize=9)
# direct labels on the first row instead of a legend box
ax.annotate("Same target", (4.13, y[0]), xytext=(4.13, y[0] + 0.62), ha="right", fontsize=9, color=ACCENT, fontweight="bold")
ax.annotate("Different target", (1.11, y[0]), xytext=(1.11, y[0] + 0.62), ha="center", fontsize=9, color=MUTED, fontweight="bold")
ax.text(5.75, y[0] + 0.62, "Ratio", ha="right", fontsize=9, color=INK2)
ax.set_ylim(-0.6, len(rows) - 0.1)
fig.text(0.01, 0.975, "Molecules at the same target share nearly four times as many safety signals",
         fontsize=11.5, fontweight="bold", va="top")
fig.text(0.01, 0.915, "1,229 molecules, 216 targets; the gap holds under every control", fontsize=9.5, color=INK2, va="top")
fig.text(0.01, 0.01, "Source: FDA Adverse Event Reporting System (openFDA), ChEMBL, RxNorm. Author's analysis, 12 Nov 2026.",
         fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 0.88))
fig.savefig("fig1.png", facecolor=SURFACE); fig.savefig("fig1.svg", facecolor=SURFACE)

# ---- Figure 2: harms most of a target's molecules carry
items = [  # target, reaction, n_with, n_members
    ("HMG-CoA reductase", "Myalgia", 5, 5),
    ("BTK", "Platelet count decreased", 6, 6),
    ("PD-L1", "Pneumonitis", 5, 6),
    ("EGFR", "Stomatitis", 12, 15),
    ("PD-1", "Hypothyroidism", 11, 15),
    ("SGLT2", "Diabetic ketoacidosis", 5, 7),
    ("BTK", "Atrial fibrillation", 4, 6),
    ("PD-1", "Myocarditis", 10, 15),
    ("EGFR", "Interstitial lung disease", 9, 15),
]
fig, ax = plt.subplots(figsize=(7.2, 3.9), dpi=300)
fig.patch.set_facecolor(SURFACE); frame(ax)
y = list(range(len(items)))[::-1]
vals = [100 * a / b for _, _, a, b in items]
ax.barh(y, vals, height=0.56, color=ACCENT, edgecolor=SURFACE, linewidth=0)
for yi, v, (_, _, a, b) in zip(y, vals, items):
    ax.text(v + 1.5, yi, f"{a} of {b}", va="center", fontsize=9, color=INK2)
ax.set_yticks(y)
ax.set_yticklabels([f"{t}  ·  {r}" for t, r, _, _ in items])
ax.set_xlim(0, 112); ax.set_xticks([0, 25, 50, 75, 100]); ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
ax.set_xlabel("Share of the target's molecules with this reaction among their top-20 signals", fontsize=9)
fig.text(0.01, 0.975, "The harm belongs to the mechanism", fontsize=11.5, fontweight="bold", va="top")
fig.text(0.01, 0.915, "Reactions carried by most molecules at the same target", fontsize=9.5, color=INK2, va="top")
fig.text(0.01, 0.01, "Source: FAERS (openFDA), ChEMBL, RxNorm. Author's analysis, 12 Nov 2026. HMG-CoA reductase group lacks atorvastatin and rosuvastatin.",
         fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0.03, 1, 0.88))
fig.savefig("fig2.png", facecolor=SURFACE); fig.savefig("fig2.svg", facecolor=SURFACE)
print("ok")
