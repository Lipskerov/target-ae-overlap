#!/usr/bin/env python3
"""Main-text figures for the revised analysis, drawn only from saved result files. -> results/figures/"""
import json, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
P = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")); T = os.path.join(P, "results/tables"); F = os.path.join(P, "results/figures"); os.makedirs(F, exist_ok=True)
J = lambda f: json.load(open(os.path.join(T, f)))
SURF, INK, INK2, MUT, GRID, ACC, ACC2 = "#ffffff", "#0b0b0b", "#52514e", "#8a8880", "#e3e2dd", "#2a78d6", "#d9822b"
plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.edgecolor": GRID, "text.color": INK, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK})
def frame(ax, axis="x"):
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID); ax.tick_params(length=0); ax.grid(axis=axis, color=GRID, linewidth=0.7); ax.set_axisbelow(True)
def save(fig, name):
    fig.savefig(os.path.join(F, name + ".png"), dpi=300, facecolor=SURF); fig.savefig(os.path.join(F, name + ".svg"), facecolor=SURF); plt.close(fig)
VAR = (("prr_clean", "PRR ranking, standard cleaning", ACC), ("ic025_clean2", "IC025 ranking, extended cleaning", ACC2))

# Fig 1: molecule-averaged difference by extraction, both variants
ext = (("b", "Primary + secondary suspect"), ("c", "Primary suspect only"), ("d", "Suspect, health-professional reports"), ("e", "Sole suspect"), ("balltargets", "Suspect, all ChEMBL targets"))
fig, ax = plt.subplots(figsize=(7.2, 3.9)); frame(ax); y0 = len(ext)
for k, (e, lab) in enumerate(ext):
    for j, (v, vl, col) in enumerate(VAR):
        x = J(f"primary_R2_{e}_{v}.json")["shared_signals"]; yy = y0 - k + (0.16 if j == 0 else -0.16)
        ax.plot(x["ci95"], [yy, yy], color=col, linewidth=2.2, solid_capstyle="round"); ax.scatter([x["delta"]], [yy], s=30, color=col, zorder=3, label=vl if k == 0 else None)
        ax.text(x["ci95"][1] + 0.05, yy, f"{x['delta']:.2f}", va="center", fontsize=7.5, color=INK2)
ax.axvline(0, color=INK2, linewidth=0.8); ax.set_yticks([y0 - k for k in range(len(ext))]); ax.set_yticklabels([l for _, l in ext]); ax.set_xlim(-0.1, 2.9)
ax.set_xlabel("Difference in shared profile signals, same-target minus different-target (95% CI)", fontsize=8.5); ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
fig.tight_layout(); save(fig, "fig1_difference_by_extraction")

# Fig 2: same-target minus different-target within class strata, with CIs
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
rows = (("atc4", "same_class", "Same ATC level-4 group"), ("atc4", "different_class", "Different ATC level-4 group"), ("epc", "same_class", "Same FDA class"), ("epc", "different_class", "Different FDA class"),
        ("indication", "same_class", "Shared indication"), ("indication", "different_class", "No shared indication"))
for ax, (v, vl, col) in zip(axes, VAR):
    frame(ax); t = J(f"target_vs_class_R2_{v}.json")["within_class_contrasts"]; ys = list(range(len(rows)))[::-1]
    for yy, (k, c, lab) in zip(ys, rows):
        d = t[k][c]; ax.plot(d["ci95"], [yy, yy], color=col, linewidth=2.2, solid_capstyle="round"); ax.scatter([d["difference"]], [yy], s=30, color=col, zorder=3)
        ax.text(d["ci95"][1] + 0.04, yy, f"{d['difference']:.2f}  (n = {d['same_target_pairs']:,})", va="center", fontsize=7, color=INK2)
    ax.axvline(0, color=INK2, linewidth=0.8); ax.set_yticks(ys); ax.set_yticklabels([r[2] for r in rows]); ax.set_title(vl, fontsize=8.5, loc="left"); ax.set_xlim(-0.1, 3.3)
    ax.set_xlabel("Same-target minus different-target (95% CI)", fontsize=8)
fig.tight_layout(); save(fig, "fig2_target_within_class_strata")

# Fig 3: adjusted same-target coefficient across models
fig, ax = plt.subplots(figsize=(6.4, 2.6)); frame(ax)
mods = (("model_unadjusted", "Unadjusted"), ("model_class_only", "Class (ATC 4/3/1, FDA class)"), ("model_full", "Class + indication + volume + year gap + co-mention"))
for k, (mk, lab) in enumerate(mods):
    for j, (v, vl, col) in enumerate(VAR):
        d = J(f"target_vs_class_R2_{v}.json")[mk]["same_target"]; yy = len(mods) - k + (0.15 if j == 0 else -0.15)
        ax.plot(d["ci95"], [yy, yy], color=col, linewidth=2.2, solid_capstyle="round"); ax.scatter([d["beta"]], [yy], s=30, color=col, zorder=3, label=vl if k == 0 else None)
        ax.text(d["ci95"][1] + 0.03, yy, f"{d['beta']:.2f}", va="center", fontsize=7.5, color=INK2)
ax.axvline(0, color=INK2, linewidth=0.8); ax.set_yticks([len(mods) - k for k in range(len(mods))]); ax.set_yticklabels([l for _, l in mods]); ax.set_xlabel("Same-target coefficient, all pairs (95% CI)", fontsize=8.5); ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); save(fig, "fig3_adjusted_coefficient")

# Fig 4: co-mention exclusion, as-of analysis and trial analysis, with intervals
fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.9))
ax = axes[0]; frame(ax)
for j, (v, vl, col) in enumerate(VAR):
    c = J(f"comention_summary_{v}.json")
    for k, key in enumerate(("difference_before", "difference_after")):
        d = c[key]; yy = 1 - k + (0.15 if j == 0 else -0.15); ax.plot(d["ci95"], [yy, yy], color=col, linewidth=2.2, solid_capstyle="round"); ax.scatter([d["estimate"]], [yy], s=28, color=col, zorder=3)
ax.set_yticks([1, 0]); ax.set_yticklabels(["All reports", "Co-mention\nreports removed"]); ax.axvline(0, color=INK2, linewidth=0.8); ax.set_title("Co-mention", fontsize=8.5, loc="left"); ax.set_xlabel("Same − different target", fontsize=8)
ax = axes[1]; frame(ax)
for j, (v, vl, col) in enumerate(VAR):
    a = J(f"asof_summary_{v}.json")["clean"]
    for k, key in enumerate(("eA", "eB", "eC")):
        d = a[key]; yy = 2 - k + (0.15 if j == 0 else -0.15); ax.plot(d["ci95"], [yy, yy], color=col, linewidth=2.2, solid_capstyle="round"); ax.scatter([d["estimate"]], [yy], s=28, color=col, zorder=3)
ax.set_yticks([2, 1, 0]); ax.set_yticklabels(["vs volume-matched\nfirst molecules", "vs other\nsuccessors", "vs class-matched\nfirst molecules"], fontsize=7.5); ax.axvline(0, color=INK2, linewidth=0.8)
ax.set_title("As-of-date analysis", fontsize=8.5, loc="left"); ax.set_xlabel("First molecule − control", fontsize=8)
ax = axes[2]; frame(ax)
for k, (lab, nm) in enumerate((("R2_ctgov_all", "All pairs"), ("R2_ctgov_no_shared_trial", "No shared\ntrial"), ("R2_ctgov_no_shared_condition", "No shared\ncondition"))):
    d = J(f"primary_{lab}.json")["shared_signals"]; yy = 2 - k; ax.plot(d["ci95"], [yy, yy], color=MUT, linewidth=2.2, solid_capstyle="round"); ax.scatter([d["delta"]], [yy], s=28, color=MUT, zorder=3)
ax.set_yticks([2, 1, 0]); ax.set_yticklabels(["All pairs", "No shared\ntrial", "No shared\ncondition"], fontsize=7.5); ax.axvline(0, color=INK2, linewidth=0.8); ax.set_title("Trial-reported events", fontsize=8.5, loc="left"); ax.set_xlabel("Same − different target", fontsize=8)
fig.tight_layout(); save(fig, "fig4_comention_asof_trials"); print(sorted(os.listdir(F)))
