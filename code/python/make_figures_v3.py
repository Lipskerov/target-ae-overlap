#!/usr/bin/env python3
"""F1 — figures v3 (Fig 1, 4, 5, 6B/C; Suppl S3, S4).

Reads ONLY JSON/CSV result files in results/tables/ (Fig 1 flow numbers are hard-coded,
cited to results/RESULTS_SUMMARY.md). Writes PNG (300 dpi) + SVG to results/figures/v3/
and a manifest (results/figures/v3/figures_v3_manifest.json) with every plotted value,
its source file, and the figures that were skipped and why.

Style: journal, Okabe-Ito colorblind-safe palette, PRR and IC025 side by side,
95% CIs as bars, no significance stars. Post hoc analyses are labelled as such.

Run:
  python3 code/python/make_figures_v3.py   # from the repository root
"""
import csv
import json
import os
import sys
import traceback

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TAB = os.path.join(ROOT, "results", "tables")
OUT = os.path.join(ROOT, "results", "figures", "v3")
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito
OI = {"blue": "#0072B2", "verm": "#D55E00", "green": "#009E73", "yellow": "#F0E442",
      "sky": "#56B4E9", "orange": "#E69F00", "purple": "#CC79A7", "black": "#000000",
      "grey": "#999999"}
RANK = [("prr_clean", "PRR, standard cleaning"), ("ic025_clean2", "IC025, extended cleaning (amendment 2)")]
RANK_COL = {"prr_clean": OI["blue"], "ic025_clean2": OI["verm"]}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "svg.fonttype": "none", "savefig.dpi": 300,
})

MANIFEST = {"script": "code/python/make_figures_v3.py", "out_dir": "results/figures/v3/",
            "made": {}, "skipped": {}}


def log(*a):
    print(*a, flush=True)


def load_json(name):
    p = os.path.join(TAB, name)
    if not os.path.exists(p):
        raise FileNotFoundError(f"missing input: results/tables/{name}")
    with open(p) as f:
        return json.load(f)


def save(fig, stem, sources, values):
    paths = []
    for ext in ("png", "svg"):
        p = os.path.join(OUT, f"{stem}.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight")
        paths.append(os.path.relpath(p, ROOT))
    plt.close(fig)
    MANIFEST["made"][stem] = {"files": paths, "sources": sources, "values": values}
    log(f"[made] {stem}: {paths}")


def skip(stem, why):
    MANIFEST["skipped"][stem] = why
    log(f"[SKIP] {stem}: {why}")


def fmt(b, ci):
    return f"{b:.2f} ({ci[0]:.2f} to {ci[1]:.2f})"


def forest_axis(ax, rows, xlim, color, marker="o", show_labels=True, ref=None, ref_label=None):
    """rows: list of dicts {label, beta, ci, n, kind} or {'header': text}. y from top."""
    y = 0
    ys, yt, ylab = [], [], []
    for r in rows:
        if "header" in r:
            y += 0.4
            if show_labels:
                import matplotlib.transforms as mtrans
                ax.text(-0.01, -y, r["header"], fontsize=6.8, fontweight="bold", va="center", ha="right",
                        transform=mtrans.blended_transform_factory(ax.transAxes, ax.transData), clip_on=False)
            y += 1
            continue
        b, lo, hi = r["beta"], r["ci"][0], r["ci"][1]
        m = r.get("marker", marker)
        c = r.get("color", color)
        lo_c, hi_c = max(lo, xlim[0]), min(hi, xlim[1])
        ax.plot([lo_c, hi_c], [-y, -y], color=c, lw=1.0, solid_capstyle="butt")
        if lo < xlim[0]:
            ax.annotate("", xy=(xlim[0], -y), xytext=(xlim[0] + 0.04 * (xlim[1] - xlim[0]), -y),
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=0.8, mutation_scale=6))
            ax.text(xlim[0], -y + 0.38, f"lower {lo:.2f}", fontsize=5, color=c, ha="left", va="bottom")
        if hi > xlim[1]:
            ax.annotate("", xy=(xlim[1], -y), xytext=(xlim[1] - 0.04 * (xlim[1] - xlim[0]), -y),
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=0.8, mutation_scale=6))
            ptxt = f"point {b:.2f}; " if b > xlim[1] else ""
            ax.text(xlim[1], -y + 0.38, f"{ptxt}upper {hi:.2f}", fontsize=5, color=c, ha="right", va="bottom")
        if xlim[0] <= b <= xlim[1]:
            ax.plot(b, -y, marker=m, color=c, ms=r.get("ms", 3.6), mec=c,
                    mfc=r.get("mfc", c), mew=0.8, zorder=3)
        ys.append(-y)
        yt.append(r["label"])
        ylab.append(r)
        y += 1
    ax.set_xlim(*xlim)
    ax.set_ylim(-y + 0.4, 0.6)
    ax.axvline(0, color="black", lw=0.6)
    if ref is not None:
        ax.axvline(ref, color=OI["grey"], lw=0.8, ls="--")
        if ref_label:
            ax.text(ref, 0.75, ref_label, fontsize=5.5, color="#555555", ha="center", va="bottom")
    ax.set_yticks(ys)
    ax.set_yticklabels(yt if show_labels else [""] * len(yt))
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    return ys


# --------------------------------------------------------------------------- Fig 1
def fig1():
    stem = "fig1_design_flow"
    # Flow numbers hard-coded from results/RESULTS_SUMMARY.md (generated by
    # code/python/make_results_summary_r2.py): §1 data flow, §3 cohort, §4 / §5-8 pairs,
    # §9 ClinicalTrials.gov parallel analysis.
    faers = [
        ("FAERS quarterly files 2004Q1-2026Q3", "91 files; 25,233,828 DEMO rows"),          # RS §1
        ("Case-level deduplication", "21,054,151 cases"),                                     # RS §1
        ("Field-based duplicate removal", "20,571,050 cases (483,101 removed)"),              # RS §1
        ("Molecules with \u2265100 suspect reports", "2,104 molecules"),                          # RS §3
        ("Cohort: molecules at human protein targets", "1,039 molecules; 180 targets"),       # RS §3
        ("Molecule pairs", "539,241 pairs; 7,923 same-target"),                               # RS §4, §5-8
    ]
    trials = [
        ("ClinicalTrials.gov result arms", "94,239 arms"),                                    # RS §9
        ("Monotherapy arms (rule v2)", "18,770 arms"),                                        # RS §9
        ("Molecules meeting minimum data", "626 molecules"),                                  # RS §9
        ("Trial cohort", "357 molecules; 89 targets"),                                        # RS §9
        ("Molecule pairs", "1,187 same-target; 62,359 different-target"),                     # RS §9
    ]
    fig = plt.figure(figsize=(7.2, 4.6))
    # Panel A: schematic
    axA = fig.add_axes([0.0, 0.0, 0.47, 1.0])
    axA.set_xlim(0, 10)
    axA.set_ylim(0, 10)
    axA.axis("off")
    axA.text(0.1, 9.7, "A", fontsize=10, fontweight="bold", va="top")
    axA.text(0.7, 9.65, "Design: do same-target molecules share adverse-event signals?", fontsize=7,
             va="top", fontweight="bold")

    def box(ax, x, y, w, h, text, fc, ec="black", fs=6.3, bold=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.15",
                                    fc=fc, ec=ec, lw=0.7))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", wrap=True)

    def arrow(ax, x1, y1, x2, y2, color="black", ls="-"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8,
                                     color=color, lw=0.8, ls=ls))

    # molecules and target
    box(axA, 0.4, 7.2, 2.4, 1.2, "Molecule A\n(top-20 FAERS\nsignals)", "#DCEBF6")
    box(axA, 0.4, 4.4, 2.4, 1.2, "Molecule B\n(top-20 FAERS\nsignals)", "#DCEBF6")
    box(axA, 4.0, 5.8, 2.0, 1.2, "Shared target\n(ChEMBL)", "#FBE3D2")
    arrow(axA, 2.8, 7.6, 4.0, 6.7)
    arrow(axA, 2.8, 5.2, 4.0, 6.1)
    box(axA, 7.0, 5.8, 2.7, 1.2, "Shared signals\nper pair\n(0-20)", "#D9F0E7", bold=True)
    arrow(axA, 6.0, 6.4, 7.0, 6.4)
    axA.text(6.5, 6.6, "pair", fontsize=5.5, ha="center", va="bottom")
    axA.text(5.0, 4.8, "Estimand: same-target minus\ndifferent-target pairs (Δ);\n"
             "adjusted same-target coefficient (OLS)", fontsize=5.8, ha="center", va="top",
             style="italic")
    # confounders
    conf = ["Pharmacological class (ATC-4/3/1, FDA EPC)", "Shared indication",
            "Report volume (log sum, log difference)", "Approval-year gap",
            "Co-mention in the same report"]
    axA.add_patch(FancyBboxPatch((0.4, 0.35), 9.3, 2.75, boxstyle="round,pad=0.02,rounding_size=0.15",
                                 fc="#F2F2F2", ec=OI["grey"], lw=0.7, ls="--"))
    axA.text(0.6, 2.9, "Alternative explanations adjusted for (model_full)", fontsize=6.3,
             fontweight="bold", va="top")
    for i, c in enumerate(conf):
        axA.text(0.8, 2.35 - i * 0.42, "• " + c, fontsize=6.0, va="top")
    arrow(axA, 8.35, 3.15, 8.35, 5.75, color=OI["grey"], ls="--")
    axA.text(0.6, 3.55, "Effect modifiers (post hoc): drug type, target family, therapeutic area",
             fontsize=5.8, style="italic", va="center")

    # Panel B: flows
    axB = fig.add_axes([0.5, 0.0, 0.5, 1.0])
    axB.set_xlim(0, 10)
    axB.set_ylim(0, 10)
    axB.axis("off")
    axB.text(0.1, 9.7, "B", fontsize=10, fontweight="bold", va="top")
    axB.text(0.7, 9.65, "Data flow", fontsize=7, va="top", fontweight="bold")

    def flow(x, w, items, title, fc, top=9.0, h=1.05, gap=0.42):
        axB.text(x + w / 2, top, title, fontsize=6.5, fontweight="bold", ha="center", va="top")
        y = top - 0.35
        prev = None
        for t, n in items:
            y0 = y - h
            axB.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                         fc=fc, ec="black", lw=0.6))
            axB.text(x + w / 2, y0 + h * 0.66, t, fontsize=5.5, ha="center", va="center")
            axB.text(x + w / 2, y0 + h * 0.28, n, fontsize=5.8, ha="center", va="center",
                     fontweight="bold")
            if prev is not None:
                axB.add_patch(FancyArrowPatch((x + w / 2, prev), (x + w / 2, y), arrowstyle="-|>",
                                              mutation_scale=7, color="black", lw=0.6))
            prev = y0
            y = y0 - gap
        return y

    flow(0.2, 4.8, faers, "FAERS (spontaneous reports)", "#DCEBF6")
    flow(5.2, 4.7, trials, "ClinicalTrials.gov (parallel analysis)", "#FBE3D2")
    save(fig, stem, ["results/RESULTS_SUMMARY.md (hard-coded flow numbers, §1, §3, §4, §9)"],
         {"faers_flow": faers, "trial_flow": trials})


# --------------------------------------------------------------------------- Fig 4
FAM_ORDER = ["kinase", "other membrane receptor", "GPCR", "ion channel", "nuclear receptor",
             "other enzyme", "transporter", "protease", "secreted protein", "surface antigen",
             "structural protein", "multi-family"]
AREA_ORDER = ["neoplasms", "nervous system", "mental disorders", "cardiovascular",
              "metabolic/nutritional", "urogenital", "respiratory", "skin", "immune"]


def fig4():
    stem = "fig4_forest_adjusted_modifiers"
    tgt = load_json("effect_modification_target.json")
    ind = load_json("effect_modification_indication.json")
    adj = {lab: load_json(f"effect_modification_adjusted_{lab}.json") for lab, _ in RANK}
    joint = {}
    joint_note = None
    for lab, _ in RANK:
        p = os.path.join(TAB, f"modifier_joint_{lab}.json")
        if os.path.exists(p):
            j = load_json(f"modifier_joint_{lab}.json")
            if j.get("repro_check", {}).get("status") == "PASS":
                joint[lab] = j
            else:
                joint_note = f"modifier_joint_{lab}.json repro_check not PASS; joint block omitted"
        else:
            joint_note = f"modifier_joint_{lab}.json missing; joint block omitted"
    if len(joint) < 2:
        joint = {}
        MANIFEST["skipped"]["fig4_joint_block"] = joint_note
        log("[SKIP] fig4 joint block:", joint_note)

    values = {}
    rows_by = {}
    for lab, _ in RANK:
        rows = [{"header": "Target family (molecule label)"}]
        for f in FAM_ORDER:
            v = tgt[lab]["family"].get(f)
            if v is None or "adjusted_same_target" not in v:
                continue
            a = v["adjusted_same_target"]
            if a.get("beta") is None:
                MANIFEST["skipped"][f"fig4_row_{lab}_{f}"] = f"no adjusted estimate stored: {a.get('note')}"
                continue
            rows.append({"label": f"{f} (n={v['same_target_pairs']:,})", "beta": a["beta"],
                         "ci": a["ci95"]})
        rows.append({"header": "Therapeutic area (trial-based label)"})
        for ar in AREA_ORDER:
            v = ind[lab]["trials"].get(ar)
            if v is None or "adjusted_same_target" not in v:
                MANIFEST["skipped"][f"fig4_row_{lab}_{ar}"] = "no adjusted estimate stored"
                continue
            a = v["adjusted_same_target"]
            rows.append({"label": f"{ar} (n={v['same_target_pairs']:,})", "beta": a["beta"],
                         "ci": a["ci95"]})
        rows.append({"header": "Drug type"})
        for key, nm in [("biologic_pairs", "biologic"), ("small_molecule_pairs", "small molecule"),
                        ("peptide_pairs", "peptide")]:
            s = adj[lab]["strata"][key]
            a = s["full"]["same_target"]
            rows.append({"label": f"{nm} (n={s['same_target_pairs']:,})", "beta": a["beta"],
                         "ci": a["ci95"]})
        if joint:
            j = joint[lab]
            rows.append({"header": "Joint model, same target × modifier (WP1)"})
            for mod, nm in [("bio_pair", "× biologic pair"), ("kinase_pair", "× kinase pair"),
                            ("onc_pair", "× neoplasm pair"), ("new_pair", "× approved ≥2000 pair")]:
                a = j["joint_model"][f"same_target_x_{mod}"]
                n = j["cell_counts"][mod]["modifier_1"]["same_target_pairs"]
                rows.append({"label": f"{nm} (n={n:,})", "beta": a["beta"], "ci": a["ci95"],
                             "marker": "D", "mfc": "white", "ms": 3.4})
            a = j["joint_model"]["same_target"]
            rows.append({"label": "same target, none of the four", "beta": a["beta"], "ci": a["ci95"],
                         "marker": "D", "mfc": "white", "ms": 3.4})
        rows_by[lab] = rows
        values[lab] = [{k: r[k] for k in ("label", "beta", "ci")} for r in rows if "header" not in r]
        values[lab + "_reference_whole_cohort"] = adj[lab]["paper_model_full_same_target"]

    n_rows = len(rows_by["prr_clean"])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.155 * n_rows + 1.1), sharey=False)
    xlims = {"prr_clean": (-1.5, 2.6), "ic025_clean2": (-2.0, 4.5)}
    for k, (lab, title) in enumerate(RANK):
        ax = axes[k]
        ref = adj[lab]["paper_model_full_same_target"]["beta"]
        forest_axis(ax, rows_by[lab], xlims[lab], RANK_COL[lab], show_labels=(k == 0),
                    ref=ref, ref_label=f"whole cohort {ref:.3f}")
        ax.set_title(title, color=RANK_COL[lab], pad=12)
        ax.set_xlabel("Adjusted same-target coefficient\n(extra shared signals per pair, 95% CI)")
    if joint:
        axes[1].text(1.0, -0.12 - 0.0, "Open diamonds (joint model): interaction = increment over the same-target "
                     "term; not a stratum coefficient", transform=axes[1].transAxes, fontsize=5.3,
                     ha="right", va="top", color="#444444")
    fig.text(0.0, -0.02 if joint else 0.0,
             "Post hoc, exploratory. Strata: both molecules carry the label; model_full within the stratum "
             "(molecule bootstrap, B = 500). n = same-target pairs. Arrows: CI beyond axis (value shown).",
             fontsize=5.3, color="#444444", ha="left", va="top")
    fig.subplots_adjust(wspace=0.08)
    srcs = ["results/tables/effect_modification_target.json",
            "results/tables/effect_modification_indication.json",
            "results/tables/effect_modification_adjusted_prr_clean.json",
            "results/tables/effect_modification_adjusted_ic025_clean2.json"]
    if joint:
        srcs += ["results/tables/modifier_joint_prr_clean.json", "results/tables/modifier_joint_ic025_clean2.json"]
    save(fig, stem, srcs, values)
    for lab, _ in RANK:
        log(f"  fig4 {lab}:")
        for r in values[lab]:
            log(f"    {r['label']:<45s} {fmt(r['beta'], r['ci'])}")


# --------------------------------------------------------------------------- Fig 5
def fig5():
    stem = "fig5_biologic_robustness"
    adj = {lab: load_json(f"effect_modification_adjusted_{lab}.json") for lab, _ in RANK}
    tgt = load_json("effect_modification_target.json")
    em = {lab: load_json(f"effect_modification_{lab}.json") for lab, _ in RANK}
    chk = load_json("effect_modification_checks.json")
    values = {}
    MOD_COL = {"biologic": OI["purple"], "small molecule": OI["green"], "peptide": OI["orange"],
               "all pairs": OI["black"]}

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.25, 0.55], hspace=0.62, wspace=0.08)

    # A: ladder
    strata = [("all_pairs_refit", "all pairs"), ("small_molecule_pairs", "small molecule"),
              ("biologic_pairs", "biologic"), ("peptide_pairs", "peptide")]
    steps = [("unadjusted", "unadjusted"), ("class_only", "class only"), ("full", "full model")]
    for k, (lab, title) in enumerate(RANK):
        ax = fig.add_subplot(gs[0, k])
        vals = []
        for si, (sk, sn) in enumerate(strata):
            s = adj[lab]["strata"][sk]
            for ti, (tk, tn) in enumerate(steps):
                if tk not in s or "same_target" not in s[tk]:
                    continue
                a = s[tk]["same_target"]
                x = ti + (si - 1.5) * 0.18
                ax.plot([x, x], a["ci95"], color=MOD_COL[sn], lw=1.0)
                ax.plot(x, a["beta"], "o", color=MOD_COL[sn], ms=3.4)
                vals.append({"stratum": sn, "step": tn, "beta": a["beta"], "ci": a["ci95"],
                             "same_target_pairs": s["same_target_pairs"]})
        values[f"A_{lab}"] = vals
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(range(3))
        ax.set_xticklabels([t for _, t in steps])
        ax.set_xlim(-0.5, 2.5)
        ax.set_title(title, color=RANK_COL[lab])
        if k == 0:
            ax.set_ylabel("Same-target coefficient\n(95% CI)")
            ax.text(-0.2, 1.08, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")
        if k == 0:
            hs = [plt.Line2D([], [], color=MOD_COL[sn], marker="o", ms=3.4, lw=1.0,
                             label=f"{sn} ({adj[lab]['strata'][sk]['same_target_pairs']:,})")
                  for sk, sn in strata]
            ax.legend(handles=hs, frameon=False, loc="lower left", ncol=2, title="stratum (same-target pairs)",
                      title_fontsize=6, fontsize=5.8)

    # B: robustness (unadjusted)
    for k, (lab, title) in enumerate(RANK):
        ax = fig.add_subplot(gs[1, k])
        a = adj[lab]["strata"]
        rows = [{"header": "Δ, molecule-averaged (B = 2000)"}]
        mod = em[lab]["modality"]
        st = tgt[lab]["single_target_only"]
        ck = chk[lab]
        for nm, key in [("biologic", "biologic"), ("small molecule", "small_molecule")]:
            v = mod[key]
            rows.append({"label": f"all molecules: {nm} (n={v['same_target_pairs']:,})", "beta": v["delta"],
                         "ci": v["ci95"], "color": MOD_COL[nm]})
        for nm, key in [("biologic", "biologic"), ("small molecule", "small_molecule")]:
            v = st[key]
            rows.append({"label": f"single-target only: {nm} (n={v['same_target_pairs']:,})",
                         "beta": v["delta"], "ci": v["ci95"], "color": MOD_COL[nm]})
        v = ck["bio_excl_top3_targets"]
        rows.append({"label": f"biologic without {', '.join(v['targets_dropped'])} (n={v['same_target_pairs']:,})",
                     "beta": v["delta"], "ci": v["ci95"], "color": MOD_COL["biologic"]})
        rows.append({"header": "Pair-mean difference, same − different target"})
        for ck_key, cname in [("no_shared_indication", "no shared indication"),
                              ("no_comention", "never co-mentioned"),
                              ("different_atc4", "different ATC-4")]:
            for nm, sk in [("biologic", "biologic_pairs"), ("small molecule", "small_molecule_pairs")]:
                v = a[sk]["within_class_contrasts"][ck_key]
                rows.append({"label": f"{cname}: {nm} (n={v['same_target_pairs']:,})",
                             "beta": v["difference"], "ci": v["ci95"], "color": MOD_COL[nm]})
        xl = (-0.5, 3.6) if lab == "prr_clean" else (-0.5, 6.0)
        forest_axis(ax, rows, xl, OI["black"], show_labels=(k == 0))
        ax.set_xlabel("Extra shared signals, same vs different target (95% CI)")
        if k == 0:
            ax.text(-1.15, 1.04, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")
        values[f"B_{lab}"] = [{kk: r[kk] for kk in ("label", "beta", "ci")} for r in rows if "header" not in r]

    # C: oncology by modality
    for k, (lab, title) in enumerate(RANK):
        ax = fig.add_subplot(gs[2, k])
        nb = tgt[lab]["neoplasms_by_modality"]
        rows = []
        for nm, key in [("biologic", "biologic"), ("small molecule", "small_molecule")]:
            v = nb[key]
            rows.append({"label": f"neoplasms: {nm} (n={v['same_target_pairs']:,})", "beta": v["delta"],
                         "ci": v["ci95"], "color": MOD_COL[nm]})
        xl = (-0.5, 3.6) if lab == "prr_clean" else (-0.5, 6.0)
        forest_axis(ax, rows, xl, OI["black"], show_labels=(k == 0))
        ax.set_xlabel("Δ within trial-labelled neoplasm molecules (95% CI)")
        if k == 0:
            ax.text(-1.15, 1.12, "C", transform=ax.transAxes, fontsize=10, fontweight="bold")
        values[f"C_{lab}"] = [{kk: r[kk] for kk in ("label", "beta", "ci")} for r in rows]
    fig.text(0.0, 0.035, "Post hoc, exploratory. A: model_full within stratum (B = 500). B, C: unadjusted. "
             "n = same-target pairs.",
             fontsize=5.3, color="#444444", ha="left", va="top")
    save(fig, stem, ["results/tables/effect_modification_adjusted_{prr_clean,ic025_clean2}.json",
                     "results/tables/effect_modification_{prr_clean,ic025_clean2}.json",
                     "results/tables/effect_modification_target.json",
                     "results/tables/effect_modification_checks.json"], values)
    for key, vals in values.items():
        log(f"  fig5 {key}:")
        for r in vals:
            lab = r.get("label") or f"{r['stratum']} / {r['step']}"
            log(f"    {lab:<62s} {fmt(r['beta'], r['ci'])}")


# --------------------------------------------------------------------------- Fig 6
def fig6b():
    stem = "fig6b_positive_control"
    pcs = {lab: load_json(f"positive_control_R2_{lab}.json") for lab, _ in RANK}
    fig, ax = plt.subplots(figsize=(3.4, 1.7))
    values = {}
    for i, (lab, title) in enumerate(RANK):
        p = pcs[lab]
        y = -i
        ax.plot(p["ci95"], [y, y], color=RANK_COL[lab], lw=1.2)
        ax.plot(p["difference"], y, "o", color=RANK_COL[lab], ms=4)
        ax.text(p["ci95"][1] + 0.05, y,
                f"own {p['mean_own_overlap']:.3f} vs other {p['mean_other_overlap']:.3f}; "
                f"own > other {p['pct_molecules_own_gt_other']:.1f}%",
                fontsize=5.5, va="center", color="#333333")
        values[lab] = p
    ax.set_yticks([0, -1])
    ax.set_yticklabels(["PRR, standard", "IC025, extended"])
    for t, (lab, _) in zip(ax.get_yticklabels(), RANK):
        t.set_color(RANK_COL[lab])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlim(-0.1, 3.6)
    ax.set_ylim(-1.6, 0.6)
    ax.set_xlabel("Own minus other molecules' trial profile\n(shared signals with FAERS profile, 95% CI)")
    ax.set_title(f"B  Positive control: {pcs['prr_clean']['molecules_in_both']} molecules in both cohorts",
                 loc="left", fontweight="bold")
    save(fig, stem, ["results/tables/positive_control_R2_prr_clean.json",
                     "results/tables/positive_control_R2_ic025_clean2.json"], values)
    for lab, p in values.items():
        log(f"  fig6b {lab}: diff {fmt(p['difference'], p['ci95'])}; own {p['mean_own_overlap']}; "
            f"other {p['mean_other_overlap']}; own>other {p['pct_molecules_own_gt_other']}%")


def fig6c():
    stem = "fig6c_trial_modifiers"
    d = load_json("modifier_trials.json")
    if d.get("repro_check", {}).get("status") != "PASS":
        raise RuntimeError("modifier_trials.json repro_check not PASS")
    A = d["A_stratum_delta"]
    B = d["B_adjusted_ols"]
    sets = [("all_pairs", "all pairs", OI["black"], "o"),
            ("excluding_shared_condition_pairs", "excl. shared-condition pairs", OI["sky"], "s")]
    groups = [("modality", [("biologic", "biologic"), ("small_molecule", "small molecule"),
                            ("difference_biologic_minus_small_molecule", "biologic − small molecule")]),
              ("family", [("kinase", "kinase"), ("GPCR", "GPCR"), ("other", "other families"),
                          ("difference_kinase_minus_GPCR", "kinase − GPCR")]),
              ("area", [("neoplasms", "neoplasms"), ("other", "other areas"),
                        ("difference_neoplasms_minus_other", "neoplasms − other")])]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [1.35, 1.0],
                                                                          "wspace": 0.75})
    values = {"A": [], "B": []}
    y = 0
    yt, yl = [], []
    xl = (-1.6, 2.6)
    for g, items in groups:
        for key, nm in items:
            for si, (sk, sn, col, mk) in enumerate(sets):
                v = A[sk][g][key]
                b = v.get("delta", v.get("diff"))
                ci = v["ci95"]
                n = v.get("same_target_pairs")
                yy = -y - si * 0.32
                low_n = "note" in v
                ax1.plot([max(ci[0], xl[0]), min(ci[1], xl[1])], [yy, yy], color=col, lw=1.0,
                         ls=":" if low_n else "-")
                ax1.plot(b, yy, mk, color=col, ms=3.2, mfc="white" if low_n else col, mec=col)
                if n is not None:
                    ax1.text(xl[1], yy, f"n={n}", fontsize=5, va="center", ha="right", color=col)
                values["A"].append({"pair_set": sn, "group": g, "label": nm, "estimate": b, "ci": ci,
                                    "same_target_pairs": n, "note": v.get("note")})
            yt.append(-y - 0.16)
            yl.append(nm)
            y += 1
        y += 0.5
    ax1.set_yticks(yt)
    ax1.set_yticklabels(yl)
    ax1.tick_params(axis="y", length=0)
    ax1.spines["left"].set_visible(False)
    ax1.axvline(0, color="black", lw=0.6)
    ax1.set_xlim(*xl)
    ax1.set_xlabel("Δ in trial-reported AE profiles (95% CI)")
    ax1.set_title("C  Modifiers in trial data (WP2): stratum Δ", loc="left", fontweight="bold")
    hs = [plt.Line2D([], [], color=c, marker=m, ms=3.2, lw=1.0, label=n) for _, n, c, m in sets]
    hs.append(plt.Line2D([], [], color=OI["grey"], marker="o", mfc="white", ms=3.2, lw=1.0, ls=":",
                         label="<10 same-target pairs, not interpreted"))
    ax1.legend(handles=hs, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=2)

    # adjusted
    terms = [("same_target", "same target"), ("same_target_x_bio_pair", "× biologic pair"),
             ("same_target_x_kinase_pair", "× kinase pair"), ("same_target_x_onc_pair", "× neoplasm pair")]
    cells = B["cells"]
    for i, (t, nm) in enumerate(terms):
        a = B["model"][t]
        ax2.plot(a["ci95"], [-i, -i], color=OI["black"], lw=1.0)
        ax2.plot(a["beta"], -i, "D", color=OI["black"], mfc="white" if i else OI["black"], ms=3.4)
        mod = t.replace("same_target_x_", "")
        n = cells[mod]["same_target=1,%s=1" % mod] if mod in cells else B["same_target_pairs"]
        values["B"].append({"term": t, "beta": a["beta"], "ci": a["ci95"], "same_target_pairs_in_cell": n})
    ax2.set_yticks([-i for i in range(len(terms))])
    ax2.set_yticklabels([f"{nm} (n={v['same_target_pairs_in_cell']:,})" for (t, nm), v in zip(terms, values["B"])])
    ax2.tick_params(axis="y", length=0)
    ax2.spines["left"].set_visible(False)
    ax2.axvline(0, color="black", lw=0.6)
    ax2.set_xlim(-1.6, 2.0)
    ax2.set_ylim(-len(terms) + 0.4, 0.8)
    ax2.set_xlabel("Coefficient (95% CI)")
    ax2.set_title("Adjusted trial-pair OLS", loc="left", fontweight="bold")
    ax2.text(0.0, -0.3, f"{B['pairs']:,} pairs; class, EPC, shared trial condition adjusted.\n"
             "No report-volume or co-mention covariate exists for trials.\n"
             "Interaction = increment over the same-target term.",
             transform=ax2.transAxes, fontsize=5.3, color="#444444", va="top")
    fig.text(0.0, -0.27, f"Post hoc (amendment 3, WP2). Trial cohort {d['n_molecules']} molecules; "
             f"whole-cohort Δ {A['all_pairs']['whole_cohort_delta']} (all pairs), "
             f"{A['excluding_shared_condition_pairs']['whole_cohort_delta']} (excl. shared condition). "
             "Single trial-profile ranking.", fontsize=5.3, color="#444444")
    save(fig, stem, ["results/tables/modifier_trials.json"], values)
    for r in values["A"]:
        log(f"  fig6c A {r['pair_set']:<30s} {r['label']:<28s} {fmt(r['estimate'], r['ci'])} n={r['same_target_pairs']}"
            + (f" [{r['note']}]" if r["note"] else ""))
    for r in values["B"]:
        log(f"  fig6c B {r['term']:<28s} {fmt(r['beta'], r['ci'])} n={r['same_target_pairs_in_cell']}")


# --------------------------------------------------------------------------- S3
FAM_COL = {"GPCR": OI["blue"], "kinase": OI["verm"], "ion channel": OI["green"],
           "nuclear receptor": OI["orange"], "other enzyme": OI["sky"],
           "other membrane receptor": OI["purple"], "other": OI["grey"]}


def figS3():
    stem = "figS3_target_waterfall"
    data = {}
    for lab, _ in RANK:
        p = os.path.join(TAB, f"target_level_R2_{lab}.csv")
        if not os.path.exists(p):
            raise FileNotFoundError(f"missing input: results/tables/target_level_R2_{lab}.csv")
        with open(p) as f:
            data[lab] = list(csv.DictReader(f))
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.8), sharex=True)
    values = {}
    for k, (lab, title) in enumerate(RANK):
        ax = axes[k]
        rows = sorted(data[lab], key=lambda r: -float(r["contrast"]))
        cs = [float(r["contrast"]) for r in rows]
        cols = [FAM_COL.get(r["family"], FAM_COL["other"]) for r in rows]
        ax.bar(range(len(cs)), cs, color=cols, width=0.85, lw=0)
        mean = sum(cs) / len(cs)
        npos = sum(c > 0 for c in cs)
        ax.axhline(0, color="black", lw=0.6)
        ax.axhline(mean, color="black", lw=0.7, ls="--")
        ax.text(len(cs) - 1, mean, f"mean {mean:.4f}", fontsize=5.5, ha="right", va="bottom")
        ax.set_ylabel("Per-target contrast\n(shared signals)")
        ax.set_title(f"{title}: {len(cs)} targets, {npos} positive", color=RANK_COL[lab], loc="left")
        ax.set_xlim(-1, len(cs))
        # label top 5 targets
        top5 = "; ".join(f"{rows[i]['target_symbol']} {cs[i]:.2f}" for i in range(5))
        ax.text(0.03, 0.93, f"Top 5 (ranks 1-5): {top5}", transform=ax.transAxes, fontsize=5.5,
                ha="left", va="top")
        ax.set_ylim(min(cs) - 0.5, max(cs) * 1.12)
        values[lab] = {"targets": len(cs), "mean_contrast": round(mean, 4), "n_positive": npos,
                       "min": min(cs), "max": max(cs), "top5": top5}
    axes[1].set_xlabel("Targets with ≥3 cohort members, ranked by contrast")
    hs = [plt.Rectangle((0, 0), 1, 1, color=c, label=f) for f, c in FAM_COL.items()]
    axes[0].legend(handles=hs, frameon=False, ncol=4, loc="upper right", bbox_to_anchor=(1.0, 0.86), fontsize=5.8,
                   title="target family", title_fontsize=6)
    fig.text(0.0, 0.0, "Contrast = mean overlap among partner members minus mean overlap with non-partners "
             "(pairs_all.py target_level). Amendment 3, WP5.", fontsize=5.3, color="#444444", va="top")
    fig.tight_layout()
    save(fig, stem, ["results/tables/target_level_R2_prr_clean.csv",
                     "results/tables/target_level_R2_ic025_clean2.csv"], values)
    for lab, v in values.items():
        log(f"  figS3 {lab}: {v}")


# --------------------------------------------------------------------------- S4
GROUPS = [("overall", "All same-target pairs"), ("kinase", "Kinase pairs"), ("GPCR", "GPCR pairs"),
          ("ion_channel", "Ion-channel pairs"), ("biologic", "Biologic pairs"), ("neoplasm", "Neoplasm pairs")]


def figS4():
    for lab, title in RANK:
        stem = f"figS4_shared_terms_{lab}"
        try:
            d = load_json(f"shared_terms_by_group_{lab}.json")
            rc = d.get("reproduction_check", {}).get("status")
            if rc != "PASS":
                raise RuntimeError(f"shared_terms_by_group_{lab}.json reproduction_check={rc}")
            fig, axes = plt.subplots(3, 2, figsize=(8.6, 11.5))
            values = {}
            for ax, (g, gname) in zip(axes.flat, GROUPS):
                v = d["shared_terms"][g]
                top = v["top25"]
                terms = [t[0] for t in top][::-1]
                cnts = [t[1] for t in top][::-1]
                ax.barh(range(len(cnts)), cnts, color=RANK_COL[lab], height=0.72, lw=0)
                ax.set_yticks(range(len(cnts)))
                ax.set_yticklabels(terms, fontsize=4.9)
                ax.tick_params(axis="y", length=0)
                ax.set_ylim(-0.6, len(cnts) - 0.4)
                ax.set_title(f"{gname} ({v['same_target_pairs']:,} pairs; {v['pairs_sharing_none']:,} share none)",
                             loc="left", fontsize=6.5)
                ax.set_xlabel("same-target pairs sharing the term", fontsize=6)
                values[g] = {"same_target_pairs": v["same_target_pairs"],
                             "pairs_sharing_none": v["pairs_sharing_none"], "top25": top}
            fig.suptitle(f"Most frequently shared reactions among same-target pairs, {title}",
                         x=0.0, ha="left", fontsize=7.5, fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.985))
            fig.subplots_adjust(wspace=1.45, hspace=0.32)
            save(fig, stem, [f"results/tables/shared_terms_by_group_{lab}.json"], values)
            for g, v in values.items():
                log(f"  figS4 {lab} {g}: pairs {v['same_target_pairs']}, none {v['pairs_sharing_none']}, "
                    f"top1 {v['top25'][0]}")
        except Exception as e:  # noqa: BLE001
            skip(stem, f"{type(e).__name__}: {e}")


# --------------------------------------------------------------------------- main
def main():
    log("F1 make_figures_v3 — reads results/tables/ only; repro_check: not applicable (figures).")
    # Report upstream repro statuses that feed these figures
    ups = {}
    for name, path in [("WP1 PRR", "modifier_joint_prr_clean.json"), ("WP1 IC025", "modifier_joint_ic025_clean2.json"),
                       ("WP2", "modifier_trials.json"), ("WP5 PRR", "shared_terms_by_group_prr_clean.json"),
                       ("WP5 IC025", "shared_terms_by_group_ic025_clean2.json")]:
        try:
            d = load_json(path)
            rc = d.get("repro_check") or d.get("reproduction_check")
            ups[name] = rc
        except Exception as e:  # noqa: BLE001
            ups[name] = f"not available: {e}"
        log(f"  upstream {name}: {ups[name]}")
    MANIFEST["upstream_repro_checks"] = ups
    for fn, stem in [(fig1, "fig1_design_flow"), (fig4, "fig4_forest_adjusted_modifiers"),
                     (fig5, "fig5_biologic_robustness"), (fig6b, "fig6b_positive_control"),
                     (fig6c, "fig6c_trial_modifiers"), (figS3, "figS3_target_waterfall")]:
        try:
            fn()
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            skip(stem, f"{type(e).__name__}: {e}")
    figS4()
    with open(os.path.join(OUT, "figures_v3_manifest.json"), "w") as f:
        json.dump(MANIFEST, f, indent=1, default=str)
    log(f"made {len(MANIFEST['made'])}: {sorted(MANIFEST['made'])}")
    log(f"skipped {len(MANIFEST['skipped'])}: {MANIFEST['skipped']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
