#!/usr/bin/env python3
"""Main figures 2, 3, 6A and supplementary S1 in the v3 style (re-panelled from the saved results; replaces the
round-2 fig1-fig4 layout of make_paper_figures.py, which is kept unchanged). Reads only results/tables/*.json.
Writes results/figures/v3/fig2_overall_and_adjustment, fig3_what_explains, fig6a_trials, figS1_all_extractions."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import make_figures_v3 as F
plt, OI, RANK, RANK_COL, J = F.plt, F.OI, F.RANK, F.RANK_COL, F.load_json

def two_col(rows_fn, stem, xl, xlabel, note, sources, h):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, h), gridspec_kw=dict(wspace=0.08)); vals = {}
    for k, (lab, title) in enumerate(RANK):
        rows = rows_fn(lab)
        F.forest_axis(axes[k], rows, xl[lab], RANK_COL[lab], show_labels=(k == 0))
        axes[k].set_title(title, color=RANK_COL[lab], pad=8); axes[k].set_xlabel(xlabel)
        vals[lab] = [{"label": r["label"], "beta": r["beta"], "ci": r["ci"]} for r in rows if "header" not in r]
    fig.text(0.5, -0.04 - 0.25 / h, note, ha="center", fontsize=5.8, color="#444444")   # lower for short figures
    F.save(fig, stem, sources, vals)

def fig2():
    def rows(lab):
        p = J(f"primary_R2_b_{lab}.json")["shared_signals"]; t = J(f"target_vs_class_R2_{lab}.json")
        r = [{"header": "A  Overall difference (primary extraction b)"},
             dict(label=f"same-target minus different-target Δ ({J(f'primary_R2_b_{lab}.json')['same_target_pairs']:,} same-target pairs)", beta=p["delta"], ci=p["ci95"]),
             {"header": "B  Same-target coefficient, all 539,241 pairs"}]
        for mk, nm in (("model_unadjusted", "unadjusted"), ("model_class_only", "+ class (ATC 4/3/1, FDA EPC)"),
                       ("model_full", "+ indication, volume, year gap, co-mention")):
            d = t[mk]["same_target"]; r.append(dict(label=nm, beta=d["beta"], ci=d["ci95"], marker="D", mfc="white" if mk != "model_full" else None))
        return r
    two_col(rows, "fig2_overall_and_adjustment", {"prr_clean": (-0.1, 1.2), "ic025_clean2": (-0.1, 2.1)},
            "Extra shared top-20 signals per pair (95% CI)",
            "A: molecule-averaged Δ, molecule bootstrap. B: OLS on all pairs, molecule bootstrap B = 500. PRR = SAP v1 ranking; IC025 = amendment 2.",
            ["primary_R2_b_*.json", "target_vs_class_R2_*.json"], 2.6)

def fig3():
    def rows(lab):
        t = J(f"target_vs_class_R2_{lab}.json")["within_class_contrasts"]; c = J(f"comention_summary_{lab}.json")
        a = J(f"asof_summary_{lab}.json")["clean"]
        r = [{"header": "A  Within class and indication strata (same − different target)"}]
        for k, cl, nm in (("atc4", "same_class", "same ATC level-4 group"), ("atc4", "different_class", "different ATC level-4 group"),
                          ("epc", "same_class", "same FDA class"), ("epc", "different_class", "different FDA class"),
                          ("indication", "same_class", "shared indication"), ("indication", "different_class", "no shared indication")):
            d = t[k][cl]; r.append(dict(label=f"{nm} ({d['same_target_pairs']:,})", beta=d["difference"], ci=d["ci95"]))
        r.append({"header": "B  Co-reporting"})
        for key, nm in (("difference_before", "all reports"), ("difference_after", "reports naming both molecules removed")):
            d = c[key]; r.append(dict(label=nm, beta=d["estimate"], ci=d["ci95"]))
        r.append({"header": "C  Hindsight: first molecule's as-of profile vs successor"})
        for key, nm in (("eA", "vs volume-matched first molecules"), ("eB", "vs other successors"), ("eC", "vs class-matched first molecules")):
            d = a[key]; r.append(dict(label=nm, beta=d["estimate"], ci=d["ci95"]))
        return r
    two_col(rows, "fig3_what_explains", {"prr_clean": (-0.2, 2.0), "ic025_clean2": (-0.2, 3.4)},
            "Difference in shared signals (95% CI)",
            "A: pairs within / outside a shared class or indication (n = same-target pairs). B: amendment 2. C: SAP as-of analysis with deviation (cluster sign-flip CIs).",
            ["target_vs_class_R2_*.json", "comention_summary_*.json", "asof_summary_*.json"], 3.8)

def fig6a():
    fig, ax = plt.subplots(figsize=(4.6, 1.5)); rows = []
    for lab, nm in (("R2_ctgov_all", "all pairs"), ("R2_ctgov_no_shared_trial", "pairs sharing a trial removed"), ("R2_ctgov_no_shared_condition", "pairs sharing a trial condition removed")):
        o = J(f"primary_{lab}.json"); d = o["shared_signals"]
        rows.append(dict(label=f"{nm} ({o['same_target_pairs']:,})", beta=d["delta"], ci=d["ci95"], color=OI["black"]))
    F.forest_axis(ax, rows, (-0.1, 1.7), OI["black"])
    ax.set_xlabel("Δ in trial-reported AE profiles (95% CI)"); ax.set_title("A  Trial-reported adverse events (357 molecules)", loc="left")
    F.save(fig, "fig6a_trials", ["primary_R2_ctgov_*.json"], [{"label": r["label"], "beta": r["beta"], "ci": r["ci"]} for r in rows])

def figS1():
    ext = (("b", "primary + secondary suspect"), ("c", "primary suspect only"), ("d", "health-professional reports"), ("e", "sole suspect"), ("balltargets", "all ChEMBL targets"))
    def rows(lab):
        r = [{"header": "Extraction"}]
        for e, nm in ext:
            d = J(f"primary_R2_{e}_{lab}.json")["shared_signals"]; r.append(dict(label=nm, beta=d["delta"], ci=d["ci95"]))
        r.append({"header": "Ranking variants (extraction b)"})
        for v, nm in ((f"R2_b_{'ror025' if lab == 'prr_clean' else 'ror025'}_clean", "ROR025, standard cleaning"), ("R2_b_ror025_clean2", "ROR025, extended cleaning")):
            d = J(f"primary_{v}.json")["shared_signals"]; r.append(dict(label=nm, beta=d["delta"], ci=d["ci95"], marker="D", mfc="white"))
        return r
    two_col(rows, "figS1_all_extractions", {"prr_clean": (-0.1, 2.6), "ic025_clean2": (-0.1, 2.6)},
            "Δ, same-target minus different-target (95% CI)", "All extractions and target sets; ROR025 rows are identical in both columns (ranking-independent of PRR/IC025).",
            ["primary_R2_*.json"], 3.0)

for f in (fig2, fig3, fig6a, figS1):
    try: f()
    except Exception as e: F.skip(f.__name__, repr(e))
import json
mp = os.path.join(F.OUT, "figures_v3_manifest.json"); m = json.load(open(mp)) if os.path.exists(mp) else {"made": {}, "skipped": {}}
m["made"].update(F.MANIFEST["made"]); m["skipped"].update(F.MANIFEST["skipped"]); json.dump(m, open(mp, "w"), indent=1)
