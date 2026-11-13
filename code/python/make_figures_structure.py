#!/usr/bin/env python3
"""Supplementary figures for the biologic structure analysis (S9-S11) and the reporting-artefact test (S12).
Reads only results/tables/biologic_structure.json and era_reporting.json; style and forest helper from make_figures_v3.py.
Writes results/figures/v3/figS9..figS12 (PNG 300 dpi + SVG) and adds them to figures_v3_manifest.json."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import make_figures_v3 as F          # importing only defines helpers and style (its figures run under __main__)
plt, OI, RANK, RANK_COL = F.plt, F.OI, F.RANK, F.RANK_COL
S = F.load_json("biologic_structure.json")
GREY = OI["grey"]

def crow(lab, key, label, kind="contrast"):
    c = S[lab]["contrasts"][key]
    small = min(c["pairs_1"], c["pairs_0"]) < 10
    r = dict(label=f"{label} ({c['pairs_1']:,} vs {c['pairs_0']:,})", beta=c["diff"], ci=c["ci95"])
    if small: r.update(color=GREY, mfc="white", label=r["label"] + " *")
    return r

ROWS = [("header", "Generic structure, different-target biologic pairs"),
        ("H2a_isotype_difftarget_same_vs_diff", "same vs different IgG isotype"),
        ("H3a_origin_difftarget_same_vs_diff", "same vs different origin"),
        ("H4a_format_difftarget_same_vs_diff", "same vs different format"),
        ("H5b_difftarget_same_vs_diff_class", "same vs different target class"),
        ("header", "Same-target biologic pairs, structure differs"),
        ("H2b_isotype_sametarget_same_vs_diff", "same vs different isotype"),
        ("H3b_origin_sametarget_same_vs_diff", "same vs different origin"),
        ("H4b_format_sametarget_same_vs_diff", "same vs different format"),
        ("H5a_sametarget_membrane_vs_secreted", "membrane vs secreted target"),
        ("header", "Mechanism-specific cells"),
        ("H1a_ADC_diff_antigen_same_vs_diff_payload_class", "ADCs, different antigen: same vs different payload"),
        ("H1b_ADC_x_smallmol_payload_target_vs_no_shared_target", "ADC x small molecule at payload target vs none"),
        ("H4_desc_CD3E_pairs_TCE_TCE_vs_TCE_monospecific", "CD3E: T-cell engager pair vs engager x anti-CD3"),
        ("H4_desc_CD19_pairs_CART_CART_vs_CART_antibody", "CD19: CAR-T pair vs CAR-T x antibody"),
        ("H4_desc_difftarget_CART_CART_vs_other_difftarget", "CAR-T pairs, different targets vs other"),
        ("H5_desc_PD1_x_PDL1_vs_other_difftarget_antibody", "PD-1 x PD-L1 vs other different-target antibodies")]

def figS9():
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 5.6), sharey=False, gridspec_kw=dict(wspace=0.08))
    vals = {}
    for k, (lab, title) in enumerate(RANK):
        rows = []
        for key, label in ROWS:
            rows.append({"header": label} if key == "header" else crow(lab, key, label))
        add = S["post_spec_addendum"][lab]
        rows.append({"header": "Biologic same-target Delta, products removed (post-spec addendum)"})
        for nm, label in (("biologic_all", "all biologics"), ("excl_CAR_T_TCE", "without CAR-T and T-cell engagers"),
                          ("excl_CAR_T_TCE_ADC", "without CAR-T, engagers and ADCs")):
            rows.append(dict(label=f"{label} ({add[nm]['same_target_pairs']} pairs)", beta=add[nm]["delta"], ci=add[nm]["ci95"], marker="D", mfc="white"))
        xl = (-2.2, 6.0) if lab == "prr_clean" else (-3.5, 10.0)
        F.forest_axis(axes[k], rows, xl, RANK_COL[lab], show_labels=(k == 0))
        axes[k].set_title(title, color=RANK_COL[lab], pad=8)
        axes[k].set_xlabel("Difference in shared signals per pair (95% CI)\n(Delta for the last block)")
        vals[lab] = [{"label": r["label"], "beta": r["beta"], "ci": r["ci"]} for r in rows if "header" not in r]
    fig.text(0.5, -0.04, "Post hoc structure analysis (spec fixed before results, sha256 6bad6214...). Pair counts: group 1 vs group 0.\n"
             "* fewer than 10 pairs in a group: shown, not interpreted. Molecule-level bootstrap, B = 2000.", ha="center", fontsize=5.8, color="#444444")
    F.save(fig, "figS9_biologic_structure_contrasts", ["results/tables/biologic_structure.json"], vals)

def figS10():
    a = S["annotation"]
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.2), gridspec_kw=dict(wspace=0.55, hspace=0.45))
    for ax, key, title in zip(axes.flat, ("format", "isotype_used", "origin_used", "payload"), ("Format", "IgG isotype", "Origin", "ADC payload")):
        d = sorted(a[key].items(), key=lambda x: -x[1])
        ax.barh([k.replace("unknown/not used", "unknown").replace("_", " ") for k, _ in d][::-1], [v for _, v in d][::-1], color=OI["purple"], height=0.6)
        for i, (_, v) in enumerate(d[::-1]): ax.text(v, i, f" {v}", va="center", fontsize=5.5)
        ax.set_title(title); ax.tick_params(axis="y", length=0); ax.set_xticks([]); ax.set_xlim(0, max(v for _, v in d) * 1.18)
        ax.spines["bottom"].set_visible(False)
    fig.text(0.5, -0.02, f"{a['n_biologics']} cohort biologics; {a['thera_matched']} matched to Thera-SAbDab (downloaded 2026-11-13). "
             f"Origin: Thera-SAbDab vs INN stem agree {a['origin_qc_thera_vs_inn']['agree']}/{a['origin_qc_thera_vs_inn']['n']}. "
             "Unmatched = unknown, never imputed.", ha="center", fontsize=5.8, color="#444444")
    F.save(fig, "figS10_biologic_structure_annotation", ["results/tables/biologic_structure.json"], {k: a[k] for k in ("format", "isotype_used", "origin_used", "payload")})

def figS11():
    cells = [("H4_desc_CD3E_pairs_TCE_TCE_vs_TCE_monospecific", "CD3E T-cell engager pairs"),
             ("H4_desc_difftarget_CART_CART_vs_other_difftarget", "CAR-T pairs, different targets *"),
             ("H1b_ADC_x_smallmol_payload_target_vs_no_shared_target", "ADC x small molecule at payload target"),
             ("H5_desc_PD1_x_PDL1_vs_other_difftarget_antibody", "PD-1 x PD-L1 antibody pairs")]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 4.6), gridspec_kw=dict(wspace=1.25, hspace=0.6)); vals = {}
    for ax, (key, title) in zip(axes.flat, cells):
        c = S["prr_clean"]["contrasts"][key]; t = c.get("top_shared_terms_group_1", [])[:8]
        ax.barh([x[0].capitalize()[:52] for x in t][::-1], [x[1] for x in t][::-1], color=OI["blue"], height=0.6)
        ax.set_title(f"{title}\n({c['pairs_1']} pairs)", fontsize=6.3); ax.tick_params(axis="y", length=0, labelsize=5.3)
        ax.set_xlabel("pairs sharing the term", fontsize=5.8); vals[key] = t
    fig.text(0.5, -0.02, "PRR / standard cleaning profiles. Counts of pairs sharing each reaction; no inference. * fewer than 10 pairs.", ha="center", fontsize=5.8, color="#444444")
    F.save(fig, "figS11_mechanism_cell_shared_terms", ["results/tables/biologic_structure.json"], vals)

def figS12():
    E = F.load_json("era_reporting.json")["variants"]
    names = {"R1": "reports >= 3 y after approval", "R2": "FAERS era (>= 2014-09)", "R3": "health professionals only",
             "R4": "consumers only", "R5": "non-US reports only"}
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw=dict(wspace=0.08)); vals = {}
    for k, (lab, title) in enumerate(RANK):
        rows = []
        for term, head in (("same_target_x_bio_pair", "x biologic pair (joint model)"), ("same_target_x_new_pair", "x approved >= 2000 pair (joint model)")):
            rows.append({"header": head})
            for v in ("R1", "R2", "R3", "R4", "R5"):
                x = E[v][lab]["joint"][term]; ref = E[v][lab + "_reference_same_molecules"]["joint"][term]
                rows.append(dict(label=f"{names[v]} (n={E[v]['info']['eligible_molecules']})", beta=x["beta"], ci=x["ci95"]))
                rows.append(dict(label="   full-data reference, same molecules", beta=ref["beta"], ci=ref["ci95"], color=GREY, marker="D", mfc="white"))
        xl = (-0.6, 2.2) if lab == "prr_clean" else (-0.6, 3.6)
        F.forest_axis(axes[k], rows, xl, RANK_COL[lab], show_labels=(k == 0))
        axes[k].set_title(title, color=RANK_COL[lab], pad=8); axes[k].set_xlabel("Interaction coefficient (95% CI)")
        vals[lab] = [{"label": r["label"], "beta": r["beta"], "ci": r["ci"]} for r in rows if "header" not in r]
    fig.text(0.5, -0.06, "Amendment 5 (fixed before running). Profiles rebuilt per reporting variant; full-data rebuild reproduces the stored profiles exactly.\n"
             "Reference = full-data profiles on the same eligible molecules. Joint model A4.1b design, molecule bootstrap B = 500.", ha="center", fontsize=5.8, color="#444444")
    F.save(fig, "figS12_reporting_artefact_test", ["results/tables/era_reporting.json"], vals)

if __name__ == "__main__":
    for f in (figS9, figS10, figS11, figS12):
        try: f()
        except Exception as e: F.skip(f.__name__, repr(e))
    mp = os.path.join(F.OUT, "figures_v3_manifest.json"); m = json.load(open(mp)) if os.path.exists(mp) else {"made": {}, "skipped": {}}
    m["made"].update(F.MANIFEST["made"]); m["skipped"].update(F.MANIFEST["skipped"]); json.dump(m, open(mp, "w"), indent=1)
