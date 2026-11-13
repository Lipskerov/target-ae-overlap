#!/usr/bin/env python3
"""Post hoc results (2026-11-13) -> the manuscript's single number source.

1. Regenerates results/RESULTS_SUMMARY.md with make_results_summary_r2.py (sections 1-9, unchanged), then appends
   section 10 "Post hoc analyses" generated here from the saved JSON files.
2. Writes results/POST_HOC_COMPARISONS.md: one row per post hoc comparison (estimate and 95% CI per ranking,
   same-target pairs, source file, specification). This is the supplementary list of every comparison run.
Reads only results/tables/*.json. Every number is copied from those files.
"""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
T = os.path.join(PROJ, "results/tables"); J = lambda f: json.load(open(os.path.join(T, f)))
L = ("prr_clean", "ic025_clean2")
rows = []   # (block, spec, comparison, kind, {label: (est, ci)}, n_same, source)

def est(d):
    if d is None: return None
    for k in ("delta", "diff", "beta", "difference"):
        if k in d and d[k] is not None: return (d[k], d.get("ci95"))
    return None
def add(block, spec, comp, kind, vals, n, src):
    if any(v is not None and v[0] is not None for v in vals.values()): rows.append((block, spec, comp, kind, vals, n, src))  # descriptive sub-tables without a single estimate are skipped
def both(f, getter):
    out = {}
    for lab in L:
        try: out[lab] = est(getter(J(f.format(lab=lab))))
        except (KeyError, TypeError, FileNotFoundError): out[lab] = None
    return out

# ---- E1: one-at-a-time strata (effect_modification_*)
PH = "post hoc (E1)"
em = {lab: J(f"effect_modification_{lab}.json") for lab in L}
for fac, ref in (("modality", "small_molecule"), ("route", "oral_only"), ("atc1", "A"), ("era", "pre2000")):
    for lv, r in em["prr_clean"][fac].items():
        if "__all_partners" in lv: continue
        add(f"1 {fac}", PH, f"{fac} = {lv}: Δ within stratum", "Δ", {lab: est(em[lab][fac].get(lv)) for lab in L}, r["same_target_pairs"], "effect_modification_*.json")
    for lv in em["prr_clean"].get(f"{fac}_differences_vs_{ref}", {}):
        add(f"1 {fac}", PH, f"{fac}: {lv} − {ref}", "difference of Δ", {lab: est(em[lab][f"{fac}_differences_vs_{ref}"][lv]) for lab in L}, None, "effect_modification_*.json")
# ---- E1: adjusted within modality + pooled interaction
ad = {lab: J(f"effect_modification_adjusted_{lab}.json") for lab in L}
for st, r in ad["prr_clean"]["strata"].items():
    for model in ("unadjusted", "class_only", "full"):
        if model in r:
            add("2 adjusted by modality", PH, f"{st}: same-target coefficient, {model}", "coefficient", {lab: est(ad[lab]["strata"][st].get(model, {}).get("same_target")) for lab in L}, r["same_target_pairs"], "effect_modification_adjusted_*.json")
for term in ("same_target", "target_x_bio", "target_x_mixed"):
    add("2 adjusted by modality", PH, f"pooled interaction model: {term}", "coefficient", {lab: est(ad[lab]["pooled_interaction"].get(term)) for lab in L}, None, "effect_modification_adjusted_*.json")
# ---- E2: therapeutic area
ind = J("effect_modification_indication.json")
for area, r in ind["prr_clean"]["trials"].items():
    add("3 therapeutic area", "post hoc (E2)", f"area = {area}: Δ", "Δ", {lab: est(ind[lab]["trials"].get(area)) for lab in L}, r["same_target_pairs"], "effect_modification_indication.json")
    if "adjusted_same_target" in r:
        add("3 therapeutic area", "post hoc (E2)", f"area = {area}: adjusted same-target coefficient", "coefficient", {lab: est(ind[lab]["trials"][area].get("adjusted_same_target")) for lab in L}, r["same_target_pairs"], "effect_modification_indication.json")
# ---- E3/E4: target family, oncology by modality, single target
tg = J("effect_modification_target.json")
for fam, r in tg["prr_clean"]["family"].items():
    add("4 target family", "post hoc (E3)", f"family = {fam}: Δ", "Δ", {lab: est(tg[lab]["family"].get(fam)) for lab in L}, r["same_target_pairs"], "effect_modification_target.json")
    if "adjusted_same_target" in r and r["adjusted_same_target"].get("beta") is not None:
        add("4 target family", "post hoc (E3)", f"family = {fam}: adjusted same-target coefficient", "coefficient", {lab: est(tg[lab]["family"][fam].get("adjusted_same_target")) for lab in L}, r["same_target_pairs"], "effect_modification_target.json")
for mo in ("small_molecule", "biologic"):
    add("4 target family", "post hoc (E3)", f"neoplasms, {mo}: Δ", "Δ", {lab: est(tg[lab]["neoplasms_by_modality"][mo]) for lab in L}, tg["prr_clean"]["neoplasms_by_modality"][mo]["same_target_pairs"], "effect_modification_target.json")
for k in ("small_molecule", "biologic", "difference"):
    add("4 target family", "post hoc (E4)", f"single-target molecules only: {k}", "Δ" if k != "difference" else "difference of Δ", {lab: est(tg[lab]["single_target_only"][k]) for lab in L}, tg["prr_clean"]["single_target_only"].get(k, {}).get("same_target_pairs"), "effect_modification_target.json")
# ---- A3 WP1 joint + A4.2
for f, spec in (("modifier_joint_{lab}.json", "amendment 3 WP1"), ("modifier_joint_a4_{lab}.json", "amendment 4 A4.2")):
    for term in ("same_target", "same_target_x_bio_pair", "same_target_x_kinase_pair", "same_target_x_onc_pair", "same_target_x_new_pair"):
        add("5 joint modifier model", spec, f"joint model: {term}", "coefficient", both(f, lambda o, t=term: o["joint_model"][t]), None, f.replace("{lab}", "*"))
# ---- A3 WP2 trials (single trial profile)
tr = J("modifier_trials.json")
for pset in ("all_pairs", "excluding_shared_condition_pairs"):
    for fac, d in tr["A_stratum_delta"][pset].items():
        if not isinstance(d, dict): continue
        for k, v in d.items():
            n = v.get("same_target_pairs") if isinstance(v, dict) else None
            add("6 trial-reported AEs", "amendment 3 WP2", f"{pset}: {fac} {k}", "Δ" if "difference" not in k else "difference of Δ", {"trial": est(v)}, n, "modifier_trials.json")
for term in ("same_target", "same_target_x_bio_pair", "same_target_x_kinase_pair", "same_target_x_onc_pair"):
    add("6 trial-reported AEs", "amendment 3 WP2", f"adjusted trial-pair model: {term}", "coefficient", {"trial": est(tr["B_adjusted_ols"]["model"].get(term))}, None, "modifier_trials.json")
# ---- A3 WP3 / A4.3
for f, spec in (("action_direction_{lab}.json", "amendment 3 WP3"), ("action_direction_a4_{lab}.json", "amendment 4 A4.3")):
    for scope in ("all_targets", "gpcr_nuclear_receptor"):
        d0 = J(f.format(lab="prr_clean"))[scope]
        for key in [k for k in d0 if isinstance(d0[k], dict) and est(d0[k])]:
            add("7 opposite action", spec, f"{scope}: {key}", "difference", both(f, lambda o, s=scope, k=key: o[s][k]), sum(d0["pair_counts"].values()), f.replace("{lab}", "*"))
# ---- A3 WP4 grid
for g in J("sensitivity_grid_R2.json")["grid"]:
    lab = "prr_clean" if str(g.get("ranking", "")).startswith("prr") else "ic025_clean2"
    add("8 profile length / floor", "amendment 3 WP4", f"K = {g['K']}, floor {g['floor']}", "Δ", {lab: (g["delta"], [g["ci_lo"], g["ci_hi"]])}, g.get("same_pairs"), "sensitivity_grid_R2.json")
# ---- A3 WP5 / A4.4
mb = J("target_level_molboot_a4.json")
add("9 per-target contrast", "amendment 3 WP5 / A4.4", "mean per-target contrast (molecule bootstrap CI)", "mean", {lab: (mb[lab]["point"], mb[lab]["ci95_molecule_bootstrap"]) for lab in L}, None, "target_level_molboot_a4.json")
# ---- A5 reporting variants
er = J("era_reporting.json")["variants"]
for v, vo in er.items():
    if v == "R0": continue
    for term in ("same_target_x_bio_pair", "same_target_x_new_pair"):
        add("10 reporting artefact", "amendment 5", f"{v}: {term}", "coefficient", {lab: est(vo[lab]["joint"].get(term)) for lab in L}, vo["prr_clean"].get("joint_same_target_pairs"), "era_reporting.json")
        add("10 reporting artefact", "amendment 5", f"{v}: {term} — full-data reference, same molecules", "coefficient", {lab: est(vo[lab + "_reference_same_molecules"]["joint"].get(term)) for lab in L}, None, "era_reporting.json")
    add("10 reporting artefact", "amendment 5", f"{v}: biologic − small molecule Δ", "difference of Δ", {lab: est(vo[lab]["biologic_minus_small_molecule"]) for lab in L}, None, "era_reporting.json")
bs = J("biologic_subtypes.json")
for st, r in bs["prr_clean"].items():
    add("11 biologic subtypes", "amendment 5 B-part", f"{st}: Δ", "Δ", {lab: est(bs[lab].get(st)) for lab in L}, r["same_target_pairs"], "biologic_subtypes.json")
# ---- structure test
sx = J("biologic_structure.json")
for key, c in sx["prr_clean"]["contrasts"].items():
    add("12 biologic structure", "structure SPEC", key, "difference", {lab: est(sx[lab]["contrasts"].get(key)) for lab in L}, f"{c.get('pairs_1')} vs {c.get('pairs_0')}", "biologic_structure.json")
for k, r in sx["post_spec_addendum"]["prr_clean"].items():
    add("12 biologic structure", "structure post-spec addendum", f"biologic Δ: {k}", "Δ", {lab: est(sx["post_spec_addendum"][lab].get(k)) for lab in L}, r.get("same_target_pairs"), "biologic_structure.json")

def f(v):
    if not v or v[0] is None: return "—"
    e, ci = v
    if e != e: return "n/a"
    return f"{e:.3f}" + (f" ({ci[0]:.3f} to {ci[1]:.3f})" if ci and ci[0] == ci[0] else "")
lines = ["# Every post hoc comparison (2026-11-13)", "",
         "Generated by `code/python/make_post_hoc_summary.py` from `results/tables/*.json`; do not edit by hand.",
         "No multiplicity correction was applied. Strata with fewer than 10 same-target pairs are listed but not interpreted.",
         "Specifications: E1–E4 `LIMITATIONS_AND_DEVIATIONS.md` (b2); amendments 3–5 `protocol/AMENDMENT_3–5.md`; structure `private/structure/SPEC.md`.", "",
         "| # | block | specification | comparison | estimate | PRR (95% CI) | IC025 (95% CI) | trial (95% CI) | same-target pairs | source |", "|---|---|---|---|---|---|---|---|---|---|"]
for i, (b, s, c, k, v, n, src) in enumerate(rows, 1):
    lines.append(f"| {i} | {b} | {s} | {c} | {k} | {f(v.get('prr_clean'))} | {f(v.get('ic025_clean2'))} | {f(v.get('trial'))} | {n if n is not None else ''} | `{src}` |")
lines += ["", f"**Total rows: {len(rows)}** (each row = one comparison, reported under one or both rankings or in trial data)."]
open(os.path.join(PROJ, "results/POST_HOC_COMPARISONS.md"), "w").write("\n".join(lines) + "\n")

# ---- RESULTS_SUMMARY: regenerate sections 1-9 then append section 10
subprocess.run([sys.executable, os.path.join(HERE, "make_results_summary_r2.py")], check=True, cwd=PROJ)
blocks = {}
for r in rows: blocks.setdefault(r[0], []).append(r)
out = ["", "## 10. Post hoc analyses, 2026-11-13 (all post hoc relative to SAP v1; specifications fixed before each analysis where stated)", "",
       "Generated by `code/python/make_post_hoc_summary.py`. The full list with sources is in `results/POST_HOC_COMPARISONS.md`. "
       "Interpretation and context: `results/EFFECT_MODIFICATION.md` §1–12.", ""]
for b, rs in blocks.items():
    out += [f"### 10.{b}", "", "| comparison | spec | PRR | IC025 | trial | n same-target pairs |", "|---|---|---|---|---|---|"]
    out += [f"| {c} | {s} | {f(v.get('prr_clean'))} | {f(v.get('ic025_clean2'))} | {f(v.get('trial'))} | {n if n is not None else ''} |" for (_, s, c, k, v, n, _) in rs]
    out.append("")
with open(os.path.join(PROJ, "results/RESULTS_SUMMARY.md"), "a") as fh: fh.write("\n".join(out))
print(f"{len(rows)} post hoc rows; RESULTS_SUMMARY.md regenerated with section 10")
