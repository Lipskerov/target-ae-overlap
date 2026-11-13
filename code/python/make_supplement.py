#!/usr/bin/env python3
"""Supplementary data tables for the revised analysis. Copies and joins saved files only."""
import json, shutil, os
import pandas as pd
P = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")); os.chdir(P)
S = "private/supplement/tables"; b = "data/derived/faers_q/b"; kw = dict(dtype=str, keep_default_na=False); T = "results/tables"
tot = pd.read_csv(f"{b}/molecule_totals.csv", **kw); coh = pd.read_csv(f"{b}/D2_cohort_molecules.csv", **kw); mt = pd.read_csv(f"{b}/D4_molecule_target.csv", **kw)
sig = pd.read_csv(f"{b}/signals.csv.gz", dtype={"ingredient_rxcui": str}, keep_default_na=False)
for c in ("clean", "clean2"): sig[c] = sig[c].astype(str).str.lower().isin(["true", "1"])
nm = dict(zip(coh.ingredient_rxcui, coh.name))
s1 = (coh.merge(tot[["ingredient_rxcui", "reports"]], on="ingredient_rxcui", how="left")
      .merge(mt.groupby("ingredient_rxcui").target_symbol.apply(lambda x: "|".join(sorted(x))).rename("targets"), on="ingredient_rxcui")
      .merge(sig[sig.clean].groupby("ingredient_rxcui").size().rename("signals_standard_cleaning"), on="ingredient_rxcui", how="left")
      .merge(sig[sig.clean2].groupby("ingredient_rxcui").size().rename("signals_extended_cleaning"), on="ingredient_rxcui", how="left"))
s1.to_csv(f"{S}/S1_cohort_molecules.csv", index=False)
s2 = mt.assign(name=mt.ingredient_rxcui.map(nm)).groupby("target_symbol").agg(n_molecules=("ingredient_rxcui", "nunique"), molecules=("name", lambda x: "|".join(sorted(set(x))))).reset_index()
s2[s2.n_molecules >= 3].to_csv(f"{S}/S2_cohort_targets.csv", index=False)
for k, f in (("S3a_profiles_top20_prr_standard_cleaning", "profile_clean_prr_k20.csv"), ("S3b_profiles_top20_ic025_extended_cleaning", "profile_clean2_ic025_k20.csv")):
    p = pd.read_csv(f"{b}/{f}", dtype={"ingredient_rxcui": str}, keep_default_na=False); p.insert(1, "name", p.ingredient_rxcui.map(nm)); p.to_csv(f"{S}/{k}.csv", index=False)
shutil.copy(f"{b}/signals.csv.gz", f"{S}/S4_all_signals_extraction_b.csv.gz")
rows = []
import glob
for fn in sorted(glob.glob(f"{T}/primary_R2_*.json")):
    r = json.load(open(fn))
    for m in ("shared_signals", "jaccard"):
        if m in r: x = r[m]; rows.append(dict(analysis=os.path.basename(fn)[11:-5], measure=m, molecules=r["n_molecules"], same_target_pairs=r["same_target_pairs"], different_target_pairs=r["different_target_pairs"], pair_mean_same=x["pair_mean_same"], pair_mean_different=x["pair_mean_different"], delta=x["delta"], ci_lo=x["ci95"][0], ci_hi=x["ci95"][1], p_perm=x["p_perm_one_sided"], permutations=r["B_perm"], bootstrap=r["B_boot"]))
pd.DataFrame(rows).to_csv(f"{S}/S5_all_primary_and_sensitivity_runs.csv", index=False)
for v in ("prr_clean", "ic025_clean2"):
    shutil.copy(f"{T}/pairs_all_R2_{v}.csv.gz", f"{S}/S6_all_pairs_{v}.csv.gz"); shutil.copy(f"{T}/comention_pairs_{v}.csv.gz", f"{S}/S7_comention_pairs_{v}.csv.gz"); shutil.copy(f"{T}/asof_sets_{v}.csv", f"{S}/S8_asof_matched_sets_{v}.csv")
    json.dump(json.load(open(f"{T}/target_vs_class_R2_{v}.json")), open(f"{S}/S9_target_vs_class_{v}.json", "w"), indent=1)
c = "data/derived/ctgov/profile"; shutil.copy(f"{c}/molecule_totals.csv", f"{S}/S10_trial_molecule_totals.csv"); shutil.copy(f"{c}/signals.csv.gz", f"{S}/S11_trial_signals.csv.gz")
shutil.copy("data/derived/ctgov/arm_attribution.csv.gz", f"{S}/S12_trial_arm_attribution_v2.csv.gz")
pd.concat([pd.read_csv(f"private/reviews/P3_arm_audit_{k}_verdicts.csv", **kw).assign(audit=l) for k, l in (("tier1", "rule v1 tier 1"), ("tier2", "rule v1 tier 2"), ("v2", "rule v2 (final)"))]).to_csv(f"{S}/S13_arm_audit_verdicts.csv", index=False)
shutil.copy("data/raw/faers_quarterly/MANIFEST.tsv", f"{S}/S14_faers_quarterly_files_manifest.tsv")
json.dump(dict(report_flow=json.load(open("data/derived/faers_q/report_flow.json")), case_characteristics=json.load(open("data/derived/faers_q/case_characteristics.json"))), open(f"{S}/S15_report_flow_and_case_characteristics.json", "w"), indent=1)
shutil.copy("data/derived/faers_q/name_map.csv.gz", f"{S}/S16_drug_string_to_ingredient_map.csv.gz")
open(f"{S}/README.md", "w").write("""# Supplementary data tables (revised analysis)

Read CSV files with missing-value detection off (pandas `keep_default_na=False`): the target symbol `NA` is neuraminidase.

| File | Content |
|---|---|
| S1_cohort_molecules.csv | Cohort molecules (extraction b): RxNorm ingredient, name, suspect reports, human protein targets, signal counts under both cleaning rules |
| S2_cohort_targets.csv | Cohort targets (three or more member molecules) with members |
| S3a_profiles_top20_prr_standard_cleaning.csv | Top-20 profile per molecule ranked by PRR, with report count |
| S3b_profiles_top20_ic025_extended_cleaning.csv | Top-20 profile per molecule ranked by IC025 after extended cleaning, with report count |
| S4_all_signals_extraction_b.csv.gz | Every signal of every cohort molecule: a, PRR, chi-square, ROR lower limit, IC025, cleaning flags, ranks |
| S5_all_primary_and_sensitivity_runs.csv | Every primary and sensitivity run: molecules, pairs, means, difference, interval, p, resample counts |
| S6_all_pairs_*.csv.gz | One row per cohort pair with overlap, same-target flag and all covariates |
| S7_comention_pairs_*.csv.gz | Pair-level co-mention counts and overlap before and after exclusion |
| S8_asof_matched_sets_*.csv | As-of matched sets with first molecule, successor, date, overlaps and the three controls |
| S9_target_vs_class_*.json | Stratified contrasts, adjusted models, target-level analysis, top shared terms |
| S10_trial_molecule_totals.csv, S11_trial_signals.csv.gz | Trial analysis: trials, subjects at risk, trial-signals per molecule |
| S12_trial_arm_attribution_v2.csv.gz | Label of every results arm under the final attribution rule |
| S13_arm_audit_verdicts.csv | The three 100-arm audits with verdicts and notes |
| S14_faers_quarterly_files_manifest.tsv | The 91 FDA quarterly files: URL, bytes, sha256, row counts |
| S15_report_flow_and_case_characteristics.json | Report flow and case characteristics |
| S16_drug_string_to_ingredient_map.csv.gz | Every drug string with ≥50 rows (and every mapped string): ingredient(s), mapping source, row count |
""")
print("S1", len(s1), "| S2", int((s2.n_molecules >= 3).sum()))
