#!/usr/bin/env python3
"""Amendment 3, WP3 (post hoc): same target, opposite action. Seed index 12.

Direction of ChEMBL action_type (R3) per the amendment table; R3 normalized_name -> ingredient_rxcui via R2;
restricted to cohort molecule x target rows in faers_q/b/D4. A molecule with both directions at a target is excluded
at that target. A same-target pair is classified over its shared targets where both molecules have a defined
direction: same (same at every such target), opposite (opposite at every such target), mixed (otherwise).

Reported per ranking (prr_clean, ic025_clean2):
  - pair counts per class; mean shared signals same vs opposite; difference with molecule-bootstrap CI (B=2000)
  - adjusted: OLS on same-target same/opposite pairs, overlap ~ opposite + model_full covariates except same_target
    (ea.fit, molecule bootstrap B=500)
  - both again within GPCR + nuclear receptor shared targets: classification over the shared directional targets whose
    family (data/derived/target_family.csv) is GPCR or nuclear receptor; pairs with none such are out of this subset.
Reproduction: pairs table same_target sum = 7,923 and equals the D4-derived same-target set; the count of opposite-
direction molecule-pair x target combinations is reported next to the reviewer's 856 (not forced).
Writes results/tables/action_direction_<label>.json
"""
import csv, json, os, sys
from collections import defaultdict, Counter
from itertools import combinations
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", "..")); sys.path.insert(0, HERE)
import effect_modification_adjusted as ea
ea.B = 500
T = os.path.join(PROJ, "results/tables"); CD = os.path.join(PROJ, "data/derived/faers_q/b"); RAW = os.path.join(PROJ, "data/raw/trialytic_db")
B_DIFF = 2000
SEED = np.random.SeedSequence(20261113).spawn(15)[12]
REVIEWER_OPPOSITE = 856

ACT = {"AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR", "OPENER"}
INH = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "INVERSE AGONIST", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
       "ALLOSTERIC ANTAGONIST", "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DEGRADER"}
def direction(a):
    a = (a or "").strip().upper()
    return "act" if a in ACT else "inh" if a in INH else None

# ---- cohort molecule x target (D4)
d4 = defaultdict(set); tg = defaultdict(set)
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")):
    d4[r["ingredient_rxcui"]].add(r["target_symbol"]); tg[r["target_symbol"]].add(r["ingredient_rxcui"])
cohort = set(d4)
fam = {r["target_symbol"]: r["family"] for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/target_family.csv")))}
GN = {"GPCR", "nuclear receptor"}

# ---- R3 -> rxcui via R2
r2 = defaultdict(set)
for r in csv.DictReader(open(f"{RAW}/R2_name_to_ingredient_map.csv")):
    if r["ingredient_rxcui"]: r2[r["normalized_name"]].add(r["ingredient_rxcui"])
dirs = defaultdict(set)          # (rxcui, target) -> set of directions (cohort rows)
dirs_all = defaultdict(set)      # unrestricted (any mapped molecule, any target) for documentation
r3_rows = r3_joined = r3_cohort = 0
for r in csv.DictReader(open(f"{RAW}/R3_chembl_target_drugs.csv")):
    r3_rows += 1; rx = r2.get(r["normalized_name"], set())
    if rx: r3_joined += 1
    d = direction(r["action_type"])
    for m in rx:
        if d: dirs_all[(m, r["target_symbol"])].add(d)
        if m in cohort and r["target_symbol"] in d4[m]:
            r3_cohort += 1
            if d: dirs[(m, r["target_symbol"])].add(d)
cohort_mt = sum(len(v) for v in d4.values())
mol_dir = {k: next(iter(v)) for k, v in dirs.items() if len(v) == 1}
conflict = sum(1 for v in dirs.values() if len(v) > 1)

def opp_combos(dd, targets_of):
    """count (molecule pair, target) combos with defined direction for both, by same/opposite."""
    c = Counter()
    for t, ms in targets_of.items():
        ms = sorted(m for m in ms if (m, t) in dd)
        for a, b in combinations(ms, 2): c["opposite" if dd[(a, t)] != dd[(b, t)] else "same"] += 1
    return c
combos = opp_combos(mol_dir, tg)
md_all = {k: next(iter(v)) for k, v in dirs_all.items() if len(v) == 1}
tg_all = defaultdict(set)
for (m, t) in md_all: tg_all[t].add(m)
combos_all = opp_combos(md_all, tg_all)

def classify(i, j, fam_filter=None):
    shared = d4[i] & d4[j]
    ds = [mol_dir[(i, t)] == mol_dir[(j, t)] for t in shared
          if (i, t) in mol_dir and (j, t) in mol_dir and (fam_filter is None or fam.get(t) in fam_filter)]
    if not ds: return None
    return "same" if all(ds) else "opposite" if not any(ds) else "mixed"

COV = [c for c in ea.COLS if c != "same_target"]

def boot_diff(Q, seed):
    mols = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(mols)}
    ii, jj = Q.rx_i.map(ix).values, Q.rx_j.map(ix).values; y = Q.overlap.values.astype(float); op = Q.opposite.values == 1
    rng = np.random.default_rng(seed); nm = len(mols); bs = np.empty(B_DIFF)
    for b in range(B_DIFF):
        w = rng.multinomial(nm, np.full(nm, 1 / nm)).astype(float); pw = w[ii] * w[jj]
        s1, s0 = pw[op].sum(), pw[~op].sum()
        bs[b] = ((pw[op] * y[op]).sum() / s1 - (pw[~op] * y[~op]).sum() / s0) if s1 > 0 and s0 > 0 else np.nan
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return dict(molecules=nm, B=B_DIFF, n_boot_undefined=int(np.isnan(bs).sum()), ci95=[round(float(lo), 4), round(float(hi), 4)])

def analyse(P, cls_col, seeds):
    out = dict(pair_counts={k: int(v) for k, v in P[cls_col].value_counts().items()})
    for k in ("same", "opposite", "mixed"):
        s = P[cls_col] == k; out[f"mean_overlap_{k}"] = round(float(P.overlap[s].mean()), 4) if s.any() else None
    Q = P[P[cls_col].isin(["same", "opposite"])].reset_index(drop=True); Q["opposite"] = (Q[cls_col] == "opposite").astype(int)
    if Q.opposite.sum() < 10 or (1 - Q.opposite).sum() < 10:
        out["note"] = "fewer than 10 pairs in a class: difference not measured"; return out
    out["difference_opposite_minus_same"] = dict(diff=round(float(Q.overlap[Q.opposite == 1].mean() - Q.overlap[Q.opposite == 0].mean()), 4),
                                                 **boot_diff(Q, seeds[0]))
    mols = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(mols)}; Q["ii"], Q["jj"] = Q.rx_i.map(ix), Q.rx_j.map(ix)
    res, _, names = ea.fit(Q, ["opposite"] + COV, mols, seeds[1])
    out["adjusted_ols"] = dict(n_pairs=int(len(Q)), molecules=len(mols), B=ea.B, dropped_constant=[c for c in ["opposite"] + COV if c not in names], coefficients=res)
    return out

if __name__ == "__main__":
    print(f"R3 rows {r3_rows}; joined to rxcui via R2 {r3_joined}; rows on cohort molecule x D4 target {r3_cohort}")
    print(f"cohort molecule x target combos (D4) {cohort_mt}; with any direction {len(dirs)}; single direction {len(mol_dir)}; conflicting (excluded) {conflict}")
    print(f"REPRO opposite-direction molecule-pair x target combos (cohort, D4) = {combos['opposite']} (same {combos['same']}) vs reviewer 856 -> diff {combos['opposite'] - REVIEWER_OPPOSITE}")
    print(f"     unrestricted (all R2-mapped R3 molecules/targets): opposite {combos_all['opposite']} same {combos_all['same']}")
    kids = SEED.spawn(8)
    for li, label in enumerate(("prr_clean", "ic025_clean2")):
        P = pd.read_csv(f"{T}/pairs_all_R2_{label}.csv.gz", dtype={"rx_i": str, "rx_j": str}); P["atc_missing"] = 1 - P.atc_both
        n_st = int(P.same_target.sum())
        st_d4 = np.array([bool(d4[a] & d4[b]) for a, b in zip(P.rx_i, P.rx_j)])
        repro_ok = n_st == 7923 and bool((st_d4 == (P.same_target.values == 1)).all())
        print(f"{label} REPRO same-target pairs in table = {n_st} (expected 7923); D4-derived set identical: {bool((st_d4 == (P.same_target.values == 1)).all())} -> {'PASS' if repro_ok else 'FAIL'}")
        if not repro_ok:
            json.dump(dict(label=label, repro_check="FAIL", same_target_pairs=n_st), open(f"{T}/action_direction_{label}.json", "w"), indent=1); sys.exit(1)
        S = P[P.same_target == 1].reset_index(drop=True)
        S["cls_all"] = [classify(a, b) for a, b in zip(S.rx_i, S.rx_j)]
        S["cls_gn"] = [classify(a, b, GN) for a, b in zip(S.rx_i, S.rx_j)]
        out = dict(label=label, amendment="AMENDMENT_3 WP3", seed="SeedSequence(20261113).spawn(15)[12]", post_hoc=True,
                   repro=dict(same_target_pairs=n_st, expected=7923, d4_set_identical=True,
                              opposite_combos_cohort=combos["opposite"], same_combos_cohort=combos["same"], reviewer_opposite_combos=REVIEWER_OPPOSITE,
                              difference_vs_reviewer=combos["opposite"] - REVIEWER_OPPOSITE,
                              opposite_combos_unrestricted=combos_all["opposite"], same_combos_unrestricted=combos_all["same"],
                              note="difference documented, not forced; reviewer mapping not available here"),
                   join=dict(r3_rows=r3_rows, r3_rows_joined_rxcui=r3_joined, r3_rows_on_cohort_d4=r3_cohort, cohort_molecule_target=cohort_mt,
                             molecule_target_with_direction=len(dirs), single_direction=len(mol_dir), conflicting_excluded=conflict),
                   same_target_pairs_without_defined_direction=int(S.cls_all.isna().sum()),
                   all_targets=analyse(S[S.cls_all.notna()].reset_index(drop=True), "cls_all", kids[4 * li:4 * li + 2]),
                   gpcr_nuclear_receptor=analyse(S[S.cls_gn.notna()].reset_index(drop=True), "cls_gn", kids[4 * li + 2:4 * li + 4]))
        for k in ("all_targets", "gpcr_nuclear_receptor"):
            r = out[k]; d = r.get("difference_opposite_minus_same", {}); a = r.get("adjusted_ols", {}).get("coefficients", {}).get("opposite", {})
            print(f"{label} {k}: counts {r['pair_counts']} | mean same {r['mean_overlap_same']} opposite {r['mean_overlap_opposite']} mixed {r['mean_overlap_mixed']}"
                  f" | diff {d.get('diff')} CI {d.get('ci95')} | adjusted {a.get('beta')} CI {a.get('ci95')}", flush=True)
        json.dump(out, open(f"{T}/action_direction_{label}.json", "w"), indent=1)
    print("done")
