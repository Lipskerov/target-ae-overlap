#!/usr/bin/env python3
"""Post-hoc (not in SAP v1): (1) Delta by target family (ChEMBL protein classification; molecule label = family of its
targets, 'multi-family' if several); (2) modality x therapeutic area (trial-based) cross-check; (3) annotation-specificity
check: biologic vs small-molecule Delta restricted to molecules with exactly one cohort target.
Writes results/tables/effect_modification_target.json. Needs data/derived/target_family.csv (fetch_target_family.py),
results/tables/molecule_types.csv and indication_types.csv."""
import csv, json, os, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import effect_modification as em, effect_modification_adjusted as ea
T = os.path.join(em.PROJ, "results/tables")
fam = {r["target_symbol"]: r["family"] for r in csv.DictReader(open(os.path.join(em.PROJ, "data/derived/target_family.csv")))}
mtg = defaultdict(set)
for r in csv.DictReader(open(f"{em.CD}/D4_molecule_target.csv")): mtg[r["ingredient_rxcui"]].add(r["target_symbol"])
famlab = {m: (next(iter({fam[t] for t in mtg[m]})) if len({fam[t] for t in mtg[m]}) == 1 else "multi-family") for m in em.mols}
ntg = np.array([len(mtg[m]) for m in em.mols])
area = pd.read_csv(f"{T}/indication_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").area_trials
mod = em.mod
def adj(label, labels, lv):
    P = pd.read_csv(f"{T}/pairs_all_R2_{label}.csv.gz", dtype={"rx_i": str, "rx_j": str}); P["atc_missing"] = 1 - P.atc_both
    sel = (P.rx_i.map(labels) == lv) & (P.rx_j.map(labels) == lv); Q = P[sel.values].reset_index(drop=True)
    ms = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(ms)}; Q["ii"], Q["jj"] = Q.rx_i.map(ix), Q.rx_j.map(ix)
    return ea.fit(Q, ea.COLS, ms, 3)[0].get("same_target", {"beta": None, "ci95": None, "note": "same_target constant in stratum"})
out = dict(note="post hoc; not in SAP v1", family_counts=Counter(famlab.values()),
           single_target_share=dict(biologic=round(float((ntg[mod == "biologic"] == 1).mean()), 3), small_molecule=round(float((ntg[mod == "small_molecule"] == 1).mean()), 3)),
           modality_x_area=pd.crosstab(pd.Series([area[m] for m in em.mols], name="area"), pd.Series(mod, name="modality")).to_dict())
for label, prof in (("prr_clean", f"{em.CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{em.CD}/profile_clean2_ic025_k20.csv")):
    O = em.overlap(prof); rng = np.random.default_rng(np.random.SeedSequence(20261113).spawn(4)[3]); res = {"family": {}}
    arr = np.array([famlab[m] for m in em.mols])
    for f, k in Counter(arr).most_common():
        if k < 10: continue
        r, _ = em.stratum(O, arr == f, rng)
        if r["same_target_pairs"] >= 30: r["adjusted_same_target"] = adj(label, famlab, f)
        res["family"][f] = r; print(label, "family", f"{f:24s}", k, r["same_target_pairs"], r["delta"], r["ci95"], r.get("adjusted_same_target", {}).get("beta"), r.get("adjusted_same_target", {}).get("ci95"), flush=True)
    ar = np.array([area[m] for m in em.mols]); res["neoplasms_by_modality"] = {}
    for mo in ("small_molecule", "biologic"):
        r, _ = em.stratum(O, (ar == "neoplasms") & (mod == mo), rng); res["neoplasms_by_modality"][mo] = r
        print(label, "neoplasms", mo, r["n_molecules"], r["same_target_pairs"], r["delta"], r["ci95"], flush=True)
    res["single_target_only"] = {}; bs = {}
    for mo in ("small_molecule", "biologic"):
        r, b = em.stratum(O, (ntg == 1) & (mod == mo), rng); res["single_target_only"][mo] = r; bs[mo] = b
        print(label, "single-target", mo, r["n_molecules"], r["same_target_pairs"], r["delta"], r["ci95"], flush=True)
    d = bs["biologic"] - bs["small_molecule"]; lo, hi = np.nanpercentile(d, [2.5, 97.5])
    res["single_target_only"]["difference"] = dict(diff=round(res["single_target_only"]["biologic"]["delta"] - res["single_target_only"]["small_molecule"]["delta"], 3), ci95=[round(float(lo), 3), round(float(hi), 3)])
    print(label, "single-target bio-sm diff", res["single_target_only"]["difference"], flush=True)
    out[label] = res
json.dump(out, open(f"{T}/effect_modification_target.json", "w"), indent=1, default=dict)
