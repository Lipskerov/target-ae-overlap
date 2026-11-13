# AMENDMENT_4 A4.1b + A4.2 sensitivity (copy of modifier_joint_model.py; kinase_pair = kinase-family AND small molecule; era_missing dropped)
#!/usr/bin/env python3
"""Amendment 3, WP1 — joint modifier model (seed index 10). Post hoc relative to SAP v1.

OLS over all pairs (results/tables/pairs_all_R2_<label>.csv.gz), outcome = shared signals (overlap):
  model_full covariates (pairs_all.py) + bio_pair, kinase_pair, onc_pair, new_pair + area_missing, era_missing
  + same_target x each of the 4 modifiers.
Molecule bootstrap B=500 (effect_modification_adjusted.fit). Leave-one-modifier-out refits drop the modifier, its
same_target interaction and (area/era) its missing indicator. Same-target pair counts per modifier cell.
Reproduction check first: model_full without added terms, same_target = 0.143 (PRR) / 0.299 (IC025) to 3 decimals.

Labels (amendment 3): modality = results/tables/molecule_types.csv; family = molecule rule of
effect_modification_target.py (single family across its D4 targets, else 'multi-family'); area =
indication_types.csv area_trials (unknown/mixed = missing); approval year = R10 earliest_orig_ap, else R11
first_approval_year (none = missing).
Writes results/tables/modifier_joint_<label>.json.
"""
import csv, json, os, sys
from collections import defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import effect_modification_adjusted as ea   # main is guarded; import only reads molecule_types.csv

PROJ = ea.PROJ; T = ea.T
CD = os.path.join(PROJ, "data/derived/faers_q/b"); RAW = os.path.join(PROJ, "data/raw/trialytic_db")
SEED = np.random.SeedSequence(20261113).spawn(15)[10]
REPRO = {"prr_clean": 0.143, "ic025_clean2": 0.299}
MODS = ["bio_pair", "kinase_pair", "onc_pair", "new_pair"]
MISS = {"onc_pair": "area_missing"}  # A4.1b: era_missing dropped (aliased with year_missing)

# ---- molecule labels
mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]
mset = set(mols)
modality = ea.mt
fam = {r["target_symbol"]: r["family"] for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/target_family.csv")))}
mtg = defaultdict(set)
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): mtg[r["ingredient_rxcui"]].add(r["target_symbol"])
famlab = {m: (next(iter({fam[t] for t in mtg[m]})) if len({fam[t] for t in mtg[m]}) == 1 else "multi-family") for m in mols}
area = pd.read_csv(f"{T}/indication_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").area_trials
yr = {}
for r in csv.DictReader(open(f"{RAW}/R11_fda_class_all.csv")):
    if r["ingredient_rxcui"] in mset and r["first_approval_year"]: yr[r["ingredient_rxcui"]] = int(float(r["first_approval_year"]))
for r in csv.DictReader(open(f"{RAW}/R10_drug_approval_all.csv")):      # R10 takes precedence
    if r["ingredient_rxcui"] in mset and r["earliest_orig_ap"]: yr[r["ingredient_rxcui"]] = int(r["earliest_orig_ap"][:4])
# cross-check against the era column stored on 2026-11-13
era_stored = pd.read_csv(f"{T}/molecule_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").era
era_mine = {m: ("unknown" if m not in yr else ("pre2000" if yr[m] < 2000 else "2000plus")) for m in mols}
era_mismatch = int(sum(era_mine[m] != era_stored[m] for m in mols))

lab = pd.DataFrame(dict(bio=[int(modality[m] == "biologic") for m in mols],
                        kinase=[int(famlab[m] == "kinase" and modality[m] == "small_molecule") for m in mols],  # A4.2: kinase_sm
                        onc=[int(area.get(m, "unknown") == "neoplasms") for m in mols],
                        area_unk=[int(area.get(m, "unknown") in ("unknown", "mixed")) for m in mols],
                        new=[int(m in yr and yr[m] >= 2000) for m in mols],
                        era_unk=[int(m not in yr) for m in mols]), index=mols)
mol_counts = {c: int(lab[c].sum()) for c in lab.columns}

def ols_point(P, xcols):
    X = np.column_stack([np.ones(len(P))] + [P[c].values.astype(float) for c in xcols])
    return dict(zip(["intercept"] + xcols, np.linalg.lstsq(X, P.overlap.values.astype(float), rcond=None)[0]))

def run(label, seeds):
    P = pd.read_csv(f"{T}/pairs_all_R2_{label}.csv.gz", dtype={"rx_i": str, "rx_j": str}); P["atc_missing"] = 1 - P.atc_both
    assert P.rx_i.isin(mset).all() and P.rx_j.isin(mset).all()
    ms = sorted(set(P.rx_i) | set(P.rx_j)); ix = {m: k for k, m in enumerate(ms)}; P["ii"], P["jj"] = P.rx_i.map(ix), P.rx_j.map(ix)
    out = dict(label=label, amendment="AMENDMENT_3 WP1", seed="SeedSequence(20261113).spawn(15)[10]", B=ea.B,
               n_pairs=int(len(P)), n_molecules=len(ms), molecules_with_label=mol_counts, era_label_mismatch_vs_molecule_types=era_mismatch)
    # reproduction check (point estimate of model_full same_target)
    st = ols_point(P, ea.COLS)["same_target"]; ok = round(st, 3) == REPRO[label]
    out["repro_check"] = dict(expected=REPRO[label], got=round(float(st), 4), status="PASS" if ok else "FAIL")
    print(label, "REPRO model_full same_target", round(st, 4), "expected", REPRO[label], out["repro_check"]["status"], flush=True)
    if not ok: return out
    L = {m: lab.loc[P.rx_i.values].reset_index(drop=True) if m == "i" else lab.loc[P.rx_j.values].reset_index(drop=True) for m in ("i", "j")}
    for mod, c in (("bio_pair", "bio"), ("kinase_pair", "kinase"), ("onc_pair", "onc"), ("new_pair", "new")):
        P[mod] = (L["i"][c].values & L["j"][c].values).astype(int)
        P[f"same_target_x_{mod}"] = P.same_target * P[mod]
    P["area_missing"] = (L["i"].area_unk.values | L["j"].area_unk.values).astype(int)
    P["era_missing"] = (L["i"].era_unk.values | L["j"].era_unk.values).astype(int)
    # same-target pair counts per modifier cell
    cells = {}
    for mod in MODS:
        miss = P[MISS[mod]].values == 1 if mod in MISS else np.zeros(len(P), bool)
        d = {}
        for cell, msk in (("modifier_1", (P[mod].values == 1)), ("modifier_0_known", (P[mod].values == 0) & ~miss), ("missing", miss)):
            d[cell] = dict(pairs=int(msk.sum()), same_target_pairs=int((msk & (P.same_target.values == 1)).sum()))
        cells[mod] = d
    out["cell_counts"] = cells; print(label, "cells", json.dumps(cells), flush=True)
    joint = ea.COLS + MODS + ["area_missing"] + [f"same_target_x_{m}" for m in MODS]
    out["joint_model"] = ea.fit(P, joint, ms, seeds[0])[0]
    for k in joint:
        if k.startswith("same_target"): print(label, "joint", k, out["joint_model"][k], flush=True)
    out["leave_one_modifier_out"] = {}
    for n_, mod in enumerate(MODS):
        drop = {mod, f"same_target_x_{mod}"} | ({MISS[mod]} if mod in MISS else set())
        cols = [c for c in joint if c not in drop]
        r = ea.fit(P, cols, ms, seeds[1 + n_])[0]; out["leave_one_modifier_out"][f"without_{mod}"] = dict(dropped=sorted(drop), coefficients=r)
        print(label, "LOMO without", mod, {k: (v["beta"], v["ci95"]) for k, v in r.items() if k.startswith("same_target")}, flush=True)
    return out

if __name__ == "__main__":
    print("molecules with label:", mol_counts, "| era mismatch vs molecule_types.csv:", era_mismatch, flush=True)
    kids = SEED.spawn(2)   # one child per ranking; 5 grandchildren each (joint + 4 leave-one-out)
    for li, label in enumerate(("prr_clean", "ic025_clean2")):
        seeds = kids[li].spawn(5)
        o = run(label, seeds)
        json.dump(o, open(f"{T}/modifier_joint_a4_{label}.json", "w"), indent=1)
        print(label, "written", flush=True)
