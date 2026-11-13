#!/usr/bin/env python3
"""Amendment 3, WP2 (post hoc; seed index 11): do the FAERS effect modifiers appear in trial-reported adverse events?

Data: data/derived/ctgov/profile/ (profile_clean_k20, D2, D4, shared_condition_pairs).
Molecule labels (not edited): modality from results/tables/molecule_types.csv; target family by the molecule rule of
effect_modification_target.py (family of the molecule's targets, 'multi-family' if several) applied to the trial-cohort
target map D4 + data/derived/target_family.csv; area from results/tables/indication_types.csv column area_trials.
Trial-cohort molecules absent from a label file get no modality label (excluded from both modality strata, bio_pair = 0)
and area 'other' (amendment: unknown/mixed are in 'other'); they are counted in the output.

(A) Stratum Delta, primary_inference definition (molecule-averaged same- minus different-target mean shared signals),
    both molecules in the stratum; molecule bootstrap within the stratum, B = 2000; differences = difference of
    independent bootstrap draws. Pair sets: all pairs; excluding shared-condition pairs (removed from both groups, as
    primary_inference EXCLUDE_PAIRS).
(B) Trial-pair OLS over all pairs of trial-cohort molecules, outcome = shared signals; covariates same_target, same
    ATC-4/3/1, atc_missing (R9), same_epc, epc_missing (R11), shared_condition, bio_pair, kinase_pair, onc_pair and
    same_target x each. Molecule bootstrap (pair weight w_i*w_j), B = 500, via effect_modification_adjusted.fit.
    No report-volume or co-mention covariate exists for trials.
Repro check: whole trial cohort Delta = 1.2142 (results/tables/primary_R2_ctgov_all.json).
Writes results/tables/modifier_trials.json.  Run: python3 code/python/modifier_trials.py
"""
import csv, json, os, sys, time
from collections import Counter
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
os.environ.setdefault("B", "500")
import effect_modification_adjusted as ea          # fit(): OLS + molecule bootstrap (main guarded; reads molecule_types.csv only)
assert ea.B == 500, ea.B
CT = os.path.join(PROJ, "data/derived/ctgov/profile"); RAW = os.path.join(PROJ, "data/raw/trialytic_db"); T = os.path.join(PROJ, "results/tables")
B_STRAT = 2000; t0 = time.time()
seq = np.random.SeedSequence(20261113).spawn(15)[11]
seq_strata, seq_reg = seq.spawn(2)
def step(m): print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CT}/D2_cohort_molecules.csv"))]
idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
names = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{CT}/D2_cohort_molecules.csv"))}
tg = {}; mtg = {m: set() for m in mols}
for r in csv.DictReader(open(f"{CT}/D4_molecule_target.csv")):
    tg.setdefault(r["target_symbol"], set()).add(idx[r["ingredient_rxcui"]]); mtg[r["ingredient_rxcui"]].add(r["target_symbol"])
S0 = np.zeros((n, n), np.float32)
for g in tg.values():
    g = np.array(sorted(g)); S0[np.ix_(g, g)] = 1
np.fill_diagonal(S0, 0); D0 = 1 - S0; np.fill_diagonal(D0, 0)

rt = {}
for r in csv.DictReader(open(f"{CT}/profile_clean_k20.csv")):
    if r["ingredient_rxcui"] in idx: rt.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"].upper())
terms = sorted({t for s in rt.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
A = np.zeros((n, len(terms)), np.float32)
for i, s in rt.items(): A[i, [tix[t] for t in s]] = 1
O = A @ A.T; np.fill_diagonal(O, 0)

# shared-condition pairs
SC = np.zeros((n, n), np.float32)
for r in csv.DictReader(open(f"{CT}/shared_condition_pairs.csv")):
    if r["rx_i"] in idx and r["rx_j"] in idx: SC[idx[r["rx_i"]], idx[r["rx_j"]]] = SC[idx[r["rx_j"]], idx[r["rx_i"]]] = 1
np.fill_diagonal(SC, 0)

# ---- labels
mt = pd.read_csv(f"{T}/molecule_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").modality
ar = pd.read_csv(f"{T}/indication_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").area_trials
fam = {r["target_symbol"]: r["family"] for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/target_family.csv")))}
mod = np.array([mt.get(m, "unlabeled") for m in mols], object)
famlab = np.array([(next(iter({fam[t] for t in mtg[m]})) if len({fam[t] for t in mtg[m]}) == 1 else "multi-family") for m in mols], object)
fam3 = np.array([f if f in ("kinase", "GPCR") else "other" for f in famlab], object)
area_raw = np.array([ar.get(m, "unlabeled") for m in mols], object)
area2 = np.array(["neoplasms" if a == "neoplasms" else "other" for a in area_raw], object)
unl_mod = [names[m] for m in mols if m not in mt.index]; unl_area = [names[m] for m in mols if m not in ar.index]
step(f"n={n} modality={dict(Counter(mod))} family3={dict(Counter(fam3))} area2={dict(Counter(area2))}")
step(f"unlabeled modality: {unl_mod} | unlabeled area: {unl_area}")

# ---- primary Delta (primary_inference.delta, extended with a focal mask as in effect_modification.delta_w)
def delta_w(Smat, Dmat, w, focal):
    ns, nd = Smat @ w, Dmat @ w
    ok = focal & (w > 0) & (ns > 0) & (nd > 0)
    d = (Smat * O) @ w / np.where(ns > 0, ns, 1) - (Dmat * O) @ w / np.where(nd > 0, nd, 1)
    return (float((w[ok] * d[ok]).sum() / w[ok].sum()) if ok.any() else np.nan), int(ok.sum())

rep, _ = delta_w(S0, D0, np.ones(n, np.float32), np.ones(n, bool))
ref = json.load(open(f"{T}/primary_R2_ctgov_all.json"))["shared_signals"]["delta"]
repro_ok = abs(round(rep, 4) - ref) < 0.00005
print(f"REPRO whole trial cohort Delta = {rep:.4f} (reference {ref}) -> {'PASS' if repro_ok else 'FAIL'}", flush=True)
out = dict(wp="WP2", amendment="protocol/AMENDMENT_3.md", seed="SeedSequence(20261113).spawn(15)[11]; children [0]=strata, [1]=regression",
           B_stratum=B_STRAT, B_regression=ea.B, n_molecules=n,
           repro_check=dict(delta=round(rep, 4), reference=ref, status="PASS" if repro_ok else "FAIL"))
if not repro_ok:
    json.dump(out, open(f"{T}/modifier_trials.json", "w"), indent=1); sys.exit("repro_check FAIL")

out["labels"] = dict(modality=dict(Counter(mod)), family_molecule_rule=dict(Counter(famlab)), family3=dict(Counter(fam3)),
                     area_trials=dict(Counter(area_raw)), area2=dict(Counter(area2)),
                     unlabeled_modality=unl_mod, unlabeled_area=unl_area,
                     note="trial-cohort molecules absent from molecule_types.csv have no modality (excluded from modality strata, bio_pair=0); "
                          "absent from indication_types.csv -> area 'other'. Peptides are in neither modality stratum and have bio_pair=0.")
iu = np.triu_indices(n, 1)

def stratum(Sm, Dm_, mask, rng):
    P = np.outer(mask, mask).astype(np.float32); Sg, Dg = Sm * P, Dm_ * P
    est, used = delta_w(Sg, Dg, np.ones(n, np.float32), mask)
    gi = np.where(mask)[0]; boots = np.empty(B_STRAT)
    for b in range(B_STRAT):
        w = np.zeros(n, np.float32); np.add.at(w, rng.choice(gi, len(gi)), 1)
        boots[b], _ = delta_w(Sg, Dg, w, mask)
    s_u, d_u = Sg[iu] > 0, Dg[iu] > 0
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    r = dict(n_molecules=int(mask.sum()), molecules_used=used, same_target_pairs=int(s_u.sum()), different_target_pairs=int(d_u.sum()),
             delta=round(est, 4), ci95=[round(float(lo), 4), round(float(hi), 4)], boot_nan=int(np.isnan(boots).sum()),
             pair_mean_same=round(float(O[iu][s_u].mean()), 3) if s_u.any() else None,
             pair_mean_different=round(float(O[iu][d_u].mean()), 3) if d_u.any() else None)
    if r["same_target_pairs"] < 10: r["note"] = "fewer than 10 same-target pairs - not interpreted"
    return r, boots

rng = np.random.default_rng(seq_strata)
FACT = (("modality", mod, ["biologic", "small_molecule"], ("biologic", "small_molecule")),
        ("family", fam3, ["kinase", "GPCR", "other"], ("kinase", "GPCR")),
        ("area", area2, ["neoplasms", "other"], ("neoplasms", "other")))
out["A_stratum_delta"] = {}
for pset in ("all_pairs", "excluding_shared_condition_pairs"):
    if pset == "all_pairs": Sm, Dm_ = S0, D0
    else: Sm, Dm_ = S0 * (1 - SC), D0 * (1 - SC)
    whole, _ = delta_w(Sm, Dm_, np.ones(n, np.float32), np.ones(n, bool))
    res = dict(whole_cohort_delta=round(whole, 4), same_target_pairs=int((Sm[iu] > 0).sum()), different_target_pairs=int((Dm_[iu] > 0).sum()))
    step(f"{pset}: whole-cohort Delta {whole:.4f}, same-target pairs {res['same_target_pairs']}")
    for fac, arr, levels, (a, b) in FACT:
        res[fac] = {}; bs = {}
        for lv in levels:
            r, bb = stratum(Sm, Dm_, arr == lv, rng); res[fac][lv] = r; bs[lv] = bb
            step(f"{pset} {fac:8s} {lv:15s} mols={r['n_molecules']:3d} same-target pairs={r['same_target_pairs']:4d} Delta={r['delta']} CI={r['ci95']}")
        d = bs[a] - bs[b]; lo, hi = np.nanpercentile(d, [2.5, 97.5])
        res[fac][f"difference_{a}_minus_{b}"] = dict(diff=round(res[fac][a]["delta"] - res[fac][b]["delta"], 4), ci95=[round(float(lo), 4), round(float(hi), 4)])
        step(f"{pset} {fac} {a} - {b} = {res[fac][f'difference_{a}_minus_{b}']}")
    out["A_stratum_delta"][pset] = res

# ---- (B) adjusted trial-pair OLS
atc = {}
for r in csv.DictReader(open(f"{RAW}/R9_rxnorm_atc.csv")):
    if r["rxcui"] in idx and len(r["code"]) >= 5: atc.setdefault(idx[r["rxcui"]], set()).add(r["code"])
def share_mat(level):
    M = np.zeros((n, n), np.float32); groups = {}
    for i, cs in atc.items():
        for c in {c[:level] for c in cs}: groups.setdefault(c, []).append(i)
    for g in groups.values():
        g = np.array(g); M[np.ix_(g, g)] = 1
    np.fill_diagonal(M, 0); return M
has_atc = np.array([i in atc for i in range(n)]); A4, A3, A1 = share_mat(5), share_mat(4), share_mat(1)
epc = {}
for r in csv.DictReader(open(f"{RAW}/R11_fda_class_all.csv")):
    if r["ingredient_rxcui"] in idx and r["pharm_class_epc"]: epc[idx[r["ingredient_rxcui"]]] = r["pharm_class_epc"]
ev = np.array([epc.get(i, "") for i in range(n)]); has_epc = ev != ""
EPC_same = (ev[:, None] == ev[None, :]) & np.outer(has_epc, has_epc)
bio = mod == "biologic"; kin = famlab == "kinase"; onc = area2 == "neoplasms"
ii, jj = iu
P = pd.DataFrame(dict(rx_i=np.array(mols)[ii], rx_j=np.array(mols)[jj], ii=ii, jj=jj, overlap=O[iu].astype(int),
                      same_target=S0[iu].astype(int), same_atc4=A4[iu].astype(int), same_atc3=A3[iu].astype(int), same_atc1=A1[iu].astype(int),
                      atc_missing=(~np.outer(has_atc, has_atc))[iu].astype(int), same_epc=EPC_same[iu].astype(int),
                      epc_missing=(~np.outer(has_epc, has_epc))[iu].astype(int), shared_condition=SC[iu].astype(int),
                      bio_pair=(bio[ii] & bio[jj]).astype(int), kinase_pair=(kin[ii] & kin[jj]).astype(int), onc_pair=(onc[ii] & onc[jj]).astype(int)))
for mname in ("bio_pair", "kinase_pair", "onc_pair"): P[f"same_target_x_{mname}"] = P.same_target * P[mname]
BASE = ["same_target", "same_atc4", "same_atc3", "same_atc1", "atc_missing", "same_epc", "epc_missing", "shared_condition"]
MODS = ["bio_pair", "kinase_pair", "onc_pair"]; XCOLS = BASE + MODS + [f"same_target_x_{m}" for m in MODS]
cells = {m: {f"same_target={s},{m}={v}": int(((P.same_target == s) & (P[m] == v)).sum()) for s in (0, 1) for v in (0, 1)} for m in MODS}
step(f"pairs={len(P)} same-target={int(P.same_target.sum())} molecules_with_atc={int(has_atc.sum())} molecules_with_epc={int(has_epc.sum())} cells={cells}")
rs = seq_reg.spawn(2)
reg = dict(pairs=int(len(P)), same_target_pairs=int(P.same_target.sum()), molecules_with_atc=int(has_atc.sum()), molecules_with_epc=int(has_epc.sum()),
           shared_condition_pairs=int(P.shared_condition.sum()), cells=cells,
           note="No report-volume or co-mention covariate exists for trial-reported AEs; neither is in the model. "
                "kinase_pair uses the molecule family rule (single family = kinase); onc_pair = both area_trials == neoplasms; "
                "bio_pair = both modality == biologic. Shared-trial pairs are not a covariate in the amendment and are not used.")
reg["base_model_no_modifiers"], _, _ = ea.fit(P, BASE, mols, rs[0])
reg["model"], _, _ = ea.fit(P, XCOLS, mols, rs[1])
step(f"base same_target {reg['base_model_no_modifiers']['same_target']}")
for c in XCOLS: step(f"model {c:28s} {reg['model'].get(c)}")
out["B_adjusted_ols"] = reg
json.dump(out, open(f"{T}/modifier_trials.json", "w"), indent=1)
step("wrote results/tables/modifier_trials.json")
