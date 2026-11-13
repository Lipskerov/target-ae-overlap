#!/usr/bin/env python3
"""AMENDMENT_4 A4.4: molecule-bootstrap 95% CI for the WP5 mean per-target contrast (targets with >= 3 members).
Statistic as in the methods reviewer's check (private/reviews/V2/c7_wp5_ci.py): per target, weighted mean overlap among
member pairs minus members' weighted mean overlap with non-partners; mean over targets. B = 500,
seed SeedSequence(20261113).spawn(15)[14].spawn(3)[2]. Writes results/tables/target_level_molboot_a4.json"""
import csv, json, os
import numpy as np
PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
CD = os.path.join(PROJ, "data/derived/faers_q/b")
mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]; idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
tg = {}
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): tg.setdefault(r["target_symbol"], set()).add(idx[r["ingredient_rxcui"]])
T3 = {t: np.array(sorted(m)) for t, m in tg.items() if len(m) >= 3}
S = np.zeros((n, n), bool)
for m in tg.values(): m = np.array(sorted(m)); S[np.ix_(m, m)] = True
np.fill_diagonal(S, False); Df = (~S).astype(np.float32); np.fill_diagonal(Df, 0)
out = dict(amendment="AMENDMENT_4 A4.4", B=500, targets=len(T3))
for label, pf in (("prr_clean", "profile_clean_prr_k20.csv"), ("ic025_clean2", "profile_clean2_ic025_k20.csv")):
    rt = {}
    for r in csv.DictReader(open(f"{CD}/{pf}")):
        if r["ingredient_rxcui"] in idx: rt.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"])
    terms = sorted({t for s in rt.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
    A = np.zeros((n, len(terms)), np.float32)
    for i, s in rt.items(): A[i, [tix[t] for t in s]] = 1
    O = A @ A.T; np.fill_diagonal(O, 0)
    def stat(w):
        dmean = (Df * O) @ w / np.clip(Df @ w, 1e-9, None); v = []
        for m in T3.values():
            wm = w[m]
            if (wm > 0).sum() < 2: continue
            W = np.outer(wm, wm); np.fill_diagonal(W, 0)
            v.append((W * O[np.ix_(m, m)]).sum() / W.sum() - (wm * dmean[m]).sum() / wm.sum())
        return float(np.mean(v))
    pt = stat(np.ones(n, np.float32))
    rng = np.random.default_rng(np.random.SeedSequence(20261113).spawn(15)[14].spawn(3)[2])
    bs = np.array([stat(rng.multinomial(n, np.full(n, 1 / n)).astype(np.float32)) for _ in range(500)])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    out[label] = dict(point=round(pt, 4), ci95_molecule_bootstrap=[round(float(lo), 4), round(float(hi), 4)], boot_mean=round(float(bs.mean()), 4))
    print(label, out[label], flush=True)
json.dump(out, open(os.path.join(PROJ, "results/tables/target_level_molboot_a4.json"), "w"), indent=1)
