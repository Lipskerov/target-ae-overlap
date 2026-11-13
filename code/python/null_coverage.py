#!/usr/bin/env python3
"""Null-coverage check required before H1 is read (SAP.md §4, SAP_draft Q13).

Builds R null datasets by permuting which molecule owns which profile (target graph fixed), so the true Delta is 0, and
runs the same permutation test and molecule bootstrap on each. Reports the rejection rate at one-sided 0.025 and how often
the 95% bootstrap interval contains 0. The observed Delta of the real data is never computed here.

Usage: null_coverage.py <signals.csv> <label> [R=100] [B_perm=199] [B_boot=300]     env COHORT_DIR as in primary_inference.py"""
import csv, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
DER = os.environ.get("COHORT_DIR", os.path.join(PROJ, "data", "derived"))
sig_path, label = sys.argv[1], sys.argv[2]
R = int(sys.argv[3]) if len(sys.argv) > 3 else 100
BP = int(sys.argv[4]) if len(sys.argv) > 4 else 199
BB = int(sys.argv[5]) if len(sys.argv) > 5 else 300
rng = np.random.default_rng(np.random.SeedSequence(20261112).spawn(16)[9])      # seed index 9 = null-coverage simulation

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(os.path.join(DER, "D2_cohort_molecules.csv")))]
idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
tg = {}
for r in csv.DictReader(open(os.path.join(DER, "D4_molecule_target.csv"))):
    tg.setdefault(r["target_symbol"], []).append(idx[r["ingredient_rxcui"]])
S = np.zeros((n, n), dtype=np.float32)
for m in tg.values():
    m = np.array(sorted(set(m))); S[np.ix_(m, m)] = 1
np.fill_diagonal(S, 0); Dm = 1 - S; np.fill_diagonal(Dm, 0)
rx_terms = {}
for r in csv.DictReader(open(sig_path)):
    if r["ingredient_rxcui"] in idx: rx_terms.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"])
terms = sorted({t for s in rx_terms.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
A = np.zeros((n, len(terms)), dtype=np.float32)
for i, s in rx_terms.items(): A[i, [tix[t] for t in s]] = 1
O = A @ A.T; np.fill_diagonal(O, 0)
ones = np.ones(n, dtype=np.float32)
def delta(M, Sm, Dd, w=ones):
    ns, nd = Sm @ w, Dd @ w; ok = (w > 0) & (ns > 0) & (nd > 0)
    d = (Sm * M) @ w / np.where(ns > 0, ns, 1) - (Dd * M) @ w / np.where(nd > 0, nd, 1)
    return float((w[ok] * d[ok]).sum() / w[ok].sum())
rej = cover = 0; ests = []; t0 = time.time()
for r in range(R):
    p0 = rng.permutation(n); M = O[np.ix_(p0, p0)]            # null dataset: profiles reassigned to molecules
    est = delta(M, S, Dm); ests.append(est)
    ge = 0
    for _ in range(BP):
        p = rng.permutation(n); ge += delta(M, S[np.ix_(p, p)], Dm[np.ix_(p, p)]) >= est
    rej += (1 + ge) / (BP + 1) <= 0.025
    SM, DM_ = S * M, Dm * M; bs = np.empty(BB)
    for b in range(BB):
        w = rng.multinomial(n, np.full(n, 1 / n)).astype(np.float32)
        ns, nd = S @ w, Dm @ w; ok = (w > 0) & (ns > 0) & (nd > 0)
        d = SM @ w / np.where(ns > 0, ns, 1) - DM_ @ w / np.where(nd > 0, nd, 1)
        bs[b] = (w[ok] * d[ok]).sum() / w[ok].sum()
    lo, hi = np.percentile(bs, [2.5, 97.5]); cover += lo <= 0 <= hi
    if (r + 1) % 20 == 0: print(f"  {r+1}/{R} {time.time()-t0:.0f}s", flush=True)
res = dict(label=label, null_datasets=R, B_perm=BP, B_boot=BB, n_molecules=n,
           rejection_rate_at_0025=round(rej / R, 3), ci_contains_zero_rate=round(cover / R, 3),
           null_delta_mean=round(float(np.mean(ests)), 4), null_delta_sd=round(float(np.std(ests, ddof=1)), 4))
json.dump(res, open(os.path.join(PROJ, "results", "tables", f"null_coverage_{label}.json"), "w"), indent=1); print(json.dumps(res))
