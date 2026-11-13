#!/usr/bin/env python3
"""Primary estimand and inference (SAP_draft.md §1–2), run on a signal table.

Estimand: Delta = mean over cohort molecules i of
    d_i = mean shared signals with same-target partners - mean shared signals with different-target partners.
Inference: node permutation (QAP; target graph fixed, profiles permuted across molecules), one-sided p;
           molecule-level (pigeonhole) bootstrap, percentile 95% CI.
Also reports the pair-averaged means (descriptive) and the same contrast on Jaccard similarity.

Usage: primary_inference.py <signals.csv> <label> [B_perm] [B_boot]
  signals.csv columns: ingredient_rxcui, reaction_pt   (one row per molecule-signal)
Reads data/derived/D2_cohort_molecules.csv and D4_molecule_target.csv. Writes results/tables/primary_<label>.json.
"""
import csv, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
DER = os.environ.get("COHORT_DIR", os.path.join(PROJ, "data", "derived"))   # folder holding D2_cohort_molecules.csv and D4_molecule_target.csv
sig_path, label = sys.argv[1], sys.argv[2]
B_PERM = int(sys.argv[3]) if len(sys.argv) > 3 else 9999
B_BOOT = int(sys.argv[4]) if len(sys.argv) > 4 else 10000
MASTER_SEED = 20261112
REGISTRY = {"perm": 0, "boot": 1}          # fixed spawn index per analysis
ss = np.random.SeedSequence(MASTER_SEED).spawn(8)

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(os.path.join(DER, "D2_cohort_molecules.csv")))]
idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
tg = {}
for r in csv.DictReader(open(os.path.join(DER, "D4_molecule_target.csv"))):
    tg.setdefault(r["target_symbol"], []).append(idx[r["ingredient_rxcui"]])
S = np.zeros((n, n), dtype=np.float32)
for members in tg.values():
    m = np.array(sorted(set(members))); S[np.ix_(m, m)] = 1
np.fill_diagonal(S, 0)
Dm = 1 - S; np.fill_diagonal(Dm, 0)          # different-target indicator, self excluded

rx_terms = {}
for r in csv.DictReader(open(sig_path)):
    if r["ingredient_rxcui"] in idx: rx_terms.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"].upper())
terms = sorted({t for s in rx_terms.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
A = np.zeros((n, len(terms)), dtype=np.float32)
for i, s in rx_terms.items():
    A[i, [tix[t] for t in s]] = 1
k = A.sum(1)
O = A @ A.T; np.fill_diagonal(O, 0)
union = k[:, None] + k[None, :] - O
J = np.divide(O, union, out=np.zeros_like(O), where=union > 0); np.fill_diagonal(J, 0)
# Excess overlap X = O - E (SAP_draft §3.2): E_ij = sum_r q_ir q_jr with q_ir = k_i * c_r^(-ij) / T^(-ij), leave-pair-out,
# in closed form: sum_r (c_r - A_ir - A_jr)^2 = sum c^2 + k_i + k_j + 2 O_ij - 2 s_i - 2 s_j, s_i = sum_r c_r A_ir.
cr = A.sum(0).astype(np.float64); T = float(k.sum()); s_i = (A @ cr.astype(np.float32)).astype(np.float64)
kk = k.astype(np.float64); Tij = T - kk[:, None] - kk[None, :]
cap = float(kk.max() * cr.max() / max(Tij.min(), 1.0))            # if >= 1 the min(1, .) cap would bind and the closed form is not exact
E = (kk[:, None] * kk[None, :]) * ((cr ** 2).sum() + kk[:, None] + kk[None, :] + 2 * O - 2 * s_i[:, None] - 2 * s_i[None, :]) / np.where(Tij > 0, Tij, 1) ** 2
X = (O - E).astype(np.float32); np.fill_diagonal(X, 0)
MEASURES = {"shared_signals": O, "jaccard": J, "excess": X}
WANTED = os.environ.get("MEASURES", "shared_signals,jaccard").split(",")
# optional pair restriction: a CSV of ingredient pairs (rx_i, rx_j) to EXCLUDE from both groups (e.g. pairs sharing a trial condition)
if os.environ.get("EXCLUDE_PAIRS"):
    for r in csv.DictReader(open(os.environ["EXCLUDE_PAIRS"])):
        if r["rx_i"] in idx and r["rx_j"] in idx:
            i_, j_ = idx[r["rx_i"]], idx[r["rx_j"]]; S[i_, j_] = S[j_, i_] = 0; Dm[i_, j_] = Dm[j_, i_] = 0

def delta(M, Smat, Dmat, w=None):
    """Molecule-averaged contrast. w: node multiplicities (bootstrap) or None."""
    if w is None: w = np.ones(n, dtype=np.float32)
    ns, nd = Smat @ w, Dmat @ w
    ok = (w > 0) & (ns > 0) & (nd > 0)
    d = (Smat * M) @ w / np.where(ns > 0, ns, 1) - (Dmat * M) @ w / np.where(nd > 0, nd, 1)
    return float((w[ok] * d[ok]).sum() / w[ok].sum()), int(ok.sum())

res = dict(label=label, signals_file=os.path.relpath(sig_path, PROJ), n_molecules=n, n_targets=len(tg),
           molecules_with_signals=int((k > 0).sum()), mean_signals_per_molecule=round(float(k.mean()), 3),
           same_target_pairs=int(S.sum() / 2), different_target_pairs=int(Dm.sum() / 2),
           B_perm=B_PERM, B_boot=B_BOOT, master_seed=MASTER_SEED)
iu = np.triu_indices(n, 1); s_u = S[iu] > 0; d_u = Dm[iu] > 0
t0 = time.time()
res["excess_cap_ratio"] = round(cap, 3)
for name, M in ((m, MEASURES[m]) for m in WANTED):
    est, n_used = delta(M, S, Dm)
    # node permutation: permute molecule labels of the target graph (equivalent to permuting profiles)
    rng = np.random.default_rng(ss[REGISTRY["perm"]]); ge = 0; null = np.empty(B_PERM)
    for b in range(B_PERM):
        p = rng.permutation(n); Sp = S[np.ix_(p, p)]; Dp = Dm[np.ix_(p, p)]
        null[b], _ = delta(M, Sp, Dp); ge += null[b] >= est
    p_perm = (1 + ge) / (B_PERM + 1)
    # molecule-level bootstrap
    rng = np.random.default_rng(ss[REGISTRY["boot"]]); boots = np.empty(B_BOOT)
    SM, DM_ = S * M, Dm * M
    for b in range(B_BOOT):
        w = rng.multinomial(n, np.full(n, 1 / n)).astype(np.float32)
        ns, nd = S @ w, Dm @ w
        ok = (w > 0) & (ns > 0) & (nd > 0)
        d = SM @ w / np.where(ns > 0, ns, 1) - DM_ @ w / np.where(nd > 0, nd, 1)
        boots[b] = (w[ok] * d[ok]).sum() / w[ok].sum()
    lo, hi = np.percentile(boots, [2.5, 97.5])
    res[name] = dict(delta=round(est, 4), molecules_used=n_used, ci95=[round(float(lo), 4), round(float(hi), 4)],
                     boot_bias=round(float(boots.mean() - est), 4), boot_sd=round(float(boots.std(ddof=1)), 4),
                     p_perm_one_sided=p_perm, perm_null_mean=round(float(null.mean()), 4), perm_null_sd=round(float(null.std(ddof=1)), 4),
                     z_perm=round(float((est - null.mean()) / null.std(ddof=1)), 2), perm_null_max=round(float(null.max()), 4),
                     pair_mean_same=round(float(M[iu][s_u].mean()), 4), pair_mean_different=round(float(M[iu][d_u].mean()), 4))
    print(name, json.dumps(res[name]), f"{time.time()-t0:.0f}s", flush=True)
out = os.path.join(PROJ, "results", "tables", f"primary_{label}.json")
json.dump(res, open(out, "w"), indent=1); print("wrote", os.path.relpath(out, PROJ))
