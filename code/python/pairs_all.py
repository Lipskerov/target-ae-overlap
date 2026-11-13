#!/usr/bin/env python3
"""Amendment 2 (A6): target versus class, indication and co-mention on ALL cohort pairs.

For every pair of cohort molecules (extraction b): shared profile signals, same target, same ATC level 4 / 3 / 1 group,
same FDA EPC class (same / different / missing), shared indication (>= 1 of each molecule's ten most reported indications),
report volumes, approval-year gap (missing flagged), co-mention share (reports naming both molecules in any role divided by
the smaller molecule's any-role report count).
Outputs (results/tables/, suffix = profile label):
  pairs_all_<label>.csv.gz        one row per pair
  target_vs_class_<label>.json    stratified means with molecule-bootstrap CIs, within-class same-target contrasts,
                                  adjusted models (main effects; with target x ATC-4 interaction), target-level analysis
Usage: pairs_all.py <profile.csv> <label>     env COHORT_DIR (default extraction b), B (bootstrap resamples, default 1000)"""
import csv, gzip, json, os, re, sys, time
import duckdb, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly")); RAW = os.path.join(PROJ, "data/raw/trialytic_db")
CD = os.environ.get("COHORT_DIR", os.path.join(PROJ, "data/derived/faers_q/b")); prof_path, label = sys.argv[1], sys.argv[2]
B = int(os.environ.get("B", "1000")); ss = np.random.SeedSequence(20261112).spawn(40); t0 = time.time()
def step(m): print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)
mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]; idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
names = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))}
tg = {}
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): tg.setdefault(r["target_symbol"], set()).add(idx[r["ingredient_rxcui"]])
S = np.zeros((n, n), bool)
for m in tg.values():
    m = np.array(sorted(m)); S[np.ix_(m, m)] = True
np.fill_diagonal(S, False)
rt = {}
for r in csv.DictReader(open(prof_path)):
    if r["ingredient_rxcui"] in idx: rt.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"])
terms = sorted({t for s in rt.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
A = np.zeros((n, len(terms)), np.float32)
for i, s in rt.items(): A[i, [tix[t] for t in s]] = 1
O = A @ A.T; np.fill_diagonal(O, 0)
# covariates
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
has_atc = np.array([i in atc for i in range(n)]); atc_ok = np.outer(has_atc, has_atc)
A4, A3, A1 = share_mat(5), share_mat(4), share_mat(1)
epc = {}; yr = {}
for r in csv.DictReader(open(f"{RAW}/R11_fda_class_all.csv")):
    if r["ingredient_rxcui"] in idx:
        if r["pharm_class_epc"]: epc[idx[r["ingredient_rxcui"]]] = r["pharm_class_epc"]
        if r["first_approval_year"]: yr[idx[r["ingredient_rxcui"]]] = float(r["first_approval_year"])
for r in csv.DictReader(open(f"{RAW}/R10_drug_approval_all.csv")):
    if r["ingredient_rxcui"] in idx and r["earliest_orig_ap"]: yr[idx[r["ingredient_rxcui"]]] = float(r["earliest_orig_ap"][:4])
ev = np.array([epc.get(i, "") for i in range(n)]); has_epc = ev != ""
EPC_same = ((ev[:, None] == ev[None, :]) & np.outer(has_epc, has_epc)).astype(np.float32); EPC_miss = (~np.outer(has_epc, has_epc)).astype(np.float32)
yv = np.array([yr.get(i, np.nan) for i in range(n)]); YG = np.abs(yv[:, None] - yv[None, :]); Y_miss = np.isnan(YG).astype(np.float32); YG = np.where(np.isnan(YG), np.nanmedian(YG), YG).astype(np.float32)
con = duckdb.connect(os.path.join(STORE, "faers.duckdb"), read_only=True)
con.execute(f"SET threads=4; SET memory_limit='4GB'; SET temp_directory='{os.path.join(STORE, 'duck_tmp')}'; SET preserve_insertion_order=false")
con.execute("CREATE TEMP TABLE coh AS SELECT unnest(?) AS rx, unnest(?) AS k", [mols, list(range(n))])
con.execute("CREATE TEMP TABLE anyr AS SELECT DISTINCT d.pid, c.k FROM drug_rx d JOIN coh c USING (rx) WHERE d.pid IN (SELECT pid FROM reac)")
tot_any = np.zeros(n);
for k, v in con.execute("SELECT k, count(*) FROM anyr GROUP BY 1").fetchall(): tot_any[k] = v
CM = np.zeros((n, n), np.float64)
for bk in range(16):
    for i, j, v in con.execute(f"SELECT a.k, b.k, count(*) FROM anyr a JOIN anyr b ON a.pid = b.pid AND a.k < b.k WHERE a.pid % 16 = {bk} GROUP BY 1, 2").fetchall(): CM[i, j] += v; CM[j, i] += v
step("co-mention counts for all pairs")
tot_s = np.zeros(n)
for k, v in con.execute("SELECT c.k, count(DISTINCT d.pid) FROM drug_rx d JOIN coh c USING (rx) WHERE d.role_rank <= 2 AND NOT d.fdc_only AND d.pid IN (SELECT pid FROM reac) GROUP BY 1").fetchall(): tot_s[k] = v
UNK = ("PRODUCT USED FOR UNKNOWN INDICATION", "DRUG USE FOR UNKNOWN INDICATION", "UNKNOWN", "INDICATION UNKNOWN")
ind = con.execute(f"""SELECT k, indi_pt FROM (SELECT c.k, i.indi_pt, count(DISTINCT i.pid) n, row_number() OVER (PARTITION BY c.k ORDER BY count(DISTINCT i.pid) DESC, i.indi_pt) rn
                      FROM indi i JOIN coh c USING (rx) WHERE i.indi_pt NOT IN {UNK} GROUP BY 1, 2) WHERE rn <= 10""").fetchall(); con.close()
its = sorted({t for _, t in ind}); iix = {t: k for k, t in enumerate(its)}; IM = np.zeros((n, len(its)), np.float32)
for k, t in ind: IM[k, iix[t]] = 1
IND = ((IM @ IM.T) > 0).astype(np.float32); np.fill_diagonal(IND, 0)
SH = CM / np.minimum(tot_any[:, None], tot_any[None, :]).clip(min=1)
lr = np.log10(np.clip(tot_s, 1, None)); LS = lr[:, None] + lr[None, :]; LD = np.abs(lr[:, None] - lr[None, :])
iu = np.triu_indices(n, 1); ii, jj = iu
P = pd.DataFrame(dict(rx_i=np.array(mols)[ii], rx_j=np.array(mols)[jj], same_target=S[iu].astype(int), overlap=O[iu].astype(int), comention_reports=CM[iu].astype(int),
                      comention_share=SH[iu], same_atc4=A4[iu].astype(int), same_atc3=A3[iu].astype(int), same_atc1=A1[iu].astype(int), atc_both=atc_ok[iu].astype(int),
                      same_epc=EPC_same[iu].astype(int), epc_missing=EPC_miss[iu].astype(int), shared_indication=IND[iu].astype(int),
                      logrep_sum=LS[iu], logrep_absdiff=LD[iu], year_gap=YG[iu], year_missing=Y_miss[iu].astype(int)))
P.to_csv(os.path.join(PROJ, f"results/tables/pairs_all_{label}.csv.gz"), index=False, float_format="%.5f", compression="gzip"); step(f"pair table {len(P):,}")
out = dict(label=label, profile=os.path.relpath(prof_path, PROJ), molecules=n, pairs=int(len(P)), same_target_pairs=int(P.same_target.sum()), molecules_with_atc=int(has_atc.sum()), molecules_with_epc=int(has_epc.sum()), B=B)
W = np.random.default_rng(ss[20]).multinomial(n, np.full(n, 1 / n), size=B).astype(np.float64)       # shared bootstrap weights (molecule level)
y = P.overlap.values.astype(float); st = P.same_target.values == 1
def ci(v): lo, hi = np.percentile(v, [2.5, 97.5]); return [round(float(lo), 4), round(float(hi), 4)]
def wmean(mask, w=None):
    if w is None: return float(y[mask].mean()) if mask.any() else float("nan")
    pw = w[ii[mask]] * w[jj[mask]]; return float((pw * y[mask]).sum() / pw.sum()) if pw.sum() > 0 else float("nan")
# stratified means and within-class contrasts
strata = {}
for key, col in (("atc4", "same_atc4"), ("atc3", "same_atc3"), ("epc", "same_epc"), ("indication", "shared_indication")):
    okm = (P.atc_both.values == 1) if key.startswith("atc") else ((P.epc_missing.values == 0) if key == "epc" else np.ones(len(P), bool))
    cv = P[col].values == 1; d = {}
    for nm, m1, m0 in (("same_class", okm & cv & st, okm & cv & ~st), ("different_class", okm & ~cv & st, okm & ~cv & ~st)):
        est = wmean(m1) - wmean(m0); bs = [wmean(m1, w) - wmean(m0, w) for w in W]; bs = [b for b in bs if np.isfinite(b)]
        d[nm] = dict(same_target_pairs=int(m1.sum()), different_target_pairs=int(m0.sum()), mean_same_target=round(wmean(m1), 4), mean_different_target=round(wmean(m0), 4),
                     difference=round(est, 4), ci95=ci(bs), share_bootstrap_le0=round(float(np.mean(np.array(bs) <= 0)), 4))
    strata[key] = d
out["within_class_contrasts"] = strata; step("stratified contrasts")
# adjusted models on all pairs
cols = ["same_target", "same_atc4", "same_atc3", "same_atc1", "atc_missing", "same_epc", "epc_missing", "shared_indication", "logrep_sum", "logrep_absdiff", "year_gap", "year_missing", "comention_share"]
P["atc_missing"] = 1 - P.atc_both
def fit(xcols, seed_rows=None):
    X = np.column_stack([np.ones(len(P))] + [P[c].values.astype(float) for c in xcols]); beta = np.linalg.lstsq(X, y, rcond=None)[0]; bs = []
    for w in W[: min(B, 500)]:
        pw = w[ii] * w[jj]; m = pw > 0; Xw = X[m] * pw[m][:, None]
        try: bs.append(np.linalg.solve(Xw.T @ X[m], Xw.T @ y[m]))
        except np.linalg.LinAlgError: pass
    bs = np.array(bs); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    return {c: dict(beta=round(float(b), 4), ci95=[round(float(l), 4), round(float(h), 4)], share_bootstrap_le0=round(float((bs[:, k] <= 0).mean()), 4)) for k, (c, b, l, h) in enumerate(zip(["intercept"] + xcols, beta, lo, hi))}
out["model_unadjusted"] = fit(["same_target"])
out["model_class_only"] = fit(["same_target", "same_atc4", "same_atc3", "same_atc1", "atc_missing", "same_epc", "epc_missing"])
out["model_full"] = fit(cols)
P["target_x_atc4"] = P.same_target * P.same_atc4
out["model_full_with_interaction"] = fit(cols + ["target_x_atc4"]); step("models")
# target-level analysis (each target: mean overlap among its member pairs minus its members' mean overlap with non-partners)
Sf = S.astype(np.float32); Df = 1 - Sf; np.fill_diagonal(Df, 0)
dmean = (Df * O).sum(1) / Df.sum(1).clip(min=1)
dt = []
for t, m in tg.items():
    m = np.array(sorted(m))
    if len(m) >= 3:
        sub = O[np.ix_(m, m)]; dt.append(sub[np.triu_indices(len(m), 1)].mean() - dmean[m].mean())
dt = np.array(dt); rb = np.random.default_rng(ss[21]); bs = [dt[rb.integers(0, len(dt), len(dt))].mean() for _ in range(5000)]
out["target_level"] = dict(targets=int(len(dt)), mean_contrast=round(float(dt.mean()), 4), ci95=ci(bs), share_targets_positive=round(float((dt > 0).mean()), 4), median=round(float(np.median(dt)), 4))
# descriptives of the profile entries
if "a" in pd.read_csv(prof_path, nrows=1).columns:
    pa = pd.read_csv(prof_path).a; out["profile_entry_reports"] = dict(entries=int(len(pa)), pct_le3=round(100 * float((pa <= 3).mean()), 1), pct_le10=round(100 * float((pa <= 10).mean()), 1), median=float(pa.median()))
tt = (A[ii[st]] * A[jj[st]]).sum(0); top = sorted(zip(terms, tt), key=lambda x: -x[1])[:25]
out["top_shared_terms_same_target"] = [[t, int(v)] for t, v in top]; out["shared_signal_slots_same_target"] = int(tt.sum())
out["pct_same_target_pairs_sharing_none"] = round(100 * float((y[st] == 0).mean()), 1); out["pct_different_target_pairs_sharing_none"] = round(100 * float((y[~st] == 0).mean()), 1)
json.dump(out, open(os.path.join(PROJ, f"results/tables/target_vs_class_{label}.json"), "w"), indent=1)
print(json.dumps({k: out[k] for k in ("molecules", "pairs", "same_target_pairs", "within_class_contrasts", "target_level")}, indent=1))
print("full model same_target:", out["model_full"]["same_target"], "| interaction model:", out["model_full_with_interaction"]["same_target"], out["model_full_with_interaction"]["target_x_atc4"])
