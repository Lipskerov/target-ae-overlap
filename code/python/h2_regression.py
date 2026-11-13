#!/usr/bin/env python3
"""H2 (SAP.md §5, §7): dyadic regression of shared profile signals on same-target with covariates, and the
therapeutic-area 2x2. Pair set: every same-target pair + five seeded different-target pairs each (comention_pairs.csv.gz).
Covariates: same ATC level-3 group, same ATC level-1 group, same FDA EPC class, log10 report volume (sum and absolute
difference), log10 number of clean signals (sum), approval-year gap, co-mention share (larger of the two directions).
Inference: molecule-level bootstrap (B = 2,000, pair weight w_i*w_j, weighted least squares), percentile CI,
one-sided p = share of bootstrap coefficients <= 0. Writes results/tables/h2_regression.json."""
import csv, json, os, re
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
csv.field_size_limit(10 ** 9)
df = pd.read_csv(os.path.join(PROJ, "results/tables/comention_pairs.csv.gz"), dtype={"rx_i": str, "rx_j": str})
atc, epc, yr = {}, {}, {}
for r in csv.DictReader(open(os.path.join(PROJ, "data/raw/trialytic_db/R6_fda_crosswalk_cohort.csv"))):
    rx = r["ingredient_rxcui"]
    codes = set(re.findall(r"[A-Z]\d\d[A-Z]{1,2}", r.get("atc_codes") or ""))
    if codes: atc.setdefault(rx, set()).update(codes)
    if r.get("pharm_class_epc"): epc.setdefault(rx, set()).add(r["pharm_class_epc"])
    y = r.get("first_approval_year")
    if y: yr[rx] = min(int(float(y)), yr.get(rx, 9999))
tot = {r["ingredient_rxcui"]: int(r["reports"]) for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/faers_q/b/molecule_totals.csv")))}
sig = pd.read_csv(os.path.join(PROJ, "data/derived/faers_q/b/signals.csv.gz"), usecols=["ingredient_rxcui", "clean"], dtype={"ingredient_rxcui": str})
ncl = sig[sig.clean].groupby("ingredient_rxcui").size().to_dict()
def share(d, a, b, n):
    x, y = d.get(a), d.get(b)
    if not x or not y: return np.nan
    return float(bool({c[:n] for c in x} & {c[:n] for c in y})) if n else float(bool(x & y))
df["same_atc3"] = [share(atc, a, b, 4) for a, b in zip(df.rx_i, df.rx_j)]
df["same_atc1"] = [share(atc, a, b, 1) for a, b in zip(df.rx_i, df.rx_j)]
df["same_epc"] = [share(epc, a, b, 0) for a, b in zip(df.rx_i, df.rx_j)]
li, lj = np.log10(df.rx_i.map(tot).astype(float)), np.log10(df.rx_j.map(tot).astype(float))
df["logrep_sum"], df["logrep_absdiff"] = li + lj, (li - lj).abs()
df["logsig_sum"] = np.log10(df.rx_i.map(ncl).fillna(0) + 1) + np.log10(df.rx_j.map(ncl).fillna(0) + 1)
df["year_gap"] = [abs(yr[a] - yr[b]) if a in yr and b in yr else np.nan for a, b in zip(df.rx_i, df.rx_j)]
df["comention_share"] = df[["share_i", "share_j"]].max(axis=1)
out = dict(pairs_total=int(len(df)), coverage={c: int(df[c].notna().sum()) for c in ("same_atc3", "same_epc", "year_gap")})
# therapeutic-area 2x2 (SAP §7)
tab = []
for st in (1, 0):
    for sa in (1, 0):
        d = df[(df.same_target == st) & (df.same_atc3 == sa)]
        tab.append(dict(same_target=st, same_atc3=sa, pairs=int(len(d)), mean_overlap=round(float(d.overlap_before.mean()), 4) if len(d) else None,
                        mean_overlap_comention_excluded=round(float(d.overlap_after.mean()), 4) if len(d) else None))
out["area_2x2"] = tab
mols = sorted(set(df.rx_i) | set(df.rx_j)); idx = {m: k for k, m in enumerate(mols)}; n = len(mols)
def fit(d, ycol, xcols, B=2000, seed_index=12):
    d = d.dropna(subset=xcols + [ycol]); X = np.column_stack([np.ones(len(d))] + [d[c].values.astype(float) for c in xcols]); y = d[ycol].values.astype(float)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    ii, jj = d.rx_i.map(idx).values, d.rx_j.map(idx).values
    rng = np.random.default_rng(np.random.SeedSequence(20261112).spawn(32)[seed_index]); bs = []
    for _ in range(B):
        w = rng.multinomial(n, np.full(n, 1 / n)); pw = (w[ii] * w[jj]).astype(float); m = pw > 0
        sw = np.sqrt(pw[m]); 
        try: bs.append(np.linalg.lstsq(X[m] * sw[:, None], y[m] * sw, rcond=None)[0])
        except Exception: pass
    bs = np.array(bs); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    return dict(outcome=ycol, pairs=int(len(d)), same_target_pairs=int((d.same_target == 1).sum()),
                coefficients={c: dict(beta=round(float(b), 4), ci95=[round(float(l), 4), round(float(h), 4)], p_one_sided_le0=round(float((bs[:, k] <= 0).mean()), 4))
                              for k, (c, b, l, h) in enumerate(zip(["intercept"] + xcols, beta, lo, hi))})
base = ["same_target"]
full = ["same_target", "same_atc3", "same_atc1", "logrep_sum", "logrep_absdiff", "logsig_sum", "year_gap", "comention_share"]
out["model_unadjusted"] = fit(df, "overlap_before", base)
out["model_adjusted_H2"] = fit(df, "overlap_before", full)
out["model_adjusted_comention_excluded_outcome"] = fit(df, "overlap_after", full, seed_index=13)
out["model_adjusted_with_epc"] = fit(df, "overlap_before", full + ["same_epc"], seed_index=14)
json.dump(out, open(os.path.join(PROJ, "results/tables/h2_regression.json"), "w"), indent=1)
print(json.dumps(out["coverage"])); print(json.dumps(out["area_2x2"], indent=0))
for m in ("model_unadjusted", "model_adjusted_H2", "model_adjusted_comention_excluded_outcome", "model_adjusted_with_epc"):
    print(m, out[m]["pairs"], json.dumps(out[m]["coefficients"]))
