#!/usr/bin/env python3
"""Amendment 3, WP4 — profile length and evidence floor (seed index 13).

Grid: K in {10, 20, 50} x ranking in {PRR standard cleaning ('clean', sort prr), IC025 extended cleaning ('clean2',
sort ic025)} x floor in {none, a >= 10}. Extraction b. Profiles are rebuilt from data/derived/faers_q/b/signals.csv.gz
with the ranking and tie-break of faers_q_profiles.py: filter on the cleaning flag, sort by
[ingredient_rxcui asc, measure desc, a desc, reaction_pt asc], keep the first K per molecule. The floor drops
signals with a < 10 BEFORE the top-K is taken (amendment 3, WP4).

Estimate: primary Delta (molecule-averaged same-target minus different-target mean shared signals) and the same
contrast on Jaccard, computed exactly as primary_inference.py (copied, not imported: that script runs on import and
writes its own output file). 999 node permutations, 1,000 molecule bootstrap resamples.
Seed: numpy.random.SeedSequence(20261113).spawn(15)[13]; its .spawn(2) gives the permutation stream [0] and the
bootstrap stream [1].

Reproduction check (amendment 3): the rebuilt K=20 no-floor profiles must equal the stored
profile_clean_prr_k20.csv / profile_clean2_ic025_k20.csv as row sets, and give Delta 0.8718 (PRR) / 1.7401 (IC025).
If it fails, the script stops without computing the grid.

Measures are recomputed unrounded (see RECOMPUTE below) because signals.csv.gz rounds them to 4 decimals and the
rounding creates rank ties the stored profiles do not have. Reads ~/faers_quarterly/faers.duckdb read-only (env FAERS_STORE).
Writes results/tables/sensitivity_grid_R2.csv and .json. Writes nothing else.
"""
import csv, json, os, sys, time
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
CD = os.path.join(PROJ, "data/derived/faers_q/b")
OUT_CSV = os.path.join(PROJ, "results/tables/sensitivity_grid_R2.csv")
OUT_JSON = os.path.join(PROJ, "results/tables/sensitivity_grid_R2.json")
B_PERM, B_BOOT = 999, 1000
SEED_INDEX = 13
seed = np.random.SeedSequence(20261113).spawn(15)[SEED_INDEX]
perm_ss, boot_ss = seed.spawn(2)
RANKINGS = {"prr_standard": ("clean", "prr", "profile_clean_prr_k20.csv", 0.8718),
            "ic025_extended": ("clean2", "ic025", "profile_clean2_ic025_k20.csv", 1.7401)}
KS = (10, 20, 50); FLOORS = {"none": 0, "a>=10": 10}

# ---- cohort and target graph (identical to primary_inference.py with COHORT_DIR = faers_q/b)
mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(os.path.join(CD, "D2_cohort_molecules.csv")))]
idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
tg = {}
for r in csv.DictReader(open(os.path.join(CD, "D4_molecule_target.csv"))):
    tg.setdefault(r["target_symbol"], []).append(idx[r["ingredient_rxcui"]])
S = np.zeros((n, n), dtype=np.float32)
for members in tg.values():
    m = np.array(sorted(set(members))); S[np.ix_(m, m)] = 1
np.fill_diagonal(S, 0)
Dm = 1 - S; np.fill_diagonal(Dm, 0)
iu = np.triu_indices(n, 1); s_u = S[iu] > 0; d_u = Dm[iu] > 0


def matrices(prof):
    rx_terms = {}
    for rx, pt in zip(prof.ingredient_rxcui.astype(str), prof.reaction_pt.astype(str)):
        if rx in idx: rx_terms.setdefault(idx[rx], set()).add(pt.upper())
    terms = sorted({t for s in rx_terms.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
    A = np.zeros((n, len(terms)), dtype=np.float32)
    for i, s in rx_terms.items(): A[i, [tix[t] for t in s]] = 1
    k = A.sum(1)
    O = A @ A.T; np.fill_diagonal(O, 0)
    union = k[:, None] + k[None, :] - O
    J = np.divide(O, union, out=np.zeros_like(O), where=union > 0); np.fill_diagonal(J, 0)
    return k, {"shared_signals": O, "jaccard": J}


def delta(M, Smat, Dmat, w=None):
    if w is None: w = np.ones(n, dtype=np.float32)
    ns, nd = Smat @ w, Dmat @ w
    ok = (w > 0) & (ns > 0) & (nd > 0)
    d = (Smat * M) @ w / np.where(ns > 0, ns, 1) - (Dmat * M) @ w / np.where(nd > 0, nd, 1)
    return float((w[ok] * d[ok]).sum() / w[ok].sum()), int(ok.sum())


def infer(M):
    est, n_used = delta(M, S, Dm)
    rng = np.random.default_rng(perm_ss); ge = 0; null = np.empty(B_PERM)
    for b in range(B_PERM):
        p = rng.permutation(n); null[b], _ = delta(M, S[np.ix_(p, p)], Dm[np.ix_(p, p)]); ge += null[b] >= est
    rng = np.random.default_rng(boot_ss); boots = np.empty(B_BOOT); SM, DM_ = S * M, Dm * M
    for b in range(B_BOOT):
        w = rng.multinomial(n, np.full(n, 1 / n)).astype(np.float32)
        ns, nd = S @ w, Dm @ w
        ok = (w > 0) & (ns > 0) & (nd > 0)
        d = SM @ w / np.where(ns > 0, ns, 1) - DM_ @ w / np.where(nd > 0, nd, 1)
        boots[b] = (w[ok] * d[ok]).sum() / w[ok].sum()
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return dict(delta=round(est, 4), molecules_used=n_used, ci95=[round(float(lo), 4), round(float(hi), 4)],
                p_perm_one_sided=(1 + ge) / (B_PERM + 1), perm_null_mean=round(float(null.mean()), 4),
                pair_mean_same=round(float(M[iu][s_u].mean()), 4), pair_mean_different=round(float(M[iu][d_u].mean()), 4))


sig = pd.read_csv(os.path.join(CD, "signals.csv.gz"), dtype={"ingredient_rxcui": str})
# signals.csv.gz stores prr/ic025 rounded to 4 decimals, which creates ties at the K boundary that the stored profiles
# (ranked on unrounded values) do not have. Recompute both measures unrounded with the formulas of faers_q_profiles.py:
# a (stored), ab = molecule report total (molecule_totals.csv), ac = background reports with the PT and N = background
# size, both read read-only from the FAERS DuckDB with the extraction-b background definition (case_filter TRUE).
import duckdb
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly"))
con = duckdb.connect(os.path.join(STORE, "faers.duckdb"), read_only=True)
con.execute("SET threads=4; SET memory_limit='4GB'")
con.execute("CREATE OR REPLACE TEMP TABLE base AS SELECT pid FROM cases WHERE TRUE AND pid IN (SELECT pid FROM reac)")
N = con.execute("SELECT count(*) FROM base").fetchone()[0]
ptn = con.execute("SELECT r.pt, count(*) n FROM reac r SEMI JOIN base b ON b.pid = r.pid GROUP BY 1").df().set_index("pt")["n"]
con.close()
N_summary = json.load(open(os.path.join(CD, "summary.json")))["reports_in_background"]
tot = pd.read_csv(os.path.join(CD, "molecule_totals.csv"), dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui")["reports"]
a_ = sig.a.values.astype(float); ab_ = sig.ingredient_rxcui.map(tot).values.astype(float); ac_ = sig.reaction_pt.map(ptn).values.astype(float)
c_ = ac_ - a_; d_ = N - ab_ - c_
with np.errstate(divide="ignore", invalid="ignore"):
    prr_x = (a_ / ab_) / (c_ / (c_ + d_))
    ic_x = np.log2((a_ + 0.5) / (ab_ * ac_ / N + 0.5)) - 3.3 * (a_ + 0.5) ** -0.5 - 2 * (a_ + 0.5) ** -1.5
recompute = dict(N_duckdb=int(N), N_summary=int(N_summary), rows=len(sig), missing_counts=int(np.isnan(ab_).sum() + np.isnan(ac_).sum()),
                 prr_rows_off_by_more_than_rounding=int((np.abs(np.round(prr_x, 4) - sig.prr.values) > 1.01e-4).sum()),
                 ic025_rows_off_by_more_than_rounding=int((np.abs(np.round(ic_x, 4) - sig.ic025.values) > 1.01e-4).sum()))
sig["prr"], sig["ic025"] = prr_x, ic_x
print("RECOMPUTE", json.dumps(recompute), flush=True)


def build(flag, meas, K, floor):
    s = sig[sig[flag]]
    if floor: s = s[s.a >= floor]
    s = s.sort_values(["ingredient_rxcui", meas, "a", "reaction_pt"], ascending=[True, False, False, True])
    return s[s.groupby("ingredient_rxcui").cumcount() < K][["ingredient_rxcui", "reaction_pt", "a", meas]]


t0 = time.time()
print(f"cohort molecules {n}, targets {len(tg)}, same-target pairs {int(S.sum() / 2)}, signals rows {len(sig)}", flush=True)
# ---------------- reproduction check
repro = {}; repro_ok = True
for rk, (flag, meas, fname, expected) in RANKINGS.items():
    stored = pd.read_csv(os.path.join(CD, fname), dtype={"ingredient_rxcui": str})
    rebuilt = build(flag, meas, 20, 0)
    a_set = set(zip(stored.ingredient_rxcui, stored.reaction_pt)); b_set = set(zip(rebuilt.ingredient_rxcui, rebuilt.reaction_pt))
    k, Ms = matrices(rebuilt); est, _ = delta(Ms["shared_signals"], S, Dm)
    m = stored.merge(rebuilt, on=["ingredient_rxcui", "reaction_pt"], suffixes=("_s", "_r"))
    maxdiff = float(np.abs(m[meas + "_s"] - m[meas + "_r"]).max())
    ok = (a_set == b_set) and len(stored) == len(rebuilt) and round(est, 4) == expected and maxdiff < 1e-9
    repro[rk] = dict(stored_rows=len(stored), rebuilt_rows=len(rebuilt), only_in_stored=len(a_set - b_set), only_in_rebuilt=len(b_set - a_set),
                     row_sets_equal=a_set == b_set, max_abs_diff_measure_vs_stored=maxdiff, delta_rebuilt=round(est, 4), delta_expected=expected, ok=ok)
    repro_ok &= ok
    print("REPRO", rk, json.dumps(repro[rk]), flush=True)
print("REPRO_CHECK", "PASS" if repro_ok else "FAIL", flush=True)
res = dict(recompute=recompute, wp="WP4", amendment="protocol/AMENDMENT_3.md", seed="SeedSequence(20261113).spawn(15)[13].spawn(2) -> [perm, boot]",
           B_perm=B_PERM, B_boot=B_BOOT, n_molecules=n, n_targets=len(tg), same_target_pairs=int(S.sum() / 2),
           different_target_pairs=int(Dm.sum() / 2), repro_check="PASS" if repro_ok else "FAIL", repro=repro, grid=[])
if not repro_ok:
    json.dump(res, open(OUT_JSON, "w"), indent=1); print("stopped: reproduction check failed"); sys.exit(1)

# ---------------- grid
rows = []
for rk, (flag, meas, _, _) in RANKINGS.items():
    for floor_name, floor in FLOORS.items():
        for K in KS:
            prof = build(flag, meas, K, floor)
            k, Ms = matrices(prof)
            r = dict(ranking=rk, cleaning=flag, K=K, floor=floor_name, profile_rows=len(prof),
                     molecules_with_signals=int((k > 0).sum()), mean_signals_per_molecule=round(float(k.mean()), 3),
                     molecules_at_K=int((k >= K).sum()))
            for mname, M in Ms.items():
                x = infer(M)
                pre = "" if mname == "shared_signals" else "jaccard_"
                r.update({pre + "delta": x["delta"], pre + "ci_lo": x["ci95"][0], pre + "ci_hi": x["ci95"][1],
                          pre + "p_perm": x["p_perm_one_sided"], pre + "molecules_used": x["molecules_used"],
                          pre + "pair_mean_same": x["pair_mean_same"], pre + "pair_mean_different": x["pair_mean_different"]})
            rows.append(r); res["grid"].append(r)
            print(f"[{time.time() - t0:5.0f}s]", json.dumps(r), flush=True)
            pd.DataFrame(rows).to_csv(OUT_CSV, index=False); json.dump(res, open(OUT_JSON, "w"), indent=1)
print("wrote", os.path.relpath(OUT_CSV, PROJ), os.path.relpath(OUT_JSON, PROJ))
