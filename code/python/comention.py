#!/usr/bin/env python3
"""Co-mention analysis (SAP.md §6, PV review item 2), on extraction (b) of the quarterly-file database.

Pairs: every same-target pair of the cohort plus five seeded different-target pairs per same-target pair (seed index 10).
For each pair (i, j):
  comention      deduplicated reports naming both molecules in any role
  share_i/j      share of i's (j's) suspect reports that also name the other molecule
  overlap_before shared signals between the cleaned top-20 profiles
  overlap_after  the same after each molecule's profile is recomputed WITHOUT the reports that name the partner
                 (background counts unchanged; same signal rule, cleaning and K)
Writes results/tables/comention_pairs.csv.gz and comention_summary.json."""
import csv, gzip, json, os, re, sys, time
import duckdb, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", ""); sys.path.insert(0, os.path.join(ROOT, "scripts"))  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
from compute_faers_prr import BLACKLIST          # noqa: E402
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly"))
CD = os.path.join(PROJ, "data/derived/faers_q/b"); K = 20; N_DIFF = 5
MEASURE = os.environ.get("MEASURE", "prr"); FLAG = os.environ.get("FLAG", "clean"); TAG = f"{MEASURE}_{FLAG}"   # amendment 2: profile variant
ARTEFACT = re.compile(open(os.path.join(PROJ, "private/sql/04_clean.sql")).read().split("reaction_pt ~ '")[1].split("';")[0])
UNK_IND = {"PRODUCT USED FOR UNKNOWN INDICATION", "DRUG USE FOR UNKNOWN INDICATION", "UNKNOWN", "INDICATION UNKNOWN"}
BL = {b.upper() for b in BLACKLIST}
ss = np.random.SeedSequence(20261112).spawn(32); rng = np.random.default_rng(ss[10]); t0 = time.time()
def step(m): print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]; idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
tg = {}
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): tg.setdefault(r["target_symbol"], set()).add(idx[r["ingredient_rxcui"]])
S = np.zeros((n, n), bool)
for m in tg.values():
    m = np.array(sorted(m)); S[np.ix_(m, m)] = True
np.fill_diagonal(S, False)
same = [(i, j) for i in range(n) for j in np.nonzero(S[i])[0] if i < j]
diff = set()
while len(diff) < N_DIFF * len(same):
    i, j = sorted(rng.integers(0, n, 2).tolist())
    if i != j and not S[i, j]: diff.add((i, j))
pairs = [(i, j, 1) for i, j in same] + [(i, j, 0) for i, j in sorted(diff)]
step(f"pairs: same {len(same):,}, different {len(diff):,}")

con = duckdb.connect(os.path.join(STORE, "faers.duckdb"), read_only=True)
con.execute(f"SET threads=4; SET memory_limit='4GB'; SET temp_directory='{os.path.join(STORE, 'duck_tmp')}'; SET preserve_insertion_order=false")
con.execute("CREATE TEMP TABLE coh AS SELECT unnest(?) AS rx", [mols])
con.execute("CREATE TEMP TABLE base AS SELECT pid FROM cases WHERE pid IN (SELECT pid FROM reac)")
Ntot = con.execute("SELECT count(*) FROM base").fetchone()[0]
con.execute("CREATE TEMP TABLE mr AS SELECT DISTINCT d.pid, d.rx FROM drug_rx d JOIN coh USING (rx) SEMI JOIN base b ON b.pid = d.pid WHERE d.role_rank <= 2 AND NOT d.fdc_only")
con.execute("CREATE TEMP TABLE anyrole AS SELECT DISTINCT d.pid, d.rx FROM drug_rx d JOIN coh USING (rx) SEMI JOIN base b ON b.pid = d.pid")
tot = dict(con.execute("SELECT rx, count(*) FROM mr GROUP BY 1").fetchall())
ptn = dict(con.execute("SELECT r.pt, count(*) FROM reac r SEMI JOIN base b ON b.pid = r.pid GROUP BY 1").fetchall())
cnt = con.execute("SELECT m.rx, r.pt, count(*) AS a FROM mr m JOIN reac r ON r.pid = m.pid GROUP BY 1, 2 HAVING count(*) >= 3").df()
ind = con.execute("SELECT i.rx, i.indi_pt, count(DISTINCT i.pid) n FROM indi i SEMI JOIN mr m ON m.pid = i.pid AND m.rx = i.rx GROUP BY 1, 2").df()
ind = ind[~ind.indi_pt.isin(UNK_IND)]; ind["share"] = ind.n / ind.groupby("rx").n.transform("sum")
own = set(map(tuple, ind[ind.share >= 0.01][["rx", "indi_pt"]].values))
step(f"molecule x PT counts {len(cnt):,}")
pdf = pd.DataFrame([(mols[i], mols[j]) for i, j, _ in pairs], columns=["x", "y"])
dp = pd.concat([pdf, pdf.rename(columns={"x": "y", "y": "x"})[["x", "y"]]], ignore_index=True)     # directed: x = profile owner
con.execute("CREATE TEMP TABLE dp AS SELECT * FROM dp")
con.execute("CREATE TEMP TABLE up AS SELECT * FROM pdf")
cm_any = {(x, y): c for x, y, c in con.execute("SELECT u.x, u.y, count(*) FROM up u JOIN anyrole a ON a.rx = u.x JOIN anyrole b ON b.rx = u.y AND b.pid = a.pid GROUP BY 1, 2").fetchall()}
step(f"co-mention counts for {len(cm_any):,} pairs with at least one shared report")
con.execute("CREATE TEMP TABLE sub AS SELECT d.x, d.y, m.pid FROM dp d JOIN mr m ON m.rx = d.x JOIN anyrole o ON o.rx = d.y AND o.pid = m.pid")
subn = {(x, y): c for x, y, c in con.execute("SELECT x, y, count(*) FROM sub GROUP BY 1, 2").fetchall()}
subc = con.execute("SELECT s.x, s.y, r.pt, count(*) AS a FROM sub s JOIN reac r ON r.pid = s.pid GROUP BY 1, 2, 3").df()
step(f"directed exclusion rows {len(subc):,}")
con.close()

# extended cleaning (clean2) is taken from the flags in signals.csv.gz: terms flagged clean but not clean2 are dropped
DROP2 = set()
if FLAG == "clean2":
    sg = pd.read_csv(f"{CD}/signals.csv.gz", dtype={"ingredient_rxcui": str}, keep_default_na=False)
    b_ = lambda c: sg[c].astype(str).str.lower().isin(["true", "1"])
    DROP2 = set(map(tuple, sg[b_("clean") & ~b_("clean2")][["ingredient_rxcui", "reaction_pt"]].values))
# per-molecule arrays
M = {}
for rx, g in cnt.groupby("rx"):
    pts = g.pt.values; a = g.a.values.astype(float); ac = np.array([ptn[p] for p in pts], float)
    cleanmask = np.array([(p not in BL) and not ARTEFACT.search(p) and (rx, p) not in own and (rx, p) not in DROP2 for p in pts])
    M[rx] = (pts, a, ac, cleanmask, {p: k for k, p in enumerate(pts)})
def top(rx, minus=None, nsub=0):
    pts, a, ac, cm, pos = M[rx]; a = a.copy(); ab = tot[rx] - nsub
    if minus is not None:
        for p, v in minus:
            k = pos.get(p)
            if k is not None: a[k] -= v
    b, c = ab - a, ac - a; d = Ntot - ab - c
    ok = cm & (a >= 3) & (b > 0) & (c > 0) & (d > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        prr = (a / ab) / (c / (c + d))
        e = np.stack([ab * ac / Ntot, ab * (Ntot - ac) / Ntot, (Ntot - ab) * ac / Ntot, (Ntot - ab) * (Ntot - ac) / Ntot])
        chi2 = (((np.abs(np.stack([a, b, c, d]) - e) - 0.5) ** 2) / e).sum(0)
    ok &= (np.round(prr, 2) >= 2) & (np.round(chi2, 1) >= 4)
    k = np.nonzero(ok)[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        ic025 = np.log2((a + 0.5) / (e[0] + 0.5)) - 3.3 * (a + 0.5) ** -0.5 - 2 * (a + 0.5) ** -1.5
    score = ic025 if MEASURE == "ic025" else prr
    order = sorted(k, key=lambda q: (-score[q], -a[q], pts[q]))[:K]
    return {pts[q] for q in order}
base_top = {rx: top(rx) for rx in mols if rx in M}
minus = {(x, y): list(zip(g.pt.values, g.a.values)) for (x, y), g in subc.groupby(["x", "y"])}
rows = []
for i, j, s in pairs:
    x, y = mols[i], mols[j]
    if x not in M or y not in M: continue
    tx = top(x, minus.get((x, y)), subn.get((x, y), 0)) if (x, y) in subn else base_top[x]
    ty = top(y, minus.get((y, x)), subn.get((y, x), 0)) if (y, x) in subn else base_top[y]
    rows.append((x, y, s, cm_any.get((x, y), 0), subn.get((x, y), 0) / tot[x], subn.get((y, x), 0) / tot[y],
                 len(base_top[x] & base_top[y]), len(tx & ty), len(tx), len(ty)))
df = pd.DataFrame(rows, columns=["rx_i", "rx_j", "same_target", "comention_reports", "share_i", "share_j", "overlap_before", "overlap_after", "k_i_after", "k_j_after"])
df.to_csv(os.path.join(PROJ, f"results/tables/comention_pairs_{TAG}.csv.gz"), index=False, float_format="%.5f", compression="gzip")
step("pairs scored")
df["share_max"] = df[["share_i", "share_j"]].max(axis=1)
def grp(d): return dict(pairs=int(len(d)), median_share_max=round(float(d.share_max.median()), 5), pct_share_ge_5pct=round(100 * float((d.share_max >= 0.05).mean()), 1),
                        pct_no_comention=round(100 * float((d.comention_reports == 0).mean()), 1),
                        mean_overlap_before=round(float(d.overlap_before.mean()), 3), mean_overlap_after=round(float(d.overlap_after.mean()), 3))
out = dict(extraction="b", measure=MEASURE, flag=FLAG, K=K, cohort_molecules=n, different_pairs_per_same_pair=N_DIFF, same=grp(df[df.same_target == 1]), different=grp(df[df.same_target == 0]), strata=[])
for lo, hi in ((0, 0.001), (0.001, 0.01), (0.01, 0.05), (0.05, 1.01)):
    st = dict(share_from=lo, share_to=hi)
    for nm, s in (("same", 1), ("different", 0)):
        d = df[(df.same_target == s) & (df.share_max >= lo) & (df.share_max < hi)]
        st[nm] = dict(pairs=int(len(d)), before=round(float(d.overlap_before.mean()), 3) if len(d) else None, after=round(float(d.overlap_after.mean()), 3) if len(d) else None)
    out["strata"].append(st)
# molecule-cluster bootstrap of (same - different) pair-mean difference, before and after
ii = df.rx_i.map(idx).values; jj = df.rx_j.map(idx).values; s = df.same_target.values == 1
brng = np.random.default_rng(ss[11]); B = 2000; res = {"before": [], "after": []}
for _ in range(B):
    w = brng.multinomial(n, np.full(n, 1 / n)); pw = w[ii] * w[jj]
    for nm, v in (("before", df.overlap_before.values), ("after", df.overlap_after.values)):
        ws, wd = pw[s].sum(), pw[~s].sum()
        if ws > 0 and wd > 0: res[nm].append((pw[s] * v[s]).sum() / ws - (pw[~s] * v[~s]).sum() / wd)
for nm, v in (("before", df.overlap_before), ("after", df.overlap_after)):
    est = float(v[s].mean() - v[~s].mean()); lo, hi = np.percentile(res[nm], [2.5, 97.5])
    out[f"difference_{nm}"] = dict(estimate=round(est, 3), ci95=[round(float(lo), 3), round(float(hi), 3)], bootstrap="molecule-level, B=2000")
for nm in ("before", "after"):
    v = df[f"overlap_{nm}"]; out[f"ratio_{nm}"] = round(float(v[s].mean() / v[~s].mean()), 2) if v[~s].mean() > 0 else None
json.dump(out, open(os.path.join(PROJ, f"results/tables/comention_summary_{TAG}.json"), "w"), indent=1); print(json.dumps(out, indent=1))
