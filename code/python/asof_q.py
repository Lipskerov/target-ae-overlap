#!/usr/bin/env python3
"""H3 — the as-of-date test on the quarterly-file database (SAP.md §8; SAP_draft §6.3–6.4).

Matched set: a distinct pair (F, S) where F is the first-approved member of a cohort target, S a later-approved member,
d = S's approval date >= 2004-01-01, and F has >= 100 suspect reports first received before d.
F's profile is rebuilt from cases first received before d (extraction b, full reaction universe, same signal rule,
cleaning and K = 20) and compared with S's current profile.
  Control A  up to 5 first-approved molecules of other targets, approved before d, sharing no target with S or F, with
             >= 100 as-of reports, nearest to F in log10 as-of report volume; profiled at the same d.
  Control B  F's as-of profile against every cohort molecule S' approved within +-3 years of d (2004 or later) that shares
             no target with F or S.
e_A = O(F@d, S) - mean_C O(C@d, S);  e_B = O(F@d, S) - mean_S' O(F@d, S').
Inference: clusters = connected components of the F-S graph; cluster bootstrap (B = 10,000) and one-sided cluster
sign-flip test (B = 9,999). Own-indication cleaning uses the molecule's current indication list (not rebuilt as-of).
Writes results/tables/asof_sets.csv and asof_summary.json."""
import csv, json, os, re, sys, time
from collections import defaultdict
import duckdb, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", ""); sys.path.insert(0, os.path.join(ROOT, "scripts"))  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
from compute_faers_prr import BLACKLIST          # noqa: E402
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly")); CD = os.path.join(PROJ, "data/derived/faers_q/b"); K = 20
ARTEFACT = re.compile(open(os.path.join(PROJ, "private/sql/04_clean.sql")).read().split("reaction_pt ~ '")[1].split("';")[0])
UNK_IND = {"PRODUCT USED FOR UNKNOWN INDICATION", "DRUG USE FOR UNKNOWN INDICATION", "UNKNOWN", "INDICATION UNKNOWN"}
MEASURE = os.environ.get("MEASURE", "prr"); FLAG = os.environ.get("FLAG", "clean"); TAG = f"{MEASURE}_{FLAG}"
BL = {b.upper() for b in BLACKLIST}; MIN_REP = int(os.environ.get("ASOF_MIN_REPORTS", "100"))
ss = np.random.SeedSequence(20261112).spawn(32); t0 = time.time()
def step(m): print(f"[{time.time()-t0:5.0f}s] {m}", flush=True)

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]; cs = set(mols)
name = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))}
tg = defaultdict(set); rxt = defaultdict(set)
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): tg[r["target_symbol"]].add(r["ingredient_rxcui"]); rxt[r["ingredient_rxcui"]].add(r["target_symbol"])
ap = {}
for r in csv.DictReader(open(os.path.join(PROJ, "data/raw/trialytic_db/R10_drug_approval_all.csv"))):
    if r["earliest_orig_ap"] and r["ingredient_rxcui"] in cs: ap[r["ingredient_rxcui"]] = int(r["earliest_orig_ap"].replace("-", ""))
ctargets = [t for t, m in tg.items() if len(m) >= 3]
first = {}
for t in ctargets:
    dated = sorted((ap[m], m) for m in tg[t] if m in ap)
    if dated: first[t] = dated[0][1]
firsts = sorted(set(first.values()))
sets = sorted({(first[t], s) for t in first for s in tg[t] if s in ap and ap[s] > ap[first[t]] and ap[s] >= 20040101})
dates = sorted({ap[s] for _, s in sets}); dpos = {d: k for k, d in enumerate(dates)}; D = np.array(dates)
step(f"cohort {len(mols)}, with approval date {len(ap)}, targets with a dated first {len(first)}, distinct firsts {len(firsts)}, candidate (F,S) sets {len(sets)}, dates {len(dates)}")
prof = {}
for kind, fn in (("clean", f"profile_{FLAG}_{MEASURE}_k20.csv"), ("raw", "profile_raw_k20.csv")):
    p = defaultdict(set)
    for r in csv.DictReader(open(f"{CD}/{fn}")): p[r["ingredient_rxcui"]].add(r["reaction_pt"])
    prof[kind] = p
DROP2 = set()
if FLAG == "clean2":
    sg = pd.read_csv(f"{CD}/signals.csv.gz", dtype={"ingredient_rxcui": str}, keep_default_na=False)
    b_ = lambda c: sg[c].astype(str).str.lower().isin(["true", "1"])
    DROP2 = set(map(tuple, sg[b_("clean") & ~b_("clean2")][["ingredient_rxcui", "reaction_pt"]].values))
atc1 = defaultdict(set)
for r in csv.DictReader(open(os.path.join(PROJ, "data/raw/trialytic_db/R9_rxnorm_atc.csv"))): atc1[r["rxcui"]].add(r["code"][:1])

con = duckdb.connect(os.path.join(STORE, "faers.duckdb"), read_only=True)
con.execute(f"SET threads=4; SET memory_limit='4GB'; SET temp_directory='{os.path.join(STORE, 'duck_tmp')}'; SET preserve_insertion_order=false")
# bin of a case = number of successor dates <= its initial receipt date; the case is in the as-of data at date index t iff bin <= t
c = con.execute("SELECT pid, init_dt FROM cases WHERE pid IN (SELECT pid FROM reac)").fetchnumpy()
cb = pd.DataFrame({"pid": c["pid"], "bin": np.searchsorted(D, c["init_dt"], side="right").astype(np.int32)}); del c
con.execute("CREATE TEMP TABLE cb AS SELECT * FROM cb"); del cb
nb = len(dates) + 1
Nbin = np.zeros(nb);
for b, v in con.execute("SELECT bin, count(*) FROM cb GROUP BY 1").fetchall(): Nbin[b] = v
Ncum = np.cumsum(Nbin)
bg = con.execute("SELECT r.pt, cb.bin, count(*) n FROM reac r JOIN cb USING (pid) GROUP BY 1, 2").df()
pts = sorted(bg.pt.unique()); pix = {p: k for k, p in enumerate(pts)}
BG = np.zeros((nb, len(pts)), dtype=np.float32); BG[bg.bin.values, bg.pt.map(pix).values] = bg.n.values; BG = np.cumsum(BG, axis=0); del bg
step(f"background snapshots {BG.shape}")
con.execute("CREATE TEMP TABLE fr AS SELECT unnest(?) AS rx", [firsts])
con.execute("CREATE TEMP TABLE mrf AS SELECT DISTINCT d.pid, d.rx, cb.bin FROM drug_rx d JOIN fr USING (rx) JOIN cb USING (pid) WHERE d.role_rank <= 2 AND NOT d.fdc_only")
vol = np.zeros((len(firsts), nb)); fpos = {f: k for k, f in enumerate(firsts)}
for rx, b, v in con.execute("SELECT rx, bin, count(*) FROM mrf GROUP BY 1, 2").fetchall(): vol[fpos[rx], b] = v
vol = np.cumsum(vol, axis=1)                                  # vol[f, t] = F's suspect reports first received before dates[t]
mc = con.execute("SELECT m.rx, m.bin, r.pt, count(*) n FROM mrf m JOIN reac r USING (pid) GROUP BY 1, 2, 3").df()
ind = con.execute("SELECT i.rx, i.indi_pt, count(DISTINCT i.pid) n FROM indi i JOIN fr USING (rx) SEMI JOIN mrf m ON m.pid = i.pid AND m.rx = i.rx GROUP BY 1, 2").df()
con.close()
ind = ind[~ind.indi_pt.isin(UNK_IND)]; ind["share"] = ind.n / ind.groupby("rx").n.transform("sum")
own = set(map(tuple, ind[ind.share >= 0.01][["rx", "indi_pt"]].values))
step(f"first-molecule day-bin x PT rows {len(mc):,}")
by_rx = {rx: g for rx, g in mc.groupby("rx")}; del mc
artef = np.array([(p in BL) or bool(ARTEFACT.search(p)) for p in pts]); blk = np.array([p in BL for p in pts])
cache = {}
def profile(rx, t):
    """(raw top-20, clean top-20) of first molecule rx from cases first received before dates[t]."""
    if (rx, t) in cache: return cache[(rx, t)]
    g = by_rx.get(rx); res = (set(), set())
    if g is not None:
        m = g.bin.values <= t; a = np.zeros(len(pts)); np.add.at(a, g.pt.map(pix).values[m], g.n.values[m])
        ab = vol[fpos[rx], t]; N = Ncum[t]; ac = BG[t].astype(float)
        b, c, d = ab - a, ac - a, N - ab - (ac - a)
        with np.errstate(divide="ignore", invalid="ignore"):
            prr = (a / ab) / (c / (c + d))
            e = np.stack([ab * ac / N, ab * (N - ac) / N, (N - ab) * ac / N, (N - ab) * (N - ac) / N])
            chi2 = (((np.abs(np.stack([a, b, c, d]) - e) - 0.5) ** 2) / e).sum(0)
        ok = (a >= 3) & (b > 0) & (c > 0) & (d > 0) & (np.round(prr, 2) >= 2) & (np.round(chi2, 1) >= 4) & ~blk
        with np.errstate(divide="ignore", invalid="ignore"):
            ic025 = np.log2((a + 0.5) / (e[0] + 0.5)) - 3.3 * (a + 0.5) ** -0.5 - 2 * (a + 0.5) ** -1.5
        score = ic025 if MEASURE == "ic025" else prr
        k = np.nonzero(ok)[0]; order = sorted(k, key=lambda q: (-prr[q], -a[q], pts[q]))
        raw = {pts[q] for q in order[:K]}
        order = sorted(k, key=lambda q: (-score[q], -a[q], pts[q]))
        clean = [pts[q] for q in order if not artef[q] and (rx, pts[q]) not in own and (rx, pts[q]) not in DROP2][:K]
        res = (raw, set(clean))
    cache[(rx, t)] = res; return res

rows = []; dropped = dict(F_under_min_reports=0, fewer_than_2_controls=0)
for F, S in sets:
    t = dpos[ap[S]]; vF = vol[fpos[F], t]
    if vF < MIN_REP: dropped["F_under_min_reports"] += 1; continue
    bad = rxt[S] | rxt[F]
    cand = [(abs(np.log10(vol[fpos[c_], t]) - np.log10(vF)), c_) for c_ in firsts
            if c_ != F and c_ != S and ap[c_] < ap[S] and not (rxt[c_] & bad) and vol[fpos[c_], t] >= MIN_REP]
    ctrl = [c_ for _, c_ in sorted(cand)[:5]]
    ctrlC = [c_ for _, c_ in sorted(x for x in cand if atc1.get(x[1], set()) & atc1.get(F, set()))[:5]]     # control C: same ATC level-1 group as F
    if len(ctrl) < 2: dropped["fewer_than_2_controls"] += 1; continue
    y0 = ap[S] // 10000
    Sp = [m for m in mols if m in ap and m not in (F, S) and ap[m] >= 20040101 and abs(ap[m] // 10000 - y0) <= 3 and not (rxt[m] & bad)]
    rec = dict(F=F, S=S, F_name=name[F], S_name=name[S], d=ap[S], F_asof_reports=int(vF), n_controls_A=len(ctrl), n_controls_B=len(Sp))
    for ki, kind in enumerate(("raw", "clean")):
        pf = profile(F, t)[ki]; ps = prof[kind].get(S, set())
        oF = len(pf & ps); oA = [len(profile(c_, t)[ki] & ps) for c_ in ctrl]; oB = [len(pf & prof[kind].get(m, set())) for m in Sp]
        if ki == 1:
            oC = [len(profile(c_, t)[ki] & ps) for c_ in ctrlC]
            rec.update(n_controls_C=len(ctrlC), meanC_clean=float(np.mean(oC)) if oC else np.nan, eC_clean=(oF - float(np.mean(oC))) if len(oC) >= 2 else np.nan)
        rec.update({f"O_F_{kind}": oF, f"meanA_{kind}": float(np.mean(oA)), f"meanB_{kind}": float(np.mean(oB)) if oB else np.nan,
                    f"eA_{kind}": oF - float(np.mean(oA)), f"eB_{kind}": (oF - float(np.mean(oB))) if oB else np.nan,
                    f"thetaA_{kind}": float(np.mean([(oF > x) + 0.5 * (oF == x) for x in oA])), f"anyA_{kind}": float(np.mean([x >= 1 for x in oA]))})
    rows.append(rec)
df = pd.DataFrame(rows); df.to_csv(os.path.join(PROJ, f"results/tables/asof_sets_{TAG}.csv"), index=False)
step(f"matched sets {len(df)}; dropped {dropped}")
# clusters = connected components of the F-S graph
par = {}
def find(x):
    while par.setdefault(x, x) != x: par[x] = par[par[x]]; x = par[x]
    return x
for F, S in zip(df.F, df.S): par[find(F)] = find(S)
df["cluster"] = [find(F) for F in df.F]; cl = df.cluster.unique(); gi = df.groupby("cluster").indices
sizes = df.groupby("cluster").size().sort_values(ascending=False)
big = sizes.index[0]; cm_ = df.groupby("cluster").eA_clean.mean()
out_extra = dict(largest_cluster_sets=int(sizes.iloc[0]), clusters_with_1_set=int((sizes == 1).sum()), eA_clean_without_largest_cluster=round(float(df[df.cluster != big].eA_clean.mean()), 4),
                 eA_clean_unweighted_cluster_mean=round(float(cm_.mean()), 4), sets_with_control_C=int(df.eC_clean.notna().sum()))
out = dict(measure=MEASURE, flag=FLAG, **out_extra, min_reports=MIN_REP, candidate_sets=len(sets), matched_sets=int(len(df)), dropped=dropped, clusters=int(len(cl)), distinct_F=int(df.F.nunique()), distinct_S=int(df.S.nunique()))
rb = np.random.default_rng(ss[9]); rp = np.random.default_rng(ss[10])
for kind in ("clean", "raw"):
    o = dict(mean_O_F=round(float(df[f"O_F_{kind}"].mean()), 3), mean_control_A=round(float(df[f"meanA_{kind}"].mean()), 3), mean_control_B=round(float(df[f"meanB_{kind}"].mean()), 3),
             share_sets_F_ge1=round(float((df[f"O_F_{kind}"] >= 1).mean()), 3), share_control_A_ge1=round(float(df[f"anyA_{kind}"].mean()), 3))
    for stat in (("eA", "eB", "thetaA", "eC") if kind == "clean" else ("eA", "eB", "thetaA")):
        v = df[f"{stat}_{kind}"].values; ok = ~np.isnan(v); est = float(np.nanmean(v))
        csum = np.array([np.nansum(v[gi[c_]]) for c_ in cl]); cn = np.array([np.sum(ok[gi[c_]]) for c_ in cl])
        bs = []
        for _ in range(10000):
            pick = rb.integers(0, len(cl), len(cl)); bs.append(csum[pick].sum() / max(cn[pick].sum(), 1))
        lo, hi = np.percentile(bs, [2.5, 97.5]); o[stat] = dict(estimate=round(est, 4), ci95=[round(float(lo), 4), round(float(hi), 4)])
        if stat != "thetaA":
            ge = sum((rp.choice([-1, 1], len(cl)) * csum).sum() / cn.sum() >= est for _ in range(9999)); o[stat]["p_signflip_one_sided"] = (1 + ge) / 10000
    out[kind] = o
json.dump(out, open(os.path.join(PROJ, f"results/tables/asof_summary_{TAG}.json"), "w"), indent=1); print(json.dumps(out, indent=1))
