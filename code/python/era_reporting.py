#!/usr/bin/env python3
"""AMENDMENT_5: is the biologic / newer-drug same-target excess a reporting artefact?

Rebuilds extraction-b profiles for the fixed 1,039-molecule cohort under reporting variants (R0 full = reproduction
check; R1 late >= approval+3 y; R2 FAERS era >= 2014-09-01; R3 health professionals; R4 consumers; R5 non-US) with the
signal/cleaning definitions of faers_q_profiles.py (its definitions block is executed unchanged), then estimates
overall Delta, Delta within biologics / small molecules, and the WP1 joint model (A4.1b design) on each variant, against a
like-for-like full-data reference on the same eligible molecules.
Env: FAERS_STORE, TRIALYTIC_ROOT (BLACKLIST), REF_DIR (effect_modification import), VARIANTS (default R0,R1,R2,R3,R4,R5).
Writes data/derived/faers_q/era/<variant>/profile_*.csv and results/tables/era_reporting.json"""
import csv, json, os, sys, time
from collections import defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))

# ---- definitions block of faers_q_profiles.py, executed unchanged (everything before the DuckDB connection)
src_path = os.path.join(HERE, "faers_q_profiles.py"); src = open(src_path).read()
NS = {"__file__": src_path, "__name__": "faers_q_profiles_defs"}
sys_argv = sys.argv; sys.argv = [src_path]
exec(compile(src[:src.index("con = duckdb.connect")], src_path, "exec"), NS)
sys.argv = sys_argv
duckdb, BL, ARTEFACT, ART2, stem, UNK_IND, K, MIN_REPORTS, EXT = (NS[k] for k in ("duckdb", "BL", "ARTEFACT", "ART2", "stem", "UNK_IND", "K", "MIN_REPORTS", "EXT"))
STORE = NS["STORE"]; HCP = NS["HCP"]

import effect_modification as em            # cohort order (mols), target graph S/Dm, stratum()
import effect_modification_adjusted as ea   # fit(), COLS
B = 1000
T = os.path.join(PROJ, "results/tables"); CD = em.CD; RAW = os.path.join(PROJ, "data/raw/trialytic_db")
mols = em.mols; n = len(mols); idx = em.idx

# approval year (R10 precedence over R11), as in modifier_joint_model.py
yr = {}
for r in csv.DictReader(open(f"{RAW}/R11_fda_class_all.csv")):
    if r["ingredient_rxcui"] in idx and r["first_approval_year"]: yr[r["ingredient_rxcui"]] = int(float(r["first_approval_year"]))
for r in csv.DictReader(open(f"{RAW}/R10_drug_approval_all.csv")):
    if r["ingredient_rxcui"] in idx and r["earliest_orig_ap"]: yr[r["ingredient_rxcui"]] = int(r["earliest_orig_ap"][:4])

VARIANTS = os.environ.get("VARIANTS", "R0,R1,R2,R3,R4,R5").split(",")
CASE = {"R0": "TRUE", "R1": "TRUE", "R2": "init_dt >= 20140901", "R3": f"occp_cod IN {HCP}", "R4": "occp_cod = 'CN'",
        "R5": "pid IN (SELECT pid FROM dx WHERE country <> '' AND country IS NOT NULL AND country <> 'US')"}
SEEDIDX = {"R0": 15, "R1": 16, "R2": 17, "R3": 18, "R4": 19, "R5": 20}

con = duckdb.connect(os.path.join(STORE, "faers.duckdb"), read_only=True)
con.execute(f"SET threads=4; SET memory_limit='4GB'; SET temp_directory='{os.path.join(STORE, 'duck_tmp')}'; SET preserve_insertion_order=false")
con.execute("CREATE OR REPLACE TEMP TABLE coh AS SELECT unnest(?) AS rx", [mols])
con.execute("CREATE OR REPLACE TEMP TABLE cut AS SELECT unnest(?) AS rx, unnest(?) AS cutoff",
            [list(yr), [(y + 3) * 10000 + 101 for y in yr.values()]])

def build(v):
    t0 = time.time(); out = os.path.join(PROJ, "data/derived/faers_q/era", v); os.makedirs(out, exist_ok=True)
    con.execute(f"CREATE OR REPLACE TEMP TABLE base AS SELECT pid FROM cases WHERE {CASE[v]} AND pid IN (SELECT pid FROM reac)")
    N = con.execute("SELECT count(*) FROM base").fetchone()[0]
    con.execute(f"CREATE OR REPLACE TEMP TABLE mr AS SELECT DISTINCT d.pid, d.rx FROM drug_rx d SEMI JOIN base b ON b.pid = d.pid WHERE {EXT['b']}")
    if v == "R1":   # per-molecule late window; molecules without an approval year keep all reports
        con.execute("""CREATE OR REPLACE TEMP TABLE mr AS SELECT m.pid, m.rx FROM mr m JOIN cases c USING (pid) LEFT JOIN cut USING (rx)
                       WHERE cut.cutoff IS NULL OR c.init_dt >= cut.cutoff""")
    tot = dict(con.execute("SELECT rx, count(*) FROM mr GROUP BY 1").fetchall())
    elig = sorted(rx for rx in mols if tot.get(rx, 0) >= MIN_REPORTS)
    con.execute("CREATE OR REPLACE TEMP TABLE el AS SELECT unnest(?) AS rx", [elig])
    ptn = con.execute("SELECT r.pt, count(*) n FROM reac r SEMI JOIN base b ON b.pid = r.pid GROUP BY 1").df().set_index("pt")["n"]
    cnt = pd.concat([con.execute(f"""SELECT m.rx, r.pt, count(*) AS a FROM mr m JOIN el USING (rx) JOIN reac r ON r.pid = m.pid
                         WHERE hash(m.rx) % 16 = {bk} GROUP BY 1, 2 HAVING count(*) >= 3""").df() for bk in range(16)], ignore_index=True)
    ind = con.execute("""SELECT i.rx, i.indi_pt, count(DISTINCT i.pid) n FROM indi i SEMI JOIN mr m ON m.pid = i.pid AND m.rx = i.rx
                         JOIN el USING (rx) GROUP BY 1, 2""").df()
    ind = ind[~ind.indi_pt.isin(UNK_IND)]; ind["share"] = ind.n / ind.groupby("rx").n.transform("sum")
    own = set(map(tuple, ind[ind.share >= 0.01][["rx", "indi_pt"]].values)); own2 = {}
    for rx_, t_ in ind[ind.share >= 0.001][["rx", "indi_pt"]].values: own2.setdefault(rx_, set()).add(stem(t_))
    a = cnt.a.values.astype(float); ab = cnt.rx.map(tot).values.astype(float); ac = cnt.pt.map(ptn).values.astype(float)
    b, c = ab - a, ac - a; d = N - ab - c; ok = (b > 0) & (c > 0) & (d > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        prr = (a / ab) / (c / (c + d))
        e = np.stack([ab * ac / N, ab * (N - ac) / N, (N - ab) * ac / N, (N - ab) * (N - ac) / N]); o = np.stack([a, b, c, d])
        chi2 = (((np.abs(o - e) - 0.5) ** 2) / e).sum(0)
        ic025 = np.log2((a + 0.5) / (e[0] + 0.5)) - 3.3 * (a + 0.5) ** -0.5 - 2 * (a + 0.5) ** -1.5
    cnt["prr"], cnt["chi2"], cnt["ic025"] = prr, chi2, ic025
    sig = cnt[ok & (np.round(prr, 2) >= 2) & (np.round(chi2, 1) >= 4) & ~cnt.pt.isin(BL)].copy()
    sig["artefact"] = sig.pt.map(lambda p: bool(ARTEFACT.search(p)))
    sig["own_indication"] = [(r, p) in own for r, p in zip(sig.rx, sig.pt)]
    sig["clean"] = ~(sig.artefact | sig.own_indication)
    sig["clean2"] = sig.clean & ~sig.pt.map(lambda p: bool(ART2.search(p))) & ~np.array([stem(p) in own2.get(r, ()) for r, p in zip(sig.rx, sig.pt)])
    prof = {}
    for flag, meas, lab in (("clean", "prr", "prr_clean"), ("clean2", "ic025", "ic025_clean2")):
        t_ = sig[sig[flag]].sort_values(["rx", meas, "a", "pt"], ascending=[True, False, False, True])
        t_ = t_[t_.groupby("rx").cumcount() < K]
        t_[["rx", "pt"]].rename(columns={"rx": "ingredient_rxcui", "pt": "reaction_pt"}).to_csv(f"{out}/profile_{lab}_k{K}.csv", index=False)
        prof[lab] = t_
    print(f"[{v}] N={N:,} eligible={len(elig)} signals={len(sig):,} {time.time()-t0:.0f}s", flush=True)
    return dict(background_reports=int(N), eligible_molecules=len(elig), eligible=elig, signals=int(len(sig))), out

def overlap_from(path):
    return em.overlap(path)

# ---- pair-level labels for the joint model (A4.1b design: kinase per A3 rule, no era_missing)
fam = {r["target_symbol"]: r["family"] for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/target_family.csv")))}
mtg = defaultdict(set)
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): mtg[r["ingredient_rxcui"]].add(r["target_symbol"])
famlab = {m: (next(iter({fam[t] for t in mtg[m]})) if len({fam[t] for t in mtg[m]}) == 1 else "multi-family") for m in mols}
area = pd.read_csv(f"{T}/indication_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").area_trials
lab = pd.DataFrame(dict(bio=[int(em.mod[idx[m]] == "biologic") for m in mols], kin=[int(famlab[m] == "kinase") for m in mols],
                        onc=[int(area.get(m, "unknown") == "neoplasms") for m in mols],
                        aunk=[int(area.get(m, "unknown") in ("unknown", "mixed")) for m in mols],
                        new=[int(m in yr and yr[m] >= 2000) for m in mols]), index=mols)
PAIRS = {}
def pairs(label):
    if label not in PAIRS:
        P = pd.read_csv(f"{T}/pairs_all_R2_{label}.csv.gz", dtype={"rx_i": str, "rx_j": str}); P["atc_missing"] = 1 - P.atc_both
        Li, Lj = lab.loc[P.rx_i].reset_index(drop=True), lab.loc[P.rx_j].reset_index(drop=True)
        for nm, c in (("bio_pair", "bio"), ("kinase_pair", "kin"), ("onc_pair", "onc"), ("new_pair", "new")):
            P[nm] = (Li[c].values & Lj[c].values).astype(int)
        P["area_missing"] = (Li.aunk.values | Lj.aunk.values).astype(int)
        for nm in ("bio_pair", "kinase_pair", "onc_pair", "new_pair"): P[f"same_target_x_{nm}"] = P.same_target * P[nm]
        P["ii"] = P.rx_i.map(idx); P["jj"] = P.rx_j.map(idx); PAIRS[label] = P
    return PAIRS[label]
JOINT = ea.COLS + ["bio_pair", "kinase_pair", "onc_pair", "new_pair", "area_missing"] + [f"same_target_x_{m}" for m in ("bio_pair", "kinase_pair", "onc_pair", "new_pair")]

def estimates(O, elig, seed_idx, label):
    ss = np.random.SeedSequence(20261113).spawn(30)[seed_idx]; rng = np.random.default_rng(ss.spawn(2)[0])
    em.B = B
    emask = np.zeros(n, bool); emask[[idx[m] for m in elig]] = True
    res = {}
    r, _ = em.stratum(O, emask, rng); res["overall"] = r
    bs = {}
    for mo in ("biologic", "small_molecule"):
        r, b = em.stratum(O, emask & (em.mod == mo), rng); res[mo] = r; bs[mo] = b
    dd = bs["biologic"] - bs["small_molecule"]; lo, hi = np.nanpercentile(dd, [2.5, 97.5])
    res["biologic_minus_small_molecule"] = dict(diff=round(res["biologic"]["delta"] - res["small_molecule"]["delta"], 3), ci95=[round(float(lo), 3), round(float(hi), 3)])
    P = pairs(label).copy(); P["overlap"] = O[P.ii.values, P.jj.values]
    sel = emask[P.ii.values] & emask[P.jj.values]; Q = P[sel].reset_index(drop=True)
    ms = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(ms)}; Q["ii"], Q["jj"] = Q.rx_i.map(ix), Q.rx_j.map(ix)
    ea.B = 500
    full, _, _ = ea.fit(Q, JOINT, ms, int(ss.generate_state(1)[0] % (2**31)))
    res["joint"] = {k: full[k] for k in ("same_target", "same_target_x_bio_pair", "same_target_x_new_pair", "same_target_x_kinase_pair", "same_target_x_onc_pair") if k in full}
    res["joint_pairs"] = int(len(Q)); res["joint_same_target_pairs"] = int(Q.same_target.sum())
    return res

out = dict(amendment="AMENDMENT_5", variants={})
stored = {"prr_clean": f"{CD}/profile_clean_prr_k20.csv", "ic025_clean2": f"{CD}/profile_clean2_ic025_k20.csv"}
for v in VARIANTS:
    info, vdir = build(v); vo = dict(info=dict((k, v_) for k, v_ in info.items() if k != "eligible"))
    if v == "R0":   # reproduction: rebuilt full-data profiles must equal the stored ones
        rep = {}
        for label, sp in stored.items():
            a_ = set(map(tuple, pd.read_csv(sp, dtype=str)[["ingredient_rxcui", "reaction_pt"]].values))
            b_ = set(map(tuple, pd.read_csv(f"{vdir}/profile_{label}_k{K}.csv", dtype=str).values))
            rep[label] = dict(stored_rows=len(a_), rebuilt_rows=len(b_), identical=a_ == b_, sym_diff=len(a_ ^ b_))
        vo["repro"] = rep; print("[R0] repro", rep, flush=True)
        if not all(r["identical"] for r in rep.values()):
            vo["note"] = "REPRO FAIL - variants not interpretable"; out["variants"][v] = vo
            json.dump(out, open(f"{T}/era_reporting.json", "w"), indent=1); sys.exit(1)
    for label in ("prr_clean", "ic025_clean2"):
        Ov = overlap_from(f"{vdir}/profile_{label}_k{K}.csv")
        vo[label] = estimates(Ov, info["eligible"], SEEDIDX[v], label)
        if v != "R0":   # like-for-like reference: full-data profiles on the same eligible molecules
            O0 = overlap_from(stored[label]); vo[label + "_reference_same_molecules"] = estimates(O0, info["eligible"], 15, label)
        j = vo[label]["joint"]; print(f"[{v}] {label} overall {vo[label]['overall']['delta']} | bio-sm {vo[label]['biologic_minus_small_molecule']} | x_bio {j.get('same_target_x_bio_pair')} | x_new {j.get('same_target_x_new_pair')}", flush=True)
        if v != "R0":
            jr = vo[label + "_reference_same_molecules"]["joint"]; print(f"[{v}] {label} REFERENCE x_bio {jr.get('same_target_x_bio_pair')} | x_new {jr.get('same_target_x_new_pair')}", flush=True)
    out["variants"][v] = vo
    json.dump(out, open(f"{T}/era_reporting.json", "w"), indent=1, default=str)
print("done", flush=True)
