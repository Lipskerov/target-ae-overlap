#!/usr/bin/env python3
"""R10 step 3: disproportionality signals, profiles and cohort from the deduplicated quarterly-file database (SAP.md §1–3).

For each extraction (a any role, b primary+secondary suspect, c primary suspect only, d = b restricted to
health-professional reporters) writes under data/derived/faers_q/<ext>/:
  molecule_totals.csv      rx, name, reports
  signals.csv.gz           every signal of every eligible molecule: a, PRR, chi2, ROR lower limit, IC025, flags, PRR rank
  D2_cohort_molecules.csv, D4_molecule_target.csv   (same shapes as the article freeze, so primary_inference.py runs unchanged)
  profile_<clean|raw>_k20.csv                       ingredient_rxcui, reaction_pt
  summary.json
Usage: faers_q_profiles.py [ext ...]   default: b a c d
Optional env ASOF=YYYYMMDD restricts to cases first received before that date (used by the as-of test)."""
import csv, gzip, json, math, os, re, sys, time
import duckdb, numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", "")  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from compute_faers_prr import BLACKLIST          # noqa: E402
STORE = os.environ.get("FAERS_STORE", os.path.expanduser("~/faers_quarterly"))
OUTROOT = os.path.join(PROJ, "data", "derived", "faers_q")
RAW = os.path.join(PROJ, "data", "raw", "trialytic_db")
ARTEFACT = re.compile(open(os.path.join(PROJ, "private/sql/04_clean.sql")).read().split("reaction_pt ~ '")[1].split("';")[0])
UNK_IND = {"PRODUCT USED FOR UNKNOWN INDICATION", "DRUG USE FOR UNKNOWN INDICATION", "UNKNOWN", "INDICATION UNKNOWN"}
BL = {b.upper() for b in BLACKLIST}
# Amendment 2 (A4): extended cleaning, written after internal review. Disease-course, treatment-failure, resistance and drug-level terms.
ART2 = re.compile(r"(RESISTANCE|TREATMENT FAILURE|THERAPEUTIC RESPONSE|THERAPY NON-RESPONDER|DRUG LEVEL|VIRAL LOAD|VIROLOGIC|VIRAL MUTATION|MUTATION|RECURRENT|REFRACTORY|RELAPSE|RECURRENCE|RECURRENT|NEOPLASM|BLAST CRISIS|TRANSFORMATION|DISEASE COMPLICATION|CONCOMITANT DISEASE|DISEASE RECURRENCE|REMISSION|TRANSPLANT REJECTION|GRAFT)")
QUAL = {"RECURRENT", "REFRACTORY", "STAGE", "I", "II", "III", "IV", "METASTATIC", "ADVANCED", "CHRONIC", "ACUTE", "AGGRAVATED", "PROGRESSIVE", "UNSPECIFIED", "NOS", "DISEASE", "DISORDER", "TYPE", "1", "2"}
def stem(t): return frozenset(w for w in re.split(r"[^A-Z0-9]+", t) if w and w not in QUAL)
MIN_REPORTS, K = 100, 20
HCP = ("MD", "PH", "OT", "HP", "RN")
EXT = {"e": "role_rank <= 2 AND NOT fdc_only", "a": "TRUE", "b": "role_rank <= 2 AND NOT fdc_only", "c": "role_rank = 1 AND NOT fdc_only", "d": "role_rank <= 2 AND NOT fdc_only"}

# ---- names and targets (frozen Trialytic database exports, https://trialytic.ai) -------------------------------------------------
name2rx = {r["normalized_name"]: r["ingredient_rxcui"] for r in csv.DictReader(open(os.path.join(RAW, "R2_name_to_ingredient_map.csv"))) if r["ingredient_rxcui"]}
rxname = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/D1_molecules_all.csv")))}
for n, rx in sorted(name2rx.items()): rxname.setdefault(rx, n)
saltwords = set(re.findall(r"'([a-z]+)'", open(os.path.join(PROJ, "private/sql/07_salt_fix.sql")).read().split("ARRAY[")[1].split("]")[0]))
# Amendment 2 (A5): by default only human protein targets (single protein, family, complex); TARGETS=all keeps every ChEMBL target.
TARGETS = os.environ.get("TARGETS", "human")
tinfo = {r["target_chembl_id"]: r for r in csv.DictReader(open(os.path.join(RAW, "R4_target_info.csv")))}
def target_ok(r):
    if TARGETS == "all": return True
    t = tinfo.get(r["target_chembl_id"])
    return bool(t) and t["organism"] == "Homo sapiens" and t["target_type"] in ("SINGLE PROTEIN", "PROTEIN FAMILY", "PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP")
rx_targets = {}; salt_links = 0
# Amendment 2 (A1): the RxNorm dictionary also resolves ChEMBL names (salts, brands) to ingredients.
rxn1 = {}
for r in csv.DictReader(open(os.path.join(RAW, "R8_rxnorm_name_to_ingredient.csv"))): rxn1.setdefault(r["name"], set()).add(r["ingredient_rxcui"])
for n_, v_ in rxn1.items():
    if len(v_) == 1 and n_ not in name2rx: name2rx[n_] = next(iter(v_))
for r in csv.DictReader(open(os.path.join(RAW, "R3_chembl_target_drugs.csv"))):
    if not target_ok(r): continue
    nm = r["normalized_name"]; rx = name2rx.get(nm)
    if rx is None:                                   # salt-form recovery: "<base> <counter-ion>" where the base is a known molecule
        w = nm.split()
        if len(w) >= 2 and w[-1] in saltwords and " ".join(w[:-1]) not in saltwords:
            rx = name2rx.get(" ".join(w[:-1])); salt_links += rx is not None
    if rx: rx_targets.setdefault(rx, set()).add(r["target_symbol"])

con = duckdb.connect(os.path.join(STORE, "faers.duckdb"), read_only=True)
con.execute(f"SET threads=4; SET memory_limit='4GB'; SET temp_directory='{os.path.join(STORE, 'duck_tmp')}'; SET preserve_insertion_order=false")
ASOF = os.environ.get("ASOF")
t0 = time.time()
for ext in (sys.argv[1:] or ["b", "a", "c", "d"]):
    out = os.path.join(OUTROOT, (ext if TARGETS == "human" else f"{ext}_alltargets") if not ASOF else f"{ext}_asof_{ASOF}"); os.makedirs(out, exist_ok=True)
    case_filter = "TRUE" + (f" AND occp_cod IN {HCP}" if ext == "d" else "") + (f" AND init_dt < {int(ASOF)}" if ASOF else "")
    con.execute(f"CREATE OR REPLACE TEMP TABLE base AS SELECT pid FROM cases WHERE {case_filter} AND pid IN (SELECT pid FROM reac)")
    N = con.execute("SELECT count(*) FROM base").fetchone()[0]
    con.execute(f"CREATE OR REPLACE TEMP TABLE mr AS SELECT DISTINCT d.pid, d.rx FROM drug_rx d SEMI JOIN base b ON b.pid = d.pid WHERE {EXT[ext]}")
    if ext == "e":   # sole suspect: the molecule is the only mapped ingredient with a suspect role in the report
        con.execute("""CREATE OR REPLACE TEMP TABLE mr AS SELECT m.pid, m.rx FROM mr m
                       JOIN (SELECT pid FROM drug_rx WHERE role_rank <= 2 GROUP BY pid HAVING count(DISTINCT rx) = 1) s USING (pid)""")
    tot = dict(con.execute("SELECT rx, count(*) FROM mr GROUP BY 1").fetchall())
    elig = {rx for rx, n in tot.items() if n >= MIN_REPORTS}
    con.execute("CREATE OR REPLACE TEMP TABLE el AS SELECT unnest(?) AS rx", [sorted(elig)])
    ptn = con.execute("SELECT r.pt, count(*) n FROM reac r SEMI JOIN base b ON b.pid = r.pid GROUP BY 1").df().set_index("pt")["n"]
    # aggregated in 16 disjoint molecule buckets: the single-pass group-by exceeds memory for the any-role extraction
    cnt = pd.concat([con.execute(f"""SELECT m.rx, r.pt, count(*) AS a FROM mr m JOIN el USING (rx) JOIN reac r ON r.pid = m.pid
                         WHERE hash(m.rx) % 16 = {bk} GROUP BY 1, 2 HAVING count(*) >= 3""").df() for bk in range(16)], ignore_index=True)
    ind = con.execute("""SELECT i.rx, i.indi_pt, count(DISTINCT i.pid) n FROM indi i SEMI JOIN mr m ON m.pid = i.pid AND m.rx = i.rx
                         JOIN el USING (rx) GROUP BY 1, 2""").df()
    ind = ind[~ind.indi_pt.isin(UNK_IND)]
    ind["share"] = ind.n / ind.groupby("rx").n.transform("sum")
    own = set(map(tuple, ind[ind.share >= 0.01][["rx", "indi_pt"]].values))
    own2 = {}
    for rx_, t_ in ind[ind.share >= 0.001][["rx", "indi_pt"]].values: own2.setdefault(rx_, set()).add(stem(t_))
    a = cnt.a.values.astype(float); ab = cnt.rx.map(tot).values.astype(float); ac = cnt.pt.map(ptn).values.astype(float)
    b, c = ab - a, ac - a; d = N - ab - c
    ok = (b > 0) & (c > 0) & (d > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        prr = (a / ab) / (c / (c + d))
        e = np.stack([ab * ac / N, ab * (N - ac) / N, (N - ab) * ac / N, (N - ab) * (N - ac) / N])
        o = np.stack([a, b, c, d])
        chi2 = (((np.abs(o - e) - 0.5) ** 2) / e).sum(0)
        ror_lo = np.exp(np.log((a * d) / (b * c)) - 1.96 * np.sqrt(1 / a + 1 / b + 1 / c + 1 / d))
        ic025 = np.log2((a + 0.5) / (e[0] + 0.5)) - 3.3 * (a + 0.5) ** -0.5 - 2 * (a + 0.5) ** -1.5
    cnt["prr"], cnt["chi2"], cnt["ror025"], cnt["ic025"] = prr, chi2, ror_lo, ic025
    sig = cnt[ok & (np.round(prr, 2) >= 2) & (np.round(chi2, 1) >= 4) & ~cnt.pt.isin(BL)].copy()
    sig["artefact"] = sig.pt.map(lambda p: bool(ARTEFACT.search(p)))
    sig["own_indication"] = [(r, p) in own for r, p in zip(sig.rx, sig.pt)]
    sig["clean"] = ~(sig.artefact | sig.own_indication)
    sig["clean2"] = sig.clean & ~sig.pt.map(lambda p: bool(ART2.search(p))) & ~np.array([stem(p) in own2.get(r, ()) for r, p in zip(sig.rx, sig.pt)])
    sig = sig.sort_values(["rx", "prr", "a", "pt"], ascending=[True, False, False, True])
    sig["rank_raw"] = sig.groupby("rx").cumcount() + 1
    sig["rank_clean"] = np.where(sig.clean, sig.groupby("rx").clean.cumsum(), 0)
    # cohort: eligible molecules with >=1 signal at a target that has >=3 such molecules
    have = set(sig.rx)
    tg = {}
    for rx in have:
        for t in rx_targets.get(rx, ()): tg.setdefault(t, set()).add(rx)
    ctargets = {t for t, m in tg.items() if len(m) >= 3}
    cohort = sorted({rx for t in ctargets for rx in tg[t]})
    with open(os.path.join(out, "D2_cohort_molecules.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ingredient_rxcui", "name"]); w.writerows([rx, rxname.get(rx, "")] for rx in cohort)
    with open(os.path.join(out, "D4_molecule_target.csv"), "w", newline="") as f:   # ALL targets of cohort molecules (any-shared-target rule)
        w = csv.writer(f); w.writerow(["ingredient_rxcui", "target_symbol"])
        w.writerows([rx, t] for rx in cohort for t in sorted(rx_targets[rx]))
    with open(os.path.join(out, "molecule_totals.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ingredient_rxcui", "name", "reports", "in_cohort"])
        w.writerows([rx, rxname.get(rx, ""), n, int(rx in set(cohort))] for rx, n in sorted(tot.items(), key=lambda x: -x[1]) if rx in elig)
    cs = sig[sig.rx.isin(cohort)]
    cs.rename(columns={"rx": "ingredient_rxcui", "pt": "reaction_pt"}).to_csv(os.path.join(out, "signals.csv.gz"), index=False, float_format="%.4f", compression="gzip")
    cs[cs.rank_raw <= K][["rx", "pt"]].rename(columns={"rx": "ingredient_rxcui", "pt": "reaction_pt"}).to_csv(os.path.join(out, f"profile_raw_k{K}.csv"), index=False)
    cs[cs.clean & (cs.rank_clean <= K)][["rx", "pt"]].rename(columns={"rx": "ingredient_rxcui", "pt": "reaction_pt"}).to_csv(os.path.join(out, f"profile_clean_k{K}.csv"), index=False)
    for flag in ("clean", "clean2"):
        for meas in ("prr", "ic025", "ror025"):
            t_ = cs[cs[flag]].sort_values(["rx", meas, "a", "pt"], ascending=[True, False, False, True])
            t_ = t_[t_.groupby("rx").cumcount() < K]
            t_[["rx", "pt", "a", meas]].rename(columns={"rx": "ingredient_rxcui", "pt": "reaction_pt"}).to_csv(os.path.join(out, f"profile_{flag}_{meas}_k{K}.csv"), index=False)
    nsig = cs[cs.clean].groupby("rx").size().reindex(cohort, fill_value=0)
    summ = dict(extraction=ext, asof=ASOF, targets=TARGETS, clean2_signals_cohort=int(cs.clean2.sum()), definition=EXT[ext] + (" ; health-professional reporters" if ext == "d" else ""),
                reports_in_background=int(N), distinct_pt=int(len(ptn)), ingredients_with_reports=len(tot), eligible_molecules=len(elig),
                molecules_with_signal=len(have), cohort_molecules=len(cohort), cohort_targets=len(ctargets), salt_links_recovered=int(salt_links),
                signals_cohort=int(len(cs)), clean_signals_cohort=int(cs.clean.sum()),
                removed_from_top20_pct=round(100 * (1 - cs[(cs.rank_raw <= K)].clean.mean()), 1) if len(cs) else None,
                median_clean_signals_per_molecule=float(nsig.median()) if len(cohort) else None,
                pct_molecules_under_20_clean_signals=round(100 * float((nsig < K).mean()), 1) if len(cohort) else None)
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(f"[{time.time()-t0:5.0f}s]", json.dumps(summ), flush=True)
