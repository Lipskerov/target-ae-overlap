#!/usr/bin/env python3
"""Amendment 3, WP5 (seed index 14): per-target contrasts and shared reactions. Post hoc; not in SAP v1.

Per-target contrast = pairs_all.py `target_level` definition, targets with >= 3 cohort members (extraction b):
  mean overlap among the target's member pairs  minus  mean over members of each member's mean overlap with its
  non-partners (molecules sharing no cohort target with it).
Reproduction check: mean over targets (PRR, standard cleaning) = 1.2494 (results/tables/target_vs_class_R2_prr_clean.json).
Shared reactions: among same-target pairs, number of pairs sharing each reaction; overall and within kinase / GPCR /
ion-channel / biologic / neoplasm pairs (both molecules in the group). Counts only.

Labels (read from files, never recomputed; effect_modification.py is NOT imported because its import rewrites
molecule_types.csv):
  target family   data/derived/target_family.csv (target row); molecule family = effect_modification_target.py rule
  modality        results/tables/molecule_types.csv
  area            results/tables/indication_types.csv, area_trials
Outputs: results/tables/target_level_R2_<label>.csv, results/tables/shared_terms_by_group_<label>.json
"""
import csv, hashlib, json, os, sys
from collections import Counter, defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
CD = os.path.join(PROJ, "data/derived/faers_q/b"); T = os.path.join(PROJ, "results/tables")
SEED = np.random.SeedSequence(20261113).spawn(15)[14]; NBOOT = 5000   # pairs_all.py target_level uses 5000 target resamples
REPRO = 1.2494

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
LABEL_FILES = [os.path.join(T, "molecule_types.csv"), os.path.join(T, "indication_types.csv"), os.path.join(PROJ, "data/derived/target_family.csv")]
h0 = {p: sha(p) for p in LABEL_FILES}

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]; idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
tg = {}; mtg = defaultdict(set)
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")):
    tg.setdefault(r["target_symbol"], set()).add(idx[r["ingredient_rxcui"]]); mtg[r["ingredient_rxcui"]].add(r["target_symbol"])
S = np.zeros((n, n), bool)
for m in tg.values():
    m = np.array(sorted(m)); S[np.ix_(m, m)] = True
np.fill_diagonal(S, False)

fam = {r["target_symbol"]: r["family"] for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/target_family.csv")))}
famlab = {m: (next(iter({fam[t] for t in mtg[m]})) if len({fam[t] for t in mtg[m]}) == 1 else "multi-family") for m in mols}  # effect_modification_target.py rule
modal = {r["ingredient_rxcui"]: r["modality"] for r in csv.DictReader(open(os.path.join(T, "molecule_types.csv")))}
area = {r["ingredient_rxcui"]: r["area_trials"] for r in csv.DictReader(open(os.path.join(T, "indication_types.csv")))}
missing = dict(modality=sum(m not in modal for m in mols), area=sum(m not in area for m in mols), family=sum(m not in famlab for m in mols))
print("cohort molecules", n, "| targets", len(tg), "| labels missing", missing, flush=True)

def profile(path):
    rt = {}
    for r in csv.DictReader(open(path)):
        if r["ingredient_rxcui"] in idx: rt.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"])   # as pairs_all.py
    terms = sorted({t for s in rt.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
    A = np.zeros((n, len(terms)), np.float32)
    for i, s in rt.items(): A[i, [tix[t] for t in s]] = 1
    O = A @ A.T; np.fill_diagonal(O, 0); return A, O, terms

def ci(v): lo, hi = np.percentile(v, [2.5, 97.5]); return [round(float(lo), 4), round(float(hi), 4)]

def target_level(O):
    Df = (~S).astype(np.float32); np.fill_diagonal(Df, 0)
    dmean = (Df * O).sum(1) / Df.sum(1).clip(min=1)
    rows = []
    for t, m in tg.items():
        m = np.array(sorted(m))
        if len(m) < 3: continue
        sub = O[np.ix_(m, m)]; pm = sub[np.triu_indices(len(m), 1)].mean(); dm = dmean[m].mean()
        ms = [mols[i] for i in m]; ac = Counter(area[x] for x in ms); top = max(ac.values())
        rows.append(dict(target_symbol=t, family=fam.get(t, "NA"), n_members=len(m), n_biologic=sum(modal[x] == "biologic" for x in ms),
                         n_pairs=len(m) * (len(m) - 1) // 2, mean_overlap_partners=round(float(pm), 4), mean_overlap_nonpartners=round(float(dm), 4),
                         contrast=round(float(pm - dm), 4), dominant_area="|".join(sorted(a for a, k in ac.items() if k == top)),
                         dominant_area_n=top, members=";".join(sorted(ms))))
    return rows

def shared_terms(A, terms, groups):
    iu = np.triu_indices(n, 1); ii, jj = iu; st = S[iu]; out = {}
    for g, sel in groups.items():
        m = st & sel[ii] & sel[jj]; tt = (A[ii[m]] * A[jj[m]]).sum(0)
        top = sorted(((t, int(v)) for t, v in zip(terms, tt) if v > 0), key=lambda x: (-x[1], x[0]))[:25]
        out[g] = dict(same_target_pairs=int(m.sum()), molecules_in_group=int(sel.sum()), pairs_sharing_none=int(((A[ii[m]] * A[jj[m]]).sum(1) == 0).sum()) if m.any() else 0,
                      shared_slots=int(tt.sum()), top25=[[t, v] for t, v in top])
    return out

fl = np.array([famlab[m] for m in mols]); mo = np.array([modal[m] for m in mols]); ar = np.array([area[m] for m in mols])
GROUPS = {"overall": np.ones(n, bool), "kinase": fl == "kinase", "GPCR": fl == "GPCR", "ion_channel": fl == "ion channel",
          "biologic": mo == "biologic", "neoplasm": ar == "neoplasms"}

repro = None; summary_all = {}
for label, prof in (("prr_clean", f"{CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{CD}/profile_clean2_ic025_k20.csv")):
    A, O, terms = profile(prof); rows = target_level(O); dt = np.array([r["contrast"] for r in rows])
    # unrounded contrasts, same target order as rows
    Df = (~S).astype(np.float32); np.fill_diagonal(Df, 0); dmean = (Df * O).sum(1) / Df.sum(1).clip(min=1)
    raw = np.array([O[np.ix_(np.array(sorted(m)), np.array(sorted(m)))][np.triu_indices(len(m), 1)].mean() - dmean[np.array(sorted(m))].mean() for m in tg.values() if len(m) >= 3])
    if label == "prr_clean":
        ok = abs(round(float(raw.mean()), 4) - REPRO) <= 0.0001
        repro = dict(expected=REPRO, observed=round(float(raw.mean()), 4), targets=int(len(raw)), status="PASS" if ok else "FAIL")
        print("REPRO CHECK:", repro, flush=True)
        if not ok: break
    rng = np.random.default_rng(SEED.spawn(2)[0 if label == "prr_clean" else 1])
    bs = [raw[rng.integers(0, len(raw), len(raw))].mean() for _ in range(NBOOT)]
    fam_summary = {}
    for f, k in Counter(r["family"] for r in rows).most_common():
        v = np.array([x for r, x in zip(rows, raw) if r["family"] == f], float)
        fam_summary[f] = dict(targets=int(k), mean_contrast=round(float(v.mean()), 4), median=round(float(np.median(v)), 4), share_positive=round(float((v > 0).mean()), 4))
    summ = dict(targets=int(len(raw)), mean_contrast=round(float(raw.mean()), 4), ci95_target_bootstrap=ci(bs), B=NBOOT,
                median=round(float(np.median(raw)), 4), share_targets_positive=round(float((raw > 0).mean()), 4), by_target_family_descriptive=fam_summary)
    rows.sort(key=lambda r: -r["contrast"])
    with open(os.path.join(T, f"target_level_R2_{label}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys()); w.writeheader(); w.writerows(rows)
    sh = shared_terms(A, terms, GROUPS)
    out = dict(wp="WP5", amendment="protocol/AMENDMENT_3.md", label=label, profile=os.path.relpath(prof, PROJ), seed="SeedSequence(20261113).spawn(15)[14]",
               note="post hoc; counts = number of same-target pairs in which both molecules carry the reaction in their profile; group = both molecules in group (molecule labels)",
               reproduction_check=repro, target_level=summ, shared_terms=sh)
    json.dump(out, open(os.path.join(T, f"shared_terms_by_group_{label}.json"), "w"), indent=1)
    print(f"\n== {label}: target-level", json.dumps({k: v for k, v in summ.items() if k != "by_target_family_descriptive"}))
    for f, d in fam_summary.items(): print(f"   family {f:22s} targets={d['targets']:3d} mean={d['mean_contrast']:.4f} median={d['median']:.4f} pos={d['share_positive']}")
    print("   top 10 targets:", [(r["target_symbol"], r["family"], r["n_members"], r["contrast"]) for r in rows[:10]])
    print("   bottom 5 targets:", [(r["target_symbol"], r["family"], r["n_members"], r["contrast"]) for r in rows[-5:]])
    for g, d in sh.items(): print(f"   shared[{g}] pairs={d['same_target_pairs']} none={d['pairs_sharing_none']} top5={d['top25'][:5]}", flush=True)

changed = [p for p in LABEL_FILES if sha(p) != h0[p]]
print("label files unchanged:", not changed, changed)
if repro is None or repro["status"] != "PASS": sys.exit(2)
