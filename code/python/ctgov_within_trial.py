#!/usr/bin/env python3
"""Key sensitivity for H4 (SAP_draft §8.3): within-trial, reference-arm-adjusted trial profiles.

Trials with both a monotherapy arm of the molecule and a reference arm (placebo / vehicle / sham / no treatment /
observation). Per molecule and term: Mantel-Haenszel risk ratio across those trials (exposed arms pooled within trial vs
reference arms pooled within trial); trial-signal if pooled exposed a >= 3, RR_MH >= 2 and the Mantel-Haenszel
chi-square (continuity-corrected) >= 4. Same term matching and cleaning as ctgov_profiles.py. Minimum data: >= 3 such
trials and >= 5 clean trial-signals. Writes data/derived/ctgov/profile_within_trial/."""
import csv, gzip, json, os, re, sys
from collections import defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", ""); sys.path.insert(0, os.path.join(ROOT, "scripts")); sys.path.insert(0, HERE)  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
from compute_faers_prr import BLACKLIST          # noqa: E402
RAW = os.path.join(PROJ, "data/raw"); CT = os.path.join(PROJ, "data/derived/ctgov"); OUT = os.path.join(CT, "profile_within_trial"); os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10 ** 9); K = 20
src = open(os.path.join(HERE, "ctgov_profiles.py")).read()
ns = {"__file__": os.path.join(HERE, "ctgov_profiles.py")}; exec(src.split("arm_rx = {}; arm_trial = {}")[0], ns)        # reuse VOCAB, std(), ARTEFACT, BL exactly as in the main profile script
std, ARTEFACT, BL = ns["std"], ns["ARTEFACT"], ns["BL"]
lab = {}; rxof = {}; trial = {}
for r in csv.DictReader(gzip.open(f"{CT}/arm_attribution.csv.gz", "rt", newline="")):
    if r["label"] in ("exposed", "reference"): lab[r["result_group_id"]] = r["label"]; rxof[r["result_group_id"]] = r["ingredient_rxcui"]; trial[r["result_group_id"]] = r["nct_id"]
at_risk = {r["result_group_id"]: int(r["subjects_at_risk"] or 0) for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C3_arm_at_risk.csv.gz", "rt", newline=""))}
ev = pd.read_csv(f"{RAW}/ctgov/C2_events_nonzero.csv.gz", usecols=["result_group_id", "adverse_event_term", "subjects_affected"], dtype={"result_group_id": str})
ev = ev[ev.result_group_id.isin(lab)]
tmap = {t: std(t) for t in ev.adverse_event_term.dropna().unique()}; ev["pt"] = ev.adverse_event_term.map(tmap); ev = ev[ev.pt.notna()]
g = ev.groupby(["result_group_id", "pt"], as_index=False).subjects_affected.sum(); g["n"] = g.result_group_id.map(at_risk); g = g[g.n > 0]; g["a"] = np.minimum(g.subjects_affected, g.n)
ref_n = defaultdict(float); ref_a = defaultdict(lambda: defaultdict(float)); exp_n = defaultdict(float); exp_a = defaultdict(lambda: defaultdict(float))
for arm, l in lab.items():
    n = at_risk.get(arm, 0)
    if n <= 0: continue
    if l == "reference": ref_n[trial[arm]] += n
    else: exp_n[(rxof[arm], trial[arm])] += n
for r in g.itertuples():
    if lab[r.result_group_id] == "reference": ref_a[trial[r.result_group_id]][r.pt] += r.a
    else: exp_a[(rxof[r.result_group_id], trial[r.result_group_id])][r.pt] += r.a
cond = defaultdict(set)
for f, col in (("C5_trial_conditions_mesh.csv.gz", "mesh_term"), ("C6_trial_conditions_text.csv.gz", "condition")):
    for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/{f}", "rt", newline="")): cond[r["nct_id"]].add(r[col].upper())
by_rx = defaultdict(list)
for (rx, t), n1 in exp_n.items():
    if ref_n.get(t, 0) > 0: by_rx[rx].append(t)
rows = []; ntr = {}; top_cond = {}
for rx, trials in by_rx.items():
    ntr[rx] = len(trials); own_cond = set().union(*(cond.get(t, set()) for t in trials)); cc = defaultdict(int)
    for t in trials:
        for c_ in cond.get(t, ()): cc[c_] += 1
    top_cond[rx] = {c_ for c_, _ in sorted(cc.items(), key=lambda x: (-x[1], x[0]))[:10]}
    terms = set().union(*(exp_a[(rx, t)].keys() for t in trials))
    for pt in terms:
        num = den = E = V = A = 0.0
        for t in trials:
            a, n1 = exp_a[(rx, t)].get(pt, 0.0), exp_n[(rx, t)]; c, n0 = ref_a[t].get(pt, 0.0), ref_n[t]; N = n1 + n0; m1 = a + c
            num += a * n0 / N; den += c * n1 / N; A += a; E += n1 * m1 / N
            if N > 1: V += n1 * n0 * m1 * (N - m1) / (N * N * (N - 1))
        if A < 3 or V <= 0: continue
        rr = num / den if den > 0 else float("inf"); chi = (abs(A - E) - 0.5) ** 2 / V
        if rr >= 2 and chi >= 4:
            rows.append((rx, pt, A, len(trials), rr if np.isfinite(rr) else 9999.0, round(chi, 2), int(pt not in BL and not ARTEFACT.search(pt) and pt not in own_cond)))
sig = pd.DataFrame(rows, columns=["ingredient_rxcui", "reaction_pt", "a", "trials", "rr_mh", "chi2_mh", "clean"])
ncl = sig[sig.clean == 1].groupby("ingredient_rxcui").size(); qual = [rx for rx in by_rx if ntr[rx] >= 3 and ncl.get(rx, 0) >= 5]
tgt = defaultdict(set)                                  # targets: reuse the cohort file logic through the main profile's D4 where available, else recompute
name2rx = {r["normalized_name"]: r["ingredient_rxcui"] for r in csv.DictReader(open(f"{RAW}/trialytic_db/R2_name_to_ingredient_map.csv")) if r["ingredient_rxcui"]}
saltwords = set(re.findall(r"'([a-z]+)'", open(os.path.join(PROJ, "private/sql/07_salt_fix.sql")).read().split("ARRAY[")[1].split("]")[0]))
rx_targets = defaultdict(set)
_tinfo = {r["target_chembl_id"]: r for r in csv.DictReader(open(f"{RAW}/trialytic_db/R4_target_info.csv"))}   # amendment 2 (A5): human protein targets only
def _target_ok(r):
    t = _tinfo.get(r["target_chembl_id"]); return bool(t) and t["organism"] == "Homo sapiens" and t["target_type"] in ("SINGLE PROTEIN", "PROTEIN FAMILY", "PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP")
for r in csv.DictReader(open(f"{RAW}/trialytic_db/R3_chembl_target_drugs.csv")):
    if not _target_ok(r): continue
    nm = r["normalized_name"]; rx = name2rx.get(nm)
    if rx is None:
        w = nm.split()
        if len(w) >= 2 and w[-1] in saltwords and " ".join(w[:-1]) not in saltwords: rx = name2rx.get(" ".join(w[:-1]))
    if rx: rx_targets[rx].add(r["target_symbol"])
for rx in qual:
    for t in rx_targets.get(rx, ()): tgt[t].add(rx)
ctargets = {t for t, m in tgt.items() if len(m) >= 3}; cohort = sorted({rx for t in ctargets for rx in tgt[t]})
same_pairs = len({tuple(sorted((x, y))) for t in ctargets for x in tgt[t] for y in tgt[t] if x < y})
with open(f"{OUT}/D2_cohort_molecules.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ingredient_rxcui", "name"]); w.writerows([rx, ""] for rx in cohort)
with open(f"{OUT}/D4_molecule_target.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ingredient_rxcui", "target_symbol"]); w.writerows([rx, t] for rx in cohort for t in sorted(rx_targets[rx]))
s2 = sig[sig.ingredient_rxcui.isin(set(cohort)) & (sig.clean == 1)].sort_values(["ingredient_rxcui", "rr_mh", "a", "reaction_pt"], ascending=[True, False, False, True])
sig[sig.ingredient_rxcui.isin(set(cohort))].to_csv(f"{OUT}/signals.csv.gz", index=False, compression="gzip")
s2[s2.groupby("ingredient_rxcui").cumcount() < K][["ingredient_rxcui", "reaction_pt"]].to_csv(f"{OUT}/profile_clean_k20.csv", index=False)
with open(f"{OUT}/shared_condition_pairs.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["rx_i", "rx_j"])
    for i, x in enumerate(cohort):
        for y in cohort[i + 1:]:
            if top_cond[x] & top_cond[y]: w.writerow([x, y])
summ = dict(trials_with_exposed_and_reference=len({t for ts in by_rx.values() for t in ts}), molecules_with_such_trials=len(by_rx), molecules_meeting_minimum_data=len(qual),
            cohort_molecules=len(cohort), cohort_targets=len(ctargets), same_target_pairs_within_cohort_targets=same_pairs,
            assessable=bool(len(cohort) >= 200 and len(ctargets) >= 50 and same_pairs >= 1000))
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1); print(json.dumps(summ))
