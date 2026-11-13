#!/usr/bin/env python3
"""ClinicalTrials.gov replication, step 2 (SAP_draft §8.2–8.4): molecule-level trial adverse-event profiles.

Arms: monotherapy arms from arm attribution v2 (audited 93/100). Terms: upper-cased, US->UK spelling, kept only if they
equal a FAERS preferred term (exact-match dictionary). Count per arm and term = serious + other subjects affected, capped
at the arm's subjects at risk. Per molecule: pooled incidence against the pooled incidence of all other molecules'
monotherapy arms, excluding every arm of a trial in which the molecule is linked. Trial-signal: a >= 3, RR >= 2, Yates
chi-square >= 4. Cleaning: FAERS blacklist + artefact pattern + the conditions of the molecule's own trials.
Profile: 20 trial-signals with the highest RR. Minimum data: >= 3 trials, >= 300 subjects at risk, >= 5 clean trial-signals.
Writes data/derived/ctgov/profile/ (D2, D4, profile_clean_k20.csv, signals.csv.gz, shared_condition_pairs.csv, summary.json)."""
import csv, gzip, json, os, re, sys
from collections import defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
ROOT = os.environ.get("TRIALYTIC_ROOT", ""); sys.path.insert(0, os.path.join(ROOT, "scripts"))  # Trialytic repository (https://trialytic.ai): proprietary modules, not included; see README
from compute_faers_prr import BLACKLIST          # noqa: E402
RAW = os.path.join(PROJ, "data/raw"); CT = os.path.join(PROJ, "data/derived/ctgov"); OUT = os.path.join(CT, "profile"); os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10 ** 9); K = 20
ARTEFACT = re.compile(open(os.path.join(PROJ, "private/sql/04_clean.sql")).read().split("reaction_pt ~ '")[1].split("';")[0])
BL = {b.upper() for b in BLACKLIST}
# FAERS preferred-term vocabulary (exact-match dictionary, no MedDRA licence needed)
vf = os.path.join(PROJ, "data/derived/faers_q/pt_vocabulary.csv")
if not os.path.exists(vf):
    import duckdb
    con = duckdb.connect(os.path.expanduser("~/faers_quarterly/faers.duckdb"), read_only=True)
    con.execute("SET threads=2; SET memory_limit='2GB'")
    con.execute(f"COPY (SELECT pt, count(*) AS reports FROM reac GROUP BY 1 ORDER BY 1) TO '{vf}' (HEADER, DELIMITER ',')"); con.close()
VOCAB = {r["pt"] for r in csv.DictReader(open(vf))}
US_UK = [("ANEMIA", "ANAEMIA"), ("DIARRHEA", "DIARRHOEA"), ("EDEMA", "OEDEMA"), ("HEMORRHAG", "HAEMORRHAG"), ("HEMATO", "HAEMATO"), ("HEMOGLOBIN", "HAEMOGLOBIN"),
         ("HEMOLY", "HAEMOLY"), ("HEMOPTYSIS", "HAEMOPTYSIS"), ("HEMATURIA", "HAEMATURIA"), ("HEMATEMESIS", "HAEMATEMESIS"), ("LEUKOPENIA", "LEUKOPENIA"),
         ("ESOPHAG", "OESOPHAG"), ("ESTROGEN", "OESTROGEN"), ("FETAL", "FOETAL"), ("PEDIATRIC", "PAEDIATRIC"), ("TUMOR", "TUMOUR"), ("COLOR", "COLOUR"),
         ("HYPERGLYCEMIA", "HYPERGLYCAEMIA"), ("HYPOGLYCEMIA", "HYPOGLYCAEMIA"), ("GLYCEMI", "GLYCAEMI"), ("EMIA", "AEMIA"), ("DYSPNEA", "DYSPNOEA"),
         ("APNEA", "APNOEA"), ("ORTHOPNEA", "ORTHOPNOEA"), ("TACHYPNEA", "TACHYPNOEA"), ("CELIAC", "COELIAC"), ("CECAL", "CAECAL"), ("FECAL", "FAECAL"), ("FECES", "FAECES"),
         ("ISCHEMI", "ISCHAEMI"), ("LEUKEMIA", "LEUKAEMIA"), ("SEPTICEMIA", "SEPTICAEMIA"), ("BACTEREMIA", "BACTERAEMIA"), ("VIREMIA", "VIRAEMIA"), ("UREMI", "URAEMI"),
         ("PYREXIA", "PYREXIA"), ("EDEMATOUS", "OEDEMATOUS"), ("ANESTHE", "ANAESTHE"), ("GYNECO", "GYNAECO"), ("ORTHOPEDIC", "ORTHOPAEDIC"), ("HEMORRHOID", "HAEMORRHOID")]
def std(t):
    t = re.sub(r"\s+", " ", t.upper().strip()).rstrip("*.,;: ")
    if t in VOCAB: return t
    u = t
    for a, b in US_UK:
        if a in u and b not in u: u = u.replace(a, b)
    if u in VOCAB: return u
    u2 = re.sub(r"AEMIA", "AEMIA", u)
    return u2 if u2 in VOCAB else None

arm_rx = {}; arm_trial = {}
for r in csv.DictReader(gzip.open(f"{CT}/arm_attribution.csv.gz", "rt", newline="")):
    if r["label"] == "exposed": arm_rx[r["result_group_id"]] = r["ingredient_rxcui"]; arm_trial[r["result_group_id"]] = r["nct_id"]
at_risk = {r["result_group_id"]: int(r["subjects_at_risk"] or 0) for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C3_arm_at_risk.csv.gz", "rt", newline=""))}
name2rx = {r["normalized_name"]: r["ingredient_rxcui"] for r in csv.DictReader(open(f"{RAW}/trialytic_db/R2_name_to_ingredient_map.csv")) if r["ingredient_rxcui"]}
rxname = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{PROJ}/data/derived/D1_molecules_all.csv"))}
for n_, rx in sorted(name2rx.items()): rxname.setdefault(rx, n_)
rx_trials_linked = defaultdict(set)
for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C4_drug_trial_bridge.csv.gz", "rt", newline="")):
    rx = name2rx.get(r["normalized_name"])
    if rx: rx_trials_linked[rx].add(r["nct_id"])
# events on exposed arms
ev = pd.read_csv(f"{RAW}/ctgov/C2_events_nonzero.csv.gz", usecols=["nct_id", "result_group_id", "event_type", "adverse_event_term", "subjects_affected"], dtype={"result_group_id": str})
ev = ev[ev.result_group_id.isin(arm_rx)]
rows_all, subj_all = len(ev), int(ev.subjects_affected.sum())
tmap = {t: std(t) for t in ev.adverse_event_term.dropna().unique()}
ev["pt"] = ev.adverse_event_term.map(tmap); ev = ev[ev.pt.notna()]
rows_kept, subj_kept = len(ev), int(ev.subjects_affected.sum())
g = ev.groupby(["result_group_id", "pt"], as_index=False).subjects_affected.sum()
g["n"] = g.result_group_id.map(at_risk); g = g[g.n > 0]; g["a"] = np.minimum(g.subjects_affected, g.n)
g["rx"] = g.result_group_id.map(arm_rx); g["nct"] = g.result_group_id.map(arm_trial)
arms = pd.DataFrame({"arm": list(arm_rx)}); arms["rx"] = arms.arm.map(arm_rx); arms["nct"] = arms.arm.map(arm_trial); arms["n"] = arms.arm.map(at_risk).fillna(0); arms = arms[arms.n > 0]
N_all = float(arms.n.sum()); A_all = g.groupby("pt").a.sum()
n_trial = arms.groupby("nct").n.sum(); a_trial = g.groupby(["nct", "pt"]).a.sum()
n_m = arms.groupby("rx").n.sum(); trials_m = arms.groupby("rx").nct.nunique(); a_m = g.groupby(["rx", "pt"]).a.sum()
cond = defaultdict(set)
for f, col in (("C5_trial_conditions_mesh.csv.gz", "mesh_term"), ("C6_trial_conditions_text.csv.gz", "condition")):
    for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/{f}", "rt", newline="")): cond[r["nct_id"]].add(r[col].upper())
a_trial_d = defaultdict(dict)
for (nct, pt), v in a_trial.items(): a_trial_d[nct][pt] = v
sig_rows = []; own_trials = arms.groupby("rx").nct.apply(set).to_dict(); top_cond = {}
for rx, nm in n_m.items():
    excl = set(rx_trials_linked.get(rx, ())) | own_trials[rx]            # background leaves out every trial in which the molecule appears
    n0 = N_all - float(n_trial.reindex(list(excl & set(n_trial.index))).sum())
    sub = defaultdict(float)
    for t in excl:
        for pt, v in a_trial_d.get(t, {}).items(): sub[pt] += v
    own_cond = set().union(*(cond.get(t, set()) for t in own_trials[rx]))
    cc = defaultdict(int)
    for t in own_trials[rx]:
        for c_ in cond.get(t, ()): cc[c_] += 1
    top_cond[rx] = {c_ for c_, _ in sorted(cc.items(), key=lambda x: (-x[1], x[0]))[:10]}
    if rx not in a_m.index.get_level_values(0): continue
    for pt, a in a_m.loc[rx].items():
        if a < 3: continue
        a0 = float(A_all[pt]) - sub.get(pt, 0.0)
        b, c, d = nm - a, a0, n0 - a0; zero_bg = 0
        if a0 <= 0 or b <= 0 or d <= 0:
            rr, chi2 = (float("inf") if a0 <= 0 else None), None
            if a0 > 0: continue
            # term never seen in the background: a background count of 0.5 is used; the RR >= 2 and chi-square >= 4 rule still applies
            a0 = 0.5; c, d = a0, n0 - a0; zero_bg = 1
        rr = (a / nm) / (a0 / n0); Nn = nm + n0
        e = [nm * (a + c) / Nn, nm * (b + d) / Nn, n0 * (a + c) / Nn, n0 * (b + d) / Nn]
        chi2 = sum((abs(o - x) - 0.5) ** 2 / x for o, x in zip([a, b, c, d], e) if x > 0)
        if rr >= 2 and chi2 >= 4:
            clean = pt not in BL and not ARTEFACT.search(pt) and pt not in own_cond
            sig_rows.append((rx, pt, a, nm, round(rr, 3), round(chi2, 2), int(clean), zero_bg))
sig = pd.DataFrame(sig_rows, columns=["ingredient_rxcui", "reaction_pt", "a", "n_at_risk", "rr", "chi2", "clean", "zero_background"])
nclean = sig[sig.clean == 1].groupby("ingredient_rxcui").size()
qual = [rx for rx in n_m.index if trials_m[rx] >= 3 and n_m[rx] >= 300 and nclean.get(rx, 0) >= 5]
# targets (same rule as FAERS, incl. salt-form recovery)
saltwords = set(re.findall(r"'([a-z]+)'", open(os.path.join(PROJ, "private/sql/07_salt_fix.sql")).read().split("ARRAY[")[1].split("]")[0]))
rx_targets = defaultdict(set)
_tinfo = {r["target_chembl_id"]: r for r in csv.DictReader(open(f"{RAW}/trialytic_db/R4_target_info.csv"))}   # amendment 2 (A5): human protein targets only
def _target_ok(r):
    t = _tinfo.get(r["target_chembl_id"]); return bool(t) and t["organism"] == "Homo sapiens" and t["target_type"] in ("SINGLE PROTEIN", "PROTEIN FAMILY", "PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP")
for r in csv.DictReader(open(f"{RAW}/trialytic_db/R3_chembl_target_drugs.csv")):
    if not _target_ok(r): continue
    nm_ = r["normalized_name"]; rx = name2rx.get(nm_)
    if rx is None:
        w = nm_.split()
        if len(w) >= 2 and w[-1] in saltwords and " ".join(w[:-1]) not in saltwords: rx = name2rx.get(" ".join(w[:-1]))
    if rx: rx_targets[rx].add(r["target_symbol"])
tg = defaultdict(set)
for rx in qual:
    for t in rx_targets.get(rx, ()): tg[t].add(rx)
ctargets = {t for t, m in tg.items() if len(m) >= 3}; cohort = sorted({rx for t in ctargets for rx in tg[t]})
cs = set(cohort); same_pairs = len({tuple(sorted((x, y))) for t in ctargets for x in tg[t] for y in tg[t] if x < y})
with open(f"{OUT}/D2_cohort_molecules.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ingredient_rxcui", "name"]); w.writerows([rx, rxname.get(rx, "")] for rx in cohort)
with open(f"{OUT}/D4_molecule_target.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ingredient_rxcui", "target_symbol"]); w.writerows([rx, t] for rx in cohort for t in sorted(rx_targets[rx]))
s2 = sig[sig.ingredient_rxcui.isin(cs)].sort_values(["ingredient_rxcui", "rr", "a", "reaction_pt"], ascending=[True, False, False, True])
s2.to_csv(f"{OUT}/signals.csv.gz", index=False, compression="gzip")
c2 = s2[s2.clean == 1]; c2 = c2[c2.groupby("ingredient_rxcui").cumcount() < K]
c2[["ingredient_rxcui", "reaction_pt"]].to_csv(f"{OUT}/profile_clean_k20.csv", index=False)
with open(f"{OUT}/shared_condition_pairs.csv", "w", newline="") as f:      # pairs sharing >= 1 of their ten most frequent trial conditions
    w = csv.writer(f); w.writerow(["rx_i", "rx_j"]); n_sc = 0
    for i, x in enumerate(cohort):
        for y in cohort[i + 1:]:
            if top_cond[x] & top_cond[y]: w.writerow([x, y]); n_sc += 1
with open(f"{OUT}/shared_trial_pairs.csv", "w", newline="") as f:           # pairs with >= 1 trial contributing arms to both molecules
    w = csv.writer(f); w.writerow(["rx_i", "rx_j"]); n_st = 0
    for i, x in enumerate(cohort):
        for y in cohort[i + 1:]:
            if own_trials[x] & own_trials[y]: w.writerow([x, y]); n_st += 1
with open(f"{OUT}/molecule_totals.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["ingredient_rxcui", "name", "trials", "subjects_at_risk", "clean_trial_signals", "in_cohort"])
    w.writerows([rx, rxname.get(rx, ""), int(trials_m[rx]), int(n_m[rx]), int(nclean.get(rx, 0)), int(rx in cs)] for rx in n_m.index)
summ = dict(exposed_arms=len(arm_rx), arms_with_at_risk=int(len(arms)), event_rows_on_exposed_arms=rows_all, event_rows_kept_after_term_match=rows_kept,
            pct_rows_kept=round(100 * rows_kept / rows_all, 1), pct_affected_subjects_kept=round(100 * subj_kept / subj_all, 1),
            distinct_raw_terms=len(tmap), distinct_terms_matched=len({v for v in tmap.values() if v}),
            molecules_with_arms=int(len(n_m)), molecules_meeting_minimum_data=len(qual), cohort_molecules=len(cohort), cohort_targets=len(ctargets),
            same_target_pairs_within_cohort_targets=same_pairs, pairs_sharing_a_top10_condition=n_sc, pairs_sharing_a_trial=n_st, cohort_signals_with_zero_background=int(sig[sig.ingredient_rxcui.isin(cs)].zero_background.sum()), profile_entries_with_zero_background=int(c2.zero_background.sum()), zero_background_entries_whose_term_is_in_another_profile=int((c2[c2.zero_background == 1].reaction_pt.map(c2.groupby("reaction_pt").ingredient_rxcui.nunique()) > 1).sum()), profile_entries=int(len(c2)),
            assessable=bool(len(cohort) >= 200 and len(ctargets) >= 50 and same_pairs >= 1000),
            median_clean_signals=float(nclean.reindex(cohort).median()) if cohort else None)
json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
