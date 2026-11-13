#!/usr/bin/env python3
"""Post-hoc (not in SAP v1): does the same-target effect differ by therapeutic area (indication type)?

Each cohort molecule gets a dominant therapeutic area = MeSH disease tree branch (C01 infections, C04 neoplasms, ...,
F03 mental disorders), defined two independent ways:
  trials : branch counted over the molecule's ClinicalTrials.gov trials (C4 drug-trial bridge x C5 MeSH conditions);
           a trial counts once per branch; >= 3 trials required; dominant = most frequent branch; tie -> 'mixed'
  faers  : branch counted over the molecule's ten most-reported FAERS indications (D8), matched to MeSH by exact name
           (>= 2 matched indications required)
Generic branches (C23 pathological conditions/signs, C26 wounds) are ignored.
Stratum estimand = paper's Delta with both molecules in the area (effect_modification.stratum); adjusted = model_full
within the area (effect_modification_adjusted.fit). Areas with < 30 same-target pairs are reported but not modelled.
Needs REF_DIR/mesh_tree.csv (MeSH terms with tree numbers, exported from AACT ctgov.mesh_terms) plus the files used by
effect_modification.py. Writes results/tables/effect_modification_indication_<label>.json and indication_types.csv.
"""
import csv, gzip, json, os, sys
from collections import Counter, defaultdict
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import effect_modification as em                      # cohort, target graph, stratum()
import effect_modification_adjusted as ea             # fit(), COLS  (module runs its own main only under __main__)

PROJ = em.PROJ; RAW = os.path.join(PROJ, "data/raw"); T = os.path.join(PROJ, "results/tables")
AREA = {"C01": "infections", "C04": "neoplasms", "C05": "musculoskeletal", "C06": "digestive", "C07": "stomatognathic",
        "C08": "respiratory", "C09": "ENT", "C10": "nervous system", "C11": "eye", "C12": "urogenital", "C13": "urogenital",
        "C14": "cardiovascular", "C15": "blood/lymphatic", "C16": "congenital", "C17": "skin", "C18": "metabolic/nutritional",
        "C19": "endocrine", "C20": "immune", "C25": "chemically induced", "F03": "mental disorders"}
tree = defaultdict(set)
for r in csv.DictReader(open(os.path.join(em.REF, "mesh_tree.csv"))):
    a = AREA.get(r["tree_number"][:3])
    if a: tree[r["downcase_mesh_term"].strip()].add(a)

# ---- trial-based
name2rx = {r["normalized_name"]: r["ingredient_rxcui"] for r in csv.DictReader(open(f"{RAW}/trialytic_db/R2_name_to_ingredient_map.csv"))}
trials = defaultdict(set)
for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C4_drug_trial_bridge.csv.gz", "rt")):
    rx = name2rx.get(r["normalized_name"])
    if rx in em.idx: trials[rx].add(r["nct_id"])
cond = defaultdict(set)
for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C5_trial_conditions_mesh.csv.gz", "rt")):
    cond[r["nct_id"]] |= tree.get(r["mesh_term"].lower().strip(), set())
def dominant(c, nmin, n):
    if n < nmin or not c: return "unknown"
    top = c.most_common(2)
    return "mixed" if len(top) > 1 and top[0][1] == top[1][1] else top[0][0]
lab_tr, lab_fa, rows = {}, {}, []
for m in em.mols:
    c = Counter(a for t in trials[m] for a in cond.get(t, ()))
    lab_tr[m] = dominant(c, 3, len(trials[m]))
# ---- FAERS-indication-based
ind = defaultdict(list)
for r in csv.DictReader(open(os.path.join(PROJ, "data/derived/D8_top_indications.csv"))): ind[r["ingredient_rxcui"]].append(r["indication"].lower().strip())
for m in em.mols:
    hits = [tree[i] for i in ind.get(m, []) if i in tree]
    lab_fa[m] = dominant(Counter(a for h in hits for a in h), 2, len(hits))
for m in em.mols:
    rows.append(dict(ingredient_rxcui=m, area_trials=lab_tr[m], area_faers=lab_fa[m], n_trials=len(trials[m])))
pd.DataFrame(rows).to_csv(f"{T}/indication_types.csv", index=False)
agree = [(lab_tr[m], lab_fa[m]) for m in em.mols if lab_tr[m] not in ("unknown", "mixed") and lab_fa[m] not in ("unknown", "mixed")]

def adjusted_within(label, labels, area):
    P = pd.read_csv(f"{T}/pairs_all_R2_{label}.csv.gz", dtype={"rx_i": str, "rx_j": str}); P["atc_missing"] = 1 - P.atc_both
    sel = (P.rx_i.map(labels) == area) & (P.rx_j.map(labels) == area)
    Q = P[sel.values].reset_index(drop=True); mols = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(mols)}
    Q["ii"], Q["jj"] = Q.rx_i.map(ix), Q.rx_j.map(ix)
    full, _, _ = ea.fit(Q, ea.COLS, mols, 3); return full["same_target"]

out = dict(note="post hoc; not in SAP v1", molecules=len(em.mols),
           counts_trials=Counter(lab_tr.values()), counts_faers=Counter(lab_fa.values()),
           agreement=dict(both_labelled=len(agree), same_area=sum(a == b for a, b in agree)))
print(json.dumps({k: out[k] for k in ("counts_trials", "counts_faers", "agreement")}, default=dict), flush=True)
for label, prof in (("prr_clean", f"{em.CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{em.CD}/profile_clean2_ic025_k20.csv")):
    O = em.overlap(prof); rng = np.random.default_rng(np.random.SeedSequence(20261113).spawn(3)[2]); res = {}
    for defn, labels in (("trials", lab_tr), ("faers", lab_fa)):
        arr = np.array([labels[m] for m in em.mols]); res[defn] = {}
        for area, k in Counter(arr).most_common():
            mask = arr == area
            if mask.sum() < 10: continue
            r, _ = em.stratum(O, mask, rng)
            if r["same_target_pairs"] >= 30 and area not in ("unknown", "mixed"):
                r["adjusted_same_target"] = adjusted_within(label, labels, area)
            elif r["same_target_pairs"] < 10: r["note"] = "fewer than 10 same-target pairs — not interpreted"
            res[defn][area] = r
            print(label, defn, f"{area:22s}", r["n_molecules"], r["same_target_pairs"], r["delta"], r["ci95"],
                  r.get("adjusted_same_target", {}).get("beta"), r.get("adjusted_same_target", {}).get("ci95"), flush=True)
    out[label] = res
json.dump(out, open(f"{T}/effect_modification_indication.json", "w"), indent=1, default=dict)
