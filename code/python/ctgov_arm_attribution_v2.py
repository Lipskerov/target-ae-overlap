#!/usr/bin/env python3
"""Arm attribution v2 — protocol amendment after the v1 audit failed the 90% gate (tier 1: 73/100, tier 2: 24/98).

Changes, all tightening, decided from the auditors' error categories and without reading any adverse-event data:
  * tier 2 (fallback) is dropped;
  * an arm is "exposed, monotherapy" only if its title/description names exactly one mapped ingredient AND
      - contains no reference word (placebo, vehicle, sham, no treatment, observation, untreated, control, usual/standard care),
      - names no OTHER intervention of the same trial (registry intervention names and other-names, any type except
        placebo/behavioural labels), which catches code-named agents and regimens,
      - contains no combination marker (+, plus, combination, combined with, co-administered, followed by, then, regimen/
        chemotherapy/radiotherapy words) and no code-name token (letters + 3 or more digits) that is not an alias of the molecule.
Writes arm_attribution.csv.gz (v2), attribution_summary.json and a fresh 100-arm audit sample (seed index 17; the one
redraw the plan allows)."""
import csv, gzip, json, os, re, random
from collections import defaultdict, Counter
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(PROJ, "data/raw"); OUT = os.path.join(PROJ, "data/derived/ctgov")
csv.field_size_limit(10 ** 9)
name2rx = {r["normalized_name"]: r["ingredient_rxcui"] for r in csv.DictReader(open(f"{RAW}/trialytic_db/R2_name_to_ingredient_map.csv")) if r["ingredient_rxcui"]}
aliases = defaultdict(set)
for n, rx in name2rx.items():
    if len(n) >= 4: aliases[rx].add(n)
rxname = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{PROJ}/data/derived/D1_molecules_all.csv"))}
trial_rx = defaultdict(set)
for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C4_drug_trial_bridge.csv.gz", "rt", newline="")):
    rx = name2rx.get(r["normalized_name"])
    if rx: trial_rx[r["nct_id"]].add(rx)
DOSE = re.compile(r"\b\d+(\.\d+)?\s*(mg|mcg|µg|ug|g|ml|iu|units?|%|mg/kg|mg/m2)\b")
def norm(s): return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9\-\+ ]", " ", DOSE.sub(" ", s.lower()))).strip()
GENERIC = {"placebo", "drug", "tablet", "tablets", "capsule", "capsules", "injection", "infusion", "treatment", "therapy", "dose", "oral", "active",
           "standard of care", "standard care", "usual care", "control", "comparator", "study drug", "investigational product", "open label", "open-label"}
trial_iv = defaultdict(set)
for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C8_interventions.csv.gz", "rt", newline="")):
    if r["intervention_type"] in ("Behavioral", "Diagnostic Test"): continue
    n = norm(r["name"])
    if len(n) >= 4 and n not in GENERIC and "placebo" not in n: trial_iv[r["nct_id"]].add(n)
REF = re.compile(r"\b(placebo|vehicle|sham|no treatment|observation|untreated|control|usual care|standard care|standard of care|best supportive care)\b", re.I)
COMBO = re.compile(r"(\+|\bplus\b|\bcombination\b|\bcombined with\b|\bco-?administ|\bfollowed by\b|\bthen\b|\badd-?on\b|\bregimen\b|\bchemotherap|\bradiotherap|\bradiation\b|\bchemoradi|\bbackground\b|\bconcomitant\b)", re.I)
CODE = re.compile(r"(?<![a-z0-9])[a-z]{2,6}[- ]?\d{3,}[a-z]?(?![a-z0-9])")
rx_cache = {}
def rx_re(rx):
    if rx not in rx_cache:
        al = sorted(aliases.get(rx, ()), key=len, reverse=True)
        rx_cache[rx] = re.compile(r"(?<![a-z0-9])(" + "|".join(re.escape(a) for a in al) + r")(?![a-z0-9])") if al else None
    return rx_cache[rx]
arms = list(csv.DictReader(gzip.open(f"{RAW}/ctgov/C1_arms.csv.gz", "rt", newline="")))
rows = []; cnt = Counter(); why = Counter()
for a in arms:
    linked = trial_rx.get(a["nct_id"], set()); text = f'{a["title"]} {a["description"]}'.lower(); ntext = norm(text)
    found = sorted(rx for rx in linked if rx_re(rx) is not None and rx_re(rx).search(text))
    is_ref = bool(REF.search(text)); label, rx = "unattributed", ""
    if len(found) >= 2: label, rx = "combination", "|".join(found)
    elif len(found) == 1:
        m = found[0]; own = aliases[m]
        other_iv = [n for n in trial_iv.get(a["nct_id"], ()) if n in ntext and not rx_re(m).search(n) and not any(n in al for al in own)]
        codes = [c for c in CODE.findall(text) if not rx_re(m).search(c) and c.replace(" ", "").replace("-", "") not in {x.replace(" ", "").replace("-", "") for x in own}]
        if is_ref: why["reference word"] += 1; label = "excluded_ambiguous"
        elif other_iv: why["names another trial intervention"] += 1; label = "combination"
        elif COMBO.search(text): why["combination marker"] += 1; label = "combination"
        elif codes: why["code-name token"] += 1; label = "combination"
        else: label, rx = "exposed", m
    elif is_ref: label = "reference"
    cnt[label] += 1
    rows.append(dict(result_group_id=a["result_group_id"], nct_id=a["nct_id"], label=label, tier=1 if label == "exposed" else 0,
                     ingredient_rxcui=rx if label == "exposed" else "", name=rxname.get(rx, "") if label == "exposed" else "", title=a["title"]))
with gzip.open(f"{OUT}/arm_attribution.csv.gz", "wt", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
exp = [r for r in rows if r["label"] == "exposed"]; trials_rx = defaultdict(set)
for r in exp: trials_rx[r["ingredient_rxcui"]].add(r["nct_id"])
cohort = {r["ingredient_rxcui"] for r in csv.DictReader(open(f"{PROJ}/data/derived/D2_cohort_molecules.csv"))}
summ = dict(version=2, arms=len(arms), labels=dict(cnt), tier1_candidates_removed_by=dict(why), molecules_with_exposed_arm=len(trials_rx),
            molecules_with_ge3_trials=sum(len(v) >= 3 for v in trials_rx.values()),
            article_cohort_molecules_with_ge3_trials=sum(len(v) >= 3 for k, v in trials_rx.items() if k in cohort),
            trials_with_exposed_and_reference_arm=len({r["nct_id"] for r in exp} & {r["nct_id"] for r in rows if r["label"] == "reference"}))
json.dump(summ, open(f"{OUT}/attribution_summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
rng = random.Random(int(np.random.SeedSequence(20261112).spawn(32)[17].generate_state(1)[0]))
desc = {a["result_group_id"]: a["description"] for a in arms}
with open(f"{PROJ}/private/reviews/P3_arm_audit_v2_input.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["result_group_id", "nct_id", "tier", "attributed_molecule", "arm_title", "arm_description"])
    for r in rng.sample(exp, 100): w.writerow([r["result_group_id"], r["nct_id"], 1, r["name"], r["title"], desc[r["result_group_id"]][:600]])
