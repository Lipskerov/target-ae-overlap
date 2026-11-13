#!/usr/bin/env python3
"""ClinicalTrials.gov replication, step 1 (SAP_draft §8.1): attribute result arms to molecules by text.

Tier 1  the arm title/description names exactly one ingredient linked to the trial (word-boundary alias match)
Tier 2  the trial is linked to exactly one ingredient in total; every non-reference arm that names no ingredient
        counts as exposed to it
reference    matches placebo|vehicle|sham|no treatment|observation and names no ingredient
combination  names two or more ingredients (excluded from monotherapy profiles)
unattributed everything else

Deviations from the draft, both stricter: (1) "ingredient" = any RxNorm ingredient in the frozen name map, not only cohort
molecules, so an arm that also names a non-cohort drug is a combination; (2) aliases shorter than 4 characters are not used
for text matching (they match dose units and abbreviations).
Looks at no adverse-event data. Writes data/derived/ctgov/arm_attribution.csv.gz, attribution_summary.json and the
200-arm audit sample private/reviews/P3_arm_audit_sample.csv (seed index 13)."""
import csv, gzip, json, os, re, random
from collections import defaultdict, Counter
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
RAW = os.path.join(PROJ, "data/raw"); OUT = os.path.join(PROJ, "data/derived/ctgov"); os.makedirs(OUT, exist_ok=True)
csv.field_size_limit(10 ** 9)
name2rx = {r["normalized_name"]: r["ingredient_rxcui"] for r in csv.DictReader(open(f"{RAW}/trialytic_db/R2_name_to_ingredient_map.csv")) if r["ingredient_rxcui"]}
aliases = defaultdict(set)
for n, rx in name2rx.items():
    if len(n) >= 4: aliases[rx].add(n)
rxname = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{PROJ}/data/derived/D1_molecules_all.csv"))}
trial_rx = defaultdict(set); unmapped_bridge = 0
for r in csv.DictReader(gzip.open(f"{RAW}/ctgov/C4_drug_trial_bridge.csv.gz", "rt", newline="")):
    rx = name2rx.get(r["normalized_name"])
    if rx: trial_rx[r["nct_id"]].add(rx)
    else: unmapped_bridge += 1
REF = re.compile(r"\b(placebo|vehicle|sham|no treatment|observation)\b", re.I)
arms = list(csv.DictReader(gzip.open(f"{RAW}/ctgov/C1_arms.csv.gz", "rt", newline="")))
rows = []; cnt = Counter(); rx_cache = {}
def rx_re(rx):
    if rx not in rx_cache:
        al = sorted(aliases.get(rx, ()), key=len, reverse=True)
        rx_cache[rx] = re.compile(r"(?<![a-z0-9])(" + "|".join(re.escape(a) for a in al) + r")(?![a-z0-9])") if al else None
    return rx_cache[rx]
for a in arms:
    linked = trial_rx.get(a["nct_id"], set())
    text = f'{a["title"]} {a["description"]}'.lower()
    found = sorted(rx for rx in linked if rx_re(rx) is not None and rx_re(rx).search(text))
    is_ref = bool(REF.search(text))
    if len(found) == 1: label, rx, tier = "exposed", found[0], 1
    elif len(found) >= 2: label, rx, tier = "combination", "|".join(found), 0
    elif is_ref: label, rx, tier = "reference", "", 0
    elif len(linked) == 1: label, rx, tier = "exposed", next(iter(linked)), 2
    else: label, rx, tier = "unattributed", "", 0
    # a tier-1 arm that also says placebo (e.g. "drug X + placebo") stays exposed; flagged for the audit
    cnt[(label, tier)] += 1
    rows.append(dict(result_group_id=a["result_group_id"], nct_id=a["nct_id"], label=label, tier=tier, ingredient_rxcui=rx,
                     name=rxname.get(rx, ""), mentions_reference_word=int(is_ref), title=a["title"]))
with gzip.open(f"{OUT}/arm_attribution.csv.gz", "wt", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
exp = [r for r in rows if r["label"] == "exposed"]
per_rx = Counter(r["ingredient_rxcui"] for r in exp); trials_rx = defaultdict(set)
for r in exp: trials_rx[r["ingredient_rxcui"]].add(r["nct_id"])
cohort = {r["ingredient_rxcui"] for r in csv.DictReader(open(f"{PROJ}/data/derived/D2_cohort_molecules.csv"))}
summ = dict(arms=len(arms), trials=len({a["nct_id"] for a in arms}), bridge_rows_unmapped=unmapped_bridge,
            labels={f"{k[0]}{'_tier'+str(k[1]) if k[1] else ''}": v for k, v in sorted(cnt.items())},
            molecules_with_exposed_arm=len(per_rx), molecules_with_ge3_trials=sum(len(v) >= 3 for v in trials_rx.values()),
            article_cohort_molecules_with_ge3_trials=sum(len(v) >= 3 for k, v in trials_rx.items() if k in cohort),
            trials_with_exposed_and_reference_arm=len({r["nct_id"] for r in exp} & {r["nct_id"] for r in rows if r["label"] == "reference"}))
json.dump(summ, open(f"{OUT}/attribution_summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
rng = random.Random(int(np.random.SeedSequence(20261112).spawn(16)[13].generate_state(1)[0]))
t1 = [r for r in exp if r["tier"] == 1]; t2 = [r for r in exp if r["tier"] == 2]
sample = rng.sample(t1, 100) + rng.sample(t2, 100)
desc = {a["result_group_id"]: a["description"] for a in arms}
with open(f"{PROJ}/private/reviews/P3_arm_audit_sample.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["result_group_id", "nct_id", "tier", "attributed_molecule", "arm_title", "arm_description", "auditor_verdict(correct/wrong/unclear)", "auditor_note"])
    for r in sample: w.writerow([r["result_group_id"], r["nct_id"], r["tier"], r["name"], r["title"], desc[r["result_group_id"]][:600], "", ""])
print("audit sample written: 100 tier-1 + 100 tier-2 arms")
