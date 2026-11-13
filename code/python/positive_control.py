#!/usr/bin/env python3
"""Positive control for the trial analysis (SAP_draft §8.5): does a molecule's FAERS profile overlap its OWN trial profile
more than other molecules' trial profiles? If not, term mapping between the sources has failed.
Usage: positive_control.py <faers_profile.csv> <label>"""
import csv, json, os, sys
import numpy as np
P = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")); fp, label = sys.argv[1], sys.argv[2]
def load(path):
    d = {}
    for r in csv.DictReader(open(path)): d.setdefault(r["ingredient_rxcui"], set()).add(r["reaction_pt"])
    return d
F = load(fp); T = load(os.path.join(P, "data/derived/ctgov/profile/profile_clean_k20.csv")); both = sorted(set(F) & set(T))
own = np.array([len(F[m] & T[m]) for m in both]); oth = np.array([np.mean([len(F[m] & T[x]) for x in both if x != m]) for m in both])
rng = np.random.default_rng(np.random.SeedSequence(20261112).spawn(40)[25]); d = own - oth
bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(5000)]
out = dict(label=label, molecules_in_both=len(both), mean_own_overlap=round(float(own.mean()), 4), mean_other_overlap=round(float(oth.mean()), 4), difference=round(float(d.mean()), 4),
           ci95=[round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)], pct_molecules_own_gt_other=round(100 * float((own > oth).mean()), 1))
json.dump(out, open(os.path.join(P, f"results/tables/positive_control_{label}.json"), "w"), indent=1); print(json.dumps(out))
