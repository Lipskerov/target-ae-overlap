#!/usr/bin/env python3
"""Pre-specified sensitivity grid (SAP.md §2): ranking measure x K on extraction (b), plus extractions (c) and (d)
and the uncleaned profile. Builds each profile from <ext>/signals.csv.gz and runs primary_inference.py with
999 permutations and 1,000 bootstrap resamples (exploratory precision). Collects results/tables/sensitivity_grid.csv."""
import csv, json, os, subprocess, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
FQ = os.path.join(PROJ, "data/derived/faers_q"); TMP = os.path.join(FQ, "_grid"); os.makedirs(TMP, exist_ok=True)
runs = [("b", m, k, True) for m in ("prr", "ror025", "ic025") for k in (10, 20, 50)] + [("b", "prr", 20, False), ("c", "prr", 20, True), ("d", "prr", 20, True), ("c", "ic025", 20, True), ("d", "ic025", 20, True), ("a", "prr", 20, True), ("a", "ic025", 20, True), ("b", "ic025", 50, True)]
rows = []
for ext, measure, k, clean in runs:
    label = f"grid_{ext}_{measure}_k{k}_{'clean' if clean else 'raw'}"
    out = os.path.join(PROJ, "results/tables", f"primary_{label}.json")
    if not os.path.exists(out):
        s = pd.read_csv(os.path.join(FQ, ext, "signals.csv.gz"))
        if clean: s = s[s.clean]
        s = s.sort_values(["ingredient_rxcui", measure, "a", "reaction_pt"], ascending=[True, False, False, True])
        s = s[s.groupby("ingredient_rxcui").cumcount() < k]
        pf = os.path.join(TMP, label + ".csv"); s[["ingredient_rxcui", "reaction_pt"]].to_csv(pf, index=False)
        subprocess.run([sys.executable, os.path.join(HERE, "primary_inference.py"), pf, label, "999", "1000"],
                       env=dict(os.environ, COHORT_DIR=os.path.join(FQ, ext)), check=True, stdout=subprocess.DEVNULL)
    r = json.load(open(out)); x, j = r["shared_signals"], r["jaccard"]
    rows.append(dict(extraction=ext, ranking=measure, K=k, cleaned=int(clean), molecules=r["n_molecules"], same_pairs=r["same_target_pairs"],
                     pair_mean_same=x["pair_mean_same"], pair_mean_different=x["pair_mean_different"], delta=x["delta"], ci_lo=x["ci95"][0], ci_hi=x["ci95"][1],
                     p_perm=x["p_perm_one_sided"], jaccard_delta=j["delta"], jaccard_ci_lo=j["ci95"][0], jaccard_ci_hi=j["ci95"][1]))
    print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(PROJ, "results/tables/sensitivity_grid.csv"), index=False)
