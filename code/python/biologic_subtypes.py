#!/usr/bin/env python3
"""AMENDMENT_5 B-part: Delta within biologic subtypes (INN-stem rules fixed in the amendment); descriptive.
Writes results/tables/biologic_subtypes.json and lists the targets carrying biologic same-target pairs."""
import json, os, re, sys
from collections import Counter
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import effect_modification as em
T = os.path.join(em.PROJ, "results/tables")
mt = pd.read_csv(f"{T}/molecule_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui")
def subtype(name):
    n = (name or "").lower()
    if re.search(r"mab\b|mab ", n): return "monoclonal antibody"
    if n.endswith("cept"): return "fusion protein"
    if re.search(r"cel$|leucel|vec$|gene ", n): return "cell/gene therapy"
    if re.search(r"insulin|tropin|glucagon|parathyroid|follitropin|somatropin", n): return "insulin/hormone"
    if re.search(r"interferon|kin$|leukin|grastim|poetin", n): return "cytokine/interferon"
    if re.search(r"ase\b|ase alfa|factor|antithrombin|esterase|plasminogen", n): return "enzyme/blood factor"
    return "other biologic"
lab = np.array([subtype(mt.name.get(m, "")) if em.mod[em.idx[m]] == "biologic" else "not biologic" for m in em.mols])
out = dict(amendment="AMENDMENT_5 B-part", counts=Counter(lab[lab != "not biologic"]))
# targets carrying biologic-biologic same-target pairs
bio = np.where(em.mod == "biologic")[0]; tcount = Counter()
for t, g in em.tg.items():
    gb = [i for i in g if em.mod[i] == "biologic"]
    if len(gb) >= 2: tcount[t] = len(gb) * (len(gb) - 1) // 2
out["targets_by_biologic_pairs"] = tcount.most_common(15)
for label, prof in (("prr_clean", f"{em.CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{em.CD}/profile_clean2_ic025_k20.csv")):
    O = em.overlap(prof); rng = np.random.default_rng(np.random.SeedSequence(20261113).spawn(30)[21]); out[label] = {}
    for st, k in Counter(lab).most_common():
        if st == "not biologic": continue
        r, _ = em.stratum(O, lab == st, rng)
        if r["same_target_pairs"] < 10: r["note"] = "fewer than 10 same-target pairs - not interpreted"
        out[label][st] = r; print(label, f"{st:22s}", r["n_molecules"], r["same_target_pairs"], r["delta"], r["ci95"], r.get("note", ""), flush=True)
json.dump(out, open(f"{T}/biologic_subtypes.json", "w"), indent=1, default=dict)
print(out["counts"]); print(out["targets_by_biologic_pairs"])
