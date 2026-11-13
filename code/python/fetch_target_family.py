#!/usr/bin/env python3
"""Fetch ChEMBL protein classification for every cohort target (public ChEMBL REST API).
Writes data/derived/target_family.csv: target_symbol, target_chembl_id, protein_class_desc, family.
family = coarse class from the ChEMBL protein classification string (first matching rule)."""
import csv, json, os, urllib.request
from concurrent.futures import ThreadPoolExecutor
PROJ = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
API = "https://www.ebi.ac.uk/chembl/api/data"
def get(u): return json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"Accept": "application/json", "User-Agent": "research-script"}), timeout=60))
syms = {r["target_symbol"] for r in csv.DictReader(open(f"{PROJ}/data/derived/faers_q/b/D4_molecule_target.csv"))}
ids = {}
for r in csv.DictReader(open(f"{PROJ}/data/raw/trialytic_db/R4_target_info.csv")):
    if r["target_symbol"] in syms and r["organism"] == "Homo sapiens": ids.setdefault(r["target_symbol"], r["target_chembl_id"])
RULES = [("7tm", "GPCR"), ("kinase", "kinase"), ("ion channel", "ion channel"), ("nuclear receptor", "nuclear receptor"),
         ("transporter", "transporter"), ("protease", "protease"), ("cytokine", "cytokine/receptor"), ("enzyme", "other enzyme"),
         ("membrane receptor", "other membrane receptor"), ("surface antigen", "surface antigen"), ("secreted", "secreted protein"),
         ("structural", "structural protein"), ("transcription factor", "transcription factor"), ("epigenetic", "epigenetic regulator")]
def one(sym):
    try:
        t = get(f"{API}/target/{ids[sym]}.json"); descs = []
        for c in t.get("target_components", []):
            for p in get(f"{API}/target_component/{c['component_id']}.json").get("protein_classifications") or []:
                descs.append(get(f"{API}/protein_classification/{p['protein_classification_id']}.json")["protein_class_desc"])
        d = " | ".join(sorted(set(descs)))
        fam = next((f for k, f in RULES if k in d.lower()), "other" if d else "unclassified")
        return dict(target_symbol=sym, target_chembl_id=ids[sym], protein_class_desc=d, family=fam)
    except Exception as e:
        return dict(target_symbol=sym, target_chembl_id=ids.get(sym, ""), protein_class_desc=f"ERROR {e}", family="unclassified")
with ThreadPoolExecutor(6) as ex: rows = list(ex.map(one, sorted(syms)))
with open(f"{PROJ}/data/derived/target_family.csv", "w", newline="") as f:
    w = csv.DictWriter(f, rows[0].keys()); w.writeheader(); w.writerows(rows)
from collections import Counter; print(Counter(r["family"] for r in rows)); print("errors", sum("ERROR" in r["protein_class_desc"] for r in rows))
