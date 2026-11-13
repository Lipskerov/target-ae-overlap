#!/usr/bin/env python3
"""Post-hoc effect modification of the primary estimand by drug type (NOT in SAP v1 — exploratory).

Strata per cohort molecule (extraction b):
  modality : biologic (any BLA application, or FDA Purple Book proper-name match) vs small molecule
  route    : oral-only / injected-only / local-only (topical, ophthalmic, inhaled, nasal, otic, vaginal, rectal,
             transdermal, dental) / multiple route classes / unknown — from FDA NDC directory + Purple Book,
             matched by application number (R6) or by single-ingredient substance name
  atc1     : ATC level-1 anatomical group (molecules with >1 group -> 'multi')
  era      : first US approval < 2000 vs >= 2000

Estimand inside a stratum G: the paper's molecule-averaged Delta (same-target minus different-target mean shared
signals), computed only over pairs with BOTH molecules in G. CI: molecule bootstrap within G (B resamples).
Difference between two strata: difference of independent bootstrap draws.
Also: all-partner version (focal molecule in G, partners anywhere) and cross-modality same-target pairs.

Usage: effect_modification.py   (env B, default 2000). Needs ndc.csv / pb.csv: the public FDA NDC directory (application_number, route, dosage_form, substance_name) and
FDA Purple Book (bla_number, proper_name, route, dosage_form), in folder env REF_DIR. Writes results/tables/effect_modification_<label>.json + molecule_types.csv.
"""
import csv, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
CD = os.path.join(PROJ, "data/derived/faers_q/b"); RAW = os.path.join(PROJ, "data/raw/trialytic_db")
REF = os.environ["REF_DIR"]; B = int(os.environ.get("B", "2000"))
ss = np.random.SeedSequence(20261113)

mols = [r["ingredient_rxcui"] for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv"))]
idx = {m: i for i, m in enumerate(mols)}; n = len(mols)
names = {m: set() for m in mols}
for r in csv.DictReader(open(f"{CD}/D2_cohort_molecules.csv")):
    if r["name"]: names[r["ingredient_rxcui"]].add(r["name"].lower())
for r in csv.DictReader(open(f"{RAW}/R9_rxnorm_atc.csv")):
    if r["rxcui"] in names and r["tty"] == "IN": names[r["rxcui"]].add(r["str"].lower())

# ---- applications from the FDA crosswalk
apps = {m: set() for m in mols}; bla_flag = {m: False for m in mols}
for r in csv.DictReader(open(f"{RAW}/R6_fda_crosswalk_cohort.csv")):
    m = r["ingredient_rxcui"]
    if m not in apps: continue
    for col in ("application_numbers", "application_numbers_ref"):
        apps[m] |= set(re.findall(r"[A-Z]+\d+", r[col] or ""))
    if r["first_bla_approval"]: bla_flag[m] = True
    if r["ingredient_name"]: names[m].add(r["ingredient_name"].lower())

# ---- NDC + Purple Book
ndc_by_app, ndc_by_sub = {}, {}
for r in csv.DictReader(open(f"{REF}/ndc.csv")):
    rt = {x.strip() for x in (r["route"] or "").split(";") if x.strip()}
    ndc_by_app.setdefault(r["application_number"], set()).update(rt)
    sub = (r["substance_name"] or "").strip()
    if sub and ";" not in sub: ndc_by_sub.setdefault(sub, set()).update(rt)
pb_names, pb_by_bla = {}, {}
for r in csv.DictReader(open(f"{REF}/pb.csv")):
    rt = {x.strip().upper() for x in re.split(r"[;,]", r["route"] or "") if x.strip()}
    pb_by_bla.setdefault("BLA" + r["bla_number"], set()).update(rt)
    pb_names.setdefault(r["proper_name"].split("-")[0].strip() if re.search(r"-[a-z]{4}$", r["proper_name"]) else r["proper_name"], set()).update(rt)

def sub_match(name, sub):   # "pantoprazole sodium" matches "pantoprazole"
    return sub == name or sub.startswith(name + " ")

# INN stems: antibodies (-mab), fusion proteins (-cept), enzymes (-ase/-kinase/-plase), interleukins (-kin), cell and
# gene therapies (-cel, -gene, -vec), biosimilar 4-letter suffix, "alfa/beta" glycoproteins
BIO_RE = r"(mab|mab [a-z]+|cept|kinase|plase|ase alfa|kin|leucel|cel|gene [a-z]+|vec|-[a-z]{4})$|( alfa| beta| alfa-\w+| beta-\w+)$|^(insulin|interferon|peginterferon|factor |antithrombin|alpha-1|immune globulin)"
# synthetic peptides marketed under NDAs (stems -tide, -relin, -relix, -pressin, -reotide, -paratide)
PEP_RE = r"(tide|relin|relix|pressin|reotide|paratide|tocin)$"
ORAL = {"ORAL", "SUBLINGUAL", "BUCCAL"}
INJ = {"INTRAVENOUS", "SUBCUTANEOUS", "INTRAMUSCULAR", "INTRATHECAL", "INTRADERMAL", "INTRAVITREAL", "INTRA-ARTERIAL",
       "INTRAARTICULAR", "INTRA-ARTICULAR", "INTRAVESICAL", "INTRALESIONAL", "INTRACAVERNOUS", "EPIDURAL", "INTRAOCULAR",
       "INTRAPERITONEAL", "INTRAMUSCULAR; INTRAVENOUS", "PERINEURAL", "INFILTRATION", "INTRACAMERAL", "INTRATUMORAL",
       "INTRAVASCULAR", "INTRACARDIAC", "INTRAPLEURAL", "INTRACAVITARY", "SUBMUCOSAL", "INTRAOSSEOUS", "PERCUTANEOUS",
       "INTRAVENTRICULAR", "INTRASYNOVIAL", "SOFT TISSUE", "INTRAUTERINE", "SUBCONJUNCTIVAL", "RETROBULBAR",
       "INTRACORONARY", "INTRA-AMNIOTIC", "INTRASPINAL", "INTRAVENOUS DRIP", "INTRAVENOUS BOLUS", "PARENTERAL", "INJECTION"}
def route_class(rt):
    c = set()
    for x in rt:
        x = x.upper()
        if x in ORAL: c.add("oral")
        elif x in INJ or x.startswith("INTRA") or x in ("SUBCUTANEOUS", "INTRAVENOUS"): c.add("injected")
        elif x in ("", "NOT APPLICABLE"): continue
        else: c.add("local")
    return c

rows = []
mod = np.empty(n, object); rte = np.empty(n, object); atc1 = np.empty(n, object); era = np.empty(n, object)
atc = {}
for r in csv.DictReader(open(f"{RAW}/R9_rxnorm_atc.csv")):
    if r["rxcui"] in idx and len(r["code"]) >= 1: atc.setdefault(r["rxcui"], set()).add(r["code"][0])
yr = {}
for r in csv.DictReader(open(f"{RAW}/R11_fda_class_all.csv")):
    if r["ingredient_rxcui"] in idx and r["first_approval_year"]: yr[r["ingredient_rxcui"]] = int(float(r["first_approval_year"]))
for r in csv.DictReader(open(f"{RAW}/R10_drug_approval_all.csv")):
    if r["ingredient_rxcui"] in idx and r["earliest_orig_ap"]: yr[r["ingredient_rxcui"]] = int(r["earliest_orig_ap"][:4])

for m in mols:
    i = idx[m]; rt = set(); src = set()
    for a in apps[m]:
        if a in ndc_by_app: rt |= ndc_by_app[a]; src.add("ndc_app")
        if a in pb_by_bla: rt |= pb_by_bla[a]; src.add("pb_app")
    pb_hit = False
    for nm in names[m]:
        for sub, r_ in ndc_by_sub.items():
            if sub_match(nm, sub): rt |= r_; src.add("ndc_name")
        for pn, r_ in pb_names.items():
            if pn == nm or pn.startswith(nm + " "): rt |= r_; pb_hit = True; src.add("pb_name")
    is_bio = bla_flag[m] or any(a.startswith("BLA") for a in apps[m]) or pb_hit
    nm0 = " ".join(sorted(names[m]))
    if not is_bio and any(re.search(BIO_RE, x) for x in names[m]): is_bio = True; src.add("inn_suffix")
    is_pep = (not is_bio) and any(re.search(PEP_RE, x) for x in names[m]) and "defibrotide" not in names[m]  # oligonucleotide, not a peptide
    mod[i] = "biologic" if is_bio else ("peptide" if is_pep else "small_molecule")
    rc = route_class(rt)
    rte[i] = "unknown" if not rc else (next(iter(rc)) + "_only" if len(rc) == 1 else "multiple")
    a1 = atc.get(m, set()); atc1[i] = "none" if not a1 else (next(iter(a1)) if len(a1) == 1 else "multi")
    era[i] = "unknown" if m not in yr else ("pre2000" if yr[m] < 2000 else "2000plus")
    rows.append(dict(ingredient_rxcui=m, name=sorted(names[m])[0] if names[m] else "", modality=mod[i], route_class=rte[i],
                     routes=";".join(sorted(rt)), atc1=atc1[i], era=era[i], route_source=";".join(sorted(src))))
with open(os.path.join(PROJ, "results/tables/molecule_types.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, rows[0].keys()); w.writeheader(); w.writerows(rows)

# ---- target graph
S = np.zeros((n, n), np.float32)
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): pass
tg = {}
for r in csv.DictReader(open(f"{CD}/D4_molecule_target.csv")): tg.setdefault(r["target_symbol"], set()).add(idx[r["ingredient_rxcui"]])
for g in tg.values():
    g = np.array(sorted(g)); S[np.ix_(g, g)] = 1
np.fill_diagonal(S, 0); Dm = 1 - S; np.fill_diagonal(Dm, 0)

def overlap(prof):
    rt = {}
    for r in csv.DictReader(open(prof)):
        if r["ingredient_rxcui"] in idx: rt.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"].upper())
    terms = sorted({t for s in rt.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
    A = np.zeros((n, len(terms)), np.float32)
    for i, s in rt.items(): A[i, [tix[t] for t in s]] = 1
    O = A @ A.T; np.fill_diagonal(O, 0); return O

def delta_w(O, Smat, Dmat, w, focal):
    ns, nd = Smat @ w, Dmat @ w
    ok = focal & (w > 0) & (ns > 0) & (nd > 0)
    d = (Smat * O) @ w / np.where(ns > 0, ns, 1) - (Dmat * O) @ w / np.where(nd > 0, nd, 1)
    return float((w[ok] * d[ok]).sum() / w[ok].sum()) if ok.any() else np.nan, int(ok.sum())

def stratum(O, mask, rng, partners="within"):
    """partners='within': both molecules in stratum; 'all': focal in stratum, partners anywhere."""
    if partners == "within":
        P = np.outer(mask, mask).astype(np.float32); Sg, Dg = S * P, Dm * P
    else:
        Sg, Dg = S, Dm
    w0 = np.ones(n, np.float32)
    est, used = delta_w(O, Sg, Dg, w0, mask)
    gi = np.where(mask)[0]; boots = np.empty(B)
    for b in range(B):
        w = np.zeros(n, np.float32)
        if partners == "within":
            np.add.at(w, rng.choice(gi, len(gi)), 1)
        else:   # resample focal molecules from stratum, partners from whole cohort
            w = rng.multinomial(n, np.full(n, 1 / n)).astype(np.float32)
            wf = np.zeros(n, np.float32); np.add.at(wf, rng.choice(gi, len(gi)), 1)
            ns, nd = Sg @ w, Dg @ w; ok = (wf > 0) & (ns > 0) & (nd > 0)
            d = (Sg * O) @ w / np.where(ns > 0, ns, 1) - (Dg * O) @ w / np.where(nd > 0, nd, 1)
            boots[b] = (wf[ok] * d[ok]).sum() / wf[ok].sum() if ok.any() else np.nan; continue
        boots[b], _ = delta_w(O, Sg, Dg, w, mask)
    iu = np.triu_indices(n, 1); Pm = (np.outer(mask, mask) if partners == "within" else (mask[:, None] | mask[None, :]))[iu]
    s_u, d_u = (S[iu] > 0) & Pm, (Dm[iu] > 0) & Pm
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    ms, md = float(O[iu][s_u].mean()) if s_u.any() else np.nan, float(O[iu][d_u].mean())
    zero = float((O[iu][s_u] == 0).mean()) if s_u.any() else np.nan
    return dict(n_molecules=int(mask.sum()), molecules_used=used, same_target_pairs=int(s_u.sum()), different_target_pairs=int(d_u.sum()),
                delta=round(est, 3), ci95=[round(float(lo), 3), round(float(hi), 3)],
                pair_mean_same=round(ms, 3), pair_mean_different=round(md, 3),
                ratio_same_vs_diff=round(ms / md, 2) if md > 0 else None, pct_same_target_pairs_sharing_none=round(100 * zero, 1)), boots

def run(label, prof):
    O = overlap(prof); rng = np.random.default_rng(ss.spawn(1)[0]); out = dict(label=label, profile=os.path.relpath(prof, PROJ), B=B)
    for fac, arr, levels in (("modality", mod, ["small_molecule", "biologic", "peptide"]),
                             ("route", rte, ["oral_only", "injected_only", "local_only", "multiple", "unknown"]),
                             ("atc1", atc1, sorted(set(atc1))), ("era", era, ["pre2000", "2000plus", "unknown"])):
        out[fac] = {}; bs = {}
        for lv in levels:
            mask = arr == lv
            if mask.sum() < 10: continue
            res, b = stratum(O, mask, rng, "within")
            if res["same_target_pairs"] < 10: res["note"] = "fewer than 10 same-target pairs — not interpreted"
            out[fac][lv] = res; bs[lv] = b
        if fac == "modality":
            for lv in levels:
                out[fac][lv + "__all_partners"], _ = stratum(O, mod == lv, rng, "all")
        # pairwise differences vs reference (first level)
        ref = levels[0]; out[fac + "_differences_vs_" + ref] = {}
        for lv in bs:
            if lv == ref or ref not in bs: continue
            d = bs[lv] - bs[ref]; lo, hi = np.nanpercentile(d, [2.5, 97.5])
            p2 = 2 * min((d <= 0).mean(), (d >= 0).mean())
            out[fac + "_differences_vs_" + ref][lv] = dict(diff=round(out[fac][lv]["delta"] - out[fac][ref]["delta"], 3),
                                                           ci95=[round(float(lo), 3), round(float(hi), 3)], p_boot_two_sided=round(float(max(p2, 1 / B)), 4))
    # cross-modality same-target pairs
    iu = np.triu_indices(n, 1); bio = mod == "biologic"
    sm = mod == "small_molecule"
    for nm, sel in (("bio_bio", bio[:, None] & bio[None, :]), ("sm_sm", sm[:, None] & sm[None, :]), ("bio_sm", (bio[:, None] & sm[None, :]) | (sm[:, None] & bio[None, :]))):
        s = (S[iu] > 0) & sel[iu]; d = (Dm[iu] > 0) & sel[iu]
        out.setdefault("pair_types", {})[nm] = dict(same_target_pairs=int(s.sum()), mean_shared_same=round(float(O[iu][s].mean()), 3) if s.any() else None,
                                                   mean_shared_different=round(float(O[iu][d].mean()), 3))
    json.dump(out, open(os.path.join(PROJ, f"results/tables/effect_modification_{label}.json"), "w"), indent=1)
    return out

if __name__ == "__main__":
    from collections import Counter
    print("modality", Counter(mod)); print("route", Counter(rte)); print("atc1", Counter(atc1)); print("era", Counter(era))
    for label, prof in (("prr_clean", f"{CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{CD}/profile_clean2_ic025_k20.csv")):
        o = run(label, prof); print(label, "done", flush=True)
