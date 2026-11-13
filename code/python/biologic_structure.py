#!/usr/bin/env python3
"""Structural features of biologics vs the biologic same-target excess (post hoc; spec: private/structure/SPEC.md).

Features: ADC payload (INN), Fc isotype / origin / format (Thera-SAbDab + INN), target membrane-vs-secreted (UniProt keywords).
Contrasts: pair means of shared top-20 FAERS signals, different-target pairs sharing vs not sharing a feature, and
same-target pairs likewise; molecule bootstrap. H0 re-computes the biologic Delta with payload / CD3E annotations removed.

Usage: REF_DIR=<dir with ndc.csv, pb.csv> B=2000 python3 code/python/biologic_structure.py
Inputs (gitignored): data/derived/structure/dl_therasabdab/TheraSAbDab_SeqStruc_OnlineDownload.csv,
                     data/derived/structure/dl_uniprot/uniprot_bio_targets.tsv
Output: results/tables/biologic_structure.json (aggregates only) + data/derived/structure/biologic_features.csv
"""
import csv, gzip, hashlib, json, os, re, sys
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import effect_modification as em

PROJ = em.PROJ; n = em.n; mols = em.mols; idx = em.idx; B = em.B
SD = os.path.join(PROJ, "data/derived/structure")
THERA = f"{SD}/dl_therasabdab/TheraSAbDab_SeqStruc_OnlineDownload.csv"
UNIP = f"{SD}/dl_uniprot/uniprot_bio_targets.tsv"
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()

# ------------------------------------------------------------------ annotations
name = {r["ingredient_rxcui"]: r["name"] for r in csv.DictReader(open(f"{PROJ}/results/tables/molecule_types.csv"))}
tg_of = {m: set() for m in mols}
for r in csv.DictReader(open(f"{em.CD}/D4_molecule_target.csv")): tg_of[r["ingredient_rxcui"]].add(r["target_symbol"])
bio = em.mod == "biologic"; sm = em.mod == "small_molecule"

thera = {}
for r in csv.DictReader(open(THERA, encoding="utf-8-sig")):
    for nm in [r["Therapeutic"]] + r["Alternative Therapeutic Names"].split(";"):
        nm = nm.strip().lower()
        if nm: thera.setdefault(nm, r)

PAYLOAD = {"vedotin": "MMAE", "mafodotin": "MMAF", "emtansine": "DM1", "soravtansine": "DM4", "deruxtecan": "DXd",
           "tesirine": "PBD", "govitecan": "SN-38", "ozogamicin": "calicheamicin"}
PCLASS = {"MMAE": "microtubule", "MMAF": "microtubule", "DM1": "microtubule", "DM4": "microtubule",
          "DXd": "topoisomerase-I", "SN-38": "topoisomerase-I", "PBD": "DNA", "calicheamicin": "DNA"}
PAYLOAD_TARGETS = {"TUBB4B", "TOP1"}
ORIG = {"Murine": "murine", "Chimeric": "chimeric", "Humanised": "humanized", "Genetically human": "human"}

def inn_origin(nm):
    w = nm.split()[0].split("-")[0]
    for stem, o in (("ximab", "chimeric"), ("zumab", "humanized"), ("umab", "human"), ("omab", "murine")):
        if w.endswith(stem): return o
    return None

feat = {}
for m in mols:
    if not bio[idx[m]]: continue
    nm = name[m].lower(); words = nm.split()
    th = next((thera[c] for c in (nm, nm.replace("ado-", ""), nm.replace("ado-", "").split()[0], nm.split("-")[0]) if c in thera), None)
    pay = next((PAYLOAD[w] for w in words[1:] if w in PAYLOAD), None)
    fmt_t = th["Format"] if th else ""; th_targets = th["Target"] if th else ""
    if pay or "ADC" in fmt_t: fmt = "ADC"
    elif ("Bispecific" in fmt_t or "T-Cell Engager" in fmt_t) and "CD3E" in th_targets: fmt = "TCE"
    elif "Bispecific" in fmt_t: fmt = "bispecific_other"
    elif fmt_t == "Whole mAb": fmt = "whole_mAb"
    elif fmt_t in ("Fab", "scFv", "Nanobody"): fmt = "fragment"
    elif "cabtagene" in nm or nm.endswith("leucel"): fmt = "CAR_T"
    elif words[0].endswith("cept"): fmt = "Fc_fusion"
    elif th: fmt = "other_antibody"
    else: fmt = "non_antibody_protein"
    iso = None
    if th:
        v = {x.strip() for x in th["CH1 Isotype"].split(";")}
        if len(v) == 1 and next(iter(v)) in ("G1", "G2", "G3", "G4"): iso = next(iter(v))
    org, org_src = None, None
    if th:
        v = {x.strip() for x in th["Genetics (Bispecifics delimited with semicolon)"].split(";")}
        if len(v) == 1 and next(iter(v)) in ORIG: org, org_src = ORIG[next(iter(v))], "thera"
    elif fmt not in ("CAR_T", "Fc_fusion", "non_antibody_protein") or inn_origin(nm):
        if inn_origin(nm): org, org_src = inn_origin(nm), "inn"
    payt = PAYLOAD_TARGETS | ({"EEF2"} if "diftitox" in nm else set())
    ant = tg_of[m] - payt if (pay or "diftitox" in nm) else set(tg_of[m])
    feat[m] = dict(ingredient_rxcui=m, name=nm, thera_match=th["Therapeutic"] if th else "", format=fmt, payload=pay or "",
                   payload_class=PCLASS.get(pay, "") if pay else "", isotype=iso or "", origin=org or "", origin_source=org_src or "",
                   targets=";".join(sorted(tg_of[m])), antigen_targets=";".join(sorted(ant)), _ant=ant, _pt=tg_of[m] - ant)

# UniProt target class
kw = {}
for r in csv.DictReader(open(UNIP), delimiter="\t"):
    for g in r["Gene Names (primary)"].split(";"):
        kw.setdefault(g.strip(), set()).update(x.strip() for x in r["Keyword ID"].split(";"))
def tclass(t):
    k = kw.get(t)
    if k is None: return None
    if "KW-0812" in k or "KW-0336" in k: return "membrane"
    if "KW-0964" in k: return "secreted"
    return "other"
for m, f in feat.items():
    cl = {tclass(t) for t in f["_ant"]}
    f["target_class"] = next(iter(cl)) if len(cl) == 1 and None not in cl else ""

ANTIBODY_FMTS = {"ADC", "TCE", "bispecific_other", "whole_mAb", "fragment", "other_antibody"}
FMT_SET = ANTIBODY_FMTS | {"CAR_T", "Fc_fusion"}
FC_UNCONJ = {"whole_mAb", "bispecific_other", "TCE"}
for m, f in feat.items():
    if f["format"] not in FC_UNCONJ or f["format"] == "TCE" and not f["isotype"]: f["isotype_used"] = ""
    else: f["isotype_used"] = f["isotype"]
    if f["format"] not in ANTIBODY_FMTS: f["origin_used"] = ""
    else: f["origin_used"] = f["origin"]

def arr(key):
    a = np.full(n, "", object)
    for m, f in feat.items(): a[idx[m]] = f[key]
    return a
FMT, PAY, PCL, ISO, ORG, TCL = (arr(k) for k in ("format", "payload", "payload_class", "isotype_used", "origin_used", "target_class"))

# antigen-only same-target matrices
def build_S(target_sets):
    tg = {}
    for i, s in target_sets.items():
        for t in s: tg.setdefault(t, set()).add(i)
    S_ = np.zeros((n, n), np.float32)
    for g in tg.values():
        g = np.array(sorted(g)); S_[np.ix_(g, g)] = 1
    np.fill_diagonal(S_, 0); return S_
ts_ant = {idx[m]: (feat[m]["_ant"] if m in feat else tg_of[m]) for m in mols}
ts_nocd3 = {i: (s - {"CD3E"} if mols[i] in feat and feat[mols[i]]["format"] == "TCE" else s) for i, s in ts_ant.items()}
S_any = em.S; S_ant = build_S(ts_ant); S_nc = build_S(ts_nocd3)
assert np.array_equal(S_any, build_S({idx[m]: tg_of[m] for m in mols}))
def Dof(S_):
    D_ = 1 - S_; np.fill_diagonal(D_, 0); return D_

# shared indication per pair
SI = np.zeros((n, n), np.int8)
with gzip.open(f"{PROJ}/results/tables/pairs_all_R2_prr_clean.csv.gz", "rt") as fh:
    for r in csv.DictReader(fh):
        if r["shared_indication"] == "1" and r["rx_i"] in idx and r["rx_j"] in idx:
            a, b = idx[r["rx_i"]], idx[r["rx_j"]]; SI[a, b] = SI[b, a] = 1

def load_A(prof):
    rt = {}
    for r in csv.DictReader(open(prof)):
        if r["ingredient_rxcui"] in idx: rt.setdefault(idx[r["ingredient_rxcui"]], set()).add(r["reaction_pt"].upper())
    terms = sorted({t for s in rt.values() for t in s}); tix = {t: k for k, t in enumerate(terms)}
    A = np.zeros((n, len(terms)), np.float32)
    for i, s in rt.items(): A[i, [tix[t] for t in s]] = 1
    O = A @ A.T; np.fill_diagonal(O, 0); return A, np.array(terms), O

# ------------------------------------------------------------------ estimators
def r3(x): return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 3)

class Contrast:
    """Pair-mean contrast g1 - g0 over pair list (ii, jj); molecule bootstrap over `universe` (bool mask over n)."""
    def __init__(s, O, A, terms, ii, jj, g1, g0, universe, rng, label1, label0):
        s.O = O; s.ii, s.jj = ii, jj; s.g1, s.g0 = g1, g0
        y = O[ii, jj]; s.res = {}
        m1 = y[g1].mean() if g1.any() else np.nan; m0 = y[g0].mean() if g0.any() else np.nan
        u = np.where(universe)[0]; bs = np.empty(B)
        for b in range(B):
            w = np.zeros(n); w[u] = rng.multinomial(len(u), np.full(len(u), 1 / len(u)))
            pw = w[ii] * w[jj]
            s1, s0 = pw[g1].sum(), pw[g0].sum()
            bs[b] = ((pw[g1] * y[g1]).sum() / s1 if s1 > 0 else np.nan) - ((pw[g0] * y[g0]).sum() / s0 if s0 > 0 else np.nan)
        lo, hi = np.nanpercentile(bs, [2.5, 97.5]) if np.isfinite(bs).any() else (np.nan, np.nan)
        small = int(g1.sum()) < 10 or int(g0.sum()) < 10
        s.res = {"group_1": label1, "group_0": label0, "pairs_1": int(g1.sum()), "pairs_0": int(g0.sum()),
                 "molecules_1": int(len(set(ii[g1]) | set(jj[g1]))), "molecules_0": int(len(set(ii[g0]) | set(jj[g0]))),
                 "mean_1": r3(m1), "mean_0": r3(m0), "diff": r3(m1 - m0), "ci95": [r3(lo), r3(hi)],
                 "share_boot_le0": r3(float(np.nanmean(bs <= 0)))}
        if small: s.res["note"] = "fewer than 10 pairs in a group - reported, not interpreted"
        if g1.any():
            c = (A[ii[g1]] * A[jj[g1]]).sum(0); top = np.argsort(-c, kind="stable")[:10]
            s.res["top_shared_terms_group_1"] = [[str(terms[k]), int(c[k])] for k in top if c[k] > 0]

def pairs_in(mask):
    ii, jj = np.triu_indices(n, 1); k = mask[ii] & mask[jj]; return ii[k], jj[k]

def delta_bio(O, S_, w, mask):
    P = np.outer(mask, mask).astype(np.float32); D_ = Dof(S_)
    return em.delta_w(O, S_ * P, D_ * P, w.astype(np.float32), mask)

# ------------------------------------------------------------------ run
def run(label, prof, rng):
    A, terms, O = load_A(prof); out = {"profile": os.path.relpath(prof, PROJ)}
    assert np.array_equal(O, em.overlap(prof))
    # H0 attribution
    gi = np.where(bio)[0]; mats = {"S_any": S_any, "S_antigen": S_ant, "S_antigen_noCD3": S_nc}
    est = {k: delta_bio(O, v, np.ones(n), bio)[0] for k, v in mats.items()}
    bs = {k: np.empty(B) for k in mats}
    for b in range(B):
        w = np.zeros(n); np.add.at(w, rng.choice(gi, len(gi)), 1)
        for k, v in mats.items(): bs[k][b] = delta_bio(O, v, w, bio)[0]
    ref = json.load(open(f"{PROJ}/results/tables/effect_modification_{label}.json"))["modality"]["biologic"]["delta"]
    iu = np.triu_indices(n, 1); bb = bio[iu[0]] & bio[iu[1]]
    h0 = {"reproduction_reference_delta": ref, "reproduction_ok": round(est["S_any"], 3) == ref}
    for k in mats:
        lo, hi = np.nanpercentile(bs[k], [2.5, 97.5]); same = (mats[k][iu] > 0) & bb
        h0[k] = {"delta": r3(est[k]), "ci95": [r3(lo), r3(hi)], "same_target_bio_pairs": int(same.sum()),
                 "pair_mean_same": r3(O[iu][same].mean()), "pair_mean_different": r3(O[iu][~same & bb].mean())}
        if k != "S_any":
            d = bs[k] - bs["S_any"]; lo, hi = np.nanpercentile(d, [2.5, 97.5])
            h0[k]["diff_vs_S_any"] = r3(est[k] - est["S_any"]); h0[k]["diff_ci95"] = [r3(lo), r3(hi)]
    lost = (S_any[iu] > 0) & (S_ant[iu] == 0) & bb; lost2 = (S_ant[iu] > 0) & (S_nc[iu] == 0) & bb
    h0["bio_pairs_same_only_via_payload_target"] = int(lost.sum()); h0["bio_pairs_same_only_via_CD3E_of_TCE"] = int(lost2.sum())
    h0["mean_shared_payload_only_pairs"] = r3(O[iu][lost].mean()) if lost.any() else None
    h0["mean_shared_CD3E_only_pairs"] = r3(O[iu][lost2].mean()) if lost2.any() else None
    out["H0_attribution"] = h0

    res = {}
    def add(key, ii, jj, g1, g0, universe, l1, l0, primary=True):
        c = Contrast(O, A, terms, ii, jj, g1, g0, universe, rng, l1, l0).res
        if primary:
            noind = SI[ii, jj] == 0
            c["no_shared_indication"] = {k: v for k, v in Contrast(O, A, terms, ii, jj, g1 & noind, g0 & noind, universe, rng, l1, l0).res.items()
                                         if not k.startswith("top_")}
        res[key] = c

    Sant = S_ant; Sa = S_any
    # H1a ADC-ADC, no shared antigen: same vs different payload class
    adc = FMT == "ADC"; ii, jj = pairs_in(adc); nd = Sant[ii, jj] == 0
    add("H1a_ADC_diff_antigen_same_vs_diff_payload_class", ii, jj, nd & (PCL[ii] == PCL[jj]), nd & (PCL[ii] != PCL[jj]), adc,
        "same payload class", "different payload class")
    # H1b ADC x small molecule
    uni = adc | sm; ii, jj = pairs_in(uni); x = (adc[ii] & sm[jj]) | (sm[ii] & adc[jj])
    a_i = np.where(adc[ii], ii, jj); s_i = np.where(adc[ii], jj, ii)
    hit = np.array([bool(xx) and bool(feat[mols[a]]["_pt"] & tg_of[mols[s]]) for a, s, xx in zip(a_i, s_i, x)])
    add("H1b_ADC_x_smallmol_payload_target_vs_no_shared_target", ii, jj, x & hit, x & (Sa[ii, jj] == 0), uni,
        "small molecule hits ADC payload target", "no shared target")
    # H2 isotype
    m2 = ISO != ""; ii, jj = pairs_in(m2); same = ISO[ii] == ISO[jj]; st = Sa[ii, jj] > 0
    add("H2a_isotype_difftarget_same_vs_diff", ii, jj, ~st & same, ~st & ~same, m2, "same isotype", "different isotype")
    add("H2b_isotype_sametarget_same_vs_diff", ii, jj, st & same, st & ~same, m2, "same isotype", "different isotype")
    res["H2_desc_difftarget_by_isotype"] = {g: dict(pairs=int(k.sum()), mean=r3(O[ii, jj][k].mean()) if k.any() else None)
                                           for g in ("G1", "G2", "G4") for k in [~st & (ISO[ii] == g) & (ISO[jj] == g)]}
    # H3 origin
    m3 = ORG != ""; ii, jj = pairs_in(m3); same = ORG[ii] == ORG[jj]; st = Sa[ii, jj] > 0
    add("H3a_origin_difftarget_same_vs_diff", ii, jj, ~st & same, ~st & ~same, m3, "same origin", "different origin")
    add("H3b_origin_sametarget_same_vs_diff", ii, jj, st & same, st & ~same, m3, "same origin", "different origin")
    bc = ~st & (ORG[ii] == "chimeric") & (ORG[jj] == "chimeric"); bh = ~st & (ORG[ii] == "human") & (ORG[jj] == "human")
    add("H3_desc_difftarget_both_chimeric_vs_both_human", ii, jj, bc, bh, m3, "both chimeric", "both human", primary=False)
    # H4 format
    m4 = np.isin(FMT, list(FMT_SET)); ii, jj = pairs_in(m4); same = FMT[ii] == FMT[jj]; st = Sa[ii, jj] > 0
    add("H4a_format_difftarget_same_vs_diff", ii, jj, ~st & same, ~st & ~same, m4, "same format", "different format")
    add("H4b_format_sametarget_same_vs_diff", ii, jj, st & same, st & ~same, m4, "same format", "different format")
    res["H4_desc_difftarget_same_format_by_format"] = {g: dict(pairs=int(k.sum()), mean=r3(O[ii, jj][k].mean()) if k.any() else None)
                                                     for g in sorted(FMT_SET) for k in [~st & (FMT[ii] == g) & (FMT[jj] == g)]}
    tset = lambda t: np.array([t in tg_of[mols[i]] for i in range(n)])
    cd3 = tset("CD3E"); ii, jj = pairs_in(cd3 & bio); tt = (FMT[ii] == "TCE") & (FMT[jj] == "TCE"); tm = ((FMT[ii] == "TCE") ^ (FMT[jj] == "TCE"))
    add("H4_desc_CD3E_pairs_TCE_TCE_vs_TCE_monospecific", ii, jj, tt, tm, cd3 & bio, "TCE-TCE", "TCE-anti-CD3 monospecific", primary=False)
    cd19 = tset("CD19"); ii, jj = pairs_in(cd19 & bio); cc = (FMT[ii] == "CAR_T") & (FMT[jj] == "CAR_T")
    ca = (FMT[ii] == "CAR_T") ^ (FMT[jj] == "CAR_T"); aa = (FMT[ii] != "CAR_T") & (FMT[jj] != "CAR_T")
    add("H4_desc_CD19_pairs_CART_CART_vs_CART_antibody", ii, jj, cc, ca, cd19 & bio, "CAR_T-CAR_T", "CAR_T-antibody", primary=False)
    add("H4_desc_CD19_pairs_antibody_antibody_vs_CART_antibody", ii, jj, aa, ca, cd19 & bio, "antibody-antibody", "CAR_T-antibody", primary=False)
    ii, jj = pairs_in(m4); st = Sa[ii, jj] > 0
    add("H4_desc_difftarget_CART_CART_vs_other_difftarget", ii, jj, ~st & (FMT[ii] == "CAR_T") & (FMT[jj] == "CAR_T"),
        ~st & ~((FMT[ii] == "CAR_T") & (FMT[jj] == "CAR_T")), m4, "CAR_T-CAR_T different target", "other different-target", primary=False)
    # H5 target class
    m5 = (TCL == "membrane") | (TCL == "secreted"); ii, jj = pairs_in(m5); st = Sa[ii, jj] > 0
    add("H5a_sametarget_membrane_vs_secreted", ii, jj, st & (TCL[ii] == "membrane"), st & (TCL[ii] == "secreted"), m5,
        "membrane target", "secreted target")
    add("H5b_difftarget_same_vs_diff_class", ii, jj, ~st & (TCL[ii] == TCL[jj]), ~st & (TCL[ii] != TCL[jj]), m5, "same class", "different class")
    for cl in ("membrane", "secreted"):
        k1 = st & (TCL[ii] == cl); k0 = ~st & (TCL[ii] == cl) & (TCL[jj] == cl)
        res[f"H5_desc_{cl}_same_minus_diff_within_class"] = dict(pairs_same=int(k1.sum()), pairs_diff=int(k0.sum()),
            mean_same=r3(O[ii, jj][k1].mean()) if k1.any() else None, mean_diff=r3(O[ii, jj][k0].mean()) if k0.any() else None)
    ab = np.isin(FMT, list(ANTIBODY_FMTS)); ii, jj = pairs_in(ab); st = Sa[ii, jj] > 0
    pd1, pl1 = tset("PDCD1"), tset("CD274"); ax = ~st & ((pd1[ii] & pl1[jj]) | (pl1[ii] & pd1[jj]))
    add("H5_desc_PD1_x_PDL1_vs_other_difftarget_antibody", ii, jj, ax, ~st & ~ax, ab, "PD-1 x PD-L1 (different target)",
        "other different-target antibody pairs", primary=False)
    out["contrasts"] = res
    return out

if __name__ == "__main__":
    for f in feat.values():
        f["isotype_used"] = f.get("isotype_used", "")
    rows = [{k: v for k, v in f.items() if not k.startswith("_")} for f in feat.values()]
    with open(f"{SD}/biologic_features.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, rows[0].keys()); w.writeheader(); w.writerows(rows)
    ann = {"n_biologics": len(feat), "thera_matched": sum(1 for f in feat.values() if f["thera_match"]),
           "format": dict(Counter(f["format"] for f in feat.values())), "payload": dict(Counter(f["payload"] for f in feat.values() if f["payload"])),
           "isotype_used": dict(Counter(f["isotype_used"] or "unknown/not used" for f in feat.values())),
           "origin_used": dict(Counter(f["origin_used"] or "unknown/not used" for f in feat.values())),
           "origin_source": dict(Counter(f["origin_source"] for f in feat.values() if f["origin_used"])),
           "target_class": dict(Counter(f["target_class"] or "unknown" for f in feat.values()))}
    agree = [(f["origin"], inn_origin(f["name"])) for f in feat.values() if f["origin_source"] == "thera" and inn_origin(f["name"])]
    ann["origin_qc_thera_vs_inn"] = {"n": len(agree), "agree": sum(a == b for a, b in agree),
                                     "disagreements": Counter(f"{a}|inn:{b}" for a, b in agree if a != b)}
    out = {"spec": "private/structure/SPEC.md", "spec_sha256": sha(f"{PROJ}/private/structure/SPEC.md"), "B": B,
           "seeds": "SeedSequence(20261113).spawn(40)[30] PRR, [31] IC025",
           "sources": {"therasabdab": {"url": "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv",
                                       "downloaded": "2026-11-13", "sha256": sha(THERA)},
                       "uniprot": {"url": "https://rest.uniprot.org/uniprotkb/search (query in dl_uniprot/query_url.txt)",
                                   "downloaded": "2026-11-13", "sha256": sha(UNIP)}},
           "annotation": ann}
    ch = np.random.SeedSequence(20261113).spawn(40)
    for k, (label, prof) in enumerate((("prr_clean", f"{em.CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{em.CD}/profile_clean2_ic025_k20.csv"))):
        out[label] = run(label, prof, np.random.default_rng(ch[30 + k])); print(label, "done", flush=True)
    # POST-SPEC ADDENDUM (added after seeing results): biologic Delta (S_any) excluding cell-engaging formats.
    out["post_spec_addendum"] = {"note": "added after results were seen; not in SPEC.md", "seeds": "spawn(40)[32] PRR, [33] IC025"}
    for k, (label, prof) in enumerate((("prr_clean", f"{em.CD}/profile_clean_prr_k20.csv"), ("ic025_clean2", f"{em.CD}/profile_clean2_ic025_k20.csv"))):
        O = em.overlap(prof); rng = np.random.default_rng(ch[32 + k]); add_ = {}
        for nm_, excl in (("biologic_all", []), ("excl_CAR_T_TCE", ["CAR_T", "TCE"]), ("excl_CAR_T_TCE_ADC", ["CAR_T", "TCE", "ADC"])):
            mask = bio & ~np.isin(FMT, excl); gi = np.where(mask)[0]; e = delta_bio(O, S_any, np.ones(n), mask)[0]; bs = np.empty(B)
            for b in range(B):
                w = np.zeros(n); np.add.at(w, rng.choice(gi, len(gi)), 1); bs[b] = delta_bio(O, S_any, w, mask)[0]
            iu = np.triu_indices(n, 1); pm = mask[iu[0]] & mask[iu[1]]
            add_[nm_] = {"n_molecules": int(mask.sum()), "same_target_pairs": int(((S_any[iu] > 0) & pm).sum()), "delta": r3(e),
                         "ci95": [r3(x) for x in np.nanpercentile(bs, [2.5, 97.5])]}
        out["post_spec_addendum"][label] = add_; print("addendum", label, add_, flush=True)
    json.dump(out, open(f"{PROJ}/results/tables/biologic_structure.json", "w"), indent=1, default=lambda o: dict(o) if isinstance(o, Counter) else str(o))
    print(json.dumps(ann, indent=1, default=dict))
