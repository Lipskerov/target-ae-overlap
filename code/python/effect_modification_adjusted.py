#!/usr/bin/env python3
"""Post-hoc (not in SAP v1): the paper's H2 adjusted pair model (pairs_all.py, model_full) refitted within modality strata,
plus one pooled model with same_target x pair-type interactions.

Inputs: results/tables/pairs_all_<label>.csv.gz (covariates as in the paper) and molecule_types.csv (effect_modification.py).
Model: OLS, overlap ~ same_target + same_atc4 + same_atc3 + same_atc1 + atc_missing + same_epc + epc_missing + shared_indication
       + logrep_sum + logrep_absdiff + year_gap + year_missing + comention_share   (identical to model_full)
CI: molecule-level bootstrap (pair weight = w_i * w_j), resampling the molecules of the stratum. B=500 as in the paper.
Also: within-class contrasts (same-target minus different-target mean overlap among pairs sharing ATC-4 / sharing an indication).
Writes results/tables/effect_modification_adjusted_<label>.json
"""
import json, os
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.abspath(os.path.join(HERE, "..", ".."))
T = os.path.join(PROJ, "results/tables"); B = int(os.environ.get("B", "500"))
COLS = ["same_target", "same_atc4", "same_atc3", "same_atc1", "atc_missing", "same_epc", "epc_missing", "shared_indication",
        "logrep_sum", "logrep_absdiff", "year_gap", "year_missing", "comention_share"]
mt = pd.read_csv(f"{T}/molecule_types.csv", dtype={"ingredient_rxcui": str}).set_index("ingredient_rxcui").modality

def fit(P, xcols, mols, seed):
    """OLS on pairs P; bootstrap over molecules `mols` (index into P.ii/P.jj)."""
    y = P.overlap.values.astype(float)
    X = np.column_stack([np.ones(len(P))] + [P[c].values.astype(float) for c in xcols])
    keep = [k for k in range(X.shape[1]) if k == 0 or np.ptp(X[:, k]) > 0]          # drop constant columns in a stratum
    names = (["intercept"] + xcols); names = [names[k] for k in keep]; X = X[:, keep]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    rng = np.random.default_rng(seed); nm = len(mols); bs = []
    for _ in range(B):
        w = rng.multinomial(nm, np.full(nm, 1 / nm)).astype(float)
        pw = w[P.ii.values] * w[P.jj.values]; m = pw > 0; Xw = X[m] * pw[m][:, None]
        try: bs.append(np.linalg.lstsq(Xw.T @ X[m], Xw.T @ y[m], rcond=None)[0])
        except np.linalg.LinAlgError: pass
    bs = np.array(bs); lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
    return {c: dict(beta=round(float(b), 3), ci95=[round(float(l), 3), round(float(h), 3)], share_boot_le0=round(float((bs[:, k] <= 0).mean()), 3))
            for k, (c, b, l, h) in enumerate(zip(names, beta, lo, hi))}, bs, names

def contrasts(P, mols, seed):
    y = P.overlap.values.astype(float); st = P.same_target.values == 1; out = {}
    rng = np.random.default_rng(seed); nm = len(mols)
    W = [rng.multinomial(nm, np.full(nm, 1 / nm)).astype(float) for _ in range(B)]
    for key, mask in (("same_atc4", (P.atc_both.values == 1) & (P.same_atc4.values == 1)),
                      ("different_atc4", (P.atc_both.values == 1) & (P.same_atc4.values == 0)),
                      ("shared_indication", P.shared_indication.values == 1),
                      ("no_shared_indication", P.shared_indication.values == 0),
                      ("no_comention", P.comention_reports.values == 0)):
        m1, m0 = mask & st, mask & ~st
        if m1.sum() < 10: out[key] = dict(same_target_pairs=int(m1.sum()), note="fewer than 10 same-target pairs"); continue
        def d(w=None):
            if w is None: return y[m1].mean() - y[m0].mean()
            r = []
            for m in (m1, m0):
                pw = w[P.ii.values[m]] * w[P.jj.values[m]]; r.append((pw * y[m]).sum() / pw.sum() if pw.sum() > 0 else np.nan)
            return r[0] - r[1]
        bs = np.array([d(w) for w in W]); lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        out[key] = dict(same_target_pairs=int(m1.sum()), different_target_pairs=int(m0.sum()), mean_same=round(float(y[m1].mean()), 3),
                        mean_different=round(float(y[m0].mean()), 3), difference=round(float(d()), 3), ci95=[round(float(lo), 3), round(float(hi), 3)])
    return out

if __name__ == "__main__":
    for label in ("prr_clean", "ic025_clean2"):
        P = pd.read_csv(f"{T}/pairs_all_R2_{label}.csv.gz", dtype={"rx_i": str, "rx_j": str}); P["atc_missing"] = 1 - P.atc_both
        P["mi"], P["mj"] = P.rx_i.map(mt), P.rx_j.map(mt)
        assert P.mi.notna().all() and P.mj.notna().all()
        paper = json.load(open(f"{T}/target_vs_class_R2_{label}.json"))["model_full"]["same_target"]
        out = dict(label=label, B=B, paper_model_full_same_target=paper, strata={})
        for nm, sel in (("all_pairs_refit", np.ones(len(P), bool)),
                        ("biologic_pairs", (P.mi == "biologic") & (P.mj == "biologic")),
                        ("small_molecule_pairs", (P.mi == "small_molecule") & (P.mj == "small_molecule")),
                        ("peptide_pairs", (P.mi == "peptide") & (P.mj == "peptide"))):
            Q = P[np.asarray(sel)].reset_index(drop=True)
            mols = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(mols)}; Q["ii"], Q["jj"] = Q.rx_i.map(ix), Q.rx_j.map(ix)
            res = dict(molecules=len(mols), pairs=int(len(Q)), same_target_pairs=int(Q.same_target.sum()))
            if res["same_target_pairs"] < 10: res["note"] = "too few same-target pairs"; out["strata"][nm] = res; continue
            res["unadjusted"], _, _ = fit(Q, ["same_target"], mols, 1)
            res["class_only"], _, _ = fit(Q, ["same_target", "same_atc4", "same_atc3", "same_atc1", "atc_missing", "same_epc", "epc_missing"], mols, 2)
            res["full"], _, _ = fit(Q, COLS, mols, 3)
            res["within_class_contrasts"] = contrasts(Q, mols, 4)
            out["strata"][nm] = res
            f = res["full"]["same_target"]; u = res["unadjusted"]["same_target"]
            print(label, nm, res["molecules"], res["same_target_pairs"], "unadj", u["beta"], u["ci95"], "| full", f["beta"], f["ci95"], flush=True)
        # pooled interaction model: does the adjusted same-target coefficient differ for biologic-biologic pairs?
        Q = P.copy(); mols = sorted(set(Q.rx_i) | set(Q.rx_j)); ix = {m: k for k, m in enumerate(mols)}; Q["ii"], Q["jj"] = Q.rx_i.map(ix), Q.rx_j.map(ix)
        Q["bio_pair"] = ((Q.mi == "biologic") & (Q.mj == "biologic")).astype(int)
        Q["mixed_pair"] = (((Q.mi == "biologic") & (Q.mj != "biologic")) | ((Q.mj == "biologic") & (Q.mi != "biologic"))).astype(int)
        Q["target_x_bio"] = Q.same_target * Q.bio_pair; Q["target_x_mixed"] = Q.same_target * Q.mixed_pair
        out["pooled_interaction"], _, _ = fit(Q, COLS + ["bio_pair", "mixed_pair", "target_x_bio", "target_x_mixed"], mols, 5)
        print(label, "interaction target_x_bio", out["pooled_interaction"]["target_x_bio"], flush=True)
        json.dump(out, open(f"{T}/effect_modification_adjusted_{label}.json", "w"), indent=1)
