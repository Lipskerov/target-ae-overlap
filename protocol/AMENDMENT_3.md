# Amendment 3 — specification of further post hoc analyses

**Fixed:** 2026-11-13, before any of the analyses below was run. The sha256 of this file is recorded in the project log before execution.

**Status.** Post hoc relative to SAP v1 and amendment 2. Written after the effect-modification results of 2026-11-13 (`results/EFFECT_MODIFICATION.md` §1–8) were seen. Pre-specified relative to the results of the analyses it defines. A change after any result is seen becomes amendment 4.

**Reporting.** Point estimates with 95% CIs only; no claims of significance. Every comparison run is listed in one supplementary table. The paper's headline remains the SAP v1 primary estimate.

**Common settings:**
- Ranking: PRR with standard cleaning first, IC025 with extended cleaning second.
- Profiles: extraction b, K = 20, unless stated.
- Bootstrap: molecule level. B = 2000 for stratum Δ, B = 500 for regressions.
- Seeds: `numpy.random.SeedSequence(20261113).spawn(15)[index]`, with the index per analysis below.

**Molecule labels** (as produced on 2026-11-13; not to be edited):
- Modality: `results/tables/molecule_types.csv`.
- Target family: `data/derived/target_family.csv`, using the molecule rule in `code/python/effect_modification_target.py`.
- Therapeutic area: `results/tables/indication_types.csv`, column `area_trials`.
- Approval year: R10 `earliest_orig_ap`, else R11 `first_approval_year`.

## WP1 — joint modifier model (seed index 10)

**Data:** `results/tables/pairs_all_R2_{prr_clean,ic025_clean2}.csv.gz`, all 539,241 pairs.

**Pair modifiers.** Each is 1 if *both* molecules carry the label, else 0:

| modifier | label |
|---|---|
| `bio_pair` | modality = biologic |
| `kinase_pair` | family = kinase |
| `onc_pair` | area = neoplasms |
| `new_pair` | first approval year ≥ 2000 |

For each modifier, `<mod>_missing` = 1 if either molecule's label is unknown. Unknown means: area in {unknown, mixed}, or no approval year. Modality and family are never unknown, so they get no missing indicator.

**Model.** OLS, outcome = shared signals:
- the `model_full` covariates of `pairs_all.py`: same_target, same_atc4, same_atc3, same_atc1, atc_missing, same_epc, epc_missing, shared_indication, logrep_sum, logrep_absdiff, year_gap, year_missing, comention_share
- plus the 4 modifiers, the 2 missing indicators (area, era), and `same_target × modifier` for each of the 4.

**Reported:**
- every coefficient with its CI
- the leave-one-modifier-out refits (4)
- the count of same-target pairs in each modifier cell

**Reproduction check.** Without the added terms, `same_target` = 0.143 (PRR) and 0.299 (IC025), to 3 decimals.

**Reading.** An interaction whose CI excludes 0 in the joint model is reported as a modifier. One whose CI includes 0 is reported as "not separable from the others", with its estimate.

## WP2 — modifiers in trial-reported adverse events (seed index 11)

**Data:** `data/derived/ctgov/profile/` (profile_clean_k20, D2, D4, shared_condition_pairs, shared_trial_pairs).

**(A) Stratum Δ.** Both molecules in the stratum, primary-Δ definition. Strata:
- modality: biologic / small molecule
- family: kinase / GPCR / all other
- area: neoplasms / all other (unknown and mixed included in "other")

Differences: biologic − small molecule, kinase − GPCR, neoplasms − other.

Pair sets: all pairs, and pairs excluding shared-condition pairs.

**(B) Adjusted.** Trial-pair OLS over all pairs of trial-cohort molecules. Covariates:
- same_target
- same ATC-4/3/1 and atc_missing (R9)
- same EPC and epc_missing (R11)
- shared trial condition (shared_condition_pairs)
- the `bio_pair`, `kinase_pair`, `onc_pair` terms and their interactions with same_target

There is no report volume or co-mention in trials; the write-up states this.

**Reproduction check.** The whole trial cohort Δ = 1.2142.

**Reading.** A consistency check in an independent source. Direction plus CI, no "replication" language.

## WP3 — same target, opposite action (seed index 12)

**Direction of ChEMBL `action_type`** (R3):

| direction | action_type values |
|---|---|
| activating | AGONIST, PARTIAL AGONIST, ACTIVATOR, POSITIVE ALLOSTERIC MODULATOR, POSITIVE MODULATOR, OPENER |
| inhibiting | INHIBITOR, ANTAGONIST, BLOCKER, INVERSE AGONIST, NEGATIVE ALLOSTERIC MODULATOR, NEGATIVE MODULATOR, ALLOSTERIC ANTAGONIST, ANTISENSE INHIBITOR, RNAI INHIBITOR, DEGRADER |
| excluded | every other value, including RELEASING AGENT, DISRUPTING AGENT, BINDING AGENT, MODULATOR |

**Join.** R3 `normalized_name` → `ingredient_rxcui` via R2. Restricted to cohort molecules and cohort targets (`faers_q/b/D4`).

**Direction per molecule × target.** The set of directions over its rows. If it contains both activating and inhibiting, the molecule is excluded at that target.

**Unit.** A same-target pair with a defined direction for both molecules at the shared target:
- same direction (both activating or both inhibiting), or
- opposite direction.

A pair sharing several targets is classified as same-direction if it is same-direction at every shared target, as opposite if it is opposite at every shared target, and as "mixed" otherwise. Mixed pairs are reported separately.

**Reported:**
- the number of pairs in each class
- mean shared signals, same vs opposite, with the difference (molecule-bootstrap CI)
- the same difference adjusted: OLS on same-target pairs with an `opposite` indicator plus the `model_full` covariates except same_target
- both, repeated within GPCR + nuclear-receptor targets (target family of the shared target)

**Reproduction check.** The pair count is reported next to the 856 opposite-direction molecule-pair × target combinations counted by the methods review under its own mapping. A difference is documented, not forced.

**Reading.** Fewer shared signals for opposite-action pairs supports mechanism-driven sharing. Equal sharing points to indication or class.

## WP4 — profile length and evidence floor (seed index 13)

**Grid.** K ∈ {10, 20, 50} × ranking ∈ {PRR standard, IC025 extended} × floor ∈ {none, a ≥ 10}. Extraction b, built from `data/derived/faers_q/b/signals.csv.gz` with the same cleaning flags and ranking as the stored profiles.

**Floor.** Drop signals with a < 10 before taking the top K.

**Estimate.** Primary Δ, with 999 permutations and 1,000 bootstrap resamples (the amendment-2 setting for non-headline runs).

**Reproduction check.** K = 20 with no floor gives 0.8718 (PRR) and 1.7401 (IC025).

## WP5 — per-target contrasts and shared reactions (seed index 14)

**Per-target contrast.** The definition in `pairs_all.py` (`target_level`) for targets with ≥ 3 cohort members. One row per target, with family, number of members, number of biologic members, and the dominant area of its members.

**Reproduction check.** Mean over targets = 1.2494 (PRR).

**Shared reactions.** Among same-target pairs: the 25 most frequently shared reactions (count of pairs) within each of kinase pairs, GPCR pairs, ion-channel pairs, biologic pairs and neoplasm pairs, plus the overall top 25. Counts only.

## Decisions taken with this amendment (2026-11-13)

- **Era** stays a supplementary, unadjusted descriptor in the main text. It enters WP1 only as a modifier term. The FAERS rebuild for the era and stimulated-reporting question is **not** part of this amendment.
- **Ambiguous action types:** INVERSE AGONIST counts as inhibiting; RELEASING AGENT and DISRUPTING AGENT are excluded.
