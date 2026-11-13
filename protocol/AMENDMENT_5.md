# Amendment 5 — is the biologic / newer-drug excess a reporting artefact?

**Fixed:** 2026-11-13, before any analysis below was run. The sha256 of this file is recorded in the project log.

**Status.** Post hoc. Written after amendment 3 found that, in the joint model, the same-target excess is carried by biologic pairs (+0.94 PRR) and newer-drug pairs (+0.32), and that the biologic excess is absent in trial-reported AEs.

**Question.** Does that excess come from *how* these drugs are reported, or from the drugs themselves? Three reporting mechanisms are candidates:
- stimulated reporting in early marketing years
- solicited, consumer-heavy reporting (e.g. patient-support programmes, mostly US)
- the modern FAERS era's reporting practice

## Design

**Data.** The parsed FAERS store (`$FAERS_STORE/faers.duckdb`, the same build as the paper). No new download.

**Cohort.** Fixed to the 1,039 molecules of extraction b (`data/derived/faers_q/b/D2,D4`), with their target graph.

**Signal and cleaning definitions.** Exactly those of `faers_q_profiles.py`, extraction b: suspect roles, no FDC-only, PRR ≥ 2, a ≥ 3, χ² ≥ 4, standard and extended cleaning, K = 20, PRR / standard and IC025 / extended rankings. The definitions block is reused unchanged.

**Reporting variants:**

| id | reports used | background |
|---|---|---|
| R1 late | per molecule, only reports first received ≥ 1 January of (first US approval year + 3). Molecules without an approval year keep all reports | all cases |
| R2 modern | cases first received on or after 2014-09-01 (FAERS era) | same subset |
| R3 HCP | reporter occupation MD, PH, OT, HP, RN | same subset |
| R4 consumer | reporter occupation CN | same subset |
| R5 non-US | reporter country present and ≠ US | same subset |

**Molecule eligibility within a variant.** ≥ 100 reports, as in the paper. Pairs are restricted to molecules eligible in the variant. Counts are reported.

## Estimates (per variant × ranking)

| | estimate | detail |
|---|---|---|
| (a) | overall Δ | molecule bootstrap B = 1,000 |
| (b) | Δ within biologics and within small molecules, and their difference | B = 1,000 |
| (c) | the WP1 joint model with the A4.1b design | overlap recomputed from the variant's profiles; covariates as stored (full-data covariates, stated as a limitation); B = 500; same_target, × biologic pair, × new pair reported |

**Reference.** The full-data values (`modifier_joint_a4_*`, `EFFECT_MODIFICATION.md` §5, §9) are recomputed on the same eligible-molecule subset, so each comparison is like for like.

**Seeds.** `SeedSequence(20261113).spawn(20)`, index 15 + variant number (R1 → 16 … R5 → 20; the reference → 15). If a spawn index runs past the 20 children, use `spawn(30)`.

## Pre-stated reading

| result | reading |
|---|---|
| × biologic (or × new) falls by ≥ 50% from the like-for-like reference in R1 and in R3 or R5 | **reporting artefact** (stimulated or solicited) |
| it keeps a CI excluding 0 in R1, R3 and R5 | **not explained by these reporting mechanisms** |
| otherwise | **mixed**; report as is |

## Biologic subtypes (B-part)

**Subtypes** (from INN stems in `molecule_types.csv` names):
- monoclonal antibody (-mab, including conjugates)
- fusion protein (-cept)
- cell or gene therapy (-cel, -vec, "gene")
- insulin / hormone (insulin, -tropin, -glucagon, parathyroid, follitropin)
- cytokine / interferon (interferon, -kin, -leukin, filgrastim, epoetin)
- enzyme / blood factor (-ase, factor, alfa enzymes, antithrombin, C1 esterase)
- other

**Estimate.** Stratum Δ within each subtype (both molecules in the subtype), shown when there are ≥ 10 same-target pairs. Descriptive only.
