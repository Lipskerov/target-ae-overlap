# Statistical analysis plan v1 — binding specification

Fixed 2026-11-12, before any confirmatory result was read. The full reasoning is in `SAP_draft.md` (biostatistician) and
`../private/reviews/P1_pv_fda_expert.md` (PV/FDA review); section numbers below refer to `SAP_draft.md`. Where this file and the
draft differ, this file rules. Decisions marked **[lead]** were taken by the lead on Fedor's instruction to proceed and can
be vetoed by him before the confirmatory run; any later change is logged in `PROJECT_LOG.md` as a protocol amendment.

Not a public registration. The sha256 of this file is written to the project log before the confirmatory run.

## 1. Data
- **Source [lead]:** FDA FAERS/AERS quarterly ASCII files, 2004Q1–2026Q3 (91 files), downloaded 2026-11-12.
  URL, size and sha256 per file: `data/raw/faers_quarterly/MANIFEST.tsv`. openFDA results are kept as a sensitivity
  extraction ("article extraction"), not as the confirmatory source.
- **Deduplication:** one row per case. FAERS era (2012Q4 on): per `caseid` keep the version with the latest `fda_dt`
  (ties: highest `primaryid`). Legacy era: per `CASE` keep the highest `ISR`. A case present in both eras keeps the FAERS-era
  version. Initial receipt date of a case = earliest `fda_dt` (or `init_fda_dt`) over all its versions.
- **Drug identity:** `prod_ai` when present, otherwise `drugname`, normalised with the repository's `normalize()` and mapped to
  an RxNorm ingredient through the frozen name map (`data/raw/trialytic_db/R2_name_to_ingredient_map.csv`, exported from the Trialytic database, https://trialytic.ai; not redistributed, see `data/README.md`).
  A product string that maps to two or more ingredients is a **fixed-dose combination**: it is excluded from single-molecule
  profiles in the primary extraction and included in a sensitivity extraction.
- **Extractions (reported side by side, PV review §6 item 1):** (a) any role; **(b) primary + secondary suspect — primary**;
  (c) primary suspect only; (d) as (b), health-professional reporters only (`occp_cod` in MD, PH, OT, HP, RN).
- **Molecule eligibility:** at least 100 deduplicated reports in the extraction.

## 2. Signals and profiles
- 2×2 table per molecule × MedDRA preferred term against all deduplicated reports (every PT has a background; no top-N cut).
- Signal: a ≥ 3, PRR ≥ 2, Yates chi-square ≥ 4 (Evans 2001).
- Removed before ranking, pre-specified here although first written after results were seen (stated as such in the paper):
  the reaction blacklist in `scripts/compute_faers_prr.py`; the artefact pattern in `private/sql/04_clean.sql`;
  the molecule's own indications (any PT that makes up ≥ 1% of that molecule's indication mentions in INDI).
- **Profile:** the K = 20 signals with the highest PRR. Sensitivities: K = 10, 50, all; ranking by the lower 95% limit of the
  ROR and by IC025; uncleaned profiles.

## 3. Cohort and pairs
- Targets: ChEMBL mechanism rows as frozen in `R3_chembl_target_drugs.csv`. Cohort targets have ≥ 3 eligible molecules with ≥ 1 signal.
- **Same-target pair:** the two molecules share any ChEMBL target (draft §0.2; the article's narrower rule misclassified 6 pairs).
- Salt-form name recovery (`07_salt_fix.sql`) is applied in the cohort.

## 4. Primary estimand and inference (draft §1–2)
- Δ = mean over molecules of (mean shared profile signals with same-target partners − mean with different-target partners).
- Node permutation, B = 9,999, one-sided; molecule-level bootstrap, B = 10,000, percentile 95% CI; master seed 20261112.
  Implemented in `code/python/primary_inference.py`.
- Before H1 is read: null-coverage simulation of the bootstrap and permutation test on permuted data (draft Q13).
- Reported with it: pair-averaged means (descriptive), Jaccard contrast, expected-overlap-adjusted contrast (draft §3.2).

## 5. Confirmatory hierarchy (fixed sequence, one-sided α = 0.025; draft §9)
H1 Δ > 0 on extraction (b) → H2 same-target coefficient > 0 in the dyadic regression with shared therapeutic area,
same class, report volumes, signal counts, approval-year gap and **co-mention share** as covariates → H3 as-of test →
H4 ClinicalTrials.gov replication. Everything else is exploratory.

## 6. Co-mention (PV review item 2) [lead]
- For every same-target pair and a seeded sample of five different-target pairs per same-target pair: the number of
  deduplicated reports naming both molecules (any role), and each molecule's profile recomputed without those reports.
- Reported: co-mention share by group; Δ-type contrast on the co-mention-excluded profiles; the covariate in H2.

## 7. Therapeutic-area control (PV review item 5)
ATC level 3 where `atc_codes` exists for both molecules, otherwise FDA EPC. Reported as the 2×2 (same/different target ×
same/different area) with all four cells; the within-area contrast is the key secondary.

## 8. As-of test (draft §6, PV review item 4)
First-approved molecule's profile from cases whose initial receipt date precedes the successor's approval date, same
extraction and rules as §1–2; successor profile current. Controls A and B as in draft §6.3; cluster inference over
connected components (draft §6.4). Successor-early-years and labelled/unlabelled splits are exploratory.

## 9. ClinicalTrials.gov replication (draft §8)
Arm-level, placebo-adjusted where a placebo arm exists; minimum data per molecule and the supportive / null / contradictory
criteria exactly as proposed in draft §8.4–8.5. Arm mapping audited on a 200-arm sample before any overlap is computed.

## 10. Reporting
READUS-PV (Fusaroli et al., Drug Saf 2024;47(6):575-584) checklist completed in `private/submission/`. A negative or null
result at any step of the hierarchy is reported as such.
