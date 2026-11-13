# Limitations and deviations, revised analysis (amendment 2)

Replaces the round-1 version. Consistent with `private/manuscript/manuscript_v2.md` Sect. 2.9 and 4.7 and with `private/reviews/RESPONSE_TO_REVIEWS.md` §9; where this file names something those do not, it is marked **not in MS**.

State of each item:
- **measured**: shown by a result file, a script or a command run on the project files (source given);
- **claimed**: stated in the project log, a review or the plan and not re-checked against a file.

Plan version 1: `protocol/SAP.md`, SHA-256 `e5412a331c776037a5e07707a6d2fe24b1fce46ba8aebb7a9472f06fd796f987` (measured: re-computed 2026-11-12, equal to the project log). Not registered.

## (a) Deviations from SAP v1

| # | Plan v1 said | What was done | State |
|---|---|---|---|
| D1 | Arm mapping audited before overlap; arm-level, placebo-adjusted where a placebo arm exists; supportive/null/contradictory criteria | Rule v1 failed (tier 1 73/100, tier 2 24/98); tier 2 dropped; rule v2 written from the error categories without reading adverse-event data; fresh 100-arm re-audit 93/100 (exact 95% interval 86.1–97.1%, lower bound below 90%). Main trial analysis is arm-level pooled incidence across molecules, not placebo-adjusted; the criteria are not used | measured (`private/reviews/P3_arm_audit_*`, `RESULTS_SUMMARY.md` §9) |
| D2 | Arm audit on a 200-arm sample | 100 arms per tier for v1 (198 read in tier 2), then 100 for v2; all verdicts by language model agents | measured (verdict files) |
| D3 | Legacy era: keep highest ISR per case | Code keeps latest FDA receipt date, then highest ISR | measured (`faers_q_build.py`) |
| D4 | Fixed-dose combinations included in a sensitivity extraction | Not run, in round 1 or the revision. **not in MS** | measured (no such extraction in `run_revision.sh`) |
| D5 | openFDA results kept as a sensitivity extraction | Not reported in the revision | measured (`run_revision.sh`) |
| D6 | Sensitivities K = 10, 50, all; uncleaned profiles | Run in round 1 only. K not re-run (MS 2.9). Uncleaned run not re-run. **uncleaned: not in MS** | measured (`run_revision.sh`; `results/round1/`) |
| D7 | Null-coverage simulation before H1 | Run in round 1 only (rejection rate 0.02, coverage 100/100); not re-run | measured (round-1 file); not re-run (MS 2.9) |
| D8 | Excess-overlap contrast reported with Δ | Not computed in the revision. **not in MS 2.9** | measured |
| D9 | Fixed testing sequence H1→H4 at one-sided 0.025 | Not used; intervals are the main results (MS 2.1) | measured (MS) |
| D10 | H2: dyadic regression with a Wald test on dyadic-clustered standard errors; covariates include signal counts | Molecule-level bootstrap (500 resamples for models); sum of signal counts not a covariate in the revision. **signal-count covariate: not in MS** | measured (`pairs_all.py`) |
| D11 | Therapeutic area = ATC level 3, else EPC; 2 × 2 | Replaced by ATC 4/3/1 and EPC on all pairs (amendment 2, A6) | measured |
| D12 | Same-target = share any ChEMBL target | Share any human protein target (amendment 2, A5) | measured |
| D13 | As-of report volumes may be approximated from yearly counts | Exact counts from receipt dates | measured (`asof_q.py`) |
| D14 | Not specified | Own indications used for cleaning are current, also in as-of profiles | measured (`asof_q.py`) |
| D15 | Not specified | Reference counts not recomputed after co-mention removal | measured (`comention.py`) |
| D16 | 9,999 / 10,000 resamples | Only the two main variants on extraction (b) and the three trial pair sets; all other runs 999 / 1,000 | measured (`run_revision.sh`) |
| D17 | θ reported | Rank-based θ not computed for Δ; θ_A computed for the as-of analysis and not reported | measured |
| D18 | Extraction (a) any role reported | Not re-run in the revision (disk exhausted) | measured (`PROJECT_LOG.md`, `run_revision_resume.sh`) |
| D19 | Signals: a ≥ 3, PRR ≥ 2, chi-square ≥ 4 | PRR and chi-square rounded to 2 and 1 decimals before comparison; b, c, d must be positive | measured (`faers_q_profiles.py`) |

## (b) Amendment 2 changes (all post hoc after internal review)

Labels refer to `private/reviews/RESPONSE.md`. The review had seen the round-1 results.

| # | Change | State |
|---|---|---|
| A1 | Drug strings also mapped through RxNorm names (IN, PIN, BN, MIN, SY, TMSY); unmapped drug rows 10,052,807 → 5,249,255 (6.8%) | post hoc; measured (`report_flow.json`) |
| A2 | Field-based duplicate rule; applicable to 8,128,914 cases, 483,101 removed, 20,571,050 left | post hoc; measured (`report_flow.json`) |
| A3 | IC025 profile carried through every analysis beside PRR | post hoc; measured (result files) |
| A4 | Extended cleaning (`ART2` pattern; stem rule at 0.1%) | post hoc; measured (`faers_q_profiles.py`) |
| A5 | Human protein targets; all-ChEMBL kept as sensitivity | post hoc; measured |
| A6 | Models on all 539,241 pairs with ATC 4/3/1, EPC with missing level, shared indication, volumes, year gap, co-mention; target × ATC 4 interaction; target-level analysis | post hoc; measured (`target_vs_class_R2_*.json`) |
| A7 | Sole-suspect extraction (e). Molecule-averaged Δ after co-mention removal: not done | post hoc; measured |
| A8 | As-of control C (same ATC level 1, not level 3 as asked); cluster-size sensitivities. Successor-early-years variant: not done | post hoc; measured |
| A9 | Trial analysis renamed parallel analysis; shared-trial exclusion; positive control. Serious-only: not done | post hoc; measured |
| A10 | Case characteristics for all deduplicated cases | post hoc; measured (`case_characteristics.json`) |
| A11 | Estimation framing; testing sequence dropped | post hoc; measured (MS) |
| A12 | Figures with intervals. K = all script and K re-run: not done. Fig. 1 legend overlaps the last row | post hoc; measured (figure viewed 2026-11-12) |

## (b2) Post hoc effect-modification analyses, 2026-11-13 (not in SAP v1, not in amendment 2)

| # | What was done | State |
|---|---|---|
| E1 | Primary Δ and the all-pairs adjusted model within strata of modality (biologic, peptide, small molecule), US route, ATC-1, approval era (`effect_modification.py`, `effect_modification_adjusted.py`) | post hoc, after all amendment-2 results were seen; measured (`EFFECT_MODIFICATION.md` §1–5) |
| E2 | Within therapeutic area (dominant MeSH branch of the molecule's trials). The FAERS-indication definition failed and is not used | post hoc; measured (§6) |
| E3 | Within ChEMBL target family (protein classification fetched from the public ChEMBL API) | post hoc; measured (§7) |
| E4 | Biologic vs small molecule restricted to single-target molecules | post hoc; measured (§8) |
| E5 | About 30 stratum comparisons across two rankings; no multiplicity correction; only modality, family and area were adjusted, one at a time | limitation; measured |
| E6 | Further analyses specified before they were run: see section (b3) | — |

## (b3) Amendments 3–5 and the structure test, 2026-11-13 (each fixed and hashed before its own results; all post hoc relative to SAP v1)

| # | Analysis | Specification | State |
|---|---|---|---|
| F1 | WP1 joint modifier model (biologic, kinase-family, oncology and approved-≥2000 pairs × same target) | `protocol/AMENDMENT_3.md` (sha256 `90e34dfa…`) | measured; reproduction 0.143 / 0.299 PASS; independent re-derivation 69/69 |
| F2 | WP2 modifiers in trial-reported AEs | amendment 3 | measured; reproduction 1.2142 PASS |
| F3 | WP3 same target, opposite action (ChEMBL action_type direction map) | amendment 3; INVERSE AGONIST = inhibiting, RELEASING and DISRUPTING AGENT excluded (lead's default, not chosen by the author) | measured |
| F4 | WP4 K = 10/20/50 × evidence floor a ≥ 10 | amendment 3; replaces the SAP K-sensitivity that was not re-run (A12) | measured; cells share random draws (A4.5) |
| F5 | WP5 per-target contrasts, shared reactions by group | amendment 3 | measured |
| F6 | Corrections after the adversarial methods review: era_missing aliased with year_missing; kinase label includes 21 biologics (small-molecule-only sensitivity); WP3 join misses 146/1,434 links (D4-consistent sensitivity); WP5 CI resampled targets (molecule bootstrap added) | `protocol/AMENDMENT_4.md` (sha256 `84ce901d…`) | measured; amendment-3 results stay primary |
| F7 | Reporting-artefact test: profiles rebuilt from late-only, FAERS-era, health-professional, consumer and non-US reports; biologic subtypes | `protocol/AMENDMENT_5.md` (sha256 `839ce2f6…`) | measured; full-data rebuild reproduces the stored profiles exactly |
| F8 | Biologic structure test (isotype, origin, format, ADC payload, target class; Thera-SAbDab and UniProt downloads) | `private/structure/SPEC.md` (sha256 `6bad6214…`); one labelled post-spec addendum | measured; reproduction biologic Δ 2.043 PASS |
| F9 | Therapeutic-area labels: FAERS-indication definition failed (271/1,039 labelled) and was dropped; the trial-MeSH definition is used | `EFFECT_MODIFICATION.md` §6 | measured |
| F10 | Data freeze v1 (2026-11-13) replaces provisional v0. Analysis inputs are identical; three raw files differ only in unused columns | `data/raw/trialytic_db/PROVENANCE.md` | measured |

**Multiplicity.** Across sections (b2) and (b3) there are 212 post hoc estimates (strata, differences and coefficients; most under both rankings, the trial analyses under one). No correction was applied. They are listed in `results/POST_HOC_COMPARISONS.md`. Only estimates with CIs are reported. The headline remains the SAP v1 primary estimate.

## (c) Planned or requested analyses not run, or not re-run on the rebuilt data

Listed in MS 2.9 and RESPONSE §9:
1. Any-role extraction (a): not re-run. measured.
2. Null-coverage simulation: not re-run. measured.
3. Sensitivity to profile length K: not re-run. measured.
4. Split of shared signals into labelled and unlabelled events: never run. measured.
5. As-of variant restricted to the successor's early years: not run. measured.
6. Grouping above the preferred term (no MedDRA licence): not run. measured.
7. Rank-based effect measure θ: not computed. measured.
8. Standard errors for dyadic clustering: not computed. measured.
9. Human reading of the audited arms: not done. measured (no human verdict file).

In RESPONSE §9 but not in MS 2.9 (author tasks or further analyses):
10. Final freeze of name map and target table with RxNorm, ChEMBL and AACT releases recorded, and regeneration. claimed as not done (RESPONSE §9.2); consistent with `data/MANIFEST.tsv` being freeze v0.
11. Second-person code review; co-author; affiliation, ORCID, model versions, commercial-interest statement.
12. Trial extras: serious-only, different-sponsor, MeSH-tree condition matching, relaxed within-trial analysis.
13. Sensitivity restricted to reports from 2014Q3 onward.
14. Logistic "any shared signal" model.
15. Schotland 2021 erratum check.
16. Official READUS-PV item wording.
17. BCa interval, stratified permutation, fixed-effects and binomial models, MRQAP.
18. First-approved molecule where the true first lacks an approval date.

Not in MS 2.9 and not in RESPONSE §9 (found while writing the supplement):
19. Uncleaned-profile run (plan v1 sensitivity). measured as not re-run.
20. Fixed-dose-combination sensitivity extraction (plan v1). measured as never run.
21. Excess-overlap contrast (plan v1). measured as not computed.
22. Sum of signal counts as a covariate (plan v1, round 1). measured as dropped.
23. Trial target-level analysis and single-target subset (STAT 3). RESPONSE §6 marks them Partly; not in §9.

## (d) Limitations

**How the study was done**
- The association was seen in exploratory openFDA work, part of which appeared in a trade article, before any plan existed; plan, data build and first analysis share 12 November 2026. claimed (project log).
- Amendment 2 followed a review that had seen round-1 results, so IC025 ranking, the extended cleaning list, the human-target restriction and the all-pairs models are post hoc. measured (project log order; MS 2.1).
- IC025 with extended cleaning was added after the shared terms of round 1 were inspected; its larger estimate (Δ 1.74 against 0.87) carries no more evidential weight. measured (RS §4).
- One author with language model agents; no second person examined the code; the arm-audit verdicts were not read by a human. measured (no human verdict file).
- The name map and target table come from a provisional freeze (v0) of a platform the author owns, taken while another work stream was editing drug-identity tables; RxNorm and ChEMBL releases not recorded. claimed (project log); release fields measured as absent.
- Data and code not yet deposited. measured (marker in MS).

**Spontaneous reports**
- Counts are reports, not incidence; notoriety, class labelling and stimulated reporting cannot be separated from pharmacology; the health-professional extraction gives a similar Δ (0.80 and 1.65). measured (RS §4).
- The field-based duplicate rule could be applied to 8,128,914 of 21,054,151 cases (38.6%); terms that appear to come from copied reports remain among the most shared: under PRR VAGINAL FLATULENCE 51 pairs, LUPUS VULGARIS 40, CONGENITAL HIATUS HERNIA 34; under IC025 LUPUS VULGARIS 71. measured (RS §1, §5–8). No with/without comparison of the rule alone was made.
- Under PRR ranking 60.1% of profile entries rest on ≤ 10 reports (18.8% on ≤ 3); under IC025 18.3% (1.7%). measured (RS).
- No MedDRA grouping or recoding to one version. measured.
- Cleaning is by text pattern: it removes abuse, dependence, overdose and administration-site terms of the injection, infusion, application and implant kind, and finds few disease terms. The pattern does not contain ADMINISTRATION SITE: 22 PRR and 5 IC025 profile entries contain it. measured (pattern; command on profile files 2026-11-12).
- Class-mates are in each molecule's reference group, which lowers disproportionality of class reactions (masking). claimed (not analysed).
- Co-mention removal is pairwise, reference counts not recomputed, and discards real events in switchers; the difference fell from 0.48 to 0.28 (PRR) and 1.06 to 0.80 (IC025) while the ratio stayed at 8.41→9.13 and 6.52→6.57. measured (RS).
- 6.8% of drug rows unmapped; fixed-dose combinations recognised reliably only where `prod_ai` exists (2014Q3 onward). Unmapped share measured; 2014Q3 claimed (round-1 limitation, not re-checked).
- `fdc_only` is computed over all roles, so a suspect mention inside a combination product can count as a suspect mention if the same ingredient appears alone as a concomitant drug in that report. measured (`faers_q_build.py`); extent not measured.

**Definitions**
- Target membership linked by name; mixes primary and secondary targets; no action type; parent and prodrug count as two molecules. measured (code); extent of prodrug pairs not measured.
- ATC codes for 968 and EPC for 658 of 1,039 cohort molecules. measured (RS).
- Class is partly a consequence of target and also assigned by use, so the adjusted coefficient (0.14 and 0.30) is the association not carried by class, not an effect of target. measured (estimates); interpretation per MS 4.2.
- Least squares on a mostly-zero count gives an adjusted mean difference; no dyadic-cluster standard errors. measured.
- Same-target pairs include 3 pairs that share only a non-cohort human target. measured (D4 file, 2026-11-12).

**As-of analysis**
- Successor profile current; case content latest version; indications and extended-cleaning flags current; co-mentions not removed; stimulated reporting on the successor not excluded. measured (`asof_q.py`).
- Control C matches only ATC level 1. measured.
- One cluster holds 145 of 403 sets; e_A without it 0.92 (PRR) and 1.71 (IC025); unweighted cluster mean 1.28 and 1.99, without interval. measured (RS).

**Trial analysis**
- Arm-level pooled incidence compared across molecules; same-target molecules are studied in similar populations; only the exclusions address this. Excluding shared-condition pairs left 196 of 1,187 same-target pairs (Δ 0.57). measured (RS §9).
- 7 of 100 audited arms were wrong (4 co-administered drugs, 2 fixed-dose products, 1 device arm). claimed (agent verdicts; project log).
- Serious and other events summed; events below a trial's reporting threshold count as zero; exact term matching kept 90.7% of event rows. measured.
- Zero-background terms (measured, corrected 2026-11-12): 65 of 6,030 trial-profile entries rest on a term with no background count (background set to 0.5; the RR ≥ 2 and chi-square ≥ 4 thresholds still apply). None of these terms appears in a second molecule's profile, so none contributes to a shared signal. MS 3.7 now states this; the earlier figure of 2 counted a different quantity.
- Trial-side name and target linkage uses the platform map without RxNorm names. measured (`ctgov_profiles.py`).
- Within-trial version not assessable: 84 molecules met minimum data; cohort 32 molecules, 10 targets, 91 same-target pairs. measured (RS §9).
- Positive control: own-profile overlap larger for 28.8% (PRR) and 56.4% (IC025) of 351 molecules. measured (RS).
