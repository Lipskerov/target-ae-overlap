# Target-level overlap of adverse event signals (FAERS + ClinicalTrials.gov)

Analysis code, pre-specified protocol and amendments, aggregate results and figures for the study:

> **Overlap of Adverse Event Signals Among Marketed Drugs That Share a Molecular Target, Within and Across Drug Classes: A Disproportionality Analysis of the FDA Adverse Event Reporting System and a Parallel Analysis of Trial-Reported Adverse Events**
> Fedor Lipskerov. *Manuscript in preparation.* The article text is not in this repository.

**Question.** Do two marketed molecules that act on the same human protein target share more of their top-ranked spontaneous-report safety signals than two molecules that do not? How much of that excess remains after drug class, shared indication and co-reporting are accounted for? And for which kinds of drugs does it remain?

## Repository layout

```
.
├── protocol/
│   ├── SAP.md                 Statistical analysis plan v1 (fixed before the confirmatory analysis)
│   └── AMENDMENT_3–5.md       Later analyses, each fixed and hashed before its own results (see protocol/README.md)
├── code/
│   ├── python/                Every analysis step, from the FDA quarterly files to the figures
│   └── run_revision*.sh       End-to-end driver for the amendment-2 analysis
├── data/
│   ├── README.md              Sources, availability, column dictionary  ← start here for data
│   ├── MANIFEST.tsv           Row counts and sha256 of the frozen input tables (freeze v1, 2026-11-13)
│   └── raw/faers_quarterly/MANIFEST.tsv   URL, size and sha256 of the 91 FDA quarterly files used
├── results/
│   ├── RESULTS_SUMMARY.md     Every reported number, generated from tables/ (sections 1–9 main, 10 post hoc)
│   ├── EFFECT_MODIFICATION.md Post hoc analyses explained (§1–12)
│   ├── POST_HOC_COMPARISONS.md Every post hoc estimate (212 rows) with source and specification
│   ├── LIMITATIONS_AND_DEVIATIONS.md  Deviations from the SAP, amendments, analyses not run
│   ├── tables/                Aggregate results (JSON only)
│   └── figures/v3/            Current figures (PNG 300 dpi + SVG); figures/round2/ = superseded
├── requirements.txt · CITATION.cff · LICENSE (MIT, code) · results/LICENSE.md (CC BY 4.0)
```

## Results at a glance

**Source.** All numbers are from [`results/RESULTS_SUMMARY.md`](results/RESULTS_SUMMARY.md). Each pair shows the PRR ranking (SAP v1) first, then the IC025 ranking (amendment 2).

**Main analyses:**
- **Primary contrast.** Data: 1,039 molecules; 7,923 same-target pairs. A same-target pair shares more top-20 signals than a different-target pair: Δ **0.872** (95% CI 0.731–0.985) / 1.740 (1.523–1.906).
- **All drivers fitted together.** On 539,241 pairs, adjusted for class, indication, volume, year and co-mention, the same-target coefficient falls 0.48 → 0.25 → **0.143** (0.078–0.213) / 0.299 (0.175–0.435).

**Post hoc analyses** (exploratory; [`EFFECT_MODIFICATION.md`](results/EFFECT_MODIFICATION.md)):
- **Where the adjusted remainder sits.** A joint model ([`AMENDMENT_3`](protocol/AMENDMENT_3.md)) places it in biologic pairs (interaction +0.94, 0.46–1.48) and in pairs of drugs approved from 2000 onward (+0.32, 0.10–0.58).
- **Kinase-family and oncology.** Neither interaction is separable from the others.
- **Reporting.** The biologic interaction persists in late-only, FAERS-era, health-professional-only and non-US reports ([`AMENDMENT_5`](protocol/AMENDMENT_5.md)).
- **Antibody structure.** Isotype, origin and format show no association.
- **Trial-reported adverse events.** The overall effect is present: Δ 1.21 (0.99–1.45). The biologic modifier is not: −0.16 (−0.74 to 0.54).

**Caveat.** All analyses are hypothesis-generating. A shared signal in spontaneous reports does not establish a shared causal adverse effect.

## Figures (`results/figures/v3/`)

| file | content |
|---|---|
| `fig1_design_flow` | Design schematic and data flow |
| `fig2_overall_and_adjustment` | Primary Δ and the adjustment ladder |
| `fig3_what_explains` | Within-class and indication strata, co-reporting, as-of (hindsight) test |
| `fig4_forest_adjusted_modifiers` | Adjusted same-target coefficient by target family, therapeutic area and drug type, plus the joint model |
| `fig5_biologic_robustness` | Biologic vs small molecule: adjustment ladder, robustness checks, oncology split |
| `fig6a_trials`, `fig6b_positive_control`, `fig6c_trial_modifiers` | Trial-reported adverse events |
| `figS1`, `figS3`, `figS4_*`, `figS9`–`figS12` | All extractions and rankings; target waterfall; shared reactions by group; structure analysis (S9–S11); reporting-artefact test (S12). S2 and S5–S8 are planned and not yet drawn |

## Data sources and availability

| Source | Used for | Public? |
|---|---|---|
| FDA FAERS / legacy AERS quarterly ASCII files, 2004Q1–2026Q3 (`fis.fda.gov`) | Reports, drugs, reactions, indications, reporter, country, dates | Yes |
| ClinicalTrials.gov results via AACT (September 2026 snapshot) | Trial arms, reported adverse events, MeSH conditions | Yes |
| ChEMBL mechanisms and protein classification (REST API) | Molecule → target; target family | Yes (CC BY-SA 3.0) |
| RxNorm; WHO ATC via RxNorm | Ingredient identity, ATC groups | Yes (subject to the UMLS/RxNorm licence terms) |
| openFDA, FDA NDC directory, FDA Purple Book, Drugs@FDA | Pharmacologic class, approval dates, routes, BLA status | Yes |
| Thera-SAbDab (Oxford OPIG) | Antibody format, isotype and origin (structure analysis) | Yes |
| UniProt | Target cellular location (structure analysis) | Yes (CC BY 4.0) |
| **Trialytic database** ([trialytic.ai](https://trialytic.ai)) | Curated drug-name → RxNorm-ingredient resolution, drug ↔ trial bridge, assembled target membership and the FDA crosswalk. Frozen exports, final freeze v1 on 2026-11-13 | **No.** Proprietary curated tables, not redistributed |

**What is and is not here:**
- **Not here:** row-level inputs and intermediate tables, i.e. anything exported from or built with the Trialytic database.
- **Here:** their structure and checksums, in [`data/README.md`](data/README.md) and `data/MANIFEST.tsv`.
- **Access:** contact the author about access to the frozen snapshot for verification.

## Reproducing the analysis

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
export FAERS_STORE=~/faers_quarterly          # parsed FAERS DuckDB store
export TRIALYTIC_ROOT=/path/to/trialytic      # only for steps that import Trialytic modules (see below)
export REF_DIR=/path/to/reference_exports     # ndc.csv, pb.csv (FDA NDC directory, Purple Book), mesh_tree.csv (MeSH tree numbers)

python3 code/python/faers_quarterly_ingest.py         # 1. download + parse the FDA quarterly files (public)
sh code/run_revision.sh                                # 2. main analysis (needs the frozen inputs in data/, see data/README.md)
python3 code/python/effect_modification.py             # 3. post hoc analyses, in this order:
python3 code/python/effect_modification_adjusted.py
python3 code/python/effect_modification_indication.py
python3 code/python/fetch_target_family.py && python3 code/python/effect_modification_target.py
python3 code/python/modifier_joint_model.py            #    amendment 3 (WP1–WP5) ...
python3 code/python/modifier_trials.py && python3 code/python/action_direction.py
python3 code/python/sensitivity_grid_R2.py && python3 code/python/target_level_terms.py
python3 code/python/modifier_joint_model_a4.py && python3 code/python/action_direction_a4.py && python3 code/python/target_level_molboot_a4.py   # amendment 4
python3 code/python/era_reporting.py && python3 code/python/biologic_subtypes.py   # amendment 5
python3 code/python/biologic_structure.py              #    structure analysis (spec fixed first)
python3 code/python/make_post_hoc_summary.py           # 4. RESULTS_SUMMARY.md + POST_HOC_COMPARISONS.md
python3 code/python/make_figures_v3.py && python3 code/python/make_figures_structure.py && python3 code/python/make_figures_main_v3.py   # 5. figures
```

**Built-in checks.** Every post hoc script first reproduces a known number and stops if it fails, for example 0.143 / 0.299, 1.2142, 0.8718 / 1.7401, 1.2494. The amendment-5 rebuild reproduces the stored profiles exactly.

**Seeds.** Every resampling uses `numpy.random.SeedSequence` with a registered spawn index:
- master seed 20261112 for the main analysis;
- 20261113 for the post hoc analyses, with indices listed in the amendments.

**`TRIALYTIC_ROOT`.** Shared Trialytic modules are not included here:
- the drug-name normaliser `build_drug_fda_crosswalk.normalize`;
- the FAERS reporting-artefact term list `compute_faers_prr.BLACKLIST`;
- the openFDA client.

The steps that import them are `faers_q_build.py`, `faers_q_profiles.py`, `era_reporting.py`, `comention.py`, `asof_q.py`, `ctgov_profiles.py`, `ctgov_within_trial.py`, `full_rank_fetch.py` and `coreport_sample.py`, plus the seed and `_v1` scripts. The artefact term list and the normaliser rules are described in the article's supplementary methods.

## Analysis map

| step | script | output |
|---|---|---|
| FAERS ingest, deduplication, drug mapping | `faers_quarterly_ingest.py`, `faers_q_build.py` | `$FAERS_STORE/faers.duckdb` |
| Signals, cleaning, cohort, top-K profiles | `faers_q_profiles.py` | `data/derived/faers_q/<a–e>/` |
| Primary estimand (node permutation, molecule bootstrap) | `primary_inference.py` | `primary_*.json` |
| Target vs class on all pairs, adjusted models | `pairs_all.py` | `target_vs_class_*.json` |
| Co-mention exclusion; as-of test | `comention.py`; `asof_q.py` | `comention_summary_*`, `asof_summary_*` |
| Trial adverse-event parallel analysis | `ctgov_arm_attribution_v2.py`, `ctgov_profiles.py`, `ctgov_within_trial.py`, `positive_control.py` | `primary_R2_ctgov_*`, `positive_control_*` |
| Effect modification, one factor at a time | `effect_modification*.py` | `effect_modification_*.json` |
| Joint modifier model; trial modifiers; opposite action; K/floor grid; per-target terms | `modifier_joint_model.py`, `modifier_trials.py`, `action_direction.py`, `sensitivity_grid_R2.py`, `target_level_terms.py` (+ `_a4` sensitivities) | `modifier_*`, `action_direction_*`, `sensitivity_grid_R2.*`, `shared_terms_by_group_*` |
| Reporting-artefact test; biologic subtypes; structure | `era_reporting.py`, `biologic_subtypes.py`, `biologic_structure.py` | `era_reporting.json`, `biologic_subtypes.json`, `biologic_structure.json` |
| Number source and figures | `make_results_summary_r2.py` + `make_post_hoc_summary.py`; `make_figures_*.py` | `RESULTS_SUMMARY.md`, `POST_HOC_COMPARISONS.md`, `figures/v3/` |

**Older script versions:**
- `_v1.py` files are the round-1 versions, kept so changes after internal review stay traceable.
- `h2_regression.py`, `asof_faers.py`, `make_paper_figures.py` and `make_figures_article.py` are superseded and kept for the record.

**Unpublished working files.** Paths under `private/` that are named in code or documents (manuscript, internal reviews, specifications) point to unpublished working files and are not part of this repository.

## Citation

See [`CITATION.cff`](CITATION.cff). The article reference will be added on publication.

## Licence

**Code:** MIT ([`LICENSE`](LICENSE)).
**Results, figures, protocol:** CC BY 4.0 ([`results/LICENSE.md`](results/LICENSE.md)).
**Third-party source data:** under their own terms.
