# Target-level overlap of adverse event signals (FAERS + ClinicalTrials.gov)

Analysis code, pre-specified protocol and aggregate results for the study:

> **Overlap of Adverse Event Signals Among Marketed Drugs That Share a Molecular Target, Within and Across Drug Classes: A Disproportionality Analysis of the FDA Adverse Event Reporting System and a Parallel Analysis of Trial-Reported Adverse Events**
> Fedor Lipskerov. *Manuscript in preparation.* The article text is not in this repository.

**Question.** Do two marketed molecules that act on the same human protein target share more of their top-ranked spontaneous-report safety signals than two molecules that do not? And how much of that excess is left after accounting for drug class, shared indication and co-reporting?

**Data.** FDA Adverse Event Reporting System quarterly files from 2004Q1 to 2026Q3. A parallel analysis uses adverse events reported in the results of ClinicalTrials.gov trials.

## Repository layout

```
.
├── protocol/            Statistical analysis plan (SAP v1, fixed before the confirmatory analysis)
├── code/
│   ├── python/          Every analysis step, from the FAERS quarterly files to figures
│   ├── run_revision.sh  End-to-end driver (amendment 2 analysis)
│   └── run_revision_resume.sh
├── data/
│   ├── README.md        Data sources, data dictionary, availability   ← start here
│   ├── MANIFEST.tsv     Row counts and sha256 of the frozen input tables
│   └── raw/faers_quarterly/MANIFEST.tsv   URL, size and sha256 of every FDA quarterly file used
├── results/
│   ├── RESULTS_SUMMARY.md               All reported numbers, generated from the files in tables/
│   ├── LIMITATIONS_AND_DEVIATIONS.md    Deviations from the SAP, post hoc amendments, analyses not run
│   ├── EFFECT_MODIFICATION.md           Post hoc: biologic vs small molecule, route, ATC group, era
│   ├── tables/          Aggregate results (JSON)
│   └── figures/         Main figures (PNG 300 dpi + SVG)
├── requirements.txt
├── CITATION.cff
└── LICENSE              MIT (code); results/ and protocol/ are CC BY 4.0
```

## Results at a glance

These numbers come from [`results/RESULTS_SUMMARY.md`](results/RESULTS_SUMMARY.md), the only source of numbers for the article.

- **Primary contrast.** 1,039 molecules, extraction b, top 20 signals per molecule. A same-target pair shares more signals than a different-target pair by **0.872** (95% CI 0.731–0.985) under PRR ranking with standard cleaning, and by **1.740** (1.523–1.906) under IC025 ranking with extended cleaning.
- **Adjusted pair model.** Run on all 539,241 pairs. It adjusts for ATC levels 4/3/1, FDA established pharmacologic class (EPC), shared indication, report volume, approval-year gap and co-mention. The same-target coefficient is **0.143** (0.078–0.213) under PRR and **0.299** (0.175–0.435) under IC025.
- **Post hoc** (see [`results/EFFECT_MODIFICATION.md`](results/EFFECT_MODIFICATION.md)). After the same adjustment, the effect is larger among biologics than among small molecules: 0.85 vs 0.11 under PRR. Route of administration adds nothing beyond modality.

All analyses are hypothesis-generating. A shared signal in spontaneous reports does not establish a shared causal adverse effect.

## Data sources and availability

| Source | Used for | Public? |
|---|---|---|
| FDA FAERS / legacy AERS quarterly ASCII files, 2004Q1–2026Q3 (`fis.fda.gov`) | Reports, drugs, reactions, indications. URLs and checksums are in `data/raw/faers_quarterly/MANIFEST.tsv` | Yes |
| ClinicalTrials.gov results via AACT (September 2026 snapshot) | Trial arms and reported adverse events | Yes |
| ChEMBL drug mechanisms | Molecule → target membership | Yes (CC BY-SA 3.0) |
| RxNorm, WHO ATC via RxNorm | Ingredient identity, ATC groups | Yes (subject to the UMLS/RxNorm licence terms) |
| openFDA, FDA NDC directory, FDA Purple Book, Drugs@FDA | Pharmacologic class, approval dates, routes, BLA status | Yes |
| **Trialytic database** ([trialytic.ai](https://trialytic.ai)) | Curated drug-name → RxNorm-ingredient resolution, drug ↔ trial bridge, assembled target membership and the FDA crosswalk. These were exported once, as a frozen snapshot, on 2026-11-12 | **No.** These are proprietary curated tables and are not redistributed |

**What is and is not in this repository:**
- **Not here:** row-level inputs and intermediate tables, i.e. anything exported from or built with the Trialytic database (per-molecule profiles, pair tables, name maps).
- **Here:** their structure (column dictionary) and checksums, in [`data/README.md`](data/README.md) and `data/MANIFEST.tsv`.
- **Access:** contact the author about access to the frozen snapshot for verification.

## Reproducing the analysis

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt

# 1. Download and parse the FDA quarterly files (public). Stores a DuckDB file in $FAERS_STORE (default ~/faers_quarterly)
python3 code/python/faers_quarterly_ingest.py

# 2. Place the frozen input tables in data/ (see data/README.md), then run the full pipeline:
export FAERS_STORE=~/faers_quarterly
export TRIALYTIC_ROOT=/path/to/trialytic     # only needed for steps importing Trialytic modules (see below)
sh code/run_revision.sh                      # logs to logs/rev_*.log; results to results/tables/

# 3. Summary and figures from saved results
python3 code/python/make_results_summary_r2.py
python3 code/python/make_paper_figures.py
```

**Environment variables.**

| variable | meaning |
|---|---|
| `FAERS_STORE` | Location of the parsed FAERS DuckDB store |
| `COHORT_DIR` | Which extraction's cohort to use (set by the driver) |
| `B` | Number of bootstrap resamples |
| `FDA_API_KEY` | Optional openFDA key, only for the openFDA fetch scripts. Never committed |
| `REF_DIR` | Folder holding `ndc.csv` and `pb.csv`: exports of the public FDA NDC directory and Purple Book. Used by `effect_modification.py` |
| `TRIALYTIC_ROOT` | Path to a Trialytic code checkout |

**`TRIALYTIC_ROOT` is needed by some steps.** Shared Trialytic modules are not included here: the drug-name normaliser (`build_drug_fda_crosswalk.normalize`), the FAERS reporting-artefact term list (`compute_faers_prr.BLACKLIST`) and the openFDA client. The steps that import them are `faers_q_build.py`, `faers_q_profiles.py`, `comention.py`, `asof_q.py`, `ctgov_profiles.py`, `ctgov_within_trial.py`, `full_rank_fetch.py`, `coreport_sample.py`, the seed script `asof_faers.py` and the `_v1` scripts. The artefact term list and the normaliser rules are described in the article's supplementary methods.

**Seeds.** Every resampling step uses a fixed seed (master seed 20261112; post hoc analyses 20261113). The output JSON files record their seed, B and input file.

## Analysis map

| step | script | output |
|---|---|---|
| FAERS ingest, deduplication, drug mapping | `faers_quarterly_ingest.py`, `faers_q_build.py` | `$FAERS_STORE/faers.duckdb` |
| Signals, cleaning, cohort, top-K profiles per extraction | `faers_q_profiles.py` | `data/derived/faers_q/<a–e>/` |
| Primary estimand (node permutation, molecule bootstrap) | `primary_inference.py` | `results/tables/primary_*.json` |
| Target vs class on all pairs, adjusted models | `pairs_all.py`, `h2_regression.py` | `results/tables/target_vs_class_*.json` |
| Co-mention exclusion | `comention.py` | `results/tables/comention_summary_*.json` |
| As-of-date (hindsight) test | `asof_q.py` | `results/tables/asof_summary_*.json` |
| Trial adverse-event replication | `ctgov_arm_attribution_v2.py`, `ctgov_profiles.py`, `ctgov_within_trial.py`, `positive_control.py` | `results/tables/primary_R2_ctgov_*.json`, `positive_control_*.json` |
| Null coverage, sensitivity grid | `null_coverage.py`, `sensitivity_grid.py` | `results/tables/` |
| Post hoc effect modification | `effect_modification.py`, `effect_modification_adjusted.py` | `results/tables/effect_modification_*.json` |
| Summary and figures | `make_results_summary_r2.py`, `make_paper_figures.py` | `results/RESULTS_SUMMARY.md`, `results/figures/` |

Paths under `private/` named in code or documents (manuscript, internal reviews, audit samples) are unpublished working files and are not part of this repository.

Files ending in `_v1.py` are the round-1 versions, kept so the record of what changed after internal review stays traceable.

## Citation

See [`CITATION.cff`](CITATION.cff). The article reference will be added on publication.

## Licence

Code: MIT ([`LICENSE`](LICENSE)). Results, figures and protocol: CC BY 4.0 ([`results/LICENSE.md`](results/LICENSE.md)). Third-party source data remain under their own terms.
