# Data

The row-level data used in this study are **not redistributed in this repository**. This file documents what was used, where it came from, and its structure, so that the analysis can be audited and rebuilt.

## 1. Public sources (obtain directly)

| Source | How to obtain | Used as |
|---|---|---|
| FDA AERS/FAERS quarterly ASCII files, 2004Q1–2026Q3 (91 files) | `code/python/faers_quarterly_ingest.py` downloads them from `https://fis.fda.gov/content/Exports/`. URL, size, sha256 and parsed row counts of every file are in `raw/faers_quarterly/MANIFEST.tsv` | Reports, suspect drugs, reactions (MedDRA PT), indications, case characteristics |
| ClinicalTrials.gov results, AACT September 2026 snapshot (`https://aact.ctti-clinicaltrials.org`) | AACT download | Trial arms, reported adverse events, conditions, interventions (`raw/ctgov/C1–C8`) |
| FDA NDC directory, FDA Purple Book | FDA downloads | US routes of administration and BLA status (post hoc effect-modification analysis; `REF_DIR`) |
| openFDA (`https://open.fda.gov`) | API | Report-count profiles for the stored-profile analysis (`raw/openfda_full/`) |

## 2. Trialytic database exports (not public)

Several inputs were exported once from the **Trialytic database** ([trialytic.ai](https://trialytic.ai)), a curated clinical-trial and drug-safety database built from AACT, openFDA, ChEMBL and RxNorm. The frozen snapshot was taken on 2026-11-12 and is marked provisional in the project record. The tables hold:
- curated drug-name → RxNorm ingredient resolution
- drug ↔ trial bridges
- ChEMBL-derived target membership
- an FDA class and approval crosswalk

These curated tables are proprietary and are **not** included. Their row counts and sha256 checksums are in [`MANIFEST.tsv`](MANIFEST.tsv). Access to the frozen snapshot for verification can be requested from the author.

Everything downstream of these exports is in `code/python/`:
- the FAERS deduplication and mapping
- signal detection
- profiles, pair tables and all statistics

## 3. Data dictionary (column names only)

`R*` = Trialytic database exports. `C*` = ClinicalTrials.gov extract assembled via the Trialytic drug ↔ trial bridge. `D*` and `derived/*` = tables built by this project's code from the above. All are local only (gitignored).

| file | columns |
|---|---|
| `data/raw/trialytic_db/R10_drug_approval_all.csv` | ingredient_rxcui, name, earliest_orig_ap, earliest_any_ap, n_applications |
| `data/raw/trialytic_db/R11_fda_class_all.csv` | ingredient_rxcui, pharm_class_epc, pharm_class_moa, first_approval_year |
| `data/raw/trialytic_db/R1_faers_signals_canonical.csv` | normalized_name, reaction_pt, report_count, rank, pct_of_drug, prr, ror, chi2, is_signal, signal_rank, updated_at, ingredient_rxcui, source_drug |
| `data/raw/trialytic_db/R2_name_to_ingredient_map.csv` | normalized_name, ingredient_rxcui, source_drug |
| `data/raw/trialytic_db/R3_chembl_target_drugs.csv` | target_symbol, target_name, target_chembl_id, normalized_name, action_type, mechanism_of_action, max_phase, n_trials |
| `data/raw/trialytic_db/R4_target_info.csv` | target_chembl_id, target_symbol, target_name, organism, target_type, synonyms, loaded_at |
| `data/raw/trialytic_db/R5_faers_summary_canonical.csv` | normalized_name, total_reports, n_serious, n_death, n_hospitalization, serious_fraction, top5_reactions, top_indications, updated_at, top_signals, n_trials, ingredient_rxcui |
| `data/raw/trialytic_db/R6_fda_crosswalk_cohort.csv` | normalized_name, query_name, n_trial_uses, generic_name, brand_names, unii, rxcui, application_numbers, pharm_class_moa, pharm_class_epc, approval_status, marketing_status, first_approval_year, faers_n_reports, faers_top_reactions, n_recalls, worst_recall_class, in_shortage, openfda_matched, match_method, updated_at, ingredient_rxcui, ingredient_name, atc_codes, atc_classes, identity_rxcui, ingredient_unii, approval_status_ref, first_single_approval, first_any_approval, first_bla_approval, application_numbers_ref, brand_names_ref, marketed, has_undated_applications, identity_source |
| `data/raw/trialytic_db/R7_drug_approval_cohort.csv` | ingredient_rxcui, name, name_source, first_ph3_start, n_ph3_trials, any_industry, n_applications, n_ap_submissions, earliest_orig_ap, earliest_any_ap, latest_ap, ap_after_first_ph3, resolved, loaded_at |
| `data/raw/trialytic_db/R8_rxnorm_name_to_ingredient.csv` | name, tty, rxcui, ingredient_rxcui, ingredient_name |
| `data/raw/trialytic_db/R9_rxnorm_atc.csv` | rxcui, tty, code, str |
| `data/raw/ctgov/C1_arms.csv.gz` | result_group_id, nct_id, ctgov_group_code, title, description |
| `data/raw/ctgov/C2_events_nonzero.csv.gz` | nct_id, result_group_id, event_type, organ_system, adverse_event_term, subjects_affected, subjects_at_risk, frequency_threshold, vocab |
| `data/raw/ctgov/C3_arm_at_risk.csv.gz` | nct_id, result_group_id, subjects_at_risk, frequency_threshold, event_rows |
| `data/raw/ctgov/C4_drug_trial_bridge.csv.gz` | nct_id, normalized_name, match_type |
| `data/raw/ctgov/C5_trial_conditions_mesh.csv.gz` | nct_id, mesh_term |
| `data/raw/ctgov/C6_trial_conditions_text.csv.gz` | nct_id, condition |
| `data/raw/ctgov/C7_trial_meta.csv.gz` | nct_id, phase, overall_status, study_type, lead_sponsor_name, lead_sponsor_class, start_date, enrollment, allocation, masking |
| `data/raw/ctgov/C8_interventions.csv.gz` | nct_id, intervention_type, name, source |
| `data/derived/D10_pairs.csv.gz` | a, b, same_target, n_shared, shares_ind, n_shared_clean |
| `data/derived/D1_molecules_all.csv` | ingredient_rxcui, name |
| `data/derived/D2_cohort_molecules.csv` | ingredient_rxcui, name |
| `data/derived/D3_cohort_targets.csv` | target_symbol |
| `data/derived/D4_molecule_target.csv` | ingredient_rxcui, target_symbol |
| `data/derived/D5_top20_signals_raw.csv` | ingredient_rxcui, reaction_pt |
| `data/derived/D6_top20_signals_clean.csv` | ingredient_rxcui, reaction_pt |
| `data/derived/D7_artefact_terms.csv` | reaction_pt |
| `data/derived/D8_top_indications.csv` | ingredient_rxcui, indication |
| `data/derived/D9_class_and_year.csv` | ingredient_rxcui, first_approval_year, pharm_class_epc |
| `data/derived/faers_q/b/D2_cohort_molecules.csv` | ingredient_rxcui, name |
| `data/derived/faers_q/b/D4_molecule_target.csv` | ingredient_rxcui, target_symbol |
| `data/derived/faers_q/b/molecule_totals.csv` | ingredient_rxcui, name, reports, in_cohort |
| `data/derived/faers_q/b/profile_clean2_ic025_k20.csv` | ingredient_rxcui, reaction_pt, a, ic025 |
| `data/derived/faers_q/b/profile_clean2_prr_k20.csv` | ingredient_rxcui, reaction_pt, a, prr |
| `data/derived/faers_q/b/profile_clean2_ror025_k20.csv` | ingredient_rxcui, reaction_pt, a, ror025 |
| `data/derived/faers_q/b/profile_clean_ic025_k20.csv` | ingredient_rxcui, reaction_pt, a, ic025 |
| `data/derived/faers_q/b/profile_clean_k20.csv` | ingredient_rxcui, reaction_pt |
| `data/derived/faers_q/b/profile_clean_prr_k20.csv` | ingredient_rxcui, reaction_pt, a, prr |
| `data/derived/faers_q/b/profile_clean_ror025_k20.csv` | ingredient_rxcui, reaction_pt, a, ror025 |
| `data/derived/faers_q/b/profile_raw_k20.csv` | ingredient_rxcui, reaction_pt |
| `data/derived/faers_q/b/signals.csv.gz` | ingredient_rxcui, reaction_pt, a, prr, chi2, ror025, ic025, artefact, own_indication, clean, clean2, rank_raw, rank_clean |
| `data/derived/ctgov/profile/D2_cohort_molecules.csv` | ingredient_rxcui, name |
| `data/derived/ctgov/profile/D4_molecule_target.csv` | ingredient_rxcui, target_symbol |
| `data/derived/ctgov/profile/molecule_totals.csv` | ingredient_rxcui, name, trials, subjects_at_risk, clean_trial_signals, in_cohort |
| `data/derived/ctgov/profile/profile_clean_k20.csv` | ingredient_rxcui, reaction_pt |
| `data/derived/ctgov/profile/shared_condition_pairs.csv` | rx_i, rx_j |
| `data/derived/ctgov/profile/shared_trial_pairs.csv` | rx_i, rx_j |
| `data/derived/ctgov/profile/signals.csv.gz` | ingredient_rxcui, reaction_pt, a, n_at_risk, rr, chi2, clean, zero_background |

## 4. Identifiers and vocabularies

- **Molecules:** RxNorm ingredient RxCUI (`ingredient_rxcui`).
- **Reactions:** MedDRA Preferred Terms as they appear in FAERS. MedDRA itself is licensed by ICH and is not distributed here.
- **Targets:** HGNC gene symbols of human single-protein ChEMBL targets.
- **Classes:** WHO ATC (levels 1, 3, 4) and FDA Established Pharmacologic Class (EPC).
