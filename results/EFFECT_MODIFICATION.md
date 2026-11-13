# Effect modification by drug type — post-hoc, exploratory

**Status:** run 2026-11-13. **Not in SAP v1.** These are post-hoc subgroup analyses, so they generate hypotheses and confirm nothing. If they go into the paper, they must be labelled as a deviation in `LIMITATIONS_AND_DEVIATIONS.md`.

**Code:** `code/python/effect_modification.py`, run with B = 2000 and seed 20261113.
**Outputs:**
- `results/tables/effect_modification_{prr_clean,ic025_clean2}.json`
- `effect_modification_checks.json`
- `molecule_types.csv`, the label for each molecule
- logs in `logs/effect_modification*.log`

**Reproduction check:** run over the whole cohort with no stratum, the code gives Δ = 0.872 (PRR) and 1.74 (IC025). These match `primary_R2_b_prr_clean.json` (0.8718) and `primary_R2_b_ic025_clean2.json` (1.7401).

## Method

- **Estimand.** This is the paper's primary Δ, applied inside one stratum. For each molecule, take its mean number of shared top-20 signals with same-target partners and subtract its mean with different-target partners. Average across molecules. Both molecules in a pair must be in the stratum.
- **CI.** Percentile 95% interval from a molecule bootstrap inside the stratum.
- **Stratum differences.** Computed as the difference of independent bootstrap draws. The p value is a two-sided bootstrap p.
- **Profiles.** Extraction b, top 20 signals, under two rankings: PRR with standard cleaning (the main ranking) and IC025 with extended cleaning (clean2).

### How molecules were labelled (1,039 cohort molecules)

| factor | rule | source |
|---|---|---|
| Modality: **biologic** (151) | Any BLA application, or a name match in the FDA Purple Book, or a biologic INN stem: -mab, -cept, -kinase/-plase, -kin, -cel, -vec, the biosimilar 4-letter suffix, alfa/beta, insulin/interferon/factor | R6 crosswalk, `ref_purple_book`. 9 molecules were added by stem alone: 6 non-US PD-1 antibodies (camrelizumab, sintilimab, spartalizumab, serplulimab, zimberelimab, nimotuzumab), streptokinase, lutropin alfa, factor VIII |
| Modality: **peptide** (23) | Synthetic peptides marketed under an NDA, identified by the stems -tide, -relin, -relix, -pressin, -reotide, -paratide, -tocin. Examples: semaglutide, tirzepatide, octreotide, leuprolide, teriparatide | Kept separate because they are neither classic small molecules nor BLA biologics. Defibrotide, an oligonucleotide, was excluded |
| Modality: **small molecule** (865) | Everything else | |
| Route | Every US route listed in the FDA NDC directory or the Purple Book for the molecule's applications or single-ingredient products. Routes were grouped as **oral** (oral, sublingual, buccal), **injected** (IV, SC, IM, intra-\*, percutaneous) or **local** (topical, ophthalmic, inhaled, nasal, otic, vaginal, rectal, transdermal) | `ref_ndc_product`, `ref_purple_book` |
| Route groups | Oral-only 421 · injected-only 234 · local-only 64 · multiple route classes 152 · **unknown 168**. The unknowns are mostly non-US or withdrawn molecules with no NDC product | |
| ATC-1 | Anatomical main group. Molecules in more than one group are labelled "multi" (90); molecules with none are "none" (71) | R9 |
| Era | First US approval before 2000 (441), 2000 or later (450), unknown (148) | R10/R11 |

**Modality and route overlap almost completely.** 139 of 151 biologics are injected-only. Only 80 small molecules are injected-only.

## Results

### 1. Biologic vs small molecule: the clearest modifier

| stratum | molecules | same-target pairs | Δ PRR (95% CI) | Δ IC025 (95% CI) | same-target pairs sharing no signal, PRR / IC025 |
|---|---|---|---|---|---|
| Small molecule | 865 | 7,189 | **0.67** (0.54–0.78) | **1.42** (1.21–1.58) | 78.2% / 60.2% |
| Biologic | 151 | 381 | **2.04** (1.53–2.57) | **3.72** (2.93–4.54) | 43.8% / 18.6% |
| Peptide | 23 | 47 | 1.27 (0.36–2.47) | 2.30 (0.62–4.40) | 53.2% / 38.3% |
| Biologic − small molecule | | | **+1.37** (0.86–1.94), p ≈ 0.0005 | **+2.31** (1.51–3.19), p ≈ 0.0005 | |
| Peptide − small molecule | | | +0.60 (−0.31 to 1.83), p = 0.27 | +0.88 (−0.77 to 2.96), p = 0.43 | |

- **Δ for biologics is about 3 times the small-molecule Δ under PRR and about 2.6 times under IC025.** The different-target baseline is the same for both (0.074 vs 0.070 shared signals under PRR), so the gap comes from the same-target side.
- **Partners from any modality.** Letting partners come from anywhere in the cohort barely changes the result: biologic 1.94, small molecule 0.68 (PRR).
- **Mixed pairs.** There are 259 same-target pairs of one biologic and one small molecule. They share 0.83 signals on average (PRR), which lies between biologic–biologic pairs (1.83) and small-molecule pairs (0.46).
- **Not carried by crowded targets.** The three targets with the most biologics are PDCD1 (13), IFNAR1 (10) and CD3E (10). With their molecules removed, the biologic Δ is 2.33 (1.69–3.03) under PRR and 4.05 (3.11–5.13) under IC025.

### 2. Route: the route effect is the biologic effect

| stratum (all modalities) | molecules | same-target pairs | Δ PRR (95% CI) | Δ IC025 (95% CI) |
|---|---|---|---|---|
| Oral only | 421 | 1,711 | 0.93 (0.70–1.13) | 1.87 (1.49–2.17) |
| Injected only | 234 | 491 | 1.69 (1.32–2.09) | 3.20 (2.60–3.82) |
| Local only | 64 | 130 | 1.41 (0.56–2.20) | 1.23 (0.37–2.15) |
| Multiple route classes | 152 | 342 | 0.52 (0.26–0.77) | 1.33 (0.78–1.89) |
| Unknown | 168 | 272 | 0.85 (0.37–1.41) | 1.49 (0.84–2.11) |
| Injected − oral | | | **+0.76** (0.34–1.23), p = 0.001 | **+1.33** (0.67–2.11), p ≈ 0.0005 |

The difference disappears **inside small molecules**, where modality cannot explain it:

| | Δ PRR (95% CI) | Δ IC025 (95% CI) |
|---|---|---|
| Small molecule, oral only (421) | 0.93 (0.70–1.14) | 1.87 (1.48–2.17) |
| Small molecule, injected only (80; 98 same-target pairs) | 1.01 (0.45–1.51) | 2.09 (1.12–2.91) |
| Injected − oral, small molecules | **+0.08** (−0.51 to 0.65) | **+0.22** (−0.75 to 1.12) |
| Injected-only biologic − injected-only small molecule | **+1.08** (0.29–1.87) | **+1.75** (0.65–3.09) |

- **Route on its own does not modify the effect.** The crude injected-vs-oral gap comes from the fact that biologics are injected.
- **Multi-route molecules had the smallest Δ under PRR**: 0.52, versus 0.93 for oral-only, a difference of −0.41 (−0.75 to −0.05). This did not reach significance under IC025 (−0.54, −1.13 to 0.16). These molecules (morphine, for example) have high background overlap with other drugs (0.16 shared signals with different-target partners).

### 3. Approval era

| | Δ PRR | Δ IC025 |
|---|---|---|
| Approved before 2000 (441) | 0.70 (0.47–0.84) | 1.40 (1.00–1.64) |
| Approved 2000 or later (450) | 1.20 (0.94–1.44) | 2.36 (1.98–2.71) |
| Difference | **+0.50** (0.23–0.84) | **+0.96** (0.53–1.51) |

Newer molecules show a larger Δ. Era overlaps with modality (most biologics were approved after 2000) and with reporting practice (the hindsight and stimulated-reporting concern raised in the paper). Neither was disentangled here.

### 4. ATC anatomical group (reference: A, alimentary/metabolic; Δ PRR 1.33, IC025 2.61)

| group | molecules | Δ PRR | Δ IC025 |
|---|---|---|---|
| L antineoplastic/immunomodulating | 257 | 1.16 (0.90–1.45) | 2.42 (1.98–2.87) |
| N nervous system | 209 | 0.81 (0.49–1.05) | 1.34 (0.95–1.68) |
| C cardiovascular | 101 | **0.26 (−0.01 to 0.58)** | **0.92 (0.30–1.44)** |
| R respiratory | 30 | **0.29 (−0.07 to 0.81)** | **0.38 (−0.07 to 1.09)** |
| Multi-group | 90 | 0.42 (0.11–0.87) | 0.87 (0.23–1.49) |
| G, H, B, M, S, none | 15–71 | 0.95–2.51, wide CIs | 1.17–3.09, wide CIs |
| D dermatological | 15 | 9 same-target pairs, not interpreted | |

- **Cardiovascular and respiratory drugs show almost no same-target excess.** They differ from A by −1.07 (C) and −1.04 (R) under PRR, and by −1.68 and −2.23 under IC025, all p ≤ 0.01. In these classes, drugs at different targets already share many signals (C: 0.22 shared with different-target partners, about 3 times the cohort average of 0.07).
- **Nervous-system drugs are lower than A under IC025** (−1.26, p = 0.012) but not under PRR (p = 0.09).
- 13 groups were tested, so these comparisons carry no multiplicity correction.

## 5. Adjusted model within biologics (added 2026-11-13)

**Code:** `code/python/effect_modification_adjusted.py`
**Output:** `results/tables/effect_modification_adjusted_{prr_clean,ic025_clean2}.json`
**Log:** `logs/effect_modification_adjusted.log`

**Model.** This is the paper's H2 `model_full`, unchanged. It is an OLS regression on all pairs of the stratum: shared signals ~ same target + same ATC-4/3/1 + ATC missing + same EPC + EPC missing + shared indication + log report volume (sum and difference) + approval-year gap + year missing + co-mention share. 95% CIs come from a molecule bootstrap (B = 500) inside the stratum.

**Reproduction check:** refitted on all pairs, the same-target coefficient is 0.143 (0.077–0.213). The paper has 0.143 (0.078–0.213).

### Same-target coefficient, in extra shared signals per pair

| stratum | molecules | same-target pairs | unadjusted | class only | **fully adjusted** |
|---|---|---|---|---|---|
| **PRR** all pairs | 1,039 | 7,923 | 0.48 (0.41–0.57) | 0.25 (0.19–0.33) | **0.14** (0.08–0.21) |
| PRR small molecules | 865 | 7,189 | 0.39 (0.32–0.47) | 0.19 (0.13–0.26) | **0.11** (0.05–0.16) |
| PRR **biologics** | 151 | 381 | 1.76 (1.24–2.47) | 1.12 (0.60–1.74) | **0.85** (0.33–1.37) |
| PRR peptides | 23 | 47 | 1.13 (0.34–2.31) | — | −0.07 (−1.22 to 1.05) |
| **IC025** all pairs | 1,039 | 7,923 | 1.05 (0.90–1.24) | 0.52 (0.40–0.66) | **0.30** (0.19–0.44) |
| IC025 small molecules | 865 | 7,189 | 0.87 (0.75–1.03) | 0.40 (0.29–0.53) | **0.24** (0.14–0.35) |
| IC025 **biologics** | 151 | 381 | 3.58 (2.72–4.81) | 2.01 (1.18–2.99) | **1.33** (0.57–2.14) |
| IC025 peptides | 23 | 47 | 1.74 (0.73–3.44) | — | 0.17 (−1.41 to 2.32) |

**Pooled interaction model.** This fits all pairs in one model: model_full plus biologic-pair and mixed-pair indicators, each with an interaction with same target.

| term | PRR | IC025 |
|---|---|---|
| same target (small-molecule or peptide pairs) | 0.09 (0.03–0.15) | 0.19 (0.08–0.31) |
| **same target × biologic–biologic pair** | **+1.10** (0.59–1.68) | **+2.22** (1.35–3.15) |
| same target × biologic–small-molecule pair | +0.49 (0.09–0.94) | +1.04 (0.40–1.79) |

### Pairs that share no indication, and pairs never named together in a report (unadjusted)

| | biologics, PRR | small molecules, PRR | biologics, IC025 | small molecules, IC025 |
|---|---|---|---|---|
| No shared indication | 0.68 (0.21–1.23), 93 pairs | 0.07 (0.04–0.12) | 1.57 (0.82–2.52) | 0.24 (0.16–0.33) |
| Never co-mentioned | 0.73 (0.30–1.22), 109 pairs | 0.12 (0.07–0.19) | 1.59 (0.95–2.36) | 0.35 (0.23–0.51) |
| Different ATC-4 | 0.83 (0.29–1.81), 101 pairs | 0.25 (0.18–0.33) | 1.93 (0.92–3.41) | 0.59 (0.47–0.75) |

### What the adjustment shows

- **Adjustment shrinks the biologic effect by about half to two-thirds, but it does not remove it.** The same-target coefficient for biologics falls from 1.76 to 0.85 under PRR and from 3.58 to 1.33 under IC025. All CIs stay above 0.
- **Biologics stay far above small molecules after full adjustment:** 0.85 vs 0.11 under PRR (about 8 times) and 1.33 vs 0.24 under IC025 (about 5.5 times). The pooled interaction model gives the same answer, with an interaction CI that excludes 0 under both rankings.
- **Not explained by shared indication or co-reporting.** Biologic pairs that share no indication, or that are never named in the same report, still have about 0.7 extra shared signals under PRR. That is 6 to 10 times the small-molecule value.
- **Peptides lose their excess after adjustment.** But there are only 47 same-target pairs, so the CIs are wide.
- **Not tested here:** how annotation differs between the two groups. A biologic usually has one specific protein target, while a small molecule's ChEMBL target list is broader. This alone could produce the modifier, as an artefact of annotation rather than biology.
- **Status:** post hoc, not in SAP v1.

## 6. Therapeutic area (indication type) — added 2026-11-13

**Code:** `code/python/effect_modification_indication.py`
**Output:** `results/tables/effect_modification_indication.json`
**Labels:** `indication_types.csv` (local only)

**Definition.** Each molecule's dominant MeSH disease branch across its ClinicalTrials.gov trials. A trial counts once per branch. A molecule needs at least 3 trials; a tie is labelled "mixed". Both molecules of a pair must be in the area.

A second definition, from FAERS top-10 indications matched to MeSH by exact name, **failed**. Only 271 of 1,039 molecules received a label, and no neoplasm label matched at all (MedDRA "… cancer" terms vs MeSH "… neoplasms"). It agreed with the trial-based label for 130 of 188 molecules labelled by both. It is not interpreted.

| area (trial-based) | molecules | same-target pairs | Δ PRR | Δ IC025 | **adjusted** PRR | **adjusted** IC025 |
|---|---|---|---|---|---|---|
| Neoplasms | 169 | 452 | 1.21 (0.83–1.69) | 2.47 (1.92–3.11) | **0.54 (0.24–0.95)** | **1.27 (0.73–1.86)** |
| Nervous system | 92 | 228 | 0.81 (0.36–1.44) | 1.33 (0.64–2.16) | −0.14 (−0.49 to 0.24) | −0.17 (−0.75 to 0.46) |
| Mental disorders | 76 | 430 | 0.83 (0.31–1.25) | 1.37 (0.70–1.95) | −0.06 (−0.28 to 0.12) | 0.25 (−0.12 to 0.60) |
| Cardiovascular | 65 | 97 | 0.60 (0.02–1.12) | 1.46 (0.47–2.32) | 0.26 (−0.24 to 1.01) | 0.59 (−0.35 to 1.93) |
| Metabolic/nutritional | 29 | 45 | 1.92 (0.98–3.04) | 4.93 (3.04–7.45) | 4.57 (0.05–8.34), unstable | 7.18 (0.32–12.53), unstable |
| Unknown (fewer than 3 trials) / mixed | 268 / 126 | 983 / 139 | 0.76 / 0.82 | 1.28 / 1.41 | — | — |

Smaller areas (urogenital, respiratory, skin, immune, blood, eye, digestive; 17–60 pairs) have wide CIs; see the JSON.

**Oncology is not just biologics.** Inside neoplasms, small molecules (129 molecules, 274 pairs) have Δ 1.27 (0.81–1.84) and biologics (34 molecules, 69 pairs) have Δ 1.28 (0.61–2.20) under PRR. Under IC025: 2.57 vs 2.88.

## 7. Target family — added 2026-11-13

**Code:** `code/python/fetch_target_family.py` (public ChEMBL REST API, protein classification of all 263 cohort targets, 0 errors) and `effect_modification_target.py`
**Output:** `results/tables/effect_modification_target.json`

**Definition.** A molecule's label is the family of its targets, or "multi-family" if its targets fall in more than one family.

| family | molecules | same-target pairs | Δ PRR | Δ IC025 | **adjusted** PRR | **adjusted** IC025 |
|---|---|---|---|---|---|---|
| GPCR | 292 | 2,468 | 0.74 (0.52–0.96) | 1.39 (1.06–1.72) | 0.02 (−0.07 to 0.11) | 0.14 (−0.03 to 0.31) |
| Ion channel | 115 | 1,145 | 0.43 (0.26–0.66) | 0.76 (0.49–1.08) | 0.01 (−0.09 to 0.15) | 0.13 (−0.02 to 0.33) |
| Nuclear receptor | 95 | 723 | 0.46 (0.23–0.77) | 0.76 (0.42–1.17) | 0.07 (−0.10 to 0.30) | 0.02 (−0.26 to 0.31) |
| Other enzyme | 154 | 802 | 0.89 (0.56–1.18) | 1.80 (1.26–2.29) | 0.10 (−0.10 to 0.33) | −0.09 (−0.46 to 0.30) |
| Transporter | 59 | 476 | 0.66 (0.26–1.28) | 1.70 (0.87–2.72) | 0.08 (−0.09 to 0.29) | 0.30 (−0.04 to 0.71) |
| **Kinase** | 108 | 555 | 1.08 (0.72–1.52) | 2.07 (1.50–2.71) | **0.46 (0.21–0.76)** | **1.05 (0.59–1.59)** |
| **Other membrane receptor** | 42 | 92 | 2.04 (1.06–3.03) | 3.11 (1.86–4.34) | **0.81 (0.10–2.04)** | **0.96 (0.06–2.34)** |
| Secreted protein | 39 | 82 | 2.11 (0.97–3.52) | 4.18 (2.35–6.09) | 0.46 (−1.44 to 3.27) | 1.02 (−1.37 to 4.07) |
| Protease | 39 | 95 | 1.31 (0.31–2.37) | 2.21 (1.07–3.22) | 0.61 (−0.62 to 5.51) | 0.83 (−1.56 to 5.29) |
| Surface antigen | 20 | 99 | 1.53 (0.70–2.76) | 3.44 (1.75–5.43) | 0.17 (−1.16 to 1.10) | 0.45 (−1.40 to 2.26) |

## 8. Annotation-specificity check — added 2026-11-13

Restricted to molecules with exactly one cohort target. 83.4% of biologics and 76.4% of small molecules have one target.

| | Δ PRR | Δ IC025 |
|---|---|---|
| Small molecule, single target (661 molecules, 3,387 pairs) | 0.79 (0.60–0.91) | 1.62 (1.31–1.79) |
| Biologic, single target (126 molecules, 260 pairs) | 2.30 (1.61–2.98) | 4.00 (3.04–5.01) |
| Difference | **+1.51 (0.86–2.25)** | **+2.38 (1.50–3.51)** |

The biologic modifier remains when every molecule has exactly one annotated target. Annotation *breadth* therefore does not explain it.

Two things are not tested:
- **Primary vs secondary mechanism.** A small molecule's single ChEMBL target may still not be its main mechanism.
- **Mutual adjustment.** Modality, target family, area and era have not been fitted in one model, and kinase, neoplasm and biologic labels overlap.

## 9. Amendment 3 results (specified before running; `protocol/AMENDMENT_3.md`) and amendment-4 sensitivities

**How these were run.** Workflow `wf_e263c12e-222`.
- All five reproduction checks passed: 0.143 / 0.299, 1.2142, 7,923 pairs, 0.8718 / 1.7401, and 1.2494.
- An independent re-derivation from the specification alone matched WP1 and WP3 on 69 of 69 numbers (`private/reviews/V1/V1_report.md`).
- The adversarial review found 4 confirmed issues and 1 plausible one (`private/reviews/V2_findings.md`). They are handled by `protocol/AMENDMENT_4.md`.

### WP1 — all modifiers in one model (`modifier_joint_*.json`; A4: `modifier_joint_a4_*.json`)

Same-target coefficient and its interactions. Each interaction is added on top of the base coefficient.

| term | PRR (A3) | IC025 (A3) | PRR, A4.2 kinase = small molecules only | IC025, A4.2 |
|---|---|---|---|---|
| same_target (pairs in no modifier category) | 0.054 (−0.016 to 0.115) | 0.136 (0.013–0.262) | 0.054 (−0.016 to 0.115) | 0.139 (0.015–0.265) |
| × biologic pair (381 same-target pairs) | **0.939 (0.461–1.481)** | **1.885 (0.994–2.865)** | **0.922 (0.456–1.477)** | **1.848 (0.966–2.859)** |
| × approved ≥ 2000 pair | **0.323 (0.099–0.577)** | **0.566 (0.192–0.959)** | **0.362 (0.156–0.615)** | **0.653 (0.286–1.012)** |
| × kinase-family pair (555) | 0.192 (−0.196 to 0.629) | 0.364 (−0.319 to 1.093) | 0.074 (−0.253 to 0.515) | 0.082 (−0.523 to 0.754) |
| × oncology pair | 0.130 (−0.231 to 0.550) | 0.522 (−0.031 to 1.169) | 0.179 (−0.161 to 0.570) | 0.633 (0.079–1.226) |

**Reading.** Fitted together:
- **Biologic and newer-drug** interactions keep CIs that exclude 0 under both rankings and both label versions.
- **Kinase** is not separable from the others.
- **Oncology** is separable only under IC025 in the A4.2 version.
- Once biologics and newer drugs are accounted for, the remaining same-target coefficient is near 0 under PRR.

`year_missing` and `era_missing` are aliased (A4.1); neither is quoted on its own.

### WP2 — modifiers in trial-reported adverse events (`modifier_trials.json`)

The trial cohort has 357 molecules and 1,187 same-target pairs.

| | Δ |
|---|---|
| Biologic (89 same-target pairs) | 1.04 (0.53–1.68) |
| Small molecule | 1.20 (0.95–1.49) |
| Biologic − small molecule | **−0.16 (−0.74 to 0.54)** |
| Kinase − GPCR | −0.39 (−1.11 to 0.41) |
| Neoplasms − other | +0.15 (−0.37 to 0.81) |

Adjusted (class and shared condition):
- same_target 1.09 (0.71–1.53)
- × biologic −0.59 (−1.33 to 0.28)
- × kinase −0.50 (−1.13 to 0.24)
- × oncology 0.32 (−0.34 to 1.00)

**Reading.** **None of the FAERS modifiers shows up in trial data.** The biologic point estimate goes in the opposite direction, though the CIs are wide. The overall same-target effect is present in trials.

### WP3 — same target, opposite action (`action_direction_*.json`; A4.3: `action_direction_a4_*.json`)

| | pairs: same / opposite / mixed | unadjusted opposite − same, PRR | adjusted, PRR | adjusted, IC025 |
|---|---|---|---|---|
| A3, all targets | 4,495 / 790 / 57 | −0.54 (−0.68 to −0.40) | −0.05 (−0.17 to 0.09) | −0.27 (−0.48 to −0.04) |
| A3, GPCR + nuclear receptor | 1,729 / 703 / 57 | −0.56 (−0.77 to −0.38) | 0.02 (−0.14 to 0.18) | −0.10 (−0.37 to 0.15) |
| A4.3 (D4-consistent join), all targets | 6,227 / 1,130 / 74 | −0.45 (−0.57 to −0.33) | −0.03 (−0.14 to 0.09) | −0.17 (−0.35 to 0.03) |

**Reading.** Opposite-action pairs share far fewer signals unadjusted. After adjustment for class, indication and co-mention, the gap is about 0. It stays slightly negative under IC025 in the primary join only. The unadjusted difference is therefore mostly class or indication: agonist and antagonist of one target sit in different classes and are used in different diseases.

### WP4 — profile length and evidence floor (`sensitivity_grid_R2.csv`)

Δ is > 0 with a CI excluding 0 in all 12 cells.

| | K = 10 | K = 20 | K = 50 |
|---|---|---|---|
| PRR, no floor | 0.33 (0.27–0.38) | 0.87 (0.74–0.99) | 2.84 (2.50–3.13) |
| PRR, floor a ≥ 10 | 0.47 | 1.14 (0.98–1.28) | 3.43 |
| IC025, no floor | 0.80 | 1.74 | 4.83 |
| IC025, floor a ≥ 10 | 0.73 | 1.52 (1.34–1.68) | 4.02 |

Cells share random draws (A4.5).

### WP5 — per-target contrasts and shared reactions (`target_level_R2_*.csv`, `shared_terms_by_group_*.json`, `target_level_molboot_a4.json`)

**Mean per-target contrast:**

| | point | CI resampling targets (`pairs_all.py` setting) | CI molecule bootstrap (A4.4) |
|---|---|---|---|
| PRR | 1.2494 | 1.03–1.48 | **0.95–1.38** |
| IC025 | 2.2558 | 1.95–2.60 | **1.92–2.42** |

88.9% (PRR) and 96.1% (IC025) of the 180 targets are positive.

**Most-shared reactions (PRR)** — mechanism-typical:
- **kinase-family:** DERMATITIS ACNEIFORM, PARONYCHIA, PALMAR-PLANTAR ERYTHRODYSAESTHESIA. Also insulin terms, because the "kinase-family" label includes INSR biologics (A4.2).
- **GPCR:** BLOOD PROLACTIN ABNORMAL, AKATHISIA, ANTICHOLINERGIC SYNDROME.
- **Ion channel:** HYPERTHERMIA MALIGNANT, LOCAL ANAESTHETIC SYSTEMIC TOXICITY.
- **Biologic:** CYTOKINE RELEASE SYNDROME, IMMUNE EFFECTOR CELL-ASSOCIATED NEUROTOXICITY SYNDROME, IMMUNE-MEDIATED LUNG DISEASE.
- **Neoplasm:** DERMATITIS ACNEIFORM, PARONYCHIA, SEROUS RETINOPATHY.

Counts only; whether these are on the label was not assessed.

## 10. Is the biologic / newer-drug excess a reporting artefact? (`protocol/AMENDMENT_5.md`, fixed before running)

**Code:** `code/python/era_reporting.py`
**Output:** `results/tables/era_reporting.json`
**Log:** `logs/era_reporting.log`

**Method.** Extraction-b profiles were rebuilt under five reporting variants for the fixed 1,039-molecule cohort.
- **Reproduction:** rebuilding on the full data (R0) reproduces the stored profiles exactly: 20,643 / 20,630 rows, symmetric difference 0.
- **Reference:** each variant is compared with full-data profiles on the *same* eligible molecules.

Joint-model interactions (A4.1b design), PRR (IC025 in brackets):

| variant | eligible molecules | × biologic pair: variant | × biologic: reference | × newer-drug pair: variant | × newer: reference |
|---|---|---|---|---|---|
| R1 reports ≥ 3 y after approval | 1,002 | **1.06 (0.48–1.74)** [1.98] | 0.91 [1.92] | 0.33 (0.10–0.59) [0.56] | 0.35 [0.60] |
| R2 FAERS era (≥ 2014-09) | 1,004 | **0.92 (0.43–1.63)** [1.72] | 0.94 [1.90] | 0.33 (0.08–0.57) [0.65] | 0.34 [0.60] |
| R3 health professionals only | 998 | **0.88 (0.36–1.56)** [1.74] | 0.95 [1.90] | 0.34 (0.08–0.59) [0.57] | 0.35 [0.61] |
| R4 consumers only | 826 | 0.48 (−0.04 to 0.99) [1.32 (0.24–2.35)] | 0.91 [1.84] | 0.53 (0.22–0.85) [0.57] | 0.35 [0.61] |
| R5 non-US reports only | 937 | **0.87 (0.39–1.48)** [1.78] | 0.96 [1.93] | 0.26 (0.02–0.55) [0.56] | 0.35 [0.62] |

Overall Δ (PRR) per variant: 0.85, 0.86, 0.82, 1.08, 0.83. Full data: 0.87.

**Pre-stated reading** (amendment 5): **not explained by these reporting mechanisms.**
- The biologic interaction does not fall when early-marketing reports are removed (R1). It does not fall under health-professional-only reporting (R3) or non-US-only reporting (R5).
- Its CI excludes 0 in R1, R2, R3 and R5 under both rankings.
- The newer-drug interaction behaves the same way.
- The only attenuation is in consumer-only reports (R4): 0.48 vs 0.91 under PRR, with a CI touching 0. That is the opposite of what solicited consumer reporting would predict.

**Consequence.** The FAERS-specific biologic excess is not stimulated reporting, not the modern-era reporting regime, not consumer or solicited reporting, and not US-specific reporting. That it is absent in trial-reported AEs (§9 WP2) remains unexplained. Two explanations are open and untested:
- **The trial biologic set is small:** 89 same-target pairs.
- **Trials record different events:** common, short-term events rather than the rare, mechanism-specific events where biologic pairs overlap in FAERS.

## 11. Biologic subtypes (amendment 5 B-part; `biologic_subtypes.json`)

| subtype | molecules | same-target pairs | Δ PRR | Δ IC025 |
|---|---|---|---|---|
| Monoclonal antibody | 85 | 217 | 1.98 (1.20–2.83) | 3.61 (2.40–4.65) |
| Insulin / hormone | 14 | 25 | 3.67 (2.34–5.41) | 7.24 (5.50–9.46) |
| Cytokine / interferon | 17 | 44 | 1.47 (0.62–2.53) | 2.81 (1.57–3.98) |
| Other biologic | 18 | 11 | 2.05 (0.00–5.20) | 3.17 (−0.05 to 7.63) |
| Enzyme, cell/gene, fusion | 9 / 6 / 2 | 4 / 7 / 0 | not interpreted (< 10 pairs) | |

Biologic same-target pairs concentrate at PDCD1 (78), IFNAR1 (45), CD3E (45), MS4A1 (28) and CD19 (28). Insulin analogues are near-identical molecules, so their high Δ is expected.

## 12. Does antibody or biologic *structure* explain the biologic excess? (structure agent; spec fixed first)

**Spec:** `private/structure/SPEC.md`, sha256 `6bad6214…`, written 14:23 before any result (output 14:28).
**Code:** `code/python/biologic_structure.py`
**Output:** `results/tables/biologic_structure.json`

**Data:**
- Thera-SAbDab: format, isotype and origin; 86 of 151 biologics matched.
- UniProt keywords for 64 targets.
- INN stems for ADC payload and CAR-T.

Downloads with sha256 are in `data/derived/structure/`.

**Reproduction:** biologic Δ 2.043 (PRR) / 3.721 (IC025), identical to §1.

| test | PRR | IC025 | reading |
|---|---|---|---|
| Different-target pairs, same vs different isotype (1,049 vs 1,197) | +0.02 (−0.08 to 0.15) | −0.03 (−0.29 to 0.30) | null |
| Different-target pairs, same vs different origin (chimeric / humanized / human) | +0.06 (−0.04 to 0.24) | +0.03 (−0.15 to 0.30) | null |
| Different-target pairs, same vs different format | +0.04 (−0.04 to 0.14) | −0.01 (−0.23 to 0.24) | null |
| Different-target pairs, same vs different target structural class | +0.03 (−0.03 to 0.07) | +0.22 (0.09–0.35) | not met under both rankings |
| ADCs with different antigens, same vs different payload class (11 vs 16 pairs) | +0.45 (−0.50 to 2.80) | +0.37 (−1.11 to 2.81) | small; ocular-toxicity terms |
| ADC × small molecule hitting the payload's target vs no shared target | +0.21 (−0.04 to 0.55) | +0.96 (0.16–1.89) | febrile neutropenia, neuropathy |
| CD3E pairs: T-cell engager × T-cell engager vs × anti-CD3 monospecific (28 vs 16) | +1.87 (0.49–3.64) | +4.05 (2.58–5.49) | CRS, ICANS: format effect within one target |
| CAR-T pairs with different targets (8 pairs; not interpreted) | 4.5 vs 0.12 | | CRS, ICANS, HLH |
| PD-1 × PD-L1 pairs vs other different-target antibody pairs (39 vs 3,390) | +0.80 (0.17–1.92) | +3.12 (1.41–5.00) | pathway, not structure: immune-mediated AEs |
| **Biologic Δ excluding CAR-T, T-cell engagers and ADCs** (post-spec addendum; 268 pairs) | **2.15 (1.57–2.80)** | **3.83 (3.01–4.90)** | unchanged from 2.04 / 3.72 |

**Reading.** Generic antibody structure explains none of the biologic excess: isotype, origin and format are all null. The excess survives removing payload-bearing and cell-engaging products.

Structure and format do produce **mechanism-typical sharing in small, specific cells**:
- **CD3 engagers and CAR-T:** CRS and ICANS.
- **ADC payloads:** neuropathy, cytopenias, ocular toxicity.
- **The PD-1/PD-L1 pathway** (a different target on the same pathway): immune-mediated AEs.

The last of these points to **pathway** rather than structure as the next lens.

**Not testable here:**
- Fc engineering: there is no constant-region field.
- Payload contrasts: there are only 8 ADCs in the cohort.
- Immunogenicity: this needs anti-drug-antibody (ADA) rates from labels.

## Summary

1. **Biologics are the clearest modifier, and the difference survives the paper's full adjustment** (section 5: 0.85 vs 0.11 under PRR). For biologics, same-target sharing is about 3 times that of small molecules, under both rankings, and it is not carried by the PD-1, interferon or CD3 clusters.
2. **Route explains nothing beyond modality.** Inside small molecules, injected and oral are the same (+0.08, CI crosses 0).
3. **Newer drugs show a larger Δ**, and **cardiovascular and respiratory drugs show almost none**.

## Limitations of this analysis

- **Post hoc.** There are about 25 stratum comparisons with no multiplicity correction. Only the biologic and era differences hold under both rankings with p < 0.001.
- **Adjustment was run only for modality (section 5).** The adjusted model was not run for route, era or ATC groups. Those results are unadjusted.
- **How biologics are made.** A biologic's target is usually a single protein bound with high specificity, while a small molecule's ChEMBL "target" can be one of several. This could explain the modifier by itself, as a difference in target annotation rather than in biology. Not tested.
- **Labels.** Modality comes from BLA/Purple Book plus INN stems. The 9 stem-only biologics were checked by name. The other labels were not checked by hand. Route is unknown for 168 molecules (16%). Route means routes *ever marketed in the US*, not the route behind each FAERS report.
- **Small strata.** Peptides (23), local-only (64) and the small ATC groups give wide intervals.
