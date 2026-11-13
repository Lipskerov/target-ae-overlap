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
