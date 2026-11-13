# Amendment 4 — corrections and sensitivity analyses after the methods review of amendment 3

**Fixed:** 2026-11-13, after the amendment-3 results and the adversarial review (`private/reviews/V2_findings.md`, F1–F5) were seen. This amendment is post hoc and adds no new question. It only corrects or bounds the amendment-3 analyses, which stay as run and as the primary amendment-3 results.

| # | Review finding | Action |
|---|---|---|
| A4.1 | F1: in WP1, `era_missing` is identical to `year_missing`, so the two are aliased | Report the WP1 coefficients of `year_missing` and `era_missing` as one aliased term; never quote either alone. Sensitivity A4.1b refits the joint model without `era_missing` (the other coefficients are expected to be unchanged) |
| A4.2 | F4: the WP1/WP5 "kinase" label includes 21 biologics (insulins at INSR/IGF1R, receptor antibodies, growth factors) | Sensitivity: the WP1 joint model with `kinase_pair` replaced by `kinase_sm_pair` = both molecules kinase-family **and** small molecule. All else as WP1 (A4.1b design, seed index 10). Text says "kinase-family target", never "kinase inhibitor", unless it uses A4.2 |
| A4.3 | F2: the WP3 direction join (R2 only) misses 146 of 1,434 D4 molecule × target links | Sensitivity: WP3 rerun with the D4-consistent join. That is R2 + R8 names resolving to one ingredient (unless already in R2) + salt-form recovery + human-protein target filter, exactly as in `faers_q_profiles.py`. All else as WP3 (seed index 12) |
| A4.4 | F3: the WP5 per-target CI resamples targets, not molecules | Report both: the target-resampling CI (pairs_all.py setting, labelled), and a molecule-bootstrap CI (B = 500, seed index 14) |
| A4.5 | F5: the WP4 grid cells share random draws | Documented only. Point estimates and per-cell CIs are valid; CIs across cells are not independent. No rerun |

These are reported in the supplement next to the amendment-3 results. No headline is chosen between them after seeing them: amendment 3 stays primary for WP1, WP3 and WP5 as specified.
