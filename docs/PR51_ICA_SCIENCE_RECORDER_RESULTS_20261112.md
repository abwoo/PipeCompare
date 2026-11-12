# PR #51 — ICA scientific recorder: completed cloud evidence and interpretation

**Date:** 2026-11-12. **Scope:** previously examined ERP CORE P3 participants 11, 15, 21; ICA seeds 0–5. **No use of the sealed 14 new subjects; no merger or default-recommender change.**

## GitHub verification

- [Real MATLAB six-seed scientific readout, run 43184910782](https://github.com/abwoo/PipeCompare/actions/runs/43184910782): all three subjects and collector **SUCCESS**.
- [Isolated MATLAB scientific recorder regression, run 43185657836](https://github.com/abwoo/PipeCompare/actions/runs/43185657836): **SUCCESS**.
- [Archived-only decomposition/ERP descriptive science, run 43186211965](https://github.com/abwoo/PipeCompare/actions/runs/43186211965): **SUCCESS**.
- The independent collector confirms three records, no original-run identity mismatches and no required scientific-recording failures; original recommendations and objective values match the immutable six-start control run 43169562528.
- Scientific implementation and protocol are described in [PR51_ICA_SCIENCE_RECORDER_PROTOCOL_20261112.md](PR51_ICA_SCIENCE_RECORDER_PROTOCOL_20261112.md). The complete validated per-seed data are preserved in the original runs' artifacts.

## What was captured

The optional, read-only recorder retains each candidate's actual ICA mixing/topographic matrices and channel labels, ICLabel probability distributions, removed component identities, complete Pz epoch timebase and target/standard waveforms, original retained-trial IDs, trial scores and weights, reference blink flags, and physical/weighted effective support. Its P300-window waveform integral is checked against the unchanged original estimate. A separate validator computes an IID trial-resampling bootstrap standard error from the preserved scores/weights, and compares scalp maps using sign-invariant absolute correlation with optimal component matching. Component indices without matching are **not** treated as physical component identity.

**Interpretation note:** the output called `bootstrapSE` estimates the standard error of the gain-corrected P300 difference, not ICA decomposition stability and not the variance of the entire waveform. The source trial bootstrap is IID, may understate error with correlated trials, and reuses the same development participants. These values are descriptive and cannot be used as a post hoc unbiased training/validation improvement claim.

## Six-start descriptive outcomes

The table reports *all six seeds 0–5*, not solely the unseen subset 3–5 from the prior report. The numerical readouts are from run 43186211965.

| Subject | Recipe | Mean bootstrap SE (μV) | Six-seed P300 range (μV) | Median absolute matched scalp-map correlation | Matched removed-IC Jaccard (median) | Full-waveform RMS change vs seed 0 (median μV) |
|---|---:|---:|---:|---:|---:|---:|
| 11 | Original 12 | 1.0584 | 0.000000 | 0.9598 | 1.0000 | 0.000000 |
| 11 | Minimax 25 | 1.0097 | 0.000000 | 0.9598 | 1.0000 | 0.000000 |
| 15 | Original 5 | 1.1714 | 0.007143 | 0.9960 | 0.6667 | 0.020923 |
| 15 | Minimax 2 | 1.1132 | 0.016632 | 0.9960 | 1.0000 | 0.015569 |
| 15 | Stability-first 6 | 1.1359 | 0.011148 | 0.9960 | 0.6667 | 0.022982 |
| 21 | Original 25 | 1.6713 | 0.115262 | 0.9997 | 1.0000 | 0.012920 |
| 21 | Minimax 4 | 1.5595 | 0.035693 | 0.9997 | 1.0000 | 0.012920 |

### Scientific reading

1. High overall topographic correlation does **not** imply identical artifact-removal decisions: subject 15 has a median matched removed-IC Jaccard of only 2/3 for original recipe 5 and stability-first recipe 6 despite median absolute map correlation 0.996. Thus decomposition geometry alone is not a sufficient recommendation target.
2. Subject 21's six-seed original P300 range (0.1153 μV) exceeds that of candidate 4 (0.0357 μV), while the Pz full-waveform RMS comparison summaries happen to match across the two recipes. A scalar endpoint and a full-waveform endpoint probe different aspects of reproducibility; neither demonstrates ground-truth neural accuracy.
3. Mean bootstrap SE is lower for the proposed alternatives in these three development participants, but this is an *in-sample descriptive comparison*: trials, subjects, recipe selection history and six starts are not independent confirmation evidence. None of these differences authorizes a positive superiority claim or fresh holdout launch.

## Decision gate

**Data recorder complete; algorithm promotion not justified.** No scientific data source or original 15/20 **PARTIAL** conclusion has changed. Continue separate development on a pre-specified criterion that explicitly combines eligible trial support, ICA component-removal agreement, full waveform and bootstrap precision without selecting based on the already inspected six-seed outcomes; add per-component classification and waveform sensitivity diagnostics and require no meaningful deterioration versus unchanged Standard before freezing code and opening the untouched 14-person cohort.

**Do not** merge PR #51, loosen original five-effective-trial floor, substitute ICLabel confidence for actual trial data, or assert that ICA map-matching proves neuron-level component identity.
