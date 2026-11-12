# PR51 ICA Robust: P300 effect-shift source-trial forensic audit

**Date:** 2026-11-12. **Research status:** exploratory development-only, **NOT** a new frozen/held-out P3 validation. Do not merge PR #51.

## Immutable evidence

The source run https://github.com/abwoo/PipeCompare/actions/runs/43162539734 produced three real ERP CORE development subjects (11,15,21) with every algorithmic Rank-run feasibility constraint and original >=5 HMM-unflagged effective trial support intact. Original selector is unchanged.

Forensic run https://github.com/abwoo/PipeCompare/actions/runs/43164885120 downloaded precisely the original real-participant scoring-audit artifacts from that run (not raw held-out EEG), validated **8/8** independent decomposition tests and preserved the full `real_ica_trialweight_effect_decomposition.json` artifact. Reproduced published candidate condition means from source trial scores and weights before calculating differences.

| Subject | Original ID → robust ID | Original to robust calibrated P300 difference (µV) | Same per-trial P300 scores? | Source of raw effect change | Paired bootstrap 95% CI on calibrated change (µV) |
|---|---|---|---|---|---|
| 11 | 12 → 25 | 2.99147 → 2.93485 | No | Original P3 trial-score change −0.08211 µV; no trial-weight change | [−0.65494, 0.57015] |
| 15 | 5 → 2 | −0.85802 → −0.26381 | **Yes, exactly equal for all aligned source trial P300 scores** | **Trial weighting** (raw Δ +0.58715 µV), NOT trial-score change | [−0.00145, 1.18110] |
| 21 | 25 → 4 | 1.62904 → 1.62904 (seed0) | Yes | No weighting or raw-score change at seed0 | [0,0] |

Bootstrap: 10,000 paired trial-within-condition resamples, fixed independent development audit seed 511515, reapplying each recipe's original weights/gain; no threshold or selection of bootstrap result was used. CI crossing zero **does not establish equivalence**; no claim of biological/clinical importance, ground-truth artifact correction or a successful scientific replacement. Individual within-subject trial-bootstrap intervals are not independent cross-participant confidence intervals.

For subject15, candidate 2 uses no `weight_trials` in its recorded pipeline, whereas candidate5 includes `weight_trials`; the source reference P3 scores at seed0 are byte-level equal. Thus the previously alarming P300 effect attenuation arises from changing the statistical trial weighting rule, not from removing a different set of source P300 score values in that baseline comparison. This **does not establish full ERP waveform equivalence** or downstream suitability. Candidate2 also uses a different lowpass setting (20 Hz rather than 40 Hz), requiring independent shape/waveform QC.

## Scientific decision

The current common-feasible minimax algorithm passes numeric and scientific gates but does **not** preserve a fixed estimation family across alternative candidates. Consequently, minimizing worst-case uncertainty may change the P300 effect estimand (weighted vs unweighted), even when per-trial P300 score values are identical. It is premature to promote robust ID2 as better than original ID5 on subject15 merely based on a reduced numerical objective.

Next prospectively test a distinct **estimation-family-constrained robustness** candidate on development data: do not compare weighted vs unweighted P300 recipes in the same robust selection; retain source-gated minmax and abstention rules, report all three test subjects and evaluate independent ICA seed-group stability. This is a *new method-development proposal*, **not** a retrospective revision of the 20-person study or any held-out 14-person results.

Held-out 14-person EEG remains unexamined. Preserve original PR51 validation: 15/20 valid pairs, PARTIAL.
