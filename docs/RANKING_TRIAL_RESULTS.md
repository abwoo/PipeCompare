# One small ranking experiment: paired cloud results

Run on 2026-11-03 in [GitHub Actions](https://github.com/abwoo/PipeCompare/actions/runs/41633294024). PR #51 remains unmerged.

## Decision

**Restore/retain the original production ranking; do not promote the direct-score variant.**

Mean gain-corrected population-target MSE: original **0.483153194**, variant **0.487012814 uV²**. Variant minus original: **+0.003859620**, participant-bootstrap 95% interval **[+0.000000000, +0.011578861]**. Improved participants: **0**; worsened: **1**; unchanged: **9** (1e-9 reporting tolerance).

Mean raw MSE: original **0.523975043**, variant **0.526250991 uV²**.

## Existing statistical regression

The unchanged `test_statistics/testFactorialShrinkageChoosesBetter` was also run on both versions in [a separate cloud check](https://github.com/abwoo/PipeCompare/actions/runs/41634891238), execution `376ec268ccb7dbcef77e9e002bdb703ded19515e`. Original ranking passed: simulated mean excess selection error **1.57%**. The direct-score variant failed: **2.38%**, with the logged unrounded variant value 0.023775597453442. This test uses known-truth simulated candidate families, rather than the ten public records. Its existing performance assertion was retained unchanged. A paired-data improvement alone would not justify keeping a change that regresses this existing statistical check.

## Frozen change and decision rule

The single production experiment makes `Rank.shrink` return the existing estimated total-error score directly, rather than applying factorial empirical-Bayes smoothing. The same operation is used by recommendation, crossfit and selection-frequency calculations. No candidate/trial is selected using the injected truth. All filtering, ICA, adaptive cutoff/EOG checks, constraints and two extra ICA starts remain enabled.

The rule was committed before outcomes in [the trial plan](https://github.com/abwoo/PipeCompare/blob/d2da597e2448aa040a6198b5e91cfd37dfe2f313/benchmarks/ranking_trial_plan.json): retain only if all ten pairs pass provenance/completeness/replay checks, primary mean MSE decreases and raw mean MSE does not increase; otherwise restore the original. Exactly one variant was tested. This is a development screen, not a statistical superiority or independent confirmation claim.

## Every planned participant

| Participant | Original MSE | Direct-score MSE | Change, variant minus original |
|---|---:|---:|---:|
| 1 | 1.393875 | 1.393875 | +0.000000 |
| 2 | 0.017796 | 0.056392 | +0.038596 |
| 3 | 0.240221 | 0.240221 | +0.000000 |
| 4 | 1.750526 | 1.750526 | +0.000000 |
| 5 | 0.256920 | 0.256920 | +0.000000 |
| 6 | 0.416618 | 0.416618 | +0.000000 |
| 7 | 0.012027 | 0.012027 | +0.000000 |
| 8 | 0.083436 | 0.083436 | +0.000000 |
| 9 | 0.043610 | 0.043610 | +0.000000 |
| 10 | 0.616503 | 0.616503 | +0.000000 |

## Execution and checks

Each worker first runs the original algorithm and then the variant in a fresh MATLAB process on the same checksum-verified source excerpt. All ten original MSEs match the first PR51 run within 1e-7 uV². Within each pair, source, known truth, full Standard options and EEGLAB identity agree; fixed-control estimates agree within 1e-7 uV. Both native replay and generated executable-export checks pass. No planned participant is dropped.

Original production revision: `7151821edf8a5a54c3167d7f18b66c530911b568`. Trial execution/recipe: `d2da597e2448aa040a6198b5e91cfd37dfe2f313`; the recipe copies the tracked `benchmarks/variants/Rank.m` into its production path only on the trial worker. The exact variant file SHA256 is recorded in the frozen plan. Algorithm identity for the variant identifies this explicit recipe, not an untouched checkout of its production path.

The [independent collector](https://github.com/abwoo/PipeCompare/blob/d2da597e2448aa040a6198b5e91cfd37dfe2f313/benchmarks/compare_ranking_trial.py) recomputes readout errors from estimates and known truth and checks the paired identities before making the decision. Committed [summary](../benchmarks/results/ranking_trial_20261103/summary.json) and [paired readouts](../benchmarks/results/ranking_trial_20261103/records.json) retain the numerical evidence; cloud artifacts additionally contain candidate tables, calibration rows and logs.

## Meaning and limits

These are the same ten previously examined, 60-second public resting EEG backgrounds with known P3/blink injection. They are useful for development and regression comparison. Their reuse does not create an independent validation cohort. The unchanged first PR51 result remains documented in [PUBLIC_PILOT_RESULTS.md](PUBLIC_PILOT_RESULTS.md).

Removing the secondary smoothing does not satisfy the prespecified improvement rule. This experiment does not justify replacing the existing selector, nor does it solve the previously observed mismatch between estimated error and true error. Further scoring changes require a separate hypothesis and evaluation.
