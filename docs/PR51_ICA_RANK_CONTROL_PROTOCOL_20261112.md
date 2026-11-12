# PR51 rank-controlled ICA initialization stress test — prospective *analysis amendment* on development subjects

Written on 2026-11-12 after the prior six-start scientific recorder and Pareto analysis, **before any rank-controlled EEG fitting**. This is an exploratory method-development experiment. Original ERP CORE P3 study result **15/20 paired valid measurements, PARTIAL** must remain intact. PR #51 OPEN/unmerged; original ranking and scientific thresholds are not changed.

## Confirmed causal code-path problem

Both pinned EEGLAB 2026.0.0 `pop_runica.m` and PipeCompare's matching `runIca` path estimate ICA rank on **the first 3000 samples**. The ICA initialization stress test deliberately permutes *all* fitted samples for seeds >0, so this 3000-sample rank estimator can select a different PCA dimension solely from the chosen permutation. In immutable science run [43184910782](https://github.com/abwoo/PipeCompare/actions/runs/43184910782), ERP CORE subject 11 has exactly 17 fitted channels but 15 components at seed 0 versus 16 at seeds 1–5. No components are removed from original candidate IDs 12 and 25. This makes decomposition dimension an uncontrolled factor and prevents rigorous same-dimensional component-set comparison.

## Predefined technical correction (shadow-only)

In opt-in `icaRankStable=true` mode only:

1. Use the same pre-ICA dataset after all existing preprocessing and existing ICA-fit-copy high-pass and trimming, **before sample reordering**.
2. Estimate rank from the **entire** fitted sample matrix, using EEGLAB's two rank checks on all samples: MATLAB `rank(X)` and `sum(eig(cov(X',1)) > 1e-7)`, choose the minimum. Do not use any new optimized threshold or experimental P300 outcome.
3. Explicitly pass `'pca', fixedRank` to runica for **all** ICA starts, including seed 0. Never append another automatic 3000-sample PCA dimension in the alternate legacy path.
4. Then preserve the existing exact sample permutations for seeds 1–5, original EEG channel data, HMM labels, trial floor 5, all cutoff/EOG logic, unchanged `Measure.m` and `Rank.m`. No production default, no merge.

The explicit full-data rank may differ from any original run's first-3000 rank. Therefore **point estimates, recommended recipe IDs, bootstrapped SME and science eligibility are permitted to change** and must be reported; old study outputs are historical controls, never retroactively replaced. For equal underlying fit data, all starts should have the same count of fitted ICA components; if not, report failure and do not claim this isolates initialization.

## Predefined pilot and stopping rule

Launch only already-seen ERP CORE participant **11** as the engineering smoke test. Run the same six-start standard search with full science recorder, compare original frozen source identity and sample support, independently check ICA rank/fit dimensionality across six starts for candidates 12 and 25, and verify all full Pz waveforms and original internal gain/score arithmetic. Archive failed evidence even if the numerical or scientific consistency gate rejects it.

If the engineering smoke test establishes rank isolation without new implementation defects, consider extending to development participants **15/21**, then recompute the same predefined science analyses. No newly held-out 14 subjects may be opened before an objectively supported method and source+algorithm+gate protocol are frozen. This trial is not independent confirmation of improvements in participant-level ERP accuracy, which lacks known neural ground truth.
