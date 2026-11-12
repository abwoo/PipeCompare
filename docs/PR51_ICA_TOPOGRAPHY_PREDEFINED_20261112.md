# PR51 Stage 3B — Development-only ICA component identity and rejection stability

**Protocol written before the newly fitted component topographies or removal identities were observed.** All code resides on independent branch `codex/pr51-ica-topography-20261112`. This is a mechanistic diagnostic, **not** a frozen algorithm comparison, held-out participant study, or a change to the original PR #51.

## Why the extra measurement is necessary

The six-start development evidence [run 43169562528] and read-only [ICA mechanism audit run 43177240139] show that the apparent perfect ICA-start stability in subject 11 is compatible with **zero components removed** in both seed-0 recipes. For subjects 15/21, seed-0 removal counts are respectively two and three in both original and Minimax candidates. However, a count is not an identity: components may permute, flip sign, split or merge between separate ICA fits. Without topology matching and real removal identities, the six-start P300 range does not prove an ICA decomposition is reproducible.

The independent ICA source is ERP CORE P3 subjects **11, 15, 21**, 300 seconds each, same pinned source file SHA256, EEGLAB 2026.0.0 and MATLAB R2026a. Seeds are exactly 0..5, with **0,1,2 only for Minimax training**, 3,4,5 strictly for the audit. The source recipe IDs and original scientific eligibility floor remain unchanged; previous ICA initialization methodology (shuffle sample order) remains the tested perturbation.

## Technical measurement, frozen before re-fit

For **every ICA-bearing original candidate** and every ICA seed, persist a *read-only* normalized scalp-map matrix `pinv(icaweights * icasphere)`, IC channel indices, recorded removed independent component indices, total removed count, and explicit availability/error reason. If one candidate has more than one ICA or component-rejection stage, missing or inconsistent component metadata, a degenerate map, or inconsistent channel indices, mark it unavailable rather than silently treating it as stable. No ground-truth source labels can be inferred from ICLabel or scalp maps.

The original (seed-0) and training-only Minimax-chosen pipeline IDs are selected **without seed 3–5 fingerprint input**; for each candidate compare its seed-0 fingerprint separately with seeds 1,2,3,4,5. Cross-seed ICA component numbering and polarity have no inherent identity. Align maps by recorded channels, center and normalize each spatial vector, then maximize the sum of absolute Pearson spatial correlations with SciPy's one-to-one Hungarian assignment (`linear_sum_assignment`). Record matched median and minimum absolute correlation. Map removed indices through that assignment and report Jaccard overlap between the seed-0 and independent-seed removed component sets. **Both empty removal sets yield unavailable/undefined Jaccard, never 1.** Distinct channel sets, component counts, corrupted or nonfinite maps are audit failures, not similarity evidence.

Report all seed 1–5 comparisons, highlighting independent seeds 3–5, alongside original Rank scientific feasibility, calibrated P300, and HMM-unflagged effective trial counts. Report whether two distinct recipe IDs share an ICA decomposition at seed 0; do not count this as cross-start stability. Treat null results, cross-seed differences, and absent evidence explicitly.

## Interpretation and advancement gates

A reproducible topography does **not** establish source identity, neural specificity, full ERP waveform preservation, absence of eye-artifact leakage, or improved P300 point-estimate accuracy. No post-hoc topology threshold and no optimization of existing five-trial floors are permitted. No candidate moves to production from this diagnostic. The 14 untouched ERP CORE subjects `[24,25,26,27,28,29,31,32,33,34,36,37,38,39]` are not downloaded or fitted.

**Only after an explicitly pre-specified estimator or selector with convincing development improvement and scientific admissibility is frozen** may a new independent-person study be executed. This diagnostic is not that permission. Historical ERP CORE 20-person result **15/20 valid pairs, PARTIAL**, production unchanged, PR #51 OPEN and unmerged.

## Source provenance

- Historical raw six-seed run: https://github.com/abwoo/PipeCompare/actions/runs/43169562528.
- Completed no-refit source mechanism audit: https://github.com/abwoo/PipeCompare/actions/runs/43177240139.
- MATLAB source component snapshot: `+pipecompare/+eval/IcaFingerprint.m`; instrumentation in `+pipecompare/+run/Executor.m`.
- Python component matching: `benchmarks/compare_ica_topography.py`, with MATLAB and Python negative controls in `tests/test_ica_fingerprint.m` and `tests/test_ica_topography.py`.

Real-data results will be appended only after Actions completes; no run outcome is assumed by this plan.
