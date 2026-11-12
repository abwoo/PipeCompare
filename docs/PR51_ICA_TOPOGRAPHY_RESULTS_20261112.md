# PR51 Stage 3B real ICA component topology and removal stability results

**Runs:** [original frozen readout and all 3 MATLAB real EEG jobs 43178065954](https://github.com/abwoo/PipeCompare/actions/runs/43178065954), [independent JSON singleton fix and artifact re-evaluation 43179870711](https://github.com/abwoo/PipeCompare/actions/runs/43179870711), [independent fingerprint-shape forensic 43180028433](https://github.com/abwoo/PipeCompare/actions/runs/43180028433) — all **SUCCESS**. The full 6-seed test was protocol-predeclared before reading new EEG. Original source raw run SHA `bee8cf0bd24a439429c46b0075d60382b533ff9d`. A1 changes only postprocessing scalar removal-index JSON parsing; all source EEG numerical results are immutable.

## Observed ICA dimensionality

Both recipe 12 (original) and recipe 25 (Minimax) for participant 11 have seed-0 fingerprints **17 channels × 15 components** but seeds **1, 2, 3, 4, 5** have **17 channels × 16 components**. All six source fingerprints exist; the original strict one-to-one full-component Hungarian comparison correctly **abstains** from declaring the 15- vs 16-component fits identical. All these fits removed zero ICs in the baseline reference. Thus the perfect (zero) P300 amplitude range is not evidence of stable ICA decomposition. The sample-order perturbation may affect the effective rank-estimation sample and therefore ICA dimensionality; technical root cause to be isolated before generalizing.

## Scientifically interpretable components, only subjects 15 and 21

Median matched *absolute* Pearson scalp-map correlations across **independent starts 3,4,5**:

- Subject 15, either original or Minimax recipe: **0.9999733, 0.9954555, 0.9960277**.
- Subject 21, either original or Minimax: **0.9999332, 0.9998888, 0.9997358**.

Aligned removed-component Jaccard overlap versus seed 0 (by seed 3 / 4 / 5):

| Participant | Original candidate ID | Original deletion overlap | Minimax ID | Minimax deletion overlap |
|---|---:|---|---:|---|
| 11 | 12 | **undefined** (15 vs 16 components and no removals) | 25 | **undefined** |
| 15 | 5 | 1.0 / 0.667 / 0.667 | 2 | **1.0 / 1.0 / 1.0** |
| 21 | 25 | 1.0 / 1.0 / 1.0 | 4 | 1.0 / 1.0 / 1.0 |

The Minimax selector produces more consistent *removed-component identities* in subject 15, but its **unseen-seed P300 range worsens** versus the original candidate: 0.014652 versus 0.006436 μV. Subject 21 is tied on the held-out seeds. Component identity overlap is not an ERP effect-size or precision endpoint.

## Separate rank-estimation confound and decision

The pinned upstream EEGLAB 2026.0.0 `pop_runica.m` calls `getrank(double(tmpdata(:,1:min(3000,size(tmpdata,2)))))`, sampling the **first 3,000 observed columns**. In the PipeCompare sensitivity runner, `ctx.icaSeed>0` permutes columns before the native ICA call. Consequently the rank estimator sees different subsets for different seeds. This implementation may confound the scientific question "vary ICA initialization alone" with ICA *model dimensionality*. The precise 15/16 observation is direct evidence of this mechanism's outcome; a pinned rank-lock experiment, not retroactive override, is needed to test whether keeping the same rank removes the confound.

**Gate:** Do NOT claim that a Minimax fixed ID has improved P300 measurement stability or full waveform fidelity. Preserve original 15/20 PARTIAL, original default algorithm, PR #51 OPEN/unmerged. No untouched 14-person cohort until a genuinely pre-specified improved science intervention is developed, tested and frozen separately.
