# PR51 Stage 2B results — full-record source HMM, development subjects 4 and 8

**Completed:** [GitHub Actions 43179740278](https://github.com/abwoo/PipeCompare/actions/runs/43179740278) all scientific unit, two real MATLAB subjects and collector **SUCCESS**. Full evidence `pr51-later-hmm-collect-43179740278` and per-subject artifact `pr51-later-hmm-{4,8}-attempt-1`. This test was pinned in [predefined protocol](./PR51_FULL_RECORD_HMM_PREDEFINED_20261112.md) before any full-record HMM outcomes.

SHA256 verified original ERP CORE `.set`, `.fdt`, `events.tsv` for exactly the two already observed development participants; exact original 16 EEG channels, real VEOG, original response/RT and event mapping, original HMM reference, MATLAB R2026a and EEGLAB 2026.0.0. Source 300-second artifact flags reproduced **exactly**.

| Subject | Window | Correct/RT target | Correct/RT standard | HMM-unflagged target | HMM-unflagged standard | Both source conditions ≥5 |
|---|---:|---:|---:|---:|---:|---|
| 4 | Frozen 300 s | 18 | 95 | 4 | 3 | **NO** |
| 4 | Full 544 s | 30 | 159 | 6 | 4 | **NO** |
| 8 | Frozen 300 s | 22 | 95 | 5 | 20 | YES at source; original paired failure downstream |
| 8 | Full 460 s | 30 | 126 | 6 | 24 | **YES at source** |

**Decision:** Subject 4 is still source-support-ineligible under unchanged HMM reference and fixed minimum 5. It cannot be rescued by selection/weighting alone. Subject 8 has crossed only the *upper-bound source-reference* support floor with more real data; candidate-level effective HMM-unflagged support after rejection and weights, calibration and P300 waveform are not yet established. This experiment supports a separate, pre-specified **single development subject 8 full-pipeline test** before any new independent cohort.

The frozen original 20-subject ERP CORE verdict remains **15/20 valid pairs, PARTIAL**. No historical data pooling, threshold relaxation, production activation, new 14-person freeze or PR #51 merge.
