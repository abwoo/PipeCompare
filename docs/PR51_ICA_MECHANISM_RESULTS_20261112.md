# PR51 stage 3A — independent ICA mechanism and trial-support findings

*2026-11-12. Development-only study, unchanged 300 s ERP CORE P3 records 11, 15, 21. This is a read-only reanalysis of saved records; zero new EEG fits.*

## Provenance and correctness

The input artifact is the completed six-start real EEG run [43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528), SHA `5c9388562974c8305da95cb59e7249982788d589`. The unit-guarded archive reader was run in [43177240139](https://github.com/abwoo/PipeCompare/actions/runs/43177240139), SUCCESS, 6 dedicated ICA mechanism unit tests + 10 unchanged unseen-seed tests. Its `pr51-ica-mechanism-43177240139` artifact includes raw result file hashes, selected IDs, seed-0 removed components, complete six-start corrected P300, finite scientific statuses, per-seed unflagged effective support and explicit unavailable fields. A MATLAB/Python array adapter defect was discovered by the preceding failed CI run [43177102792](https://github.com/abwoo/PipeCompare/actions/runs/43177102792) and corrected before using the audit result; the failed run remains preserved.

## Actual seed-0 ICA treatment and six-start P300

| Subject | Pipeline | Recipe ID | Seed-0 ICs removed | Min effective unflagged target/standard across six ICA starts | Full six-start gain-corrected P300 range (μV) |
|---|---|---:|---:|---:|---:|
| 11 | original | 12 | 0 | 10 / 23 | 0.000000 |
| 11 | Minimax | 25 | 0 | 10 / 23 | 0.000000 |
| 15 | original | 5 | 2 | 21.842756 / 89.627934 | 0.007143 |
| 15 | Minimax | 2 | 2 | 24 / 100 | 0.016632 |
| 21 | original | 25 | 3 | 9 / 45 | 0.115261 |
| 21 | Minimax | 4 | 3 | 9 / 51 | 0.035692 |

All these recorded six-start scientific status checks pass on the original separate validation, and the exact five-trial effective support floor is not the binding constraint **for these three development cases**. Do not extrapolate that to original 20-person cases 4/8: Stage 2 found their reference HMM-unflagged support insufficient, and there is no proven threshold or regression fix. Do not lower floors, discard failures or relabel ocular proxies.

**Key confound:** Subject 11's original and Minimax pipelines both fit ICA but *removed no independent components at seed zero*. A zero P300 range under sample-order perturbation may therefore reflect unchanged processed EEG rather than reproducible removal of identified ocular/muscle sources. Baseline seed-0 count alone cannot certify what happened in all six seeds.

The complete six-start range above **includes the three training seeds** and does not represent independent-outcome confirmation. Across only unseen seeds 3–5 the original and Minimax P300 ranges were equal for 11 and 21, and Minimax was larger for 15 (see [corrected independent audit](https://github.com/abwoo/PipeCompare/blob/codex/pr51-ica-unseen-audit-20261112/docs/PR51_ICA_UNSEEN_AUDIT_A1_20261112.md)). Retain that result as the prespecified initialization-level comparison; do not retrospectively replace it with the apparently better six-start aggregate on subject 21.

## Scientific decision and missing evidence

The negative development outcome remains unchanged. The Minimax selector does not show improved *independent* initialization-level P300 stability. Current archives contain no seed-resolved ICA mixing matrices and no full cross-seed ERP waveforms. The next, separately pre-specified study will measure one-to-one sign-invariant ICA scalp-map matching and aligned removed-component overlap on these same 3 subjects, without touching the independent 14-person roster. Until those diagnostics and a genuinely useful pre-specified algorithm change are validated, **no new frozen experiment, default change, superiority claim or PR #51 merge** is justified.

Historical real P300 original protocol: **15/20 valid pairs, PARTIAL**, preserved.
