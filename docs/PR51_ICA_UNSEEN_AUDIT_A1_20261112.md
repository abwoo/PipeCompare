# Independent audit amendment A1 — unseen ICA seed scientific gate checks

Date: 2026-11-12. This document and its code branch were written **while the six-start real EEG run 43169562528 was still in progress, before any held-out seed 3–5 results were read**. No additional subject or held-out outcome is used for rule selection.

## Identified analysis defect

The original `benchmarks/validate_ica_unseen_seeds.py` reported a seed as passing when its original `Rank.status` was `feasible`, its gain-corrected P300 was finite and its effective HMM-unflagged condition supports were at least five. But it omitted explicit checks that (1) the independent seed's numerical objective was finite and (2) its original stratum matched the seed-0 baseline-selected stratum. Both requirements are already stated in the pre-specified scientific comparisons and the training-only `ica_feasibility_shadow.py` selector. A feasible status should not be used as a proxy for these separate comparability checks in an audit.

## Corrective rule (no outcome-dependent thresholds)

For each preselected candidate and seed index 3, 4, 5, `validAllThree` now requires all of:

1. Original `Rank.status == feasible`.
2. Same stratum as the original baseline-recommended recipe at seed 0.
3. Finite original numeric objective.
4. Correct candidate ID and seed identity, and finite gain-corrected P300 with positive finite gain.
5. Each condition's effective HMM-unflagged trial count at least **five**, without rounding or imputation.

Every failed condition is recorded in `failureReasons`. No candidate selection uses unseen seeds and no ranking, threshold, EEG source, sensor setup or MATLAB computation is changed.

## Execution and interpretation

The raw EEG workflow run `43169562528` remains an **unaltered historical execution** on branch `codex/pr51-ica-unseen-seeds-20261112` at SHA `5c9388562974c8305da95cb59e7249982788d589`. This corrective audit is isolated on `codex/pr51-ica-unseen-audit-20261112`; the completed raw per-subject artifacts may be re-analysed without refitting when all are present. Positive and negative results must be reported. A successful CI job does not establish scientific validity.

The original 20-subject ERP CORE result remains **15/20, PARTIAL**. No full 14-new-subject experiment, promotion, PR #51 merge or revision of original frozen gates is authorized by this amendment alone.
