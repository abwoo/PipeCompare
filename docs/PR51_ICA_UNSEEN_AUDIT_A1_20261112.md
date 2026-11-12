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

## Additional reporting-leakage fix

The original six-start `icaSeedCheck.same` flag was copied into a truncated, training-only (0–2) `shadow()` input. This could label the three-start training recommendations unstable because seeds 3–5 disagreed, even though those starts are expressly held out. The amended inspector recomputes `same` on the first three recommendations before calling `shadow()`. This fixes reporting leakage, not the selected candidate's algorithm or the scientific outcome.

## Execution and interpretation

The raw EEG workflow run `43169562528` remains an **unaltered historical execution** on branch `codex/pr51-ica-unseen-seeds-20261112` at SHA `5c9388562974c8305da95cb59e7249982788d589`. This corrective audit is isolated on `codex/pr51-ica-unseen-audit-20261112`; the completed raw per-subject artifacts may be re-analysed without refitting when all are present. Positive and negative results must be reported. A successful CI job does not establish scientific validity.

The original 20-subject ERP CORE result remains **15/20, PARTIAL**. No full 14-new-subject experiment, promotion, PR #51 merge or revision of original frozen gates is authorized by this amendment alone.

## Completed archived-artifact audit (2026-11-12)

- Immutable original six-start execution: [run 43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528), SHA `5c9388562974c8305da95cb59e7249982788d589`, all three MATLAB subjects and collector **success**.
- Independent amended Python regression plus exact archived-artifact reanalysis: [run 43171200294](https://github.com/abwoo/PipeCompare/actions/runs/43171200294), **success**. All three raw result file SHA256 hashes and original run identities are captured in `pr51-ica-unseen-a1-reanalysis-43171200294`. No EEG refit, new subject, tuning or re-selection.
- Corrected validation results: **3/3 robust training-selected recipes remain scientifically feasible on all unseen starts 3–5**. The original baseline recipes also remain feasible in all three cases. Zero corrected feasibility verdict changes; no missing study evidence.
- Measurement-level results (Pz 300–600 ms, gain-corrected target-minus-standard, range over unseen starts only):

| ERP CORE subject | Original ID | Minimax ID | Original unseen range (μV) | Minimax unseen range (μV) | Minimax − original |
| --- | ---: | ---: | ---: | ---: | ---: |
| 11 | 12 | 25 | 0.000000 | 0.000000 | 0.000000 |
| 15 | 5 | 2 | 0.006436 | 0.014652 | +0.008216 |
| 21 | 25 | 4 | 0.005828 | 0.005828 | 0.000000 |

**Scientific decision:** All-seed feasibility supported on three previously observed development subjects, but **no evidence of better unseen-start P300 measurement stability**: 0 improved, 1 worsened, 2 tied. Seed-ID constancy by construction is not a valid gain endpoint. The data do not justify promotion or an automatic transition into the 14-subject confirmatory study. Preserve the original PARTIAL 15/20 PR51 result, leave PR51 open/unmerged and production defaults unchanged. The next experimental intervention must define and verify a genuinely informative stability or measurement-quality improvement before using the untouched 14-subject roster.
