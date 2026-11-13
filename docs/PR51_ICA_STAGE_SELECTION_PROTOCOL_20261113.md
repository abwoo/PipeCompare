# PR51 P1c — frozen selection/weight vs stage-signal arithmetic decomposition

Date: 2026-11-13. Predeclared **before implementing/running P1c analysis**; three development subjects and earlier candidate-level outcomes were already seen. This is a *retrospective, descriptive diagnostic*, not an independent hypothesis test.

## Source identity
Use **only** complete pre-removal six-start science records from GitHub real-EEG [43397978997](https://github.com/abwoo/PipeCompare/actions/runs/43397978997) (all 6 jobs SUCCESS, original Rank/P300 parity), and original six-start frozen control 43169562528. Subjects {11,15,21}; candidate rosters: 11:[12,25], 15:[5,2,6], 21:[25,4]. Zero new participants, EEG fits or changes to ICA threshold/Rank/trial selection.

## Question
For the **actual pre-IC-removal Pz event-level sensor signals** (already captured, before deleting ICA components), explain arithmetic variation across ICA initialization s=1..5 relative to seed0.

Let `Y_s(i)` be the observed baseline-corrected 300–600ms *pre-removal Pz sensor mean* at original urevent ID `i`; let `F_s` be each candidate's **actual final retained event IDs and corresponding original positive weights**, from the separate verified final science record. Define `D(Y,F)` as weighted target(11) minus standard(31) Pz *at this pre-removal stage*, never the final gain-corrected ERP estimand.

With seed 0 as the **predeclared fixed reference**:
- Total stage change `T_s = D(Y_s,F_s) - D(Y_0,F_0)`.
- Frozen-trial/weight **stage-signal change** `S_s = D(Y_s,F_0) - D(Y_0,F_0)`.
- Remaining **selection/weight arithmetic effect** `Q_s = D(Y_s,F_s) - D(Y_s,F_0)`.
- Verify `T_s = S_s + Q_s` to within 1e-6 μV independently for both target and standard conditions, and verify D(Y_s,F_s) against the previously independently checked per-candidate stage sensor Pz difference from `ica_contribution_validation.json`.

Interpretation: this algebraic split is order-dependent (seed0 reference and counterfactual *weights applied to the **same measured** pre-removal events*), does not identify a causal effect of ICA or rejection, and does not establish neural ERP ground truth. It is explicitly different from the final gain-corrected P300 endpoint, since later filtering/rejection/calibration may change signals. Do not report a `Q_s` percentage when total is near zero.

## Fail-closed eligibility
Each original event ID must be unambiguous, mapped to exactly one contract condition, and carry finite sensor Pz amplitude. All IDs in both `F_0` and `F_s` must exist in **both** pre-removal stage sources, with the same experimental condition and positive finite original weights in the science snapshot. A missing common event, duplicated ID, wrong condition, missing weight, inconsistent stage Pz/ICA reconstruction or changed original source parity terminates the case with an explicit error; no zero-filling and no trial floor reduction.

For audit completeness report separately:
- total, frozen-trial signal, and selection/weight effect in μV and corresponding target/standard components;
- whether original final retained IDs and original weights match;
- original final P300 drift and full Pz waveform RMS (a *different estimand*);
- stage ICA component rematching/weak components from the prior fixed method **without** replacing it or altering the original 0.90 diagnostic.
- original Rank/strata/trial support unchanged, original frozen 15/20 PARTIAL result and 14 untouched subject holdout sealed.

## Decision
Classify a case as `selection-or-weight-explains-stage-change` **only as descriptive arithmetic** if |S_s| is close to numeric reconstruction precision; otherwise label mixed. There is **no efficacy target or threshold fitted to this cohort**, no hypothesis test and no automatic selector promotion. Even if all development drift is selection-driven, this does not prove removing components is safe or that a different threshold is superior.

**PR #51 remains OPEN/UNMERGED; the Standard candidate ranking is not modified.**
