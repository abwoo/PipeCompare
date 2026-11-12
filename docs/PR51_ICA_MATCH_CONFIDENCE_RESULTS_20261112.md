# PR51 P1–P2 completed: mapping confidence is not the observed defect

Date 2026-11-12. Source: immutable six-start real-EEG science run [43184910782](https://github.com/abwoo/PipeCompare/actions/runs/43184910782). New, successful independent archive P1/P2 run [43195854240](https://github.com/abwoo/PipeCompare/actions/runs/43195854240), exact implementation SHA `068dfc5bfc3eb3acec7bac0e6a2e4fef95641321`. No new subject, source EEG refit, rank modification, changed thresholds or merge.

## P1/P2 forensic result

For each real ICA decomposition, the auditor verifies the previously published matched scalp-map minima and IC-removal agreement, then checks one-to-one Hungarian assignment against unconstrained, forward/backward nearest-neighbour matches, alternative-matching margins and 0→1→2 assignment cycle consistency.

| Subject | Weak matched pairs in training-start comparisons 1/2 (per original candidate) | Weak with **no** strong available counterpart | Forced weak matches caused by one-to-one assignment | Cycle consistency |
|---|---:|---:|---:|---:|
| 11 | 4 | 4 | 0 | 15/15 |
| 15 | 4 | 4 | 0 | 16/16 |
| 21 | 1 | 1 | 0 | 16/16 |

Original and exploratory candidate recipes share the same subject-level ICA decomposition; repeated counts across recipes are **not** independent replication. Some weak matches are uniquely mutual nearest neighbours despite low correlation; this does not turn them into good matches.

**Conclusion:** Replacing Hungarian one-to-one matching with unconstrained nearest-neighbour matching will **not** fix the existing weak-match findings on development starts 1/2. The mismatches are intrinsic to the observed mixing maps at the unchanged `|r|>=0.90` screen. The 100% three-start assignment cycle consistency is not proof that the ICA sources are biologically identical.

## Empirical ERP changes are not component-attributable

Full original-science Pz difference waveforms and original gain-corrected P300 were independently checked against the archived summary. Full-epoch RMS change from seed 0, for selected cases:

| Subject | Candidate | Seed1 RMS (μV) | Seed2 RMS (μV) |
|---|---:|---:|---:|
| 11 | Original 12 | 0.000000 | 0.000000 |
| 15 | Original 5 | 0.023456 | 0.019373 |
| 15 | Minimax 2 | 0.017555 | 0.014345 |
| 21 | Original 25 | 0.455729 | 0.447388 |
| 21 | Minimax 4 | 0.047400 | 0.075265 |

The sharp reduction of the *seed-anchored* waveform RMS in subject 21 candidate 4, despite sharing the same low-correlation ICA issue, is a useful hypothesis about downstream recipe decisions. It **does not** establish better EEG truth, prospectively generalize, isolate the influence of any retained ICA component or authorize promotion. Bootstrap SE changes and all per-seed good/bad evidence are stored in the run artifact `pr51-match-confidence-p1p2-43195854240`.

## Confirmed missing information

The prior science recorder saved ICA mixing matrices, ICLabel/removed-component identities and finished ERP waveforms, but **not each component's pre-removal, time-locked activation**. ICA mixing geometry alone cannot recover that activation, especially after `pop_subcomp` discards component contribution. Neither component-level causation nor a scientifically justified replacement topology gate is estimable from these archived outputs alone.

The next independent research intervention must add an **opt-in, read-only pre-IC-removal component-event signature** with trial IDs and sensor Pz projection, check stage-level reconstruction residuals and ensure original Rank/ERP parity. Only then evaluate whether low-correlation *retained* components materially contribute to stage-specific target-minus-standard P300; account for downstream filtering, rejection and weighting. Do not reinterpret a pre-removal effect as final ERP effect without an explicit linear/parity bridge.

**Verdict:** P1/P2 ambiguity diagnosis completed; no algorithm recommendation promoted. Multi-metric v1 abstention remains in force. PR #51 OPEN/UNMERGED; original 15/20 PARTIAL unchanged; remaining 14 independent participants sealed.
