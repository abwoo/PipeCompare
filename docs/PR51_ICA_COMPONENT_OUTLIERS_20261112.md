# PR #51 multi-metric v1 scientific stop and matched-component forensics

**Date:** 2026-11-12. Only archived ERP CORE development subjects 11, 15 and 21. No subject-level holdout, no new ICA fits, no production rank change and no PR merge.

## Immutable experiment and strict negative result

- [Pre-specified multi-metric v1 protocol](https://github.com/abwoo/PipeCompare/blob/codex/pr51-ica-multimetric-20261112/docs/PR51_ICA_MULTIMETRIC_PROTOCOL_20261112.md) was committed **before** outcome inspection. It requires original Rank feasibility and source-effective support >= 5 each condition, trustworthy component/ERP capture, minimum matched absolute scalp-map correlation >=0.90, no regression of matched removed-IC Jaccard, full Pz waveform RMS or P300 amplitude stability versus original, and >=2% mean bootstrap SE improvement without worse maximum SE.
- [Actions 43193549225](https://github.com/abwoo/PipeCompare/actions/runs/43193549225) **SUCCESS**, archived source/science evidence complete, all original recommendations preserved. **All 3 original baselines abstain** because their *lowest* matched scalp-map correlation under training ICA starts 1 and 2 is below 0.90. No alternative was promoted. Numerical threshold unchanged.
- Earlier [first workflow attempt 43193303017](https://github.com/abwoo/PipeCompare/actions/runs/43193303017) failed due to missing pinned NumPy/SciPy dependencies for an existing Python test; fixed before the successful run. This is a GitHub CI setup repair, not a scientific improvement.

## Component-level forensic diagnosis

A distinct *read-only* outlier investigation matched each candidate's actual mixing maps across ICA starts using channel-label alignment and sign-invariant Hungarian assignment, and cross-checked the original validated per-start minima, medians and removal Jaccard. [Actions 43194160598](https://github.com/abwoo/PipeCompare/actions/runs/43194160598) unit and real archived-artifact forensics **SUCCESS**. Complete per-component matched identities, topographic correlations, removal status and top ICLabel class are retained in that run's JSON artifact.

| Subject | Original ID | Weak matched components across seed 1 and 2 (counted per comparison) | Weak components kept by both fits | Weak with removal disagreement | Weak with top ICLabel class disagreement | Training minimum correlation |
|---|---:|---:|---:|---:|---:|---:|
| 11 | 12 | 4 | 4 | 0 | 0 | 0.661017 |
| 15 | 5 | 4 | 4 | 0 | 0 | 0.579231 |
| 21 | 25 | 1 | 1 | 0 | 1 | 0.799610 |

The same weak-map counts occur for other predeclared candidates of a participant because they share the fitted ICA decomposition. This is *not* evidence of independent subject replication.

**Important:** All weakly matched components are **retained in both fits**. Therefore a tempting rule that disregards weakly matched components simply because the removal decision agrees is scientifically unjustified: their neural-signal contribution is not known. In subject 21 one weak pair also changes its top ICLabel classification. No component index is treated as a persistent biological source without matching.

A high median absolute scalp-map correlation (~0.96–0.9997 in prior aggregate summaries) can coexist with a low minimum of ~0.58–0.80; averaging this away would conceal scientifically relevant instability. The forensic record does not establish that these components alter true P300 generation, only that the ICA decomposition is not sufficiently demonstrated stable by a per-component 0.90 criterion.

## Decision and next conditions

**Hold:** Multi-metric v1 is a fully documented rejected research candidate, not a methodologically supported replacement for Standard. Do not relax the 0.90 or trial support gate retrospectively to obtain a green selection outcome. Further development needs a *separate prospective recorder and protocol* quantifying each retained low-correlation component's time-locked ERP contribution and reconstruction sensitivity, with original Rank/ERP parity checks, before any topology-aware selector is scientifically motivated. Any changes must remain on isolated development branches.

The original ERP CORE 20-subject result is **15/20 valid pairs, PARTIAL**; PR #51 remains OPEN, UNMERGED; all 14 untouched participants remain sealed. GitHub Actions success indicates reproducible execution only, not improved algorithm accuracy or merger readiness.
