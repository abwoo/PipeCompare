# PR #51 pre-IC-removal observer path-parity repair — scientific protocol

Date: 2026-11-13. This is a **bug fix to research-only scientific observability**, not a change to ICA fitting, class threshold, Rank, trial support, P300 outcome or default recommender.

## Root cause from immutable development artifacts

Original real-EEG contribution run [43197542023](https://github.com/abwoo/PipeCompare/actions/runs/43197542023) completed each participant's **original Rank/P300 parity**, but the new pre-removal observer did not finish for subjects 11 and 21. An independent GitHub *artifact-only* forensic run [43397552462](https://github.com/abwoo/PipeCompare/actions/runs/43397552462), without EEG refits, found that exactly two original seed-0 snapshots are incomplete:

- ERP CORE subject **11**, candidate **25**, seed **0**: `preRemovalContributions` is empty. Seeds 1–5 and candidate 12 across all six seeds contain complete records.
- ERP CORE subject **21**, candidate **25**, seed **0**: empty. Seeds 1–5 and candidate 4 across all six seeds contain complete records.
- Subject **15**, candidates **5, 2, 6**, all six seeds: complete.

The source is an **execution path coverage bug**, not unsupported ERP data: seed-zero candidates created after the main search by **filter-cutoff refinement and/or eye-channel regression** are executed in `fineTuneCutoffs` or `eogCheck`, which originally skipped the observer placed in `runSubtree`. Subsequent alternate ICA seeds used `rerunWithIca`, which already captured before component removal. It is invalid to "repair" absent seed-zero evidence with a zero, another seed's value or post-removal ICA geometry.

## Repair and unchanged scientific requirements

Unify the opt-in `passivePreRemovalCapture` call at the actual step *before* `icremove`, and retain its result separately in each execution path:
1. ordinary `runSubtree`,
2. shared-trunk `runParallel`,
3. original post-search `fineTuneCutoffs` rerun,
4. original `eogCheck` eye-regression variant,
5. alternate-seed `rerunWithIca`.

The optional flag remains `scienceCapture` (false by default). On observer failure, the original EEG candidate's ICA removal and Rank still execute identically, while science evidence is marked failed and the workflow must fail closed. The existing pre-removal observer must still prove `sum(all original ICA component Pz projections) + residual = actual sensor Pz` for each original urevent, including nonzero residual for rank-deficient decompositions. Final retained trial IDs/weights are joined without inventing new epochs. No component truth or final ERP causality claim.

## Acceptance gates

- The unchanged old source six-ICA-start archive (run 43169562528) must match **all** new real-EEG scientific Rank/strata/score estimates, gain, trial support, original ERP P300 and seed alignment.
- The new record must contain exactly **one complete** pre-removal contribution source for all six seeds of *all predeclared candidates* on development subjects 11/15/21. A failure in any expected record invalidates the claim; no trial floor/ICA match threshold relaxation.
- MATLAB and Python isolated tests must pass. The reusable GitHub real ERP workflow is cloud-only and its raw failure and successful artifacts must remain archived.
- Neither a fixed observer nor the previously inspected performance differences authorize an ICA selector promotion, PR #51 merge, fresh 14-person holdout or scientific accuracy claim. Original ERP CORE result **15/20 PARTIAL** unchanged.
