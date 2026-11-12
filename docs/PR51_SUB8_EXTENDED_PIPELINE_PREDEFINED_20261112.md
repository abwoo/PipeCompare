# PR51 Stage 2C — exploratory full 460-second P300 complete-pipeline test, subject 8

**Prepared after the upstream HMM support result but before any 460-second Standard or fixed-arm search outcomes.** This is development-only, not an independent confirmation cohort and not a revision to the frozen original 20-person study.

## Frozen scientific question

Original 300-second ERP CORE subject 8 had source HMM unflagged 5 target/20 standard trials and a fixed comparator with only 4/18 effective unflagged, so the original 300 s pair was scientifically invalid. The prior fully independent source HMM calculation [run 43179740278](https://github.com/abwoo/PipeCompare/actions/runs/43179740278), using the same HMM, found **6/24 source HMM unflagged** in a full **460-second** original recording. This is only a prerequisite; now test whether the full **actual** candidate pipelines (with filtering, adaptive search, ICA, rejection, weighting and reference correction) preserve at least five effective HMM-unflagged trials in **each** condition.

## Pinned experiment before its outcomes

- Cohort is exactly subject 8, already used for development. Same source .set/.fdt/events SHA256 as frozen original ERP CORE protocol, same original correct-response/RT event rules and 16 EEG + true VEOG channels, 256 Hz, original Pz target-minus-standard 300–600 ms P300 estimand.
- Duration **460 seconds**, selected solely because it is the complete recording already enumerated in SHA-pinned task-event inventory, not based on its P300 result. Preserve and report source HMM baseline result **5/20 (300 s)** and full-reference HMM **6/24 (460 s)**; verify source hashes match the separately completed HMM run.
- Preserve original Standard search behavior, all candidate options, Rank scientific gates, ICLabel, HMM, thresholds, measured gain, original `Measure.m`, no retrospective tuning. Use an isolated copy `benchmarks/pr51_extended_pilot.m` derived from `benchmarks/erpcore_pilot.m`; the original frozen pilot and protocol **are not edited**.
- Preserve original two fixed comparator recipes Minimal and Common ICA; run exactly the 460 s same-source recording and save every scientific score and trial-level support.
- Stage 2C **paired measurement eligibility** demands that **all three arms** (Standard, Minimal, Common ICA) have finite gain-corrected P300, positive finite gain, finite reported difference SE, aligned and positive-weighted per-trial scores, and **effective HMM-unflagged n≥5 in both target and standard conditions**. If any arm fails, report ineligibility; do not impute, drop failed arms or relax the floor.
- Independently recompute 1,000 within-condition trial-resampling bootstrap draws with the existing ERP CORE collector function, and check bootstrap weighted point amplitude reproduces MATLAB. Report three-arm bootstrap SE and point estimates, *only if scientific evidence is valid*. This is a descriptive one-person comparison only — no participant-level confidence interval, accuracy superiority or inference of other subjects.
- Record source and algorithm digests, Matlab/EEGLAB versions, replay and executable export error; preserve failures and run artifact; no automatic rerun with new thresholds.

## Go/no-go

A successful 460 s three-arm result would only show a technically feasible development example where longer real acquisition restored paired measurement support without altering old scientific criteria. It **does not repair the original frozen 300 s** result. Subject 4 remains a negative source case (only 6/4 unflagged at 544 s) and must stay scientifically ineligible. The frozen ERP CORE report remains **15/20 valid paired records, PARTIAL**; **PR #51 stays OPEN and unmerged**.

No new experiment with the untouched 14 subjects, no production algorithm switch, and no new general-superiority claim is authorized by this test.
