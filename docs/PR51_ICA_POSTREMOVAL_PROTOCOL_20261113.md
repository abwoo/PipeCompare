# PR51 P2: ICA removal-step sensor bridge — prospective recorder protocol, already-explored development cohort

**Declared 2026-11-13 before first code change or new real EEG execution of this observer.** This experiment addresses a precise missing stage in the existing provenance chain: a pre-`icremove` IC topography/activation capture and a final P300 ERP are available, but no actual sensor Pz record **immediately after the original `icremove` step**.

## Fixed participants, source, algorithms and boundaries

Reuse only originally inspected ERP CORE P3 development subjects 11, 15, 21, and exactly six ICA starts 0–5, original unmodified Standard plus same predeclared candidates 11:[12,25], 15:[5,2,6], 21:[25,4]. Entire 14-participant subject-level holdout stays unopened; no new exploratory participant introduced. Preserve unchanged original Rank/Measure source, trials ≥5 per condition, full ERP and bootstrap evaluations, original ICA starts, treatment recipe and calibration. No PR #51 merge.

The **sole source modification** is passive, gated `scienceCapture` research observer code inside the existing execution paths, never the actual ICA fitting/cleaning or downstream weights/score. No observer result can be supplied to ICA removal, preprocessing, ranking or candidate recommendation.

## Capture immediately after real IC removal, for each path

- After the actual original `runStep(...)` of `in.type=='icremove'` (not a rerun or synthesized prediction), take a private copy of the resulting EEG in-memory at that step, *before next filtering, referencing, epoching, rejection, or weighting step*.
- Use the exact same `Contract` event codes 11/31, reference `urevent` IDs, Pz channel, baseline -0.2..0 s and 300–600 ms windows that the validated original pre-removal observer uses, and record a per-event baseline-corrected Pz window mean in μV.
- Preserve per-event IDs and experimental conditions, Pz channel, sampling rate and original stage indicator. Missing/ambiguous IDs, absent Pz, duplicates, nonfinite sensor signal or out-of-range windows must fail the **observer**, not the original candidate. Such scientific errors fail the end-to-end validation instead of silently skipping.
- Check pre and post have the exact same time-locked original event and condition identities. No new event or trial is permitted between the immediate pre/post stages.

## Independent stage-level linear bridge

For each event, verify the **actual measured post-removal Pz** follows the pre-removal sensor Pz minus the **actually deleted** components' pre-removal Pz projected signal, using the original archived ICLabel/removed IDs and validated pre-removal `A*W` component activations. Record per-event residual, RMS and maximum difference and sign. No reconstruction assumption imposed where rank-reduction/EEGLAB implementation might differ: if the bridge fails numerical tolerance (atol ≤ 1e-4 μV, rtol ≤ 1e-6), **fail closed**, archive both data streams and diagnose, not relax tolerance after observation. Do not assert that this step-level linear relationship explains the final ERP after subsequent nonlinear operations.

For each unchanged candidate/seed, also join original final retained trial IDs and original positive weights to both measured pre/post streams. Report pre/post target-minus-standard difference (µV) on the exact same final trials and original weights, with a stage Δ computed as post minus pre. Explicitly distinguish that from the final **gain-corrected** P300 ERP estimand; never equate them without separate calibration/weight/source checks.

## Tests and control experiment

- MATLAB unit: exact Pz and event ID capture, no mutation, missing Pz/duplicate IDs failure, pre/post synthetic difference, unchanged original scientific capture regressions.
- Python scientific unit: positive/negative component subtraction, reordering, empty removed set, rank-deficient residual, missing montage/event ID, missing seed0, changed final-trial support; reject unverified source parity. First numeric tolerance is frozen here.
- Cloud real MATLAB ERP CORE rerun: reuse unchanged historical six-start control [43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528), compare original entire Rank/strata/objectives/calibrated P3 source, trial support, full ERP, SE, and new source bridge; archive failures too.
- Only after all three existing development participants pass may we run a retrospective **mechanism diagnosis** to decide whether to prioritize ICA removals, later filtering/referencing, or final trial/weight changes. No efficacy study or confirmation claim without separate independent sample.

## Fail/no-promotion disposition

If all tests pass, the scientific advance is a reproducible, event-aligned **quantitative IC-deletion-to-sensor-signal bridge**, not evidence that the multimetric selector improves scientific outcomes. Multi-metric v1 remains ABSTAIN and historical 20-subject result 15/20 PARTIAL; PR #51 OPEN/UNMERGED, 14-subject holdout sealed.
