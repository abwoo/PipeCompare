# PR51 P1b implementation checkpoint — passive ICA component event signatures

Branch: `codex/pr51-ica-contribution-observer-20261112`. This checkpoint is an implementation record, **not** a scientific performance success declaration.

## What has been implemented

- `+pipecompare/+research/IcaContribution.m`: reads actual in-memory pre-`pop_subcomp` ICA weights/sphere and event-locked Pz signal, records component × urevent P300-window mean amplitude (same research contract baseline) plus observed pre-removal Pz and reconstruction residual. Does not mutate the input EEG. Explicit rank-reduced residual is preserved.
- `+pipecompare/+run/Executor.m`: opt-in observer only under existing `scienceCapture` flag, before the actual IC removal. Observer failure becomes a science error object, never a candidate-rank/cleaning rejection. Applies also to exploratory refined ICA replays.
- `+pipecompare/+research/IcaScience.m`: appends pre-removal signatures to existing full-waveform, ICLabel and per-trial scientific records, leaving original measurements unchanged.
- `benchmarks/validate_ica_contribution.py`: requires per-event component-sum + residual = measured stage Pz, checks ICA rank, source/target condition IDs, and joins only the candidate's **actually retained original final urevents and positive weights**. Preserves the difference between stage-specific contribution and actual final P300 metric.

## Cloud checks

- [Standalone MATLAB scientific unit #43197107982](https://github.com/abwoo/PipeCompare/actions/runs/43197107982): **SUCCESS**. Covers original EEG immutability, full-rank exact projection, nonzero reduced-rank reconstruction residual, absent Pz and duplicate urevent fail-closed.
- [Full GitHub three-participant real-EEG re-run #43197542023](https://github.com/abwoo/PipeCompare/actions/runs/43197542023): triggered on the code/manifest above for the already-observed subjects 11, 15 and 21 only. The job's **live final science status must be checked**; setup and unit success alone do not establish the new recorder works on real ERP data. Its extra scientific validator explicitly preserves original six-start Rank/P300 identity against historical run 43169562528 and archives both successes and failures.

## Limitations and stop gates

The per-component Pz signature is measured *at the actual pre-removal stage*. Even when weighted using the eventual unchanged retained trial IDs, it is **not an estimate of the final counterfactual P300 if an IC is or is not removed**. Post-removal filtering, referencing, event exclusion, temporal preprocessing and nonlinear decisions may intervene. Treat the stage's residual and source provenance as first-class quality metrics, not as nuisance terms to zero.

Original scientific trial threshold, 15/20 PARTIAL verdict, multi-metric v1 ABSTAIN and PR #51 OPEN/UNMERGED do not change. The 14 untouched participants remain sealed. No new algorithm may be promoted unless the completed scientific stage results and future *genuinely independent* validation justify it.
