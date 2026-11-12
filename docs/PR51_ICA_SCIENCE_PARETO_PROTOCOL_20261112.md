# PR #51 exploratory science-aware Pareto ICA selection — analysis protocol

Date: 2026-11-12. Research branch `codex/pr51-ica-science-pareto-20261112`. This is **NOT** prospective participant-level validation. The three ERP CORE participants 11/15/21 and a subset of their prior ICA/P300 findings have already been examined. New selection rules written before running this specific multi-candidate Pareto analysis do not make the old cohort independent.

**Original study:** 15/20 matched valid ERP CORE P3 records, verdict PARTIAL, unchanged; PR #51 OPEN/unmerged. **Archived six-start scientific recorder:** [run 43184910782](https://github.com/abwoo/PipeCompare/actions/runs/43184910782), original six-start control [run 43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528). All evidence re-used in GitHub; no new EEG fit or participant.

## Correctness gates first

1. Original `Rank.status` scientific feasibility in **every training seed** 0/1/2, common original baseline-selected measurement stratum, finite unmodified original objective and gain, Pz gain-corrected target-minus-standard measurement and at least **five** HMM-unflagged effective trials per condition on each training start.
2. Trial provenance must be literal: the recorder's reference IDs, HMM flags, finite per-trial score/weight rows and recomputed physical/effective unflagged support must match the originally saved candidate/seed counts. Neither source unflagged trials nor weights are guessed; discrepancies fail closed.
3. Across training seeds only, ICA component topographies must be matchable by labeled montage and sign-invariant Hungarian matching. Missing/degenerate component maps, missing complete Pz waveform or bootstrap SE, unresolvable class/removed decisions and changes of component dimensionality render that candidate ineligible rather than assigning a favorable score.

## Four-axis, no-arbitrary-weight selection rule

On the complete candidates satisfying the original and provenance gates, compare these **training-start-only** axes:

- minimize worst (maximum) 1,000-draw trial-bootstrap P300 difference SE across starts 0–2;
- minimize maximum full **gain-corrected Pz difference-wave RMS** drift relative to start 0;
- maximize minimum matched component **median absolute topographic correlation** across starts 1/2;
- maximize minimum matched removed-component **Jaccard** across starts 1/2.

There is no post-hoc rescaling or learned weighted sum. A candidate Pareto-dominates another only if it is no worse on **all four** axes and strictly better on at least one, with a frozen numerical comparison epsilon of `1e-12` for floating-point equality. If exactly one nondominated candidate exists, emit its ID in a **shadow readout**. If zero complete candidates or multiple nondominated candidates remain, emit explicit abstention. Do not silently fall back to a convenient candidate.

The original, Minimax and stability-first recommendations remain unchanged. *No new default recommender is created by this experiment.*

## Examination of ICA seeds 3–5

Only after training choice has been computed, inspect the chosen recipe on unseen-in-selection starts 3/4/5 for all original scientific constraints, fully recorded ICA/waveforms, P300 range and bootstrap SME. New seeds are not used to pick or break a Pareto tie. Outcomes are for development analysis, **not untouched subject confirmation**: earlier studies have already exposed some seed-specific summary information.

Preserve all failures, eligibility losses, Pareto frontier members and abstentions. Any apparent improvement is hypothesis generation only; no promotion of a new method, no threshold loosening, no PR merge and **no use of the untouched 14-participant ERP CORE holdout** until genuinely prospective acceptance criteria and a new source-algorithm-protocol freeze are justified.
