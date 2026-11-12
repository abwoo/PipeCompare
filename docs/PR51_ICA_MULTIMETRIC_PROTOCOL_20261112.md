# PR #51 — multi-metric ICA decision rule: development-only locked protocol v1

Created on 2026-11-12 **before implementing the comparator**. Not a prospective validation claim: the three development subjects and aggregate six-seed science results were previously inspected. The untouched 14-person ERP CORE holdout remains sealed.

## Data and source identities

Use only *archived* raw six-start ERP CORE results from [GitHub Actions 43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528), SHA `5c9388562974c8305da95cb59e7249982788d589`, and opt-in ICA science readouts from [43184910782](https://github.com/abwoo/PipeCompare/actions/runs/43184910782). The original source P300 evaluation remains 15/20 valid matched pairs, **PARTIAL**. No EEG recomputation, no source-trial relabeling, no production Rank code changes and no PR #51 merge.

Subjects and **predeclared limited comparison set**, not a general search: subject 11 IDs [12,25], subject 15 IDs [5,2,6], subject 21 IDs [25,4]. The first ID in each list is the original seed-zero Standard recommendation. The alternative IDs came from earlier development; do not silently evaluate or claim to optimize other candidates.

Only seeds **0,1,2** enter candidate selection. Seeds **3,4,5** are *assessment-only*, but not honestly independent of the overall research program because their aggregate outcomes were previously viewed. They are an initialization stress test, **not** a new participant-level holdout.

## Fail-closed eligibility and hard constraints

An alternative can be considered only when each selection-start has:

1. The same original Rank scientific status `feasible` and baseline-reference stratum, a finite original objective, and finite positive gain with finite gain-corrected P300.
2. Both conditions' effective unflagged trial support **>=5.0**, both physical source-unflagged counts **>=5**, and `0 <= effectiveN <= physicalN`. Invalid/missing trials, identities or provenance cause abstention, never repair by imputation.
3. Successful, aligned science capture, complete ICA topographic maps, component matching, removed-component Jaccard, full Pz difference waveform, retained trial IDs and finite bootstrap P300 difference SE for every selection-start. Missing science information causes abstention; never treat missing as zero.
4. Minimum absolute sign-invariant matched topographic correlation at least **0.90** across seed-zero-to-alternative matched components in selection starts. This threshold is a conservative technical screen, **not** a validated neural component identity or medical threshold.

## Comparison relative to the unchanged Standard baseline

For the *same participant* and *same selection starts*, compare alternative A against baseline B:

- **Removal consistency:** minimum matched removed-IC Jaccard(A, seeds 1/2) **>=** minimum(B). Equality allowed.
- **Waveform fidelity:** maximum full-waveform RMS difference vs seed 0 across seeds 1/2 for A **<=** that for B (numerical comparison epsilon `1e-6 μV`). This measures repeatability, **not similarity to ground-truth evoked neural signal**.
- **P300 amplitude drift:** P300 max-minus-min over starts 0/1/2 for A **<=** that for B (epsilon `1e-6 μV`). Do not trade worse amplitude stability for better precision.
- **Precision:** mean separately bootstrap-estimated P300 difference SE across starts 0/1/2 for A **<=0.98** times B (at least **2% improvement** as a predeclared engineering threshold), and worst SE for A **<=** worst SE for B (epsilon `1e-6 μV`). This is an *engineering screen*, **not** a hypothesis test; no confidence interval or clinical superiority is inferred from three repeated fits.
- **No tie forced:** If more than one alternative passes, rank by lowest *worst original numerical Rank objective* over selection starts; if identical, abstain. If none passes, explicitly retain the baseline, with reasons. The baseline is never called a promoted algorithm.

No retrospective threshold changes if any case fails. If the original baseline lacks complete science or fails the existing Rank/support gates, the decision is **abstain** for that participant; do not select a challenger that bypasses a broken comparator.

## Independent-start stress report (not used in choice)

Report all seeds 3/4/5 for *both* original and chosen challenger, including Rank status, stratum, objective, positive gain, physical/effective support and science availability. Separately report P300 range, full-waveform drift, removal agreement and mean/worst bootstrap SE. A green GitHub workflow indicates software execution, **not** that the new selector improved physiology or is eligible for production.

## Stop / promotion gate

This developmental test can identify which candidate-level trade-offs warrant further research. To justify even *planning* prospective confirmation, require an unmodified Standard comparator, full evidence for 3/3 development subjects, **no selected subject with a worsening in any protected metric under the assessment-only starts**, plus a substantive improvement beyond numeric noise that is not confined to a single previously inspected participant. Even meeting that condition is not independent confirmation.

**Never automatically activate a new selector, promote a recipe, merge PR #51 or open the sealed 14-subject group based on this experiment alone.** Any future experiment needs its own code/protocol freeze, externally meaningful validation targets and a separate source roster that has not been repeatedly examined. Production remains unchanged.
