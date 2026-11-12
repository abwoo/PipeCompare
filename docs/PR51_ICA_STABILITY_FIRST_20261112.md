# PR #51 stage 3B — physical trial-support integrity and stability-first ICA development comparison

Date: 2026-11-12. **Exploratory / already-observed subjects only**. This is not a retrospective amendment to the original 20-person frozen ERP CORE P3 study. Original result remains **15/20 matched pairs, PARTIAL**. PR #51 must remain open and unmerged.

## Execution identity

- Physical versus weighted-effective support integrity guard, independent branch `codex/pr51-ica-support-integrity-20261112`: [GitHub Actions run 43182116672](https://github.com/abwoo/PipeCompare/actions/runs/43182116672), unit and real archived-record audit **SUCCESS**. Original MATLAB result artifacts are from immutable six-start run [43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528). No source EEG refit.
- Stability-first shadow candidate comparator, branch `codex/pr51-ica-stability-first-20261112`: [GitHub Actions run 43182463157](https://github.com/abwoo/PipeCompare/actions/runs/43182463157), unit and development-subject comparison **SUCCESS**. Complete per-seed data preserved in artifact `pr51-ica-stability-first-dev-43182463157`. No source EEG refit.
- Neither branch edits the original production MATLAB algorithm, original source cohort, original 5-trial threshold, HMM labels, ICA initialization schedule or PR #51 head.

## Stage 1: confirmed implementation/data-semantic issue

The original ICA snapshots contain both (a) the physical count of source trials unflagged by reference HMM and (b) the Kish effective sample size of positively weighted, finite, unflagged trials. These quantities were not cross-checked in the shadow recommendation audit. An effective count cannot exceed the available physical unflagged count for nonnegative weights.

A new fail-closed read-only auditor inspects candidate ID and seed alignment, numerical finiteness, nonnegative effective support, integer source counts and the universal inequality `effectiveUnflagged <= sourceUnflagged`. It explicitly checks that any 5-effective-trial claim also has at least 5 physical unflagged trials, and does not impute missing support. It does not relabel any blink or rescue any missing trial.

**Real-data result:** all 348 candidate-start snapshots pass, with **zero integrity failures**: 120 subject 11, 114 subject 15, 114 subject 21. Therefore there is **no demonstrated source-versus-effective support count bug** in these three development subjects. This result does not prove all scientific trial-support decisions are correct, and cannot increase support in previously failing subjects 4 or 8.

## Stage 2: research-only ICA selection

A new read-only exploratory alternative uses **only seeds 0/1/2** and original common-feasible candidates. Select the candidate having the smallest training-set *gain-corrected P300 mean-amplitude range*; use the maximum original objective across training starts as a deterministic secondary tie-break; abstain if still tied. Seed 3/4/5 results are observed **after** selecting the recipe, are never inputs to its selection and must pass the same scientific eligibility, baseline-stratum, finite-objective, gain and effective/physical support gates. The original Standard and earlier Minimax choices remain unchanged.

This alternative is not truly prospective: the three subjects and their unseen-seed summary results were previously examined. It can identify follow-up hypotheses but **cannot** provide a new independent validation. Minimizing a scalar ERP range may overfit training initialization or select a low-amplitude, biased waveform; waveform morphology, IC topology, discarded components and actual P300 precision are not available from scalar snapshots.

| Subject | Original ID | Minimax ID | Stability-first ID | Original unseen range (μV) | Minimax unseen range (μV) | Stability-first unseen range (μV) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 11 | 12 | 25 | 25 | 0.000000 | 0.000000 | 0.000000 |
| 15 | 5 | 2 | 6 | 0.006436 | 0.014652 | 0.010213 |
| 21 | 25 | 4 | 4 | 0.005828 | 0.005828 | 0.005828 |

All three stability-first candidates remain scientific-gate feasible across the archived starts 3–5. Subject 15's unseen P300 range improves relative to Minimax by about 0.00444 μV **but is still 0.00378 μV worse than the original baseline**. The others tie. Training selection cannot be interpreted as independent candidate-ID stability, and these sub-0.02 μV fluctuations alone do not establish better ERP estimates.

## Decision and constraints for next research iteration

**Do not promote** Minimax or stability-first to the default recommender. Do not change previously frozen thresholds, original verdict, or PR #51 status. Do **not** run the untouched 14-participant holdout just to look for a winner among already-examined methods.

Before another true frozen evaluation, capture per-seed ICA component identities/weights and topographic correlations, removed-component decisions, full Pz difference waveforms, retained trial identities and bootstrap SME with source provenance. Verify that these can be captured without changing analysis decisions. Define new selection criteria and a full source/algorithm/protocol freeze **before accessing untouched participants**. The existing 3-subject development evidence has not demonstrated a scientifically meaningful overall improvement.
