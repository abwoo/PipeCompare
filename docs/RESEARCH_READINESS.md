# Methodology audit and research readiness

Audit baseline: main commit `7151821` (merged PR #49), 23 October 2026.
All 35 compatibility jobs of that preceding PR passed. This is implementation evidence,
not an independent evaluation of scientific accuracy. This document describes the next audit;
its numerical regressions do not establish that the entire method is free of problems.

## What the method estimates

The chain is contract/ROI and candidate plan → native EEGLAB transformations → trial measures
and matched-decision signal checks → feasibility → paired score bootstrap → factorial
empirical Bayes and recommendation → optional group empirical Bayes → adoption/script replay.

A recording's adaptive preprocessing can legitimately use that recording. The evaluation target
must match the claim:

- **Denoising these observed trials:** compare the final readout with the realized clean-trial
  truth in an independent simulation. Refit the entire algorithm separately in each simulated
  recording; the truth must never enter its decisions.
- **Estimating a population ERP:** compare with the population truth and include biological
  trial variability in the error. This differs from the realized-trial denoising target.
- **Predicting new trials:** training decisions and test trials must be separated, including
  filter/overlap buffers. Resampling stored scores cannot establish this claim.
- **Generalizing group priors to new participants:** leave entire participants out of prior,
  hyperparameter and development fitting. A held-out participant's permitted adaptation must
  be specified before evaluation.
- **Raw adopted-data ERP versus diagnostic readout:** ranking calibration uses estimate/gain.
  Adoption does not divide the entire dataset by that gain. Raw and corrected errors are
  exported separately; neither target can stand in for the other.

## Findings and repairs

| Priority | Finding | Repair or remaining requirement |
|---|---|---|
| P0 | Full-data fitted HMM/ICA, candidate eligibility, templates, cutoffs and shrinkage remain fixed in the score bootstrap and split check. | Ranking and split outputs explicitly say conditional and not fully held out. End-to-end independent evaluation remains necessary. |
| P0 | Four exploratory simulations do not establish stable calibration. The original 1.21 ratio predates changes to the error definition and cannot validate this revision. | Re-run frozen code on new simulations with recorded seed ranges and record-level intervals, then test independent real data. |
| P1 | Weighted artifact penalties changed when all weights were multiplied by a constant. Reliability weights were also mixed with physical bootstrap multiplicities. | Use squared weights per physical observation and a weighted variance correction; pass bootstrap counts separately. Test scale invariance and explicit replicated-sample equivalence. |
| P1 | Group shrinkage discarded paired covariance. Shared recording noise could be counted as pipeline-specific uncertainty and erase real individual differences. | Retain paired covariance, subtract projected observation noise in the interaction moment estimate, use GLS effects and covariance-aware posterior gains. |
| P1 | Missing candidate patterns can disconnect the recording/pipeline design. The old degrees of freedom assumed one connected component. | Use the actual design rank; report components and whether interaction variance is identifiable. Do not fill missing candidates with invented estimates. |
| P1 | A nonfinite eye-channel sample could poison all moving features; median-filled HMM observations could also be treated as observed quiet evidence. | Omit missing samples in eye features, marginalize missing emissions, fit emission parameters on finite features and never report missing samples as detected events. Detection inside gaps remains unassessed. |
| P2 | One-sided multitaper PSD doubled DC and the even-length Nyquist bin. | Leave those endpoints undoubled. Ordinary positive interior bands retain their normalization. |
| P1 | RMS individual total error was described as uncertainty of a group mean. Systematic biases may share a sign across subjects. | State the distinction explicitly. A group-contrast claim still requires null and non-null group experiments. |

The first focused diagnostic contained three direct numerical counterexamples (weight scaling,
disconnected designs, spectral endpoints) and three new API/behavior checks for multiplicity
and covariance handling. They failed before the repair and passed after it. Additional
integration regressions exercise the actual Rank → Group path, weighted bootstrap and HMM
missing observations. They are not six independent demonstrations of scientific validity.

## Assumptions that code repairs do not remove

1. **Artifact association is not causal artifact bias.** The difference between flagged and clean
   trials can contain genuine neural variation related to blinking, movement or attention.
   Soft thresholding intentionally discounts weak differences and may understate true bias.
   A zero estimate means no detected contribution under the detector, not proof of no artifact.
2. **Signal injection is conditional.** A noise-free template replays the fitted decisions. It
   measures transfer of that shape under those decisions; it does not test how refitting would
   react to a different signal, participant or artifact distribution.
3. **Group empirical Bayes is a plug-in model.** Covariance, effects and interaction variance are
   estimated, with delta-method log covariance. Their uncertainty, selective candidate failure,
   non-Gaussian errors and few-participant behavior still need evaluation. Correcting covariance
   does not prove a universal improvement over independent selection.
4. **Dependence and filtering matter.** ERP trials retain the independent-trial approximation;
   consecutive spectral segments use HAC/block estimates with finite chosen lags. Long-memory,
   overlap, drift and effects shared across trials can invalidate those approximations.
5. **Interpolation changes the spatial operator over time.** Exact replay does not prove harmless
   boundary transients, preserved local ERP/topography or preserved spectra. Rank-deficient
   recordings can leave stretch detection unassessed and must be counted in coverage reports.
6. **Mean log PSD is a specific estimand.** It is not integrated band power, nor log of the
   across-segment arithmetic mean. Taper smoothing can move edge/out-of-band power into a narrow
   band. Effective resolution metadata is not evidence that this bias is negligible.
7. **Family-wise inference is conditional and approximate.** An alpha tuned in limited simulations
   cannot establish uniform coverage over arbitrary adaptive candidate families. Selection
   frequency is not a probability that a pipeline is scientifically correct.

## Minimum validation package before a strong methods claim

Freeze the code revision, candidate family, contract/ROI/time windows, comparator settings,
simulation families, exclusion rules and stopping/precision criteria before seeing validation
outcomes. Reserve new seeds and participants; do not retune on them and keep calling them held out.

| Experiment | Minimum design | Primary outcomes |
|---|---|---|
| Estimand alignment and calibration | Report raw and gain-corrected readouts against both population and realized-trial truths. Start with an independent pilot of at least 30 recordings and extend to the predeclared interval precision; 30 is not a success threshold. | Selected-pipeline RMSE/bias; actual/reported RMS ratio and record-bootstrap CI; missing recommendation rate; failed/unassessed rate. |
| Selection and shrinkage ablation | Identical candidates for raw minimum, factorial shrinkage, group shrinkage and fixed baselines. Include correlated errors, genuine interactions, unequal precision and selective failures. Oracle ranking is a simulation diagnostic only. | Paired selected error and regret with uncertainty, harmful-change rate, subgroup performance and computational cost. |
| End-to-end independent simulations | Re-run the complete algorithm on new recordings and at least one additional generator family, including heavy tails, drift, overlapping trials and neural/artifact correlation. | Calibration and selected error after all adaptive fits, not just agreement between stored objectives. |
| Group generalization | Participant-level outer splits; train priors without the held-out person. Compare common fixed pipeline and per-record independent adaptation. Include null and real group differences. | Held-out error, group-effect bias and false-positive rate; preservation of individual differences. |
| HMM event validation | Held-out annotated eye events plus synthetic weak/dense/clipped/drifting events and missing-data controls. Keep detector thresholds frozen. | Event precision/recall, sample overlap, timing error, false-positive seconds/minute and unassessed duration. |
| Bad-stretch repair | Fault durations 1/5/10/20 s; isolated/adjacent channels; reference/rank conditions; include boundaries and genuine local neural transients. | Detection coverage, repair error, boundary discontinuity, spectral/ERP bias and replay equality. |
| Spectral validation | White and 1/f backgrounds; center/edge/out-of-band peaks; several durations/bandwidths; matched-resolution references. | Bias and RMSE of the declared log-PSD target as well as variance, with resolution and DC/Nyquist controls. |
| External real data | Discovery/development data and a separate participant cohort or site. Prefer more than one ERP task; treat spectral claims separately if not validated. | Independent measurement precision and signal/condition-effect preservation; test-retest evidence where available. Reliability alone is insufficient. |

Plan sample size from precision and paired effect uncertainty, not the number of candidate rows.
For illustration only, with independent Gaussian errors and known common reported variance,
the RMS ratio has approximate relative SE 1/sqrt(2n); a 95% half-width of 10% needs roughly
193 independent errors, and 5% roughly 769. These are not guarantees for the clustered,
adaptive and heterogeneous benchmark. Its record-level pilot should determine the actual budget.

For a null false-positive rate, a few dozen repetitions are also inadequate: at p=0.05,
1,000 independent replications have a Monte Carlo SE of about 0.0069. Use intervals and
predeclared tolerances rather than interpreting one favorable rejection-rate estimate.

## Paper framing

The promising question is whether adaptive EEG pipeline selection can reduce measurement error
while preserving signal and reporting error that remains calibrated after selection. A GUI,
HMM implementation, multitaper estimator or combination of established methods alone does not
establish that contribution. Paired covariance and GLS are established statistical ideas; the
repair should not be marketed as a new general statistical invention.

A compact paper should choose one decisive claim (for example, selection-aware calibration
of mean-amplitude ERP pipelines), with a matched ablation, independent end-to-end benchmark,
external real-data example and clear failure regimes. A broader paper can add group adaptation,
multiple ERP tasks and spectral measurements after those claims are supported.

Suggested evidence/figure sequence:
1. Estimator target and the fitting/evaluation boundary.
2. Reproduced failure modes and matched ablation of selection/shrinkage.
3. Independent selected-pipeline accuracy and calibration across stress families.
4. External data, group transfer and preservation of meaningful effects.

Nature Methods explicitly accepts concise method/tool developments as Brief Communications;
the shorter format still needs a strong contribution. Nature Communications seeks substantial
advances important to specialists. Cell Reports Methods emphasizes robust, reproducible methods.
Suitability cannot be inferred from passing tests, journal branding or manuscript length.
No acceptance probability or claim of first-in-field novelty is established by this audit.

A focused EEG methods/software paper is a practical nearer-term route after independent
validation; Journal of Neural Engineering and Medical & Biological Engineering & Computing
are relevant candidates to evaluate against the completed study. A Nature/Cell/Science-family
submission should follow the strength of the resulting scientific evidence.

## Sources for positioning

- Nature Methods, content types: https://www.nature.com/nmeth/content
- Nature Communications, aims and scope: https://www.nature.com/ncomms/aims
- Cell Reports Methods, methods criteria: https://www.cell.com/cp/special-issues/call-for-papers/new-approach-methodologies
- Kessler et al. (2025), How EEG preprocessing shapes decoding performance:
  https://www.nature.com/articles/s42003-025-08464-3
- Delorme (2023), EEG is better left alone:
  https://www.nature.com/articles/s41598-023-27528-0
- Brookshire et al. (2024), Data leakage in deep learning studies of translational EEG:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11099244/

These are positioning and validation references; they do not constitute an exhaustive novelty
search, nor do decoding studies validate PipeCompare's ERP error estimator.
