# Validation scope and remaining experiments

The 0.9.5 first benchmark covered four simulated recordings. Its actual/reported RMS-error ratio
of 1.21 was exploratory. Candidate pipelines and conditions from one recording are correlated
rows; their count is not the number of independent observations. That first ratio must not be
used as a universal calibration multiplier or as proof of stability.

## What the checks establish

- Signal injection checks the effect of the fitted transformations with the observed decisions
  replayed. It does not test how those decisions would change in a new recording.
- `ranking.byStratum.crossfit` is a conditional selection diagnostic (`fullyHeldOut = false`).
  It reselects some stored threshold/component options. Bad channels, HMMs, ICA, trial weights,
  fine-tuned cutoffs, candidate eligibility, signal checks and factorial shrinkage still use all
  trials. The legacy `honest` field is retained for compatibility; its name does not imply a
  fully independent evaluation.
- `Group.own` estimates choices for the observed group. Its two-level model takes raw objectives,
  not the already-shrunk ranking column. Sharing the observed group's data is intentional;
  performance on unseen participants remains a separate question. Group uncertainty resamples
  each recording's trials, not new participants.
- Continuous mean-PSD scores use Bartlett HAC standard errors with `floor(sqrt(T))` lags and a
  moving-block bootstrap for comparison. Neither method guarantees calibration under arbitrary
  nonstationarity or very long dependence. ERP trials retain the independent-trial assumption.
- `Measure.eyeBias` returns NaN when no artifact-free retained trial remains. This is lack of
  identifiability, rather than evidence that artifact bias vanished. No detected artifact returns
  zero conditional on the detector; sensitivity and specificity still need validation.
- A bad-stretch diagnostic with `assessed = false` is different from a successful check finding
  no bad stretches. Rank-deficient inputs are not silently made estimable by a pseudoinverse.

## Calibration outputs

`Calibration.rows` uses the same per-condition noise, artifact bias and shape bias as the ranking.
`actualError` targets the population P3, matching sampling uncertainty; `conditionalError` targets
the realized clean-trial mean. Their difference is genuine trial-to-trial signal variation.
These are RMS-error diagnostics, not confidence-interval coverage probabilities.

`Calibration.summarize` first averages squared errors within each recording and then gives
records equal weight. Its 95% percentile interval resamples whole records. Selected pipelines
have their own summary. An undefined condition excludes its recording from that summary and is
counted explicitly. Fewer than 30 records are flagged exploratory; reaching 30 does not prove
stability. The benchmark separately counts records with no recommendation.

Use disjoint simulation seeds after implementation choices are frozen:

```matlab
addpath(fullfile(pwd, 'benchmarks'));
out = pipecompare_benchmark(30, tempdir, 1001);
```

The initial 30 records are a pilot. Set the acceptable relative error and confidence-interval
precision before interpreting results, then extend the independent sample as needed. Do not
retune on the validation seeds and continue calling them independent validation.

## Minimum follow-up experiments

1. Compare fixed pipelines, independent selection, factorial shrinkage and group shrinkage on
   identical candidate results. Repeat with correlated candidate errors, artifact subgroups and
   nonrandomly missing/failed candidates. Report selection regret and actual selected-pipeline
   error, not only the agreement of estimated objectives.
2. For new-participant claims, leave whole participants out of prior/hyperparameter fitting.
   If a new recording is adapted, use its training portion to adapt and a disjoint portion to
   evaluate, with temporal buffers for filters and overlapping windows. Include a null group
   difference and a real group difference; systematic same-direction bias need not average away.
3. Sweep blink prevalence through 0, 0.3, 0.7, 0.95 and 1, including weak events, clipping and
   drift. Measure event precision/recall, timing overlap and false-positive duration; check ERP
   amplitude and latency after the complete pipeline.
4. Sweep channel-fault duration through 1, 5, 10 and 20 seconds, isolated and adjacent channels,
   and already referenced/interpolated inputs. Current detection needs at least 10 one-second
   active windows. Check interpolation boundaries, spectral transients, ERP bias and script replay.
5. Compare multitaper scores against matched equal-weight reference estimates on white and 1/f
   backgrounds, center/edge/out-of-band peaks, and several durations. Report bias and RMSE as well
   as variance. Scores are log10 mean PSD in uV^2/Hz, not integrated band power. `tapers` and
   `bandScores` expose effective smoothing width; bands below the minimum 2/T width are unresolved.

The regression suite checks the reproduced implementation failures and numerical invariants.
It does not replace these independent scientific validation experiments.
