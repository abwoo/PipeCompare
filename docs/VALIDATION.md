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
  Paired covariance is retained within each recording; missing candidate patterns use the
  actual design rank. Estimated covariance, group effects and interaction variance are plug-in
  quantities, with no general guarantee for small or selectively missing groups.
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

The calibration target is the **gain-corrected readout**, `estimate/gain`. Adopting a pipeline
does not divide the entire EEG dataset by that diagnostic gain. Calibration exports both
`rawEstimate` and `correctedEstimate`, their errors against population and realized truths
(`rawError`, `rawConditionalError`, `actualError`, `conditionalError`), and the truths themselves.
A claim about the accuracy of the raw ERP must use the raw errors. Passing calibration for
the corrected readout is not evidence of calibrated uncertainty for the raw adopted dataset.

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

## Reproducible public-data execution

The manually dispatched `independent public EEG validation` GitHub workflow uses MATLAB R2026a
and EEGLAB 2026.0.0. `benchmarks/open_data_protocol.json` freezes source versions, SHA-256
checksums, participant subsets, seeds, artifact prevalence and matching rules before execution.
Six disjoint shards evaluate 30 EEGMMIDB participants (21–50), 15 video-annotated EEG-IO
participants (05–19; 00–04 are reserved for development), and six ERP CORE P3 participants
(11–16). Missing, failed and unassessed evaluations remain visible in the artifacts.

- EEG-IO measures one-to-one blink-event precision/recall, timing, false-positive exposure and
  unavailable duration. Corrupt annotation intervals are unavailable, rather than quiet labels.
  Its timestamps are quantized; their complete grid is checked against the declared 250 Hz.
- EEGMMIDB provides external real EEG for matched `pmtm` comparisons, narrow/edge/out-of-band
  spectral probes, and known-duration isolated/adjacent-channel faults. Oracle interpolation
  error and automatic detection are reported separately; a known interval is not a detection.
- The complete simple Standard search is independently refitted on each new background with
  a known injected P3 and ocular contamination, including genuine neural/blink association.
  Raw/corrected readouts and population/realized signal targets remain separate. Random sign
  and circular shift define a zero-mean background ensemble; this is a controlled injection
  experiment, not knowledge of the recording's clean human neural ERP.
- Group transfer fits pipeline effects and interaction variance without the held-out person.
  `Group.applyPrior` keeps those parameters frozen and only adapts the new recording's nuisance
  scale/posterior. Training-only candidate availability determines each fold's family. The
  new person's trials are used for its allowed adaptation; this does not establish prediction
  of unseen trials. Prior estimation uncertainty remains plug-in.
- ERP CORE supplies an external real-task execution/precision example using the authors' P3
  stimulus mapping, correct-response rule and 26-ms display correction. Six participants and
  unknown clean truth cannot establish general population ERP accuracy or group-effect validity.

CI success for this workflow means the declared experiments executed and numerical/replay
invariants held. Precision/recall, reconstruction error and calibration can still be poor.
Inspect the result artifacts; do not tune on these outcomes and continue describing the same
participants as an untouched validation set. Raw recordings are not committed or uploaded
with the metrics. Source attribution and exact downloaded-file hashes accompany every shard.

Data attribution and reuse terms:

- Schalk (2009), EEG Motor Movement/Imagery Dataset 1.0.0, PhysioNet,
  https://doi.org/10.13026/C28G6P, Open Data Commons Attribution License 1.0.
- Agarwal and Sivakumar (2019), Blink: A Fully Automated Unsupervised Algorithm for Eye-Blink
  Detection in EEG Signals, Allerton. The authors make data freely available for research
  use with citation: https://github.com/meagmohit/BLINK. Bytes are retrieved from a pinned
  public mirror; original documentation and annotation conventions define the evaluation.
- Kappenman, Farrens, Zhang, Stewart and Luck, ERP CORE, NEMAR nm000132 v1.1.1,
  https://doi.org/10.82901/nemar.nm000132, CC BY 4.0 for this distribution;
  original article: https://doi.org/10.1016/j.neuroimage.2020.117465.
