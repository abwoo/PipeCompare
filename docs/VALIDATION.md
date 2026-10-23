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


Summary recovery reuses completed Actions artifacts and does not repeat ICA. It filters exact
shard names to avoid counting an aggregate artifact as another participant cohort. A missing
recommendation or unassessed source remains visible; execution counts and calibratable counts
are separate.

For the additional independent group test, participants 21–30 train the prior and 31–50 test it.
The split was fixed before inspecting group losses. Test-record bootstrap intervals are
conditional on that fitted training prior; leave-one-participant-out results are descriptive
and do not include the uncertainty of a new training cohort. Raw and corrected readouts
against population and realized truths have separate accuracy summaries. Only a readout with
the matching uncertainty target is called calibration.

An additional known-spectrum experiment uses independent white/1-f Gaussian segments, several
durations and narrow/ordinary bands with center, edge and outside peaks. Seed 120001 and 600
segments per case are fixed in the protocol. Bias/RMSE target the unsmoothed, bin-averaged PSD;
log/Jensen bias and taper leakage are retained rather than retuned on these outcomes. The
600 segments are Monte Carlo units, not additional human participants.

The declared pilot precision target is a record-bootstrap RMS-ratio relative half-width of 10%,
with at least 30 assessable records. A narrow interval around a biased ratio is not evidence
of calibration, and this precision target cannot establish stability under distribution shifts.


The runner saves one compact score checkpoint per recording before updating the combined score
file. A regression exercises two independent records through checkpoint save/load and the
actual group-ranking path. A MATLAB Code Analyzer check rejects syntax errors in benchmark
scripts before expensive evaluation. A PR commit containing [public-validation] requests a
full public-data rerun; it starts only after the fast numerical/checkpoint suites pass. Other
PR updates do not restart the long scientific experiment.

## Recorded PR #50 public-data results

Derived results are retained in [public CSV metrics](../benchmarks/results/pr50_public_csv_metrics.json)
and [known-spectrum metrics](../benchmarks/results/pr50_spectral_metrics.json). They contain no
raw recordings. The CSV recovery used an independent Python implementation of record-balanced
aggregation; its bootstrap uses a different random stream from MATLAB, so interval endpoints
need not agree exactly. Source run IDs and code revisions are included in the JSON files.

| Evaluation | Independent units | Observed result |
|---|---|---|
| EEG-IO blink events | 15 participants, 430 annotations | Macro precision 0.813, recall 0.968, F1 0.879; participant F1 range 0.727–0.945. Unmatched predictions are blink-specific and can include other ocular activity. |
| Real-background P3 calibration | 30 planned participants; 20 assessable selected readouts | Actual/reported RMS ratio 3.703; MATLAB record-bootstrap 95% CI 1.12–5.72. Selected actual RMS 2.382 uV versus reported 0.643 uV. This is an unfavorable calibration result. |
| Multitaper reference agreement | 30 participants, 810 rows | Maximum difference from matched pmtm/periodogram 3.01e-13 log10 PSD. Ninety rows fail the declared smoothing-resolution criterion. |
| Channel-fault detection | 30 participants, eight fault cases each | Recall 0 for 1/5-second faults and 1 for these injected 10/20-second strong faults. The configured minimum is ten windows. This does not support detection of short faults. |
| Oracle stretch interpolation | Same 240 cases; original observed EEG is the target | Mean reconstruction RMSE 3.82–4.91 uV versus approximately 100 uV before repair; maximum boundary jump 27.8 uV. This oracle result is separate from detection success. |
| ERP CORE P3 execution | Six participants | All six produce recommendations; mean selected raw rare-minus-frequent amplitude 2.05 uV. Clean neural truth is unknown. |
| Known-spectrum stress | 48 cases, 600 simulated segments per case | Mean log10 PSD bias from -0.1470 to +0.00424; eight short-window narrow-band cases are unresolved. Five legacy fractional-rate endpoint counterexamples pass after repair. These segments are not human participants. |

The largest selected-readout calibration discrepancy is in the 95% blink regime: only two
of its six planned participants remain assessable, with ratio 7.17. The all-blink regime
has no assessable recommendation. Coverage is part of the result, rather than an exclusion
that makes the cohort appear larger or easier.

A post-hoc readout sensitivity calculation, leaving the selected pipelines unchanged and
using the unshrunk artifact association, gives ratio 0.974. This diagnoses sensitivity to
the bias penalty; it does not validate a replacement ranking rule or establish causal artifact
bias. These same outcomes cannot be reused to tune a rule and then certify it as independently
validated. The 20 assessable records do not meet the prespecified 30-record precision floor.

The recovered CSV source had a score-checkpoint initialization error, so this recovered artifact
contains no group-prior or native-replay validation. The checkpointed source run retains
per-candidate scores for that separate evaluation. Execution and scientific calibration
remain separate outcomes.

Earlier public-fit timings used an implicit tic/toc clock that nested native calls reset.
A saved five-record shard reports approximately 28 seconds per fit while its explicit engine
clocks report 259–999 seconds. The runner now uses its own timer handle and records that method
in provenance. Summaries mark older fit-timing totals unassessed; calibration and detection
metrics are separate from those clock readings. See the saved timer counterexample JSON.

Group transfer reports two distinct baselines: the minimum raw score within the training-supported
family, and the program's actual frozen per-record recommendation. The actual recommendation may
use candidates outside that restricted family. Both loss differences and that coverage flag are
exported; a raw-minimum comparison is not labeled as improvement over the actual recommendation.
A regression deliberately makes the selected pipeline differ from the raw minimum and verifies
that the actual selected truth rows are used, including unavailable recommendations.

Group comparisons require finite positive point scores, finite nonnegative bootstrap variances,
and finite group scores on the same candidates. A point score can exist while its paired
bootstrap covariance is unavailable (for example, resamples with no clean comparison trials).
An all-missing group row returns no choice. Coverage tables retain every held-out participant's
point/covariance counts and reason; training point coverage is separate from covariance support.
Unassessed-prior fallback cannot change a choice on matched support. Comparisons conditional on
this support and the assessed-only subset are reported separately.

The paired-bootstrap audit found three pre-repair counterexamples: all-missing draws could
certify equivalence, sequential tests dropped missing tail mass, and the factorial variance
fit conditioned on complete draws. Undefined maxima now contribute positive infinity at the
original bootstrap sample size; the factorial covariance fit requires complete paired support.
Recommendation frequencies count all draws and retain an unavailable share.

Participants 21–50 and seeds 20001–20030 are audit/development inputs after these discoveries.
Fresh confirmation fixes EEGMMIDB participants 51–80, prior training 51–60, testing 61–80,
and simulation seeds 20031–20060 before inspecting their outcomes. The earlier numerical
results in this section describe the pre-repair audit, not independent confirmation of the
revised procedure. Unchanged EEG-IO and ERP CORE component controls are reused without counting
them as extra participants. Both versioned protocols and source hashes are retained.

Record-bootstrap calibration intervals also retain positive/zero ratios as infinity.
Zero/zero resamples retain their mass at both possible ratio bounds [0, infinity];
they cannot create a narrow interval by being omitted. Summaries report both shares.
If every reported error is zero, a positive actual RMS is explicitly an infinite
point ratio; zero actual and reported RMS leave the ratio undefined. This is a
summary correction and does not change candidate selection or saved truth rows.
