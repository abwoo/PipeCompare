# Methods: what PipeCompare computes, and why

This page states the mathematics behind the evaluation and ranking, with the sources each choice
rests on. The code that implements each part is named in brackets.

Notation. A candidate pipeline is a sequence of EEGLAB operations. For one component (ROI, time
window, measure) and one condition, trial *i* gives a score *y*ᵢ (*i* = 1…*N*), e.g. the mean
amplitude over the ROI channels and the window after baseline correction.

## Overview: how a search runs

1. **Starting point.**
   - The current state is read from the EEG structure: epoched or not, sampling rate, channels and
     their locations, ICA matrices, IC flags and reference.
   - `EEG.history` is parsed in order, without deduplication; continuation lines (`...`) are joined.
   - Where the history cannot describe the data, PipeCompare says so instead of guessing.
   - The data unit (uV or V) is judged from the amplitude scale (or set with `dataUnit`); volts are
     converted on PipeCompare's copy.
   - A floating-point sampling-rate residue (common after EDF import) is rounded on the copy, and
     the change is recorded.
2. **Legal pipelines.**
   - The plan expands into all combinations of searched values, alternatives and orders. Each is
     checked against the simulated data state; excluded combinations are counted with their reason.
   - Channels interpolated or removed after an average reference (of the plan, or of the data
     before PipeCompare) leave their share of that average in every channel unless the data are
     averaged again: such a pipeline is not legal. An order search never tries it, and a fixed
     order is refused with this reason.
   - Every legal pipeline is run; above `maxLeaves` (500 by default) the search is refused with its
     size, never sampled or truncated. To bring a large search within reach: fix the values you
     are already sure of, pin steps or add `before()` rules instead of searching every order, or
     search in stages (search the early steps, adopt the result, then search the later steps from
     that dataset). Raising `maxLeaves` is possible but every pipeline really runs.
3. **Execution.**
   - Pipelines run as a prefix tree, so a shared prefix (e.g. one ICA before several IC thresholds)
     is computed once.
   - Every step is a native EEGLAB call. Its command is printed and appended to that candidate's
     `EEG.history`.
   - Each candidate is checkpointed. `resume` continues and gives the same result as an
     uninterrupted run.
   - `parallel = true` distributes independent subtrees over a pool.
4. **What is measured.**
   - The contract defines the measures: mean amplitude, peak amplitude or peak latency per
     component.
   - No condition difference or p-value is optimized. Adaptive choices still require downstream validation.
   - The data quality of each measure is its standardized measurement error (SME): analytic for
     means, bootstrapped for peaks and latencies (Luck et al., 2021). SME falls when noise is removed
     and rises when trials are lost, so it prices the rejection trade-off.
5. **Signal preservation.** A known signal is used to catch processing that removes the effect along
   with the noise.
   - A copy containing only a known signal is carried through every candidate with the same
     operations and the same decisions as the real data (bad channels, ICA, removed components,
     rejected epochs, ASR reconstructions). Native commands that decide from the data on their own
     are re-run on the copy and flagged.
   - Amplitude error, peak shift, artifactual deflection, and waveform and topography correlation
     are each checked against a limit.
6. **Ranking.**
   - Failures are reported, never ranked. Constraint violations are listed with their reasons.
     If nothing is feasible, PipeCompare says so and relaxes nothing.
   - Objective: the **gain-corrected SME**, SME divided by the factor by which the candidate scales a
     known signal in that measure (read from the signal check). Raw SME would reward a pipeline that
     shrinks signal and noise alike; SME/gain does not, and ranking it is ranking signal-to-noise
     (Zhang, Garrett & Luck, 2024). Composite over measures that share a unit, or the one measure
     you choose; the others are reported. Derivation in section 2 below.
   - Paired bootstrap over trials (matched by `urevent`); intervals of the difference from the best
     are simultaneous over all candidates (bootstrap max statistic; White, 2000; Romano & Wolf,
     2005), so a larger search does not produce more false "worse" verdicts. *Not distinguished* is
     absence of evidence, not equivalence.
   - A multiverse summary reports, for each measure and condition, the spread of its value over the
     feasible pipelines and the searched choice behind most of it (`result.robustness`; Steegen et
     al., 2016). It shows sensitivity to processing and is not used for the ranking.
   - A candidate that ends with epochs marked for rejection but not removed (e.g. an ERPLAB artifact
     detection step without a removal) carries a note: marks do not remove epochs.
   - Candidates that differ in the reference are ranked in separate strata, never against each other.
     With several strata there is no overall recommendation: each stratum has its own (marked `*`),
     and you choose the one that fits your analysis.
   - The recommendation is the candidate with the lowest objective shrunk through the pipelines'
     settings (section 3, factorial empirical Bayes) among those not distinguished from the best;
     among those shown equivalent to the best, the least aggressive (most trials kept, then least
     distortion).

## 1. Data quality: the standardized measurement error

For a mean-amplitude score the averaged ERP's score is ȳ, and its standard error is the analytic
SME (Luck et al., 2021):

  SME = s_y / √N,  s_y = sample SD of *y*₁…*y*_N.

For a peak amplitude or latency the score of the average is not a mean of trial scores, so the SME
is bootstrapped (bSME): draw *N* trials with replacement *B* times, average, take the peak, and use
the SD of the *B* values (Luck et al., 2021). [`pipecompare.eval.Measure`]

Several measures that share a unit are combined as the root mean square over measures and
conditions. Measures in different units (µV and ms) are never added: one of them must be chosen.

**Difference waves.** A lateralized component (N2pc, LRP) is scored as ERP CORE measures it,
contralateral minus ipsilateral: each trial's score is the electrode contralateral to its side
(target side, response hand) minus the other, so the SME is that of the difference wave. It can
rank pipelines differently from scoring each electrode: noise common to both hemispheres cancels
in the difference. The signal check then injects a field centred on one of the two electrodes and
checks the difference of the two. A difference between independent conditions (MMN, deviant
minus standard) needs no such scoring: its SME is √(SME₁² + SME₂²), √2 times the root mean
square over the two conditions, so the ranking is the same.

## 2. Gain-corrected SME (the noise part of the ranking objective)

**Problem with raw SME.** SME is in µV, so it is not invariant to the scale of the data. A pipeline
that multiplies signal and noise by *c* < 1 lowers SME by the factor *c* without measuring the
signal more precisely. With the default amplitude-error limit of 10 %, a pipeline that attenuates
everything by 9 % passes every constraint and gains 9 % in raw SME, which is more than typical
real differences between pipelines. Zhang, Garrett and Luck (2024) make the same point for filters:
filters reduce the signal as well as the noise, so settings must be compared by signal-to-noise
ratio, not by noise alone.

**Model.** Given its decisions (bad channels, ICA weights, removed components, rejected epochs,
ASR reconstructions), each catalog step acts linearly on the data. Write a trial as signal plus
noise, *x*ᵢ = *a*·*s* + *n*ᵢ, where *s* is the unit signal and *a* its amplitude. The score is a
linear functional φ of the processed data, so

  *y*ᵢ = φ(P*x*ᵢ) = *a*·φ(P*s*) + φ(P*n*ᵢ) = *g*·*a* + *e*ᵢ,  with  *g* = φ(P*s*) / φ(R*s*).

*g* is the pipeline's **gain** for this measure, relative to the signal re-referenced like the
candidate (operator R: a change of reference changes what is measured but is not a distortion,
and candidates with different references are ranked in separate strata anyway).

**Estimator.** ȳ/*g* estimates the same quantity *a* for every candidate in a stratum, with standard
error

  SE(ȳ/*g*) = SME / *g*.

So the objective is SME/*g*. Since *a* is common to the candidates of a stratum, minimizing SME/*g*
is maximizing the signal-to-noise ratio *g*·*a*/SME.

**Invariance.** If P′ = *c*·P with *c* > 0, then *g*′ = *c*·*g* and SME′ = *c*·SME, so
SME′/*g*′ = SME/*g*: an overall rescaling cannot win. (Test:
`test_statistics/testScalingTheDataDoesNotChangeTheObjective`.)

**Measuring g.** The matched-decision injection already processes a copy that holds only a known
signal through the same operations and decisions (section 4). The gain is read there: for a mean
amplitude, the ratio of the recovered and expected window means over the ROI; for a peak
amplitude, the ratio of the peaks; for a latency, *g* = 1, because a latency does not change
when the amplitude is scaled. A gain that is not positive (signal lost or inverted) leaves the
objective undefined, and the candidate is rejected with that reason.
[`pipecompare.eval.Injection.compare` → `gain`; `pipecompare.eval.Rank`]

**Limits.**
- For a peak amplitude the correction is first-order: max(·) is not linear. It is exact when the
  peak is dominated by the signal.
- For steps that are re-run on the signal copy instead of replayed (data-driven native commands,
  flagged in the result), *g* is the gain on a noise-free copy, not necessarily on the real data.
- The correction does not make strong attenuation acceptable. The amplitude-error,
  waveform-correlation and artifactual-deflection limits are checked first and still reject such
  pipelines.

## 3. "Not distinguished from the best": simultaneous bootstrap

The *K* feasible candidates of a stratum are compared on the same trials, so one bootstrap draw
resamples the same trials (by `urevent`) for every candidate. This is the paired bootstrap of
Efron and Tibshirani (1993). With objective estimates θ̂ₖ and bootstrap replicates θ*ₖᵇ, the
statistic is the maximum over all ordered pairs of the centred differences:

  *M*ᵇ = max over *i*, *j* of [ (θ*ᵢᵇ − θ*ⱼᵇ) − (θ̂ᵢ − θ̂ⱼ) ].

Let *c* be the (1 − α) quantile of *M*. Candidate *i* is declared worse than *j* only when
θ̂ᵢ − θ̂ⱼ > *c*. The maximum accounts for comparison across the fitted candidate pairs, subject to
the conditional bootstrap approximation and its sampling assumptions.
This is the bootstrap reality check of White (2000) and the max-statistic step of Romano and Wolf
(2005). The set kept, every candidate not shown worse than another, is a model confidence set in
the sense of Hansen, Lunde and Nason (2011).

These are conditional intervals: preprocessing, artifact detection, candidate eligibility,
signal templates and model hyperparameters were fitted using the observed data. Resampling
stored trial scores does not repeat those fits. `R.validationScope` states this and
`R.fullyHeldOut` is false. The procedure has no general finite-sample family-wise error
guarantee for arbitrary adaptive candidate families, correlated ERP trials or population shifts.

Raw differences are used instead of studentized ones. With few trials the per-pair bootstrap SD is
itself noisy, and a maximum over many pairs is then driven by underestimated SDs. In simulation
this gave 40 % false "worse" verdicts at 30 trials and 24 candidates. α = 0.02 was set by
simulation so that the family-wise rate stays at or below about 5 % for 2–24 candidates and 20–100
trials (`test_statistics`). In the current suite the rates are 5.3 % with 2 candidates, and 0 %
(30 trials) and 2.5 % (100 trials) with 24 candidates.

**Quantile and number of replicates.** *c* is the order statistic *M*₍ₖ₎ with *k* = (*B* + 1)(1 − α),
which is exact when (*B* + 1)α is an integer (Davison & Hinkley, 1997). Monte-Carlo error enters
through the number of replicates beyond *c*. That number is binomial(*B*, α), with SD √(*B*α(1−α)):

| *B* | replicates beyond *c* (mean ± SD) | relative SD |
|---|---|---|
| 200 | 4 ± 1.98 | 49 % |
| 999 (default, peak measures) | 20 ± 4.43 | 22 % |
| 1999 (default, mean measures) | 40 ± 6.26 | 16 % |

With *B* = 200 the verdict for a borderline candidate could change with the random seed. The
defaults are therefore *B* = 1999 and, for the nested peak bootstrap, *B* = 999. [`pipecompare.eval.Rank`]

**Half-samples for peak measures.** A peak measure already carries an inner bootstrap. Resampling
with replacement around it overstates how much its bSME varies. Each outer draw therefore takes a
half-sample without replacement, the classical pseudo-replication device (McCarthy, 1969), and
rescales the result to the full sample by √(*m*/*N*), since bSME ∝ 1/√*n*. In simulation the SD
of a difference was 7.0 ms against 7.4 ms across real replications; with replacement it was 11 ms.

**Equivalence.** *Not distinguished* is absence of evidence. A candidate *i* in the kept set is
called equivalent to the best *b* when the upper bound of its difference is within a margin δ of
the best objective, θ̂ᵢ − θ̂_b + *c*_K ≤ δ·θ̂_b, with δ = 0.05 by default (`equivalenceMargin`).
This is a one-sided test of equivalence (Schuirmann, 1987): the best is the lowest, so only the
upper side is open. *c*_K is the same max-statistic quantile as *c*, but over the pairs of kept
candidates only (step-down; Romano & Wolf, 2005): a candidate already shown worse would otherwise
widen every bound, and no difference among the others could ever be shown negligible. The bound
is in `diffHiKept`, the verdict in `equivalent`.

**Selection bias and cross-fitting.** The minimum of *K* noisy objectives underestimates the
objective of the candidate it selects (the winner's curse): part of its advantage is chance, and
the more candidates, the larger that part. The selection is therefore also evaluated on trials
not used for it. For each of *S* random splits (`nSplits`, default 20), the reference trials of
every condition are divided into halves *A* and *B*. The objectives on each half are computed as
the half-samples of the peak bootstrap: on the *m* trials of the half, rescaled by √(*m*/*N*),
since SME ∝ 1/√*n*. The candidate best on *A* is scored on *B*, and the candidate best on *B* on
*A*. The mean of these 2*S* scores is a **conditional selection diagnostic**, reported against
the apparent best objective as `crossfit.overstatement = honest / apparent − 1`. `honest` is a
legacy field name. Preprocessing, candidate eligibility, signal checks and factorial shrinkage
were fitted on all trials; this does not validate the complete pipeline on independent data.
The result explicitly carries `scope = 'conditional-selection'` and `fullyHeldOut = false`.
For segments of one recording (band power), consecutive segments are
dependent, so the halves are the earlier and the later half (*S* = 1).

*Nested.* A candidate also makes choices from the objective itself: a rejection limit (section
9c), an ICLabel threshold (9d) or components chosen by the SME (9g). Chosen on all trials, they
flatter its objective in the same way, and more so the more such choices it has, so they would
favour pipelines with many of them. Each such step stores, for every option it chose among (one
for a limit; one per component set for 9d and 9g), the scores of the reference trials divided by
the gain and the limit statistic *m* of each trial (`curves`). Removing components and rejecting
epochs are decisions followed by linear steps, so these stored scores are exactly what the
candidate would have with any option and limit. On the selecting half *H* the choice is made
again: per option, its best limit on *H* (a fixed limit: that one; feasibility with the trial
limits scaled to the half); the first option (fewest components) whose objective is within the
equivalence margin of the best option's; and in it the most lenient limit within the margin of
its best. The point margin stands in for the bootstrap equivalence of sections 9c–9d, so on the
half the choice is if anything more lenient. That choice is then scored on *H* (to compare
candidates) and on the other half (the held-out score), with the candidate's trial weights
(section 9f; a trial it did not keep gets the median weight; the weighted SME does not depend on
the weights' scale). Only the last such step of a candidate is made again; components chosen
before a rejection limit come first, as their options already contain every limit. A candidate is
nested only when its own choice's stored scores reproduce its scores: the trials kept are those
up to one limit and the cosine of the centred scores is at least 0.999 (otherwise a step after
the choice changed the scores). `crossfit.nested` counts the nested candidates. The other
data-driven decisions (bad channels, the ICA decomposition, trial weights and cutoffs fine-tuned
afterwards, section 9e) are still made on all trials. In simulation (`test_statistics`: one
candidate whose limit is chosen on a rejection statistic unrelated to the scores, 60 trials), the
nested check finds the limit's chance advantage, about 3 %, where holding the kept trials fixed
finds none.

In simulation (`test_statistics`), with 24 candidates of equal
true noise and 60 trials, the apparent best objective was below the true SME by more than 5 %,
while the cross-fitted value was within 5 % of it. [`pipecompare.eval.Rank.crossfit`]

**Factorial shrinkage (which one to recommend).** With 30–40 trials per condition nearly every
pipeline is in the kept set (simulation, 18 pipelines with true differences up to 10 %: median 18
of 18), and the lowest objective among them is mostly chance. The pipelines are not unrelated,
though: they are combinations of settings, one per step (high-pass × low-pass × weights in
*Standard*). With each step a factor and its settings the levels (*absent* when a pipeline skips
the step), the design *D* has an intercept and one column per level but the first of each factor,
the main effects. With *H* = *D D*⁺ and *R* = *I* − *H*, the objectives split into *H*θ̂, the main
effects, which every pipeline sharing a setting estimates together, and *r* = *R*θ̂, the
interactions plus noise. If the interactions are independent with variance τ² and Σ is the
covariance of θ̂ (from the paired bootstrap of this section), E|*r*|² = τ² tr *R* + tr(*R*Σ), and *r*
is shrunk by the positive-part James–Stein factor

  *s* = max(0, 1 − (*df* − 2)/*df* · tr(*R*Σ)/|*r*|²),  *df* = tr *R*,

(Efron & Morris, 1973): θ̃ = *H*θ̂ + *s r*, empirical Bayes. An interaction the data show (|*r*|² well
above what the noise gives) is kept; the rest is taken for noise. The model needs *df* ≥ 3 (one
factor alone leaves nothing to borrow) and leaves out a pipeline whose cutoffs were fine-tuned
afterwards (section 9e; its settings are its own). The recommendation is, among the candidates
not distinguished from the best (unbiased first, section 15), the one with the lowest θ̃, unless
candidates shown equivalent to the best keep more trials: then the least aggressive of those
(highest minimum retention, then smallest amplitude error to 0.5 %). Keeping more trials decides
only between pipelines the data show to be equally good; among those merely not told apart it
cost 5–11 % of true error in simulation. The cross-fitted check above chooses on each half with
the same shrinkage. In simulation (`test_statistics/testFactorialShrinkageChoosesBetter`; Python
in the review notes: 18 pipelines, main effects of about 4 %, noise 3–10 %, correlated as paired
estimates are) the pipeline with the lowest θ̃ had up to 60 % less excess error over the truly best
than the one with the lowest θ̂, and never more. The table has θ̃ as `shrunk`, the model (*s*, *df*,
the number of factors) is in `byStratum.factorial`. [`pipecompare.eval.Rank.factorialModel`]

**Selection frequency (stability of the recommendation).** The recommendation is a function of the
trials. In each bootstrap draw *b* of the same paired resamples, the choice is made again on the
replicates θ*ₖᵇ: the candidates within *c* of the lowest replicate (the kept set of that draw, with
the critical value *c* of the full sample), then the recommendation rule of the ranking (the
draw's objectives shrunk with the same factorial model, below; which candidates are equivalent to
the best is that of the full sample). The share of draws in which
candidate *k* is recommended is its selection frequency *f*ₖ (`chosenShare`; Σₖ *f*ₖ = 1), the
bootstrap selection frequency used for model selection (Efron, 2014). *f* of the recommended
pipeline is reported (`byStratum.stability`); below 0.7 the result names the next most often
chosen one. In simulation (`test_statistics`) a candidate with a third of the others' noise was
chosen in more than 95 % of the draws, and among 5 candidates of equal true noise with 60 trials
the one recommended on the sample was recommended in about half of them. [`pipecompare.eval.Rank`]

**ICA's random start.** runica resets its own random state (`'rndreset','no'`), so a decomposition
is reproducible but is one of many local optima; which one depends on the order in which runica
draws the samples. Seed *s* > 0 gives runica the samples of the fitting copy in the order of a
Mersenne-Twister permutation with seed *s* (the objective of infomax does not depend on the order,
only its path does); the weights are applied to the data as before. With each seed 1 … *S*
(`icaSeeds`, 2 in the simple mode) every pipeline with an ICA step is run again from the starting
data, every later step deciding again (components, thresholds, rejected epochs) with its own
signal check. All of them run under the same seed, walking the tree as the search does, so
pipelines sharing the steps before ICA share one decomposition per seed (*Standard*: one ICA per
seed in all); a pipeline added by fine-tuning the cutoffs reuses the decomposition of a pipeline
with the same steps before it. For each seed the ranking (steps 1–5, without cross-fitting) is
made again with every pipeline under that seed, and the recommendation is compared with the
search's (`icaSeeds.same`). Comparing pipelines that all saw the same decomposition is the
comparison the search makes; re-running only a few of them against the others' seed-0 values
mixed two decompositions in one ranking. This is a sensitivity check, not part of the ranking: the recommendation and
the adopted pipeline remain the search's (seed 0). [`pipecompare.run.Executor`, icaSeedCheck]

## 4. Signal preservation: matched-decision injection

A copy of the starting data that holds only a known signal is processed with the same operations
and the same data-driven decisions as the real data. Because these operations are linear given
the decisions, the processed copy shows exactly how the candidate transfers the signal. It is
compared with the expected signal: amplitude error, peak-latency shift, artifactual opposite-polarity
deflection, and waveform and topography correlation. The artifactual-deflection limit (5 %) follows
the artifactual-peak criterion of Zhang, Garrett and Luck (2024). Filters are linear time-invariant
operators fixed by their design parameters (Widmann, Schröger & Maess, 2015), so re-running them on
the copy applies the same operator.

**Injected field.** One Gaussian per component, centred in its window and as wide at half maximum
as the window (σ = window/2.355; a component's measurement window is meant to cover most of it),
but at least 50 ms wide at half maximum (σ ≥ 21.2 ms). The width matters in both directions
(simulation of EEGLAB's `pop_eegfiltnew` at 250–1000 Hz, with assumed realistic widths of the ERP
CORE components):
- Too narrow a template rings after a low-pass that leaves the real component unchanged. N170's
  40 ms window would give σ = 17 ms; at σ = 10 ms a 20 Hz low-pass rang by 7.5 % of the peak, past
  the 5 % artifact limit, while it changes a realistic N170 (σ ≈ 20 ms) by about 0.1 %.
- Too narrow a template also under-tests high-passes on broad components: their slow parts are
  what a high-pass removes. With σ = window/4, a 1 Hz high-pass passed for the LRP although it
  fails at the LRP's realistic width. With σ = window/2.355, the pass/fail verdicts of every filter
  tested equal those at the realistic widths: a 1 Hz high-pass fails for N400, P3 and LRP, and
  nothing at 0.5 Hz or below, and no low-pass of 20 Hz or above, fails.

**The recording's own waveform** (`signalTemplate` = `data`, the simple mode's default). The time
course of component *j* is instead *ḡ*(t)·τ(t), where *ḡ* is the recording's average at its ROI
before any pipeline (each condition's epochs averaged, baseline removed, conditions weighted
equally, smoothed by a Gaussian of σ = 5 ms) and τ is 1 in the window and falls to 0 as a raised
cosine over *h* = max(window, 50 ms)/2 on each side, scaled to peak 1 in the window. It is used
when (i) its split-half reliability over the stretch is at least 0.8: *r* = ⟨*a*_odd, *a*_even⟩ /
(‖*a*_odd‖ ‖*a*_even‖) for the tapered averages of the odd and even epochs (uncentred, since the
stretch starts and ends at zero), stepped up to all epochs, 2*r*/(1 + *r*) (Spearman-Brown); and
(ii) the window holds one deflection of one sign (|mean| ≥ half the largest absolute value),
which keeps the gain's denominator away from zero; and (iii) it does not look like eye
activity: for every eye channel *e* (EOG by type or by a name with EOG or EYE; without one, Fp1,
Fp2 and FPz), averaged the same way over the same stretch, not both |⟨*a*_e, *a*⟩| / (‖*a*_e‖ ‖*a*‖) >
0.8 and ‖*a*_e‖ > 2 ‖*a*‖. A blink time-locked to the stimulus (common in oddball tasks, and at the
P3's latency) would otherwise pass every other test, since it is as reliable as a component, and
the check would then test how a pipeline keeps the blink. Otherwise, and for contralateral-minus-ipsilateral
components and band power, the Gaussian above. Any shape is a valid test signal: the copy is
noise-free and every operation is linear given the decisions, so the check measures exactly how
the pipeline transfers that shape. The data's shape only makes the test closer to the real
component (its width and asymmetry, which decide what a high-pass removes); reliability keeps
estimation noise, which a low-pass would remove and so be blamed for, out of it. In a simulation
(6 µV P3 of σ = 70 ms in pink noise of SD 8 µV; standards 2 µV) the criterion passed in 58 % of
datasets with 20 epochs per condition, 80 % with 30 and 100 % with 90, and in 0–2.5 % without any P3.

Over the
scalp it is a Gaussian in the angle θ between a channel's unit position vector and the ROI
centroid, w(θ) = exp(−θ²/(2·0.5²)), normalized to mean 1 over the ROI. Channels without coordinates
(often EOG/ECG) have no defined scalp position: they get no field outside the ROI and the ROI's
mean field inside it, so their presence does not reduce the whole field to a box over the ROI. A box
is used only when no ROI channel has a position. [`pipecompare.eval.Injection`]

**Filters alone, before anything runs.** A filter is a linear time-invariant operator fixed by its
design parameters, so its effect on the injected signal does not depend on the data. Before the
search runs, the filter sequence of each pipeline (its high-pass, low-pass and line-noise steps,
in order) is applied to a short stretch of the signal copy: from the first time-locking event to
the fifth, with a margin of 3.3/*f* + 1 s on each side, *f* being the lowest high-pass edge
(1 Hz when none is lower); EEGLAB's default high-pass at *f* Hz is about 3.3/*f* s long, and the
other filters are shorter than 3.3 s. Around those events this gives what filtering the whole copy
gives (the copy is zero apart from the signal, and a filter reaches no further than its length;
`test_signal` checks agreement within 0.005 in amplitude error and artifactual deflection). The
limits of the signal check, except the topography (a filter acts the same on every channel), are
then applied. Each sequence is checked once; a pipeline whose filters alone break a limit is
excluded with that reason (status *excluded*, reason *filters alone: …*) and never run. Other
steps act on the signal copy spatially (re-reference, interpolation, IC removal) or select epochs,
so they could only mask such a distortion, not remove it; with the check off, the same pipelines
are run and excluded by the full check (`test_engine`). Not done for band power, epoched data, or
a pipeline that resamples or epochs before its last filter. [`Injection.filterCheck`]

## 5. Preconditions checked before anything runs

- **Interpolation needs electrode positions.** Spherical-spline interpolation is defined on
  positions on the sphere (Perrin et al., 1989). Without any channel location, interpolation steps
  (`badchannels` with action interpolate, `channels` interpolate, `restore`, native `pop_interp`)
  are excluded at plan time, with the reason. When only some channels lack a position, a step that
  would interpolate one of those channels fails with the channel named. EEGLAB's `eeg_interp`
  would leave such a channel unchanged without saying so, while the step would report it
  interpolated.
- **ICLabel needs channel locations.** Its features include the components' scalp maps
  (Pion-Tonachini et al., 2019), so `icremove` and native `pop_iclabel` are excluded without them.
- **Flat channels are bad.** A channel whose SD is zero or under 1% of the median SD of the
  channels tested has no signal (e.g. a disconnected input). Kurtosis and joint probability are
  undefined on a constant, and `pop_rejchan` normalises them over the channels tested, so such a
  channel is set aside before the measures and marked bad by the `badchannels` step. The data
  line names flat EEG channels before a run (first 60 s). [`pipecompare.run.Steps.flatChannels`]
- **Bad channels by their own noise (measure `resid`).** `pop_rejchan`'s kurtosis and joint
  probability are invariant to a channel's scale: kurtosis is a ratio of moments, and joint
  probability bins each channel between its own minimum and maximum (`realproba`, 1000 bins), so
  multiplying a channel by any factor changes neither. A channel that is only noisier (Gaussian
  noise of its own, as with poor contact), or that lost its gain, therefore passes both. `resid`
  measures, for each tested channel $i$ with covariance matrix $C$ over the tested channels, the SD
  of what is left of it after its least-squares prediction from all the others,
  $r_i = \sqrt{C_{ii}(1 - R_i^2)} = \sqrt{C_{ii} / (C_{ii}\,(C^{-1})_{ii})}$, computed for all
  channels from one inverse of the correlation matrix. A channel is bad when the robust z of
  $\log r_i$ (median and MAD over the channels, the MAD at least 0.1) is above the threshold. Shared
  activity (volume conduction: brain rhythms, blinks, a nearby reference) is predicted by the
  neighbours and does not count; noise of the channel's own is not predicted. Samples where any
  channel is beyond 6 robust SDs of its own (blinks, movement) are left out of $C$, unless that is
  over half of them; at most 2·10⁵ samples are used; epochs have their mean removed. The test is
  skipped (and said) when the correlation matrix is singular (smallest eigenvalue below 10⁻⁶ of the
  largest: average-referenced or interpolated channels, each predicted exactly), and needs at
  least 8 channels. In simulations (64 and 32 channels, spatially smooth sources plus sensor
  noise, strong occipital alpha, frontal blinks, a vertex reference), kurtosis + joint probability
  found 0–12% of channels 1.5 to 5 times noisier than the rest, `resid` 75–100%, with 0–0.1 false
  alarms per recording against 3–8 for kurtosis + joint probability, which flagged blink-laden
  frontal channels. One missed channel three times noisier left 52% of the epochs under a 100 µV
  limit. The simple mode and the panel detect with `kurt+prob+resid`.
  [`pipecompare.run.Steps.residualZ`]
- **Kurtosis and joint probability without the blinks.** Blinks make the frontal channels'
  values heavy-tailed, so `pop_rejchan` marks them bad (in the simulations above, with frontal
  blinks: 7.5 good channels per 64-channel recording, all frontal; with the blinks left out, none
  frontal). The measures therefore run on the copy without the samples where an eye channel
  (`Injection.eyeChannels`: EOG, else Fp1, Fp2, FPz) is beyond 5 robust SDs of its own (median
  and MAD), or in a blink found by the eye channels' hidden Markov model (section 9i; continuous
  data), with 200 ms on either side (epoched data: without those epochs), unless that is over
  half of the data. A weak blink can stay under that limit in the eye channels while the
  frontal channels show it clearly, and one blink left in is enough to make them heavy-tailed
  again (simulated: 20-30% of the blinks missed, FPz then marked bad). So the blinks found give
  the blink's scalp pattern $a$ (every channel regressed on the eye channel with the most of
  them), and the matched spatial filter $w = C^{-1}a$, $C$ the channels' covariance away from
  the blinks, gives its time course $w^\top x(t) / w^\top a$; samples of it beyond 5 robust SDs are
  blinks too (simulated: all blinks found). A noisy channel weighs little in $w$, so a bad
  frontal channel's own noise or spikes are not taken for blinks (simulated: a spiky FPz is
  still marked bad). Spiky channels and line noise are found as before. `resid` leaves out such
  samples itself. [`pipecompare.run.Steps.withoutBlinks`]
- **Bad channels are found in the analysis band (`detectLowpass`).** The copy the bad channels are
  found on is also low-passed: at 40 Hz for ERPs (the highest low-pass compared; its stop band
  starts at 50 Hz), 5 Hz above the highest band for band power, none on epoched data, data already
  low-passed at or below that edge, or an edge within 10% of the Nyquist frequency. Line noise on
  one channel, which the pipeline's low-pass removes, otherwise makes that channel bad and
  replaces it by its neighbours (the measured electrode included). [`pipecompare.simple.Presets.detectLowpass`]
- **Which bad channels are interpolated (`repair` = `'error'`).** Interpolation is linear: with the
  set *b* interpolated the data are *M*_b *X*, *M*_b EEGLAB's spherical spline, read off by
  interpolating an identity matrix with `pop_interp` (`ChannelChoice.interpMatrix`; equal to
  `pop_interp` on the data, `test_engine/testBadChannelsInterpolatedByTheirEffectOnTheError`).
  Of the steps below, up to the epoch rejection, the re-referencing (a channel map *R*), epoching
  and baseline removal (per channel, commuting with *M*_b) and the rejection are applied, when they
  are the same in every pipeline below the step (else every bad channel is interpolated); the
  filters between are replaced by the band of the detection copy, and ICA is left out, as it
  differs between pipelines and is not linear. For each *b*: *R M*_b *X* epoched, the total error
  at the rejection's best limit (`IcThreshold` scoring, sections 9c–9d), the gain from the known-signal
  copy with the same *b*. Backward from all flagged channels, the channel whose being left gives the
  lowest objective is left, one at a time; among the sets equivalent to the best (paired bootstrap,
  equivalence margin), the one interpolating the most is chosen, so a flagged channel stays as it
  is only when interpolating it measurably raises the error. More than 12 flagged channels: all
  interpolated. Not re-made on the selecting half in the cross-fitting check (section 3), unlike
  the choices inside the pipeline after ICA. [`pipecompare.eval.ChannelChoice`]
- **Event types are matched exactly.** `pop_epoch` selects events with `strmatch(…, 'exact')`, so
  'S 1' and 's 1' are different types. The contract is validated the same way and names the
  near match.

## 6. Marks are not removals

EEGLAB keeps rejection marks (`EEG.reject.*`) separate from the data; only `pop_rejepoch` removes
epochs (Delorme & Makeig, 2004). ERPLAB's artifact detection marks `EEG.reject.rejmanual` and the
epochs' flags, which ERPLAB's averager honours (Lopez-Calderon & Luck, 2014). PipeCompare's scores use
every epoch that is still in the data. A candidate that ends with epochs marked but not removed
therefore carries a note, so a marking-only step is not mistaken for a rejection.
[`pipecompare.run.Executor`, note column of the ranking]

## 7. Multiverse summary (sensitivity, not ranking)

Every allowed pipeline is run, so the plan's whole "multiverse" of processing choices is available,
not a sample of it (Steegen et al., 2016; for ERP processing, Clayson et al., 2021). For each
measure and condition, and within each stratum, the result reports the spread of the point estimate
over the feasible pipelines (min, median, max, range). It also names the searched choice that
accounts for most of that spread. For a choice with values *v* (with *n*ᵥ pipelines and mean
estimate *x̄*ᵥ), the share is

  η² = Σᵥ *n*ᵥ (*x̄*ᵥ − *x̄*)² / Σᵢ (*x*ᵢ − *x̄*)²,

the fraction of the spread that lies between that choice's values. The summary describes how much
the measured quantity depends on processing. It is not used to rank or recommend pipelines, and it
contains no comparison between conditions. [`pipecompare.eval.Rank.robustness`, `result.robustness`]

## 8. Band power (continuous data, e.g. resting state)

**Segments.** A band-power contract with `segment` *T* cuts the continuous recording into
consecutive *T*-second segments. They are marked with EEGLAB's own `eeg_regepochs(…,
'extractepochs', 'off')`, which inserts the events and their urevents. A segment is therefore
identified by its urevent in every candidate, exactly as an ERP trial is, and a segment that a
candidate rejects is missing (NaN) in that candidate. *T* must hold at least two cycles of the
lowest band edge (*T* ≥ 2/*f*₁), and band edges must lie below Nyquist.

**Score and SME.** Per segment, the score is log₁₀ of the mean PSD in the band (uV²/Hz), rather
than the integral of power over the band. The PSD is the
one-sided multitaper power spectral density (Thomson, 1982): the segment is multiplied by each of
*K* discrete prolate spheroidal (Slepian) tapers, the orthonormal windows whose power is most
concentrated within ±*W* Hz of each frequency (Slepian, 1978; computed as the eigenvectors of the
tridiagonal matrix of Percival and Walden, 1993), and the *K* periodograms are averaged, then
averaged over the band's bins and the ROI channels. The time-bandwidth product is *NW* = *WT* = 2
(*K* = 3 tapers), lowered for a band narrower than 4/*T* Hz so that the smoothing (2*W*) stays
within the band when its width permits (*NW* = 1.5 gives 2 tapers, *NW* = 1 a single one).
The minimum smoothing width is 2/*T* Hz: narrower bands are unresolved and warn during contract
validation. `tapers` and `bandScores` expose effective resolution as a second output.
One taper's periodogram has
only two degrees of freedom per bin; averaging *K* nearly independent ones lowers the spread of
the score over segments. In simulation (an alpha rhythm in white noise, 1 s segments, 8-13 Hz) the
spread was half that of a single Hann taper (`test_statistics`), and with segments that also differ
in their true alpha power, the error of the mean log power over 40 segments was 28 % lower
(0.028 instead of 0.039). Consecutive segments use a dependence-aware mean SE, described below. With
`conditions` and `epoch` instead of `segment`, each epoch's band power is scored (event-related
band power).

**No gain correction is needed.** The log is scale-free: if the data are multiplied by *c*, every
score shifts by 2·log₁₀ *c*, and their SD (the SME) is unchanged. The gain of section 2 is
therefore 1 for band power. A step that attenuates the band itself is caught by the signal check:
sinusoids one frequency step (1/T for epochs of T s) inside each band edge and at the band's
centre, over the band's ROI, are carried through the candidate with the same decisions as the real
data, so a filter that cuts into the band is seen even when the centre is untouched. Their
recovered amplitudes (least squares per epoch, all frequencies fitted together, averaged) are
compared with the expected, re-referenced field; the worst frequency of the band counts (amplitude
error, topography correlation). The test frequencies lie on the 1/T grid, so they are orthogonal
over an epoch, and contiguous bands do not share one. A band narrower than five grid steps is
tested at its centre.

**Dependent segments: moving-block bootstrap.** Consecutive segments of one recording are not
independent, because state and noise change slowly. Resampling them one by one treats them as
independent and understates uncertainty. The point SE now uses the Bartlett HAC estimate of the
mean: with weighted residuals *z*_t = *w*_t(*x*_t − μ̂), sum Σ*z*_t² plus
2 Σ_(h=1)^L (1 − h/(L+1)) Σ_t *z*_t*z*_(t+h), multiply by *n*/(*n*−1), then take the square root
and divide by Σ*w*. *L* = floor(√*T*), where *T* counts original segment positions and *n* counts
positive-weight observations. Rejected/missing positions have zero influence and retain their
time positions. The automatic limit path uses the same SE. The comparison of section 3 resamples
segments in blocks of *l* = round(*n*^{1/3}) consecutive segments and recomputes the HAC estimate
in each resampled sequence; frequency counts alone would lose its temporal order.
This is the moving-block bootstrap (Künsch, 1989), with the block-length order of Hall, Horowitz
and Jing (1995). ERP trials, which are separated by inter-trial intervals and time-locked to
different events, are resampled one by one as before.

## 9. Repairing epochs: channels interpolated within an epoch

A rejection step (`reject_threshold`, `reject_jointprob`, `reject_kurtosis`) with `interpolate` =
*n* > 0 keeps an epoch that the test fails because of at most *n* channels, as FASTER's epoch-level
channel interpolation does (Nolan, Whelan & Reilly, 2010). The decision is the test's own
per-channel marks (`EEG.reject.rejthreshE`, `rejjpE`, `rejkurtE`), taken before anything is
changed: an epoch with 1 to *n* flagged channels, all with a location, gets those channels replaced
in that epoch by EEGLAB's spherical-spline interpolation (`eeg_interp`, Perrin et al., 1989) from
the other channels of the same epoch, its marks are cleared, and it is kept; the other marked
epochs are removed with `pop_rejepoch`. An epoch in which an electrode that is measured (any
component's or band's ROI) is flagged is never repaired, but rejected: an interpolated value is a
weighted mean of its neighbours and carries less noise than a real electrode, so repairing a
measured electrode would lower its SME without measuring the signal more precisely. `eeg_interp` is linear in the data, so its weights for a
set of channels are read once, as its output for a unit impulse on each channel, and applied to
every epoch with that set: the numbers are those `eeg_interp` gives on the epochs themselves
(checked in `test_engine`). Repaired epochs are not tested again.

Interpolated values change the average reference of their epoch: the replaced values were part of
it. As for whole channels (section 5), a plan that repairs epochs after an average reference must
average again afterwards, or it is refused with the reason; the simple mode adds that second
average reference. The signal check replays the same repairs (same channels in the same epochs)
on the copy, so the change they make to the known signal is measured and limited like any other.
The 20 % limit on interpolated channels concerns whole channels and does not count repairs.
[`pipecompare.run.Steps.rejectEpochs`, `pipecompare.run.Steps.interpolateEpochs`]

## 9a. Peak-to-peak rejection

`reject_threshold` with `method` = `peaktopeak` marks an epoch when, on any tested channel, the
range (maximum minus minimum) within a moving window of `window` ms (default 200), moved in 50 ms
steps, exceeds `uv`: the moving-window peak-to-peak test of ERPLAB (Lopez-Calderon & Luck, 2014)
used by ERP CORE (Kappenman et al., 2021). Unlike an absolute threshold it is insensitive to slow
drifts and offsets within the epoch. The marks are written where `pop_eegthresh` writes its own
(`EEG.reject.rejthresh`, per channel `rejthreshE`), so removal, repair (section 9) and replay on
the signal copy are the same. [`pipecompare.run.Steps.markPeakToPeak`]

## 9b. ICA fitted without extreme stretches

The decomposition is fitted on a copy from which the stretches far noisier than the rest are left
out: 1 s windows (or epochs) whose mean log SD over the EEG channels has a robust z (median, MAD)
above 3, at most 20 % of them, worst first. Non-EEG channels do not count, so blinks stay in the
fit. The unmixing matrix is applied to all of the data. Fitting ICA on data cleaned of
non-stereotyped artifacts is standard practice, since such stretches otherwise take up components
(Delorme & Makeig, 2004; the EEGLAB tutorial's advice to reject bad data before ICA). `fitClean` = 0
turns it off. [`pipecompare.run.Steps.cleanForIca`]

## 9c. Rejection limit chosen from the data (exact optimum)

`reject_threshold` with `uv` = `auto` (the simple mode's default when no epochs are repaired)
chooses the limit in each pipeline from its data at the rejection step
[`pipecompare.eval.Threshold`].

**The problem.** Let *m*ᵢ be epoch *i*'s test value: the largest absolute value over the tested
channels and the whole epoch (`pop_eegthresh`'s test), or the largest moving-window range (section
9a). The limit *u* keeps *K*(*u*) = {*i* : *m*ᵢ ≤ *u*}. The objective is the ranking's: *J*(*u*) =
RMS over conditions *c* and measures *q* of SME꜀q(*K*(*u*)) (or the one measure named by
`objective`), with SME = SD / √*n* of the per-trial scores (mean amplitude or log band power) of
the kept trials. The constraints are the ranking's: *n*꜀(*K*(*u*)) ≥ max(`minTrials`,
⌈`minRetention` · *N*꜀⌉) in every condition.

**Exact solution.** *K*(*u*) is nested in *u* and changes only where *u* crosses a value of *m*, so
*J* and the constraints are piecewise constant, with breakpoints at the sorted distinct values
*v*₁ < … < *v*ᵤ (at most one per epoch). The minimum of *J* over all real *u* is therefore the
minimum over the sets *K*(*v*ⱼ). For each set, *n*, Σ*y* and Σ*y*² are inner products of the
scores (centred per condition, which leaves the SD unchanged and avoids cancellation) with the
indicator matrix [*m*ᵢ ≤ *v*ⱼ], so all of them cost one matrix product. This is the global optimum
over all limits, not the best of a grid; it equals a brute-force search (`test_engine`).

**Choice.** The minimizer is chosen on the same trials it is scored on, and a slightly lower SME at
a stricter limit may be chance. The limits are therefore compared as the ranking compares
pipelines (section 3): a paired bootstrap over trials (the same resampled trials for every limit;
moving blocks for segments of one recording) with the simultaneous max statistic over all pairs
gives the limits not shown worse than another. At most 60 limits enter the bootstrap: the
feasible ones spread evenly by rank over the range, the minimizer included. Among those shown
*equivalent* to the best, the highest is used (the least aggressive), so an epoch is kept unless
keeping it makes the measure measurably less precise. Equivalence is tested in a fixed sequence
(Maurer, Hothorn & Lehmacher, 1995; Westfall & Krishen, 2001): from the best limit upwards, each
limit against the best alone, equivalent when θ̂ᵢ − θ̂_b plus the (1 − α) bootstrap quantile of
its deviations is at most δ·θ̂_b (δ = `equivalenceMargin`); the first limit that fails ends the
sequence. Each test is made only when all before it passed, so the family-wise error stays α
without a correction for their number. An earlier version bounded every difference with the
critical value simultaneous over all pairs of the limits kept (*c*_K); the pairs of limits far
apart set it, and no limit but the best was ever equivalent (simulation: the best itself in every
run). The same sequence chooses the ICLabel threshold or components (towards fewer components,
section 9d/9g) and the bad channels interpolated (towards more channels, section 5). Not distinguished alone is not enough: the lenient limits keep the most variable
epochs, so their bootstrap differences are the widest and set the critical value *c*; the set not
shown worse then reaches limits that reject almost nothing even when rejecting clearly helps (in
simulation with 15% large-artifact trials and 100 trials per condition, the true SME of the limit
chosen was 0.54 with that rule, 0.33 with equivalence and 0.33 at the true best limit; on clean
data both rules gave the best limit's SME).
The value used is the number with the fewest significant digits in [*v*ⱼ, *v*ⱼ₊₁) (any number
there keeps the same epochs). It is applied by `pop_eegthresh` (or the peak-to-peak marks) and
recorded in `EEG.history`, and the signal check removes the same epochs from its copy. When no
limit is feasible, nothing is rejected and the ranking reports the missing trials.

**Scope.** The scores are those of the data at the rejection step: the result is exact when no
later step changes them (the baseline is already removed from them). The signal gain *g* is left
out of *J*: the known signal is the same in every epoch, so *g* does not depend on *u* (with
several measures of different gains, their weights in the RMS would differ slightly). Peak
measures have no per-trial score and keep fixed limits; so do repaired epochs (section 9), where
the limit also decides which channels are interpolated.

**Trials kept against all trials** (every rejection, `keptShift`). The statistic *m* includes the
measured electrodes in the window, so a low limit can reject a trial for a large real component,
and the known signal, the same in every epoch, cannot show it. For each mean-amplitude component
and condition, with *x̄*_K the mean score of the *n*_K trials kept and *x̄* that of all *n*,
*z* = s (*x̄*_K − *x̄*) / (SD(*x*) √(1/*n*_K − 1/*n*)), the standardized difference under random
rejection (s = −1 for a negative component, so *z* < 0 means less of the component). The lowest
*z* is stored; the results warn when it is below −3 (`Presets.KeptShiftZ`). Artifacts that push
the window the same way give the same sign, so this is a warning, not a constraint.

## 9d. ICLabel threshold chosen from the data (exact joint optimum)

`icremove` with `threshold` = `auto` (the simple mode's default when no epochs are repaired)
chooses the ICLabel threshold in each pipeline jointly with the rejection limit that follows
[`pipecompare.eval.IcThreshold`].

**The problem.** Let *q*ₖ be component *k*'s largest ICLabel probability among the artifact classes
of the step (a probability equal to 1 is left out: `pop_icflag` flags *t* < *p* < 1). The threshold
*t* removes R(*t*) = {*k* : *q*ₖ > *t*}. Removing R gives the data X − A_R W_R X restricted to the
kept components, X_R = A₋R W₋R X (`pop_subcomp`'s back-projection; X itself when R is empty), with
W the unmixing and A its inverse. For a set R and a rejection limit *u*, let *J*(R, *u*) be the
objective of section 9c computed on X_R, with every measure's scores divided by its signal gain
*g*_q(R) (section 2), since removing components can scale the measured signal. Constraints: those of
section 9c, and the signal check (section 4) of X_R's known-signal copy.

**Exact solution.** R(*t*) is nested in *t* and changes only where *t* crosses a value of *q*, so
there are at most *K* + 1 sets: removing nothing, and R_j = {*k* : *q*ₖ ≥ *c*_j} for the distinct
values *c*₁ > *c*₂ > … of *q* above 0.5 (`MinProbability`: a component more likely that artifact
than anything else). Epoching and baseline removal, the only steps allowed between the removal and
the rejection (the plan checks it), select samples and subtract per-channel means, so they commute
with the spatial map A₋R W₋R: each set is applied to the epoched, baseline-corrected data once,
the same for the known-signal copy (which gives *g*(R) and the signal check). For each set the
inner problem min over *u* of *J*(R, *u*) is solved exactly by section 9c (or *u* is the step's fixed
limit, or no rejection follows). The minimum over all (*t*, *u*) is therefore the minimum over at
most *K* + 1 exact inner solutions: the global optimum of the joint problem, not the best of a grid.
Each set's value equals EEGLAB's own processing (`pop_subcomp`, `pop_epoch`, `pop_rmbase`,
brute-force limits; `test_engine`).

**Choice.** The sets are compared as the ranking compares pipelines (section 3): the same paired
bootstrap (the same resampled trials for every set, each at its own best limit) with the
simultaneous max statistic gives the sets not shown worse than another. Among those also shown
equivalent to the best (as in section 9c), the one removing the fewest components (the highest
threshold) is used, so a component is kept unless keeping it makes the measure measurably less
precise. The value used is the number with the
fewest decimals in [*c*_{j+1}, *c*_j) that removes exactly that set under `pop_icflag`'s comparison
(in the class ICLabel stores, single precision). It is applied by `pop_icflag` and `pop_subcomp`
and recorded in `EEG.history`; the rejection step that follows then chooses its own limit on the
data left by section 9c (the sets are compared at their best limits; the step uses, as always,
the highest limit equivalent to the best). When no set is
feasible, the set removing nothing is used and the ranking reports why.

**Scope.** Exact when only epoching, baseline removal and an amplitude rejection that removes
epochs follow the removal; all pipelines below the step must share those steps (otherwise the
step fails with that reason). The signal check of each set uses all epochs: the rejected epochs
change the average of the known signal only where epochs overlap; the full pipeline's own check is
the one the ranking uses. Peak measures have no per-trial score, and repaired epochs are not a
removal; both keep fixed thresholds. When adopting a pipeline, the threshold the search chose is
reused; `apply` (saved scripts) chooses again on the data it is given, with its own signal copy.

## 9e. Filter cutoffs fine-tuned from the data (decisions held fixed)

**Decomposition.** A pipeline does two kinds of things to the data: it takes *decisions* from the
data (bad channels, the ICA decomposition, the components removed, the epochs rejected) and it
applies *linear operations* (interpolation, re-referencing, filtering, IC removal, epoching,
baseline removal, the window mean). Once the decisions are fixed, the score of every trial is a
fixed linear functional of the raw data. Let *D* be the data just before the first filter
(channels × samples), *P*_R = I − A_R W_R the IC-removal matrix of the removed set *R*, *h* the
impulse response of the filters (the same for every channel) and *k*ᵢ the epoching, baseline
and window weights of trial *i* at the measured electrodes *w*. Then the score is
*x*ᵢ(*h*) = *w*ᵀ *P*_R (*D* ∗ *h*) *k*ᵢ = (*w*ᵀ *P*_R *D*) ∗ *h* · *k*ᵢ, because the temporal filter and the
spatial matrix commute. So only the measured electrodes of *P*_R *D* need to be filtered to
obtain every trial's score for any cutoff: a few channels instead of the whole montage.

**Gain.** The known signal is a time course times a scalp field, *s*(*c*, *t*) = *a*(*c*) *g*(*t*), so
its score factorises into a spatial part, unchanged by the cutoffs, and the window mean of *g* ∗ *h*.
The gain of other cutoffs *h*′ is therefore *G*(*h*′) = *G*(*h*) · *G*_F(*h*′) / *G*_F(*h*), where *G*(*h*) is
the pipeline's gain from its full signal check and *G*_F the gain of the filters alone (section 1's
filter check, which also excludes cutoffs whose filters alone distort the known signal).

**Search.** For the recommended pipeline, the cutoffs of each compared filter are tried on a grid
between the lowest and the highest value compared, in steps of 1, 2 or 5 times a power of ten
(at most about 20 values per filter, the compared ones included; e.g. 0.1, 0.15, …, 1 Hz and 20,
22, …, 40 Hz, 209 pairs). With the pipeline's components *R* and rejected epochs held fixed, the
objective of the ranking (gain-corrected SME) is computed exactly for every pair. As a check of
the assumptions, the objective at the pipeline's own cutoffs must reproduce the pipeline's
(relative tolerance 10⁻³); otherwise nothing is fine-tuned. The epochs are located once in the
continuous data (the epoching of a dataset that holds each sample's index), so each pair costs a
filtering of the measured electrodes and an index operation.

**Verification.** The pair with the lowest objective is a prediction under fixed decisions. When it
differs from the pipeline's own cutoffs, the pipeline with those cutoffs is run like every other
one, from the starting data (its ICA reused), its components and epochs decided again on its own
data, with its own full signal check, and it joins the ranking. Its result is therefore exact in
the same sense as every other pipeline's, and the ranking decides whether it is better.

**Scope.** Mean amplitudes of an ERP on continuous data; after the first filter only filters, IC
removal, epoching, baseline removal and an amplitude rejection that removes epochs (no repair);
no IC removal or resampling before the filters. It is on in the simple mode, an option in the
advanced panel (*Options… > Fine-tune filter cutoffs*) and `'refineCutoffs', true` in scripts.

## 9f. Epochs weighted by their noise

**Weights.** For epoch *i* of a condition with *n* epochs, let *r*ᵢ = *x*ᵢ − (Σ_{j≠i} *x*ⱼ)/(*n* − 1) be the
epoch minus the leave-one-out mean of the others, at the measured electrodes (all of them, over
all samples of the epoch). The part of *r* that goes with the scores is taken out first: with
*k*ⱼ the kernel of score *j* (the mean over its electrodes and window) and *u*ⱼ = (*r*ᵢ · *k*ⱼ)ᵢ the
epochs' score residuals (an *n*-vector each, *U* = [*u*₁ …]), *r* ← (*I* − *U U*⁺) *r* in each condition,
which is *r* − Σ*k* (*k*ᵀΣ*k*)⁻¹ *k*ᵀ*r* with Σ the residuals' covariance over the samples. EEG noise is
slow, so without this an epoch whose score happens to lie far from the others' also has a large
residual over the epoch and would be counted less *for its score*: the weights then depend on the
scores, and the sandwich SME below, which takes the weights as given, came out a third too low
(simulation, AR(1) noise φ = 0.95 at 250 Hz, every epoch with the same noise: reported/true 0.67;
the weighted mean looked 30 % better than equal weights while it was 5 % worse). With the
projection, what is left is uncorrelated with the scores (independent for Gaussian noise of any
size per epoch) and still follows each epoch's noise level: reported/true 1.04, and with a third of
the epochs 3 times noisier, 0.95 and an error 0.67 of equal weights. Its noise is *v*ᵢ = mean(*r*ᵢ²);
the leave-one-out mean keeps the epoch from lowering its own noise. With *v*ᵢ floored at a tenth of the condition's median, the weight is
*w*ᵢ = (1/*v*ᵢ) / mean(1/*v*). If the noise of each epoch is a common process scaled by a factor σᵢ
(*e*ᵢ = σᵢ *u*ᵢ), its score has variance ∝ σᵢ² ∝ *v*ᵢ, and the weighted mean Σ *w*ᵢ *s*ᵢ / Σ *w*ᵢ with *w*ᵢ ∝ 1/σᵢ²
has the smallest variance among weighted means (Gauss–Markov; for ERPs, weighted averaging goes
back to Hoke et al., 1984). Plain averaging (all *w* = 1) and rejection (*w* ∈ {0, 1}) are special cases;
the floor bounds any weight to 10 times the typical one. The weights are a data-driven decision
like a rejection, stored by the `urevent` of each epoch's time-locking event
(`EEG.etc.pipecompare.trialWeights`), and replayed on the signal copy.

**Score and error.** The estimate is *ŝ* = Σ *w*ᵢ *s*ᵢ / Σ *w*ᵢ, and its error is the sandwich
(heteroscedasticity-consistent) standard error of a weighted mean,

  SME = √( *n*/(*n* − 1) · Σ *w*ᵢ² (*s*ᵢ − *ŝ*)² ) / Σ *w*ᵢ,

which is SD/√*n* for equal weights and conditions on the fitted weights. It is not a guarantee
for weights learned from the same scores or for dependent trials. In each bootstrap
draw a trial drawn *k* times counts *k* times with its weight (the weights are not estimated
again); peak measures average the waveforms with the weights (inner bootstrap included). In
simulation (`test_statistics`: 20 trials of noise SD 1 and 10 of SD 4) the error of the mean was
0.22 with inverse-variance weights against 0.45 with equal weights, and the mean SME matched the
spread over replications within 8 %.

**Signal check.** The known signal is the same in every epoch, so its weighted average equals its
plain average up to what differs between epochs (filter tails of neighbouring events); the check
averages the copy with the same weights. **Scope:** ERP measures (for band power, a weighted mean
of log power would favour the epochs with less power). The rejection limit chosen from the data
(9c) is chosen for equal weights. Filter cutoffs can be fine-tuned (9e) with the weights held
fixed, since the weights use the measured electrodes only.
**Downstream averaging.** EEGLAB and ERPLAB average with equal weights, so a weighted pipeline's
data averaged there do not give the ERP it was scored on. Adopting it (or `apply`) stores the
weighted average of each condition, Σ*w*_i *x*_i / Σ*w*_i over its epochs, in
`EEG.etc.pipecompare.weightedERP` (`Measure.weightedAverages`), and the results name the same
pipeline with equal weights. In simulation (`test_engine`, a P3 whose size varies with SD 50 % and
latency with SD 50 ms), the weighted estimate stays within 3 standard errors of the P3 in the
data: the weights, from the noise over the whole epoch, do not count trials with a large P3
measurably less.

## 9g. ICA components chosen by their effect on the SME

`icremove` with `threshold` = `sme` chooses the set R of components to remove by the objective
instead of by ICLabel probability [`pipecompare.eval.IcThreshold.chooseBySme`].

**The problem.** With *J*(R) = min over *u* of *J*(R, *u*) as in section 9d (the objective at the best
limit of the rejection that follows, scores divided by the signal gain, infeasible when the
signal check fails or too few trials are kept), choose R ⊆ C, where C is the set of components
whose likeliest ICLabel class is not Brain (a component ICLabel takes for brain activity is never
a candidate). ICLabel is used only to protect, not to decide.

**Path.** The 2^|C| subsets cannot all be scored, so the sets scored are nested: the first *j*
components of an order of C. For a mean amplitude the order is by what each component does to
the score alone. Its part of each trial's score *d*_k,i (below) is exact, as the score is linear
in the components. Its own trial-to-trial variance there, as a share of the score's variance
(summed over the measured components and averaged over the conditions), is divided by its share
*s*_k of the known signal at the score (from the signal copy, at least a thousandth of the largest
share). The noisiest components that carry the least signal come first. Each set is then scored
exactly in the sense of section 9d. With P = A W X the back-projection of all components
(`pop_subcomp`'s rebuild; X itself for R = ∅), removing R gives P − A_R S_R with S = W X the
activations, so adding *k* subtracts the rank-one term A_{:,k} S_{k,:} from the current data.
Epoching and baseline removal commute with it, and the known-signal copy is updated the same
way. The inner problem over *u* is solved exactly (section 9c). The path stops after 30
components (`MaxRemoved`), or when *J*(R_j) > (1 + 2*m*) min_{i≤j} *J*(R_i) with *m* the equivalence
margin (5 %): the path has passed its minimum by more than the margin covers.

*Why not forward selection* (adding each time the component that lowers *J* most, used before).
With about 60 candidates (ICLabel calling most components Other) and 30-40 trials per condition,
the greedy step picks components whose removal lowers the *sample* variance by chance, through
chance covariances with the rest of the score. The cross-fitting of section 9e re-chooses only
among the sets of the path, so it cannot see that. Simulated (fourth review, `icsme_path.py`,
35 trials per condition): the error the greedy path reported was 38 % below its true error, while
the cross-fitting showed 3 %, and its true error was 0.66 of removing nothing. With the order
above, the reported error is 0.94 of the true one and the true error is 0.54 of removing nothing,
for |C| + 1 evaluations instead of |C| · min(|C|, 30). Band power and
contralateral-minus-ipsilateral measures, whose parts are not linear in the components or not
computed, keep the forward selection.

**Time-locked components are not candidates.** The SME counts every trial-to-trial variation as
noise, real brain activity that varies from trial to trial included, so removing a brain component
lowers it too. ICLabel's Brain class protects only the components it recognises. So before the
path, for each candidate *k* and condition *c*, the component's part of each trial's score is
*d*_i = ā_k · (S̄_k,i(window) − S̄_k,i(baseline)), with ā_k the mean of its scalp map over the
ROI and S̄ the window and baseline means of its activation in trial *i* (exact for a mean
amplitude: the score is linear in the components). *k* is protected (`info.protected`) when in
some condition |*t*| = |mean *d*| / (SD *d* / √*n*) > 3 and |mean *d*| > 0.5 SME of the scores:
clearly time-locked, and large enough to change the measure. Eye components are not protected, as
blinks can be time-locked: those whose likeliest ICLabel class is Eye, and those whose activation
follows an eye channel (`Injection.eyeChannels`) far more closely than the other components' do.
That is the largest |correlation| with one over the epochs above 0.4 (`EyeCorr`) and more than 3
robust SDs (median, MAD) above the other components', as FASTER finds eye components. ICLabel
often calls a blink component Other: simulated (`evoked_guard.py`), a blink component that
ICLabel does not name passed the guard in 73 % of datasets at a 30 % blink rate. Along the path
the guard also holds for the components removed together: a component is skipped
(`info.protected`) when its part added to theirs (eye components left out) would pass it, so that
a P3 spread over several weak components is not removed piece by piece. Mean amplitudes only;
contralateral-minus-ipsilateral components and band power are not checked.

**Choice.** As in section 9d: the sets on the path are compared with the same paired bootstrap and
simultaneous max statistic, and among those shown equivalent to the best the one removing the
fewest components is used. The components are flagged in `EEG.reject.gcompreject` (the command
is recorded in `EEG.history`) and removed by `pop_subcomp`. Cost: at most min(|C|, 30) + 1
evaluations of the inner problem (forward selection: |C| · min(|C|, 30)), against at most *K* + 1
for the threshold.

**Scope.** As section 9d (per-trial scores; only epoching, baseline removal and an amplitude
rejection that removes epochs between the removal and the rejection; the same steps after it in
every pipeline below it). The objective is that of the measures chosen, so the set chosen lowers
their noise; a component that matters only for other measures or electrodes may stay. When
adopting a pipeline the components the search chose are reused; `apply` chooses again on the data
it is given.

## 9h. Eye channels regressed out (after the ranking)

[`pipecompare.run.Steps.eogRegress`; `pipecompare.run.Executor` eogCheck] With eye channels
*e* (dedicated EOG only, `Injection.dedicatedEyes`; Fp1/Fp2 are scalp EEG and are never used) and
EEG channels *t* (all scalp EEG except *e*), the step replaces *X*_t by *X*_t − *B* *X*_e, with

  *B* = *C*_te *C*_ee⁺,

the least-squares coefficients (pseudo-inverse: two eye channels can be nearly collinear). *C* is
the covariance over the samples where some eye channel is beyond 4 robust SDs or in a blink of
the hidden Markov model (section 9i) (median, 1.4826 ×
MAD) of its median, dilated by ±100 ms; with fewer than one second of such samples, over all
samples. Before *C* is formed the time-locked average of the epochs (the contract's epoch window
around each condition's events; on epoched data the average epoch) is subtracted from every
channel, so the ERP, which also reaches the eye channels, does not enter *B* and is not removed
with the blinks. *B* is fixed once fitted, so the step is linear and the signal-check copy
(section 4) replays exactly the same *B*. What is not removed (background EEG shared between the
eye channels and the scalp, which biases *B* slightly; simulation: residual blink error 0.33 µV of
2.37 µV before) is seen by the total error's blink term (section 15).

The step is tried, after the ranking, on the recommended pipeline only: it is inserted before the
pipeline's first step among IC removal, epoching, epoch rejection and weighting; the steps before
it repeat the search's decisions (the same ICA, the same interpolated channels), the steps from it
decide again. The new pipeline joins the candidates and the ranking (section 3) is made again, so
it is recommended only when its total error is lower. In the simulation of the plan
(`sims-2026-10-16c/eog_regression.py`, blinks time-locked in 40 % of 35 trials, 8 % of them reaching Pz) the RMS error of
the P3 mean amplitude at Pz fell from 2.37 to 0.92 µV. It is not combined with cutoff fine-tuning
(section 9e). On in the simple mode (`eogCheck`).

## 9i. Two-state hidden Markov models: blinks, and channels bad for part of the recording

[`pipecompare.eval.Hmm`] A sequence *y*_1…*y*_T is modelled as emitted by a hidden state
*s*_t ∈ {quiet, active}, a Markov chain with transition matrix *A* and initial probabilities π,
and *y*_t | *s*_t = *k* ~ N(μ_k, σ_k²). The parameters (*A*, π, μ, σ) are the maximum-likelihood
estimates by EM (Baum-Welch: scaled forward-backward recursions for the state probabilities γ_t
and the expected transitions ξ, then the closed-form updates; started at μ = the median and 99th
percentile of *y*, at most 30 iterations, until the log-likelihood gains less than 10⁻⁶ of
itself). The states are the Viterbi path, the exact maximizer of P(*s*_1…*s*_T | *y*) by dynamic
programming. The active state is the one with the higher mean. The model always splits a
sequence into two states, also one without any events, so a stretch of active samples counts only
when it passes a test of its own (below); everything else stays quiet. [`Hmm.segment`]
Emission log densities are normalized before exponentiation, avoiding underflow from isolated
extreme samples. A nonfinite fit returns an inactive mask with `model.valid = false` and a
reason, rather than decoding invalid parameters. Existing fixed detection rules remain in use.
At least ten finite features are required. Missing emissions have likelihood one in both states,
are excluded from emission-parameter updates and cannot be reported as detected events. Moving
eye-channel features omit missing samples and expose `inputMissingSamples`; this prevents one
missing value from disabling detection throughout an otherwise usable recording. This treatment
does not validate events inside gaps or assume that informative missingness is harmless.

**Blinks and eye movements** [`Hmm.events`; used by `Steps.withoutBlinks`, `Steps.eogRegress`,
`Measure.sampleHits`, each united with its fixed rule]. *y*: for each eye channel (or left minus
right pair), the channel minus its moving median over 1 s, a moving mean over 1/16 s, every
⌊f_s/25⌋-th sample, the robust z (median, 1.4826 × MAD) of its absolute deviation from the
median; the largest over the channels. A stretch counts when its largest *y* is at least 3.5
(`Hmm.PeakZ`) and it lasts at least 100 ms (`Hmm.MinSeconds`). Summing the evidence of all of a
blink's samples finds blinks that never pass 5 robust SDs in one sample. Simulated (Python,
`sims-2026-10-16d/eye_hmm_port.py` and a replica of the MATLAB code; noise SD 10 µV, 15 blinks
a minute of 15–200 µV, half sines of 250–400 ms): with 1/f noise, blinks of 15–30 µV found 64 %
(fixed rule 11 %), 30–50 µV 100 % (83 %); with AR(1) noise, 28 % (3 %); about 2 stretches a
minute marked in data without blinks (the test `testHmmFindsWeakBlinks` checks these).

**Channels bad for part of the recording** [`Steps.badStretches`, `Steps.interpStretches`;
badchannels option `stretches`, on in the simple mode and the panel for continuous data]. For
the tested channels not bad as a whole, with *P* the inverse of their covariance (on the
detection copy, extreme samples left out, as for `resid`), the residual of channel *i* from its
least-squares prediction by the others is (*P x*)_i / *P*_ii. In 1 s windows *w*, *L*_iw =
log SD of that residual, minus the median over the channels in that window (what they share:
movement). Each channel's sequence is centred on its 20th percentile plus the median channel's
gap between its 20th percentile and median, and scaled by the median channel's MAD (at least
0.05), so a channel noisier for most of the recording is still measured against its quiet part.
A stretch counts when it lasts at least 10 windows (`Steps.MinWindows`) and its mean is at least
3 (`Steps.StretchZ`). The channel is then interpolated in that stretch only, by EEGLAB's spherical
spline from the channels not bad there (the rows of `eeg_interp` read off an identity, as in
`repair` `'error'`); a fixed linear map once the stretches are known, so the signal-check copy
replays it exactly and a saved script repeats it (`Steps.interpStretches`). Not done on epoched
data, with fewer than 8 channels or 30 windows, or when the channels are linearly dependent.
The second output of `badStretches` reports whether detection was assessed and why it was skipped;
the step retains this as `stretchCheck`. Unassessed does not mean no bad stretches were present.
Simulated (`sims-2026-10-16d/bad_stretch_hmm.py`: 32 channels, 20 min, window log SDs with spread
0.12, Pz *k* = 1.5–3 times noisier over the last 10–50 %): every bad window found, none marked in
the rest (0.3 % of windows in recordings without a bad stretch); Pz's noise relative to a clean
recording 1.02–1.10, against 1.06–1.64 with the whole-recording test alone (interpolated
stretches taken as 1.2 times as noisy). The simulation's windows are more regular than real data.

## 10. Montage check (before a run)

Scalp potentials vary smoothly over the head, so after an average reference a channel resembles its
nearest neighbours. `pipecompare.live.Montage.check` measures this on the located EEG channels (not
EOG/ECG/..., not ear or mastoid sites; at least 8): a sample of the data (the first 600 s, or as
many epochs as fit in 600 s with each epoch's mean removed; at most 2·10⁷ values) is band-passed
1–30 Hz in the frequency domain and average-referenced, and for each channel *nn* is its mean
Pearson correlation with its 3 nearest channels (straight-line distance between unit position
vectors). With *m* the median and *s* = max(0.05, 1.4826 · MAD) over the channels, a channel is
flagged when (*nn* − *m*)/*s* < −3, or when *nn* < 0.1 and *m* ≥ 0.3 (on sparse caps the
neighbours are far apart and the absolute rule would flag good channels). Flagging is iterative:
the worst flagged channel is set aside and *nn* is computed again with neighbours chosen among the
remaining channels, so a swapped channel does not also lower its neighbours' *nn*; at most a
quarter of the channels are flagged. For each flagged channel the two channels with the highest
correlation are named; when none reaches 0.3 the channel resembles nothing and is more likely bad
than mislabelled. The check only warns; it never changes the data or blocks a run.

## 11. ICA check (after a run)

ICLabel's features include the components' scalp maps (Pion-Tonachini et al., 2019), so channel
labels that do not match the positions leave it unable to recognise components. After a run, the
recommended pipeline's ICLabel step (else the first pipeline's with one) is read: fewer than 2
components with Brain probability ≥ 0.5, or a median Other probability above 0.8, is reported as
"recognised almost none"; no component removed in any pipeline is reported too, with the montage
as the first thing to check when recognition was poor, and as expected on clean data when it was
not. ICA fitted on fewer than 20 × n² data points (n components; the data points the ICA step
had, all epochs together) is reported in any case, after the usual rule for a reliable infomax
decomposition (Onton & Makeig, 2006). [`pipecompare.simple.Presets.icaText`]

## 12. Measured electrodes that were rebuilt (after a run)

For the same reason, a whole measured electrode interpolated as a bad (or previously removed)
channel lowers the SME of the pipelines that interpolate it. No correction is sound here: the real
electrode's noise is not known, since it was bad. When the recommended pipeline interpolated a
measured electrode, the result says so and suggests measuring at other electrodes if it really is
bad. [`pipecompare.simple.Presets.roiText`]

## 13. Channels behind the rejected epochs (after a run)

Each rejection step records, per channel, in how many of the epochs it rejected that channel was
over the limit (`overLimit`, summed over the pipeline's rejection steps). When the recommended
pipeline has a channel over the limit in at least half of its rejected epochs, and in at least 3,
the result names it (up to 3 such channels) as likely bad over stretches the bad-channel test,
which scores the whole recording, does not catch. When no pipeline passed, the same counts give
`nextStep`'s list instead (channels with at least 10% of all over-limit marks).
[`pipecompare.simple.Presets.rejectText`]

## 14. One pipeline for a group of recordings

[`pipecompare.Group`] Each recording *s* = 1…*S* is compared as in sections 1–13 (the simple mode's
settings), its decisions made on its own data. For the pipelines *k* run on every recording, with
*J*_{s,k} the gain-corrected objective of *k* on *s*: *k* passes when it passes on every recording
(recordings on which nothing passes are left out); the group objective is
*J*_k = √(mean_s *J*_{s,k}²). This is an RMS of individual measurement errors in the observed
group, not the MSE of the group mean or a group contrast. Independent zero-mean noise contributes
Σ SME_s² / *S*² to a group mean's variance; systematic same-direction biases instead contribute
(mean_s bias_s)² and need not average away. The paired bootstrap replicates of each
recording (section 3, its own trials) are combined the same way, *J*_k^(b) = √(mean_s
(*J*_{s,k}^(b))²), with independent draws per recording; the confidence and equivalence sets and
the choice among those not distinguished from the best (highest minimum retention over the
recordings, then smallest largest amplitude error, then lowest objective) follow section 3. The
uncertainty is that of these recordings' trials, not of sampling other participants: the choice is
for this group's data, as the single-recording choice is for one recording's.
Each recording uses a distinct reproducible bootstrap seed (`G.bootstrapSeeds`), paired only
across candidates within it. Unique recording names retain their stream under input reordering;
duplicate names are disambiguated by occurrence.

**Each recording's own pipeline (two-level empirical Bayes).** [`pipecompare.Group.twoLevel`] Let
*Y*_{s,k} = *J*_{s,k} for the pipelines passing on *s* (others missing) and *V*_{s,k} the variance of
its bootstrap replicates. On the log scale, *y*_{s,k} = log *Y*_{s,k} = *a*_s + *c*_k + *u*_{s,k} +
ε_{s,k}: a recording effect (how noisy it is), a pipeline effect (how good it is on average), what
is particular to this pair, *u* ~ (0, τ²), and estimation noise of variance *v*_{s,k} =
*V*_{s,k}/*Y*_{s,k}² (delta method). The observation covariance Σ_s also includes off-diagonal
terms from the same paired bootstrap: Cov(J_s,k,J_s,l)/(Y_s,k Y_s,l). Shared estimation noise
can cancel in pipeline contrasts; treating its diagonal as independent noise can erase real
recording-specific differences.

Let *D* be the recording/pipeline design on the observed cells, *H* its least-squares projection,
*R* = *I* − *H*, and Σ the block-diagonal matrix of Σ_s (recordings assumed independent).
With *r* = *R y*, *df* = rank(*R*), estimate

  τ̂² = max(0, (rᵀr − tr(RΣ)) / df).

The design drops one gauge per connected component of the observed recording/pipeline graph:
*df* = cells − (active recordings + active pipelines − components). A disconnected missingness
pattern cannot be treated as one connected design. If *df* = 0, the interaction variance is
unassessed; the observed estimates are retained and the diagnostic says so.

Given τ̂², additive effects are fitted by generalized least squares under C_s = τ̂² I + Σ_s.
The conditional empirical-Bayes estimate is

  θ̂_s = exp(D_s β̂ + τ̂² C_s⁻¹ (y_s − D_s β̂)).

This reduces to scalar residual shrinkage for diagonal Σ_s. Singular covariance blocks use a
small numerical diagonal regularizer. Covariance inputs are checked for symmetry, matching
diagonals and positive semidefiniteness. The trace correction and block gains avoid counting
shared recording noise as evidence against pipeline-specific differences.

Recording *s* gets argmin_k θ̂_{s,k} over the pipelines passing on it (`G.own`, `G.ownShrunk`);
`G.share` is the mean diagonal of the posterior gain matrices. `G.twoLevelModel` reports
covariance use, residual degrees of freedom, disconnected components and validation scope.
Effects, covariance and τ̂² are estimated on these recordings: this is an empirical-Bayes
plug-in estimate, not a calibrated probability or proof of benefit on new participants.
The independent-noise simulation remains a regression; paired-noise and missing-design
regressions cover the failures above. Broader independent validation is still required.

## 15. Total error: noise and bias together

[`pipecompare.eval.Rank` step 3b; `pipecompare.eval.Measure.blinks`, `Measure.eyeBias`;
`pipecompare.eval.Injection.shapeSpreads`] The SME (section 1) measures how much a mean would vary
from one set of trials to another. It does not see an error that moves every set the same way.
The objective is therefore the root mean square error of each gain-corrected condition mean,

  *J* = √(mean over measures and conditions of (SME/*g*)² + *b*²),

with *b* = (*e* + |ȳ| *s*) / *g* the bias of that mean (columns `noise` = SME/*g* alone and `bias` =
*b* alone; *J*² = noise² + bias²). It has two parts, each computed from the data and the
pipeline's own processing.

**What blinks add (*e*).** A blink (or other eye movement) time-locked to the stimulus, as after
targets in an oddball or n-back task, falls in a share *p* of the trials and adds about *d* µV to
their score at the measured electrodes. The mean moves by *p d*; the SME grows only by the spread
*d*·√(*p*(1−*p*)) adds to the trials' spread, divided by √*N*. Simulated (40 trials per condition,
SME 0.94 µV): with *p* = 0.3 and *d* = 3 µV the SME grows 3 % while the mean moves 0.9 µV, one
SME; with *p* = 0.5, 3 % and 1.6 SME. A choice by SME alone keeps such trials (within the
equivalence margin, the least aggressive pipeline is recommended), and epoch rejection does not
test the EOG channels.

The blink trials are found once, on the data as they are before any processing: for each
mean-amplitude component, an epoch's eye activity is the largest, over the eye channels (EOG by
type or name, else Fp1, Fp2, FPz; never a measured electrode), of |window mean − baseline mean|,
the part of an eye artifact that reaches a score; it is a blink when its robust *z* (median and
MAD over all trials) exceeds 3 (`Measure.BlinkZ`). Comparing the trials with each other fails
when most of them have a blink, as in n-back tasks with a response on every trial: the median
trial then has one (simulated, fourth review `eye_bias.py`: at 50 % blink trials the bias came out
0.69 µV for a true 1.44 µV, at 70 % 0.17 µV for 2.0 µV). So in continuous data a trial is also a
blink trial when a blink found sample by sample falls in its window or baseline: an eye channel
minus its moving median over 1 s beyond 5 robust SDs over the samples (blinks fill far fewer
samples than trials), or in a blink found by the hidden Markov model of section 9i, with 50 ms on
either side; eye movements the same on *h* below
(`Measure.blinks`, `sampleHits`). Then, for each pipeline and condition, over the
trials it keeps (with their weights *w*, when weighted), with *p* the blink trials' share of the
weight and *m*_B, *m*_N the weighted mean scores of the blink and other trials,

  ȳ = *m*_N + *p* (*m*_B − *m*_N),  *e* = *p* (*m*_B − *m*_N):

what the blink trials add to the mean (a pipeline that removed the eye activity, by ICA or by
rejecting those trials, has *e* ≈ 0). *e* is itself an estimate, with standard error
*se* = *p* *s* √(1/*n*_B + 1/*n*_N) with equal weights (*s* the pooled scores' SD); with most trials blinking it is often larger than
the SME itself. Only the part beyond what chance makes counts: (|*e*| − 2 *se*)₊², soft thresholding
(`test_statistics/testWhatBlinksAddToTheMean`). Taking only the variance off, *e*² − *se*², floored at
0, counts a bias where there is none, since its positive part has a positive mean: in simulation
with 35 trials, no real effect and 70 % blink trials it added a third to the total error on
average (twice as much in one run in ten) and pulled the rejection limit towards rejecting blink
trials; soft thresholding adds about 1 %. A real bias is underestimated by up to 2 *se*, the price
of not acting on noise. In the
bootstrap (section 3) *e* is recomputed in every resample.
Weights are normalized only through ratios. With physical multiplicities c_i and combined
weights K_i = c_i w_i, q_i = K_i²/c_i; the weighted pooled variance divides the second central
moment by 1−Σq_i/(ΣK_i)², and the standard error uses q_B/(Σ_B K_i)² + q_N/(Σ_N K_i)².
Resampling counts are supplied separately from trial weights. Multiplying all weights by a
constant leaves the bias, its standard error and the ranking unchanged, while drawing a trial
twice still counts two physical observations.

A difference of the blink trials'
mean that is real brain activity (a larger P3 in the trials with a blink) is counted as bias too;
using the blink's size only to select the trials, not as a regressor, keeps frontal brain activity
that correlates with the measured electrodes from being read as eye activity.
When every retained trial has a detected artifact, the artifact-free comparison is absent:
`eyeBias` and `eyeBiasSq` are NaN, `eyeBiasAssessed` is false, and the candidate cannot pass the
ranking with a fictitious zero bias. With no detected artifact, the estimate is zero conditional
on the detection rule. See [validation scope and experiments](VALIDATION.md).

**Eye movements and muscle.** Two more kinds of trials are found the same way (robust *z* > 3 over
all trials, per component) and counted with the blinks (a trial with any of them is in B):
- *eye movements*: |window mean − baseline mean| of *h* = left − right, the difference of two eye
  channels on either side (EOG channels by position, EEGLAB's *Y* axis pointing left, |*Y*|/‖·‖ >
  0.2; else by names ending in L and R; else F7 − F8, when neither is measured). A blink, about the
  same on both sides, cancels in *h*; a horizontal movement, opposite on the two sides, doubles.
  `Measure.blinks` reports it as `blinkInfo.horizontal`;
- *muscle*: log of the largest, over the component's measured electrodes, SD of the first
  difference *x*ₜ₊₁ − *x*ₜ in the window: a high-pass that the slow ERP hardly passes and EMG
  (20–150 Hz) does.

Since *e* = *p*(*m*_B − *m*_N) = ȳ − *m*_N, it is the distance of the mean from the mean of the trials
with none of the three; trials flagged by chance carry no bias and only add to its variance, which
is taken off. `blinkInfo.kindShare` gives the share of each kind.

**What the signal's shape changes (*s*).** The gain *g* (section 2) is measured with one known
shape (section 4). The real component's shape is not known, and a filter changes the gain of a
broad or slow component more than that of a narrow one. Simulated with EEGLAB's default
high-pass design, the window mean of a P3 at 300–600 ms after 0.5 Hz: −3.3 % for the bump, −3.4 %
for a recording's own average (faded outside the window), −5.5 % for a P3 with a slow tail to
900 ms; after 1 Hz: −20 %, −20 % and −28 %. So *s* is the largest relative change of the gain over a
family of shapes: the template used and bumps as wide at half maximum as 0.5, 1 and 2 times the
window (at least 50 ms) (`Injection.ShapeFactors`), the worst case over the family, a minimax.
The bumps are whole, out to 4 SDs either side of the window's centre, even where that is
outside the epoch. Cut at the epoch's edge, the widest bump would end in a step that a high-pass
passes. In simulation (fourth review, `shape_spread.py`, 300–600 ms window, epoch −0.2–0.8 s)
the cut bumps gave *s* = 9.1 % instead of 13.2 % at 1 Hz, and 0.9 % instead of 1.4 % at 0.5 Hz.
Only the filters make the gain depend on the time course: re-referencing, interpolation and IC
removal act on the scalp field, and rejection and weighting select trials. So *s* is computed
with each pipeline's filters alone on a short copy of the signal, once per filter sequence; *s* = 0
without filters, for latencies, band power and epoched data. *b* uses |ȳ| *s*, the change of the
measured mean if the gain were off by *s*.

**Inside the pipelines.** The rejection limit and the ICLabel threshold chosen from the data
(sections 9c and 9d) and the components chosen by their effect on the noise (9g) minimize the
same total error, with *e* (the shape part depends on the filters, not on these choices, so it is
left out there); the selection check (section 3) makes them again on half of the trials with it.
Fine-tuned cutoffs (9e) use the *s* of each pair of cutoffs.

**Recommendation.** Among the pipelines not distinguished from the best, one whose bias is above
half its noise (`Rank.BiasShare`) is recommended only when none below it is among them: keeping
more trials is no reason to keep what they carry in.

**Simulation.** `test_engine/testBlinksCountInTheTotalError`: blinks after 30 % of the targets;
the estimated *e* is within 3 standard errors of what the blinks add in fact (the same data
without them), rejecting the blink trials leaves almost none, and that pipeline is recommended.
`test_engine/testEyeMovementsAndMuscleCountInTheBias`: leftward eye movements in 30 % of the targets
and muscle bursts in 15 % of all trials, measured at F3: over 80 % of each are found, under 10 % of
the other trials, and *e* is within 3 standard errors of what the eye movements add.
`test_engine/testShapeSpreadGrowsWithTheHighPass`: *s* is 0 without a filter, below 2 % at 0.1 Hz,
above 10 % at 1 Hz.

**Limits.** Without an eye channel (EOG, or a frontal pole that is not measured) *e* is not
estimated, and the result says so. Blinks in more than half of the trials defeat the robust *z*.
Eye activity that does not reach the eye channels (small saccades) is not seen. The shape family
covers widths, not every shape.

## References

- Clayson, P. E., Baldwin, S. A., Rocha, H. A., & Larson, M. J. (2021). The data-processing
  multiverse of event-related potentials (ERPs): A roadmap for the optimization and standardization
  of ERP processing and reduction pipelines. *NeuroImage, 245*, 118712.
- Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their Application*. Cambridge
  University Press.
- Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG
  dynamics including independent component analysis. *Journal of Neuroscience Methods, 134*(1),
  9–21.
- Efron, B. (2014). Estimation and accuracy after model selection. *Journal of the American
  Statistical Association*, 109(507), 991–1007.
- Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall.
- Hall, P., Horowitz, J. L., & Jing, B.-Y. (1995). On blocking rules for the bootstrap with dependent
  data. *Biometrika, 82*(3), 561–574.
- Hansen, P. R., Lunde, A., & Nason, J. M. (2011). The model confidence set. *Econometrica, 79*(2),
  453–497.
- Hoke, M., Ross, B., Wickesberg, R., & Lütkenhöner, B. (1984). Weighted averaging: theory and
  application to electric response audiometry. *Electroencephalography and Clinical
  Neurophysiology*, 57(5), 484–489.
- Kane, N., Acharya, J., Beniczky, S., Caboclo, L., Finnigan, S., Kaplan, P. W., et al. (2017). A
  revised glossary of terms most commonly used by clinical electroencephalographers and updated
  proposal for the report format of the EEG findings. *Clinical Neurophysiology Practice, 2*,
  170–185.
- Kappenman, E. S., Farrens, J. L., Zhang, W., Stewart, A. X., & Luck, S. J. (2021). ERP CORE: An
  open resource for human event-related potential research. *NeuroImage, 225*, 117465.
- Kothe, C. A. E., & Makeig, S. (2013). BCILAB: a platform for brain-computer interface development.
  *Journal of Neural Engineering, 10*(5), 056014 (artifact subspace reconstruction).
- Künsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations.
  *The Annals of Statistics, 17*(3), 1217–1241.
- Lopez-Calderon, J., & Luck, S. J. (2014). ERPLAB: an open-source toolbox for the analysis of
  event-related potentials. *Frontiers in Human Neuroscience, 8*, 213.
- Luck, S. J., Stewart, A. X., Simmons, A. M., & Rhemtulla, M. (2021). Standardized measurement
  error: A universal metric of data quality for averaged event-related potentials.
  *Psychophysiology, 58*(6), e13793.
- McCarthy, P. J. (1969). Pseudo-replication: Half samples. *Review of the International
  Statistical Institute, 37*(3), 239–264.
- Nolan, H., Whelan, R., & Reilly, R. B. (2010). FASTER: Fully Automated Statistical Thresholding for
  EEG artifact Rejection. *Journal of Neuroscience Methods, 192*(1), 152–162.
- Onton, J., & Makeig, S. (2006). Information-based modeling of event-related brain dynamics.
  *Progress in Brain Research*, 159, 99–120.
- Percival, D. B., & Walden, A. T. (1993). *Spectral Analysis for Physical Applications: Multitaper
  and Conventional Univariate Techniques*. Cambridge University Press.
- Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F. (1989). Spherical splines for scalp
  potential and current density mapping. *Electroencephalography and Clinical Neurophysiology,
  72*(2), 184–187.
- Pion-Tonachini, L., Kreutz-Delgado, K., & Makeig, S. (2019). ICLabel: An automated
  electroencephalographic independent component classifier, dataset, and website. *NeuroImage,
  198*, 181–197.
- Romano, J. P., & Wolf, M. (2005). Stepwise multiple testing as formalized data snooping.
  *Econometrica, 73*(4), 1237–1282.
- Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power
  approach for assessing the equivalence of average bioavailability. *Journal of
  Pharmacokinetics and Biopharmaceutics, 15*(6), 657–680.
- Slepian, D. (1978). Prolate spheroidal wave functions, Fourier analysis, and uncertainty — V: The
  discrete case. *Bell System Technical Journal, 57*(5), 1371–1430.
- Steegen, S., Tuerlinckx, F., Gelman, A., & Vanpaemel, W. (2016). Increasing transparency through
  a multiverse analysis. *Perspectives on Psychological Science, 11*(5), 702–712.
- Thomson, D. J. (1982). Spectrum estimation and harmonic analysis. *Proceedings of the IEEE, 70*(9),
  1055–1096.
- White, H. (2000). A reality check for data snooping. *Econometrica, 68*(5), 1097–1126.
- Widmann, A., Schröger, E., & Maess, B. (2015). Digital filter design for electrophysiological
  data – a practical approach. *Journal of Neuroscience Methods, 250*, 34–46.
- Zhang, G., Garrett, D. R., & Luck, S. J. (2024). Optimal filters for ERP research I: A general
  approach for selecting filter settings. *Psychophysiology, 61*(6), e14531.


## Missing paired bootstrap replicates

A simultaneous maximum uses the original number of resampling draws. If a candidate comparison
is unavailable in a draw, that family's maximum for the draw is positive infinity. Sequential
equivalence tests use the same rule for their pair. The upper quantile includes those infinities:
missing tail mass cannot narrow a bound or establish equivalence. With insufficient supported
tail probability, the bound is unbounded. The identity comparison with the candidate itself is zero.

The factorial empirical-Bayes covariance fit requires complete paired support for its modeled
candidates; otherwise it returns the unshrunk point scores. Selection frequencies use all draws
as their denominator and separately retain the share with no assessed recommendation. These
rules concern conditional comparisons of the fitted candidate family and do not establish
independent calibration of measurement error.

### Artifact-association variance correction (PR #50 optimization)

Let W be total trial weight, p the weighted artifact share, d the difference
between marked and other-trial means, and Q_i the squared reliability weight
per physical observation (bootstrap multiplicities are divided out). The
association estimate is p*d. Its updated delta-method variance separates
pooled within-group variation from uncertainty in p:

    Var(p*d) approximately p^2 * s_within^2 * (Q_B/W_B^2 + Q_N/W_N^2)
                          + d^2 * Var(p)
    Var(p) = [Q_B*(1-p)^2 + Q_N*p^2] / [W^2 - Q_B - Q_N]

Weighted within-group residuals and degrees of freedom supply s_within.
The old total-variance expression mixed mean separation into contrast noise
and scaled prevalence uncertainty incorrectly. Removing every between-group
contribution would also be wrong: p is estimated, not fixed. This remains a
delta approximation under common within-group variance and fixed fitted
labels/weights, not a guarantee under heteroscedasticity or learned
preprocessing. Analytic flat-group and independent-record Monte Carlo tests
cover several artifact prevalences and association strengths.

The explicit legacy-total mode reproduces the previous convention for frozen
before/after controls. New objectives store their variance convention;
historical score checkpoints without it retain legacy bootstrap behavior.
Automatic rejection/component-choice bootstrap now passes physical counts
to the bias estimator; the legacy control preserves its previous convention.
Both corrections can change adaptive choices and require whole-procedure
outcome validation, not only rescoring selected outputs.


### Unequal artifact-group noise and guarded upgrades

The optional `artifactvariance='heteroscedastic'` mode replaces the shared
within-group variance in the expression above by separate weighted variances:

    Var(p*d) approximately p^2 * (s_B^2 * Q_B/W_B^2 + s_N^2 * Q_N/W_N^2)
                          + d^2 * Var(p)

Each variance uses its own residual degrees of freedom, W_g-Q_g/W_g.
A group without residual degrees of freedom is unassessed, including a lone
clean comparison trial. Reliability-weight scaling and physical bootstrap
multiplicities retain their separate meanings. This is a Welch-type delta
approximation conditional on fitted labels and weights. It corrects the
common-variance assumption; it does not identify a causal artifact effect,
validate the soft-thresholded bias proxy, or establish better EEG recovery.
The pooled default and explicit legacy-total control remain reproducible.
A change to the default requires paired whole-procedure validation.

`Rank.run(..., struct('upgradeBaseline', baselineId), ...)` optionally
protects a designated candidate in the same recording and comparison
stratum. The baseline must be feasible. A proposed replacement is accepted
only when its score difference plus the simultaneous all-pairs bootstrap
critical value is strictly negative. Otherwise the baseline is retained.
The same guard is used in conditional bootstrap selection frequencies.
An incomplete family can yield an infinite bound and cannot justify an
upgrade. Cross-fit diagnostics remain those of the underlying unguarded
score selector and explicitly carry `upgradeGuardApplied=false`.
This guards promotion on the fitted score; it is not a guarantee about
unknown true MSE or independent generalization. Declare the baseline
without consulting ground truth. The option is deliberately explicit.

Synthetic benchmarks now preserve compact per-trial score archives outside
disposable EEG/ICA checkpoints. These contain candidate statistics, factor
levels, reference trial identities and truth/provenance metadata, so future
ranking audits can reuse fits. A replay is a conditional development
diagnostic, not a new independent validation or a rerun of adaptive
preprocessing choices.

The completed four-seed development comparison on revision
4116ea11f3d5d39e4d3d282894dcf487172bd50d remains failed at its original
admission threshold: balanced selected MSE increased 53.84%. Seed 40103
contributes 1.04705 of the 1.04901 summed paired MSE increase. Saved-score
diagnostics reject a blanket switch from the two-SE soft threshold to
positive-part bias-squared minus variance: its raw-score-minimum MSE is
1.05960, versus 0.74934 for the saved balanced recommendations. This comparison
is descriptive (the recommendation rule also uses shrinkage and equivalence),
and supplies no justification for changing the default or opening confirmation.
