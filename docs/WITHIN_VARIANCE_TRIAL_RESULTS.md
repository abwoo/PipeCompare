# Within-group variance experiment: result and retention decision

**Retain this exact version as an exploratory candidate, not as the default production algorithm.** It meets the decision rule frozen before execution: all ten pairs valid, lower mean primary MSE, and no increase in mean raw MSE. PR #51 remains open and unmerged. The tested `Measure.m` is preserved under `benchmarks/within_variant/`; `+pipecompare/` remains unchanged from main. This retention decision is not a general superiority or release claim.

[Completed cloud run](https://github.com/abwoo/PipeCompare/actions/runs/41737985970), execution `f6265967ed5538629795ad308ee0baf3660dbf20`. Ten distinct public EEG participants, the same sixty-second excerpts and known-target injection protocol as the original pilot. These participants were previously examined: this is development reuse, not independent confirmation and not a real task-ERP ground truth.

## Paired effects

All MSEs are in uV², averaged over conditions within participant and then equally over ten participants. Lower is better.

| Metric | Original Standard | Within-group variant |
|---|---:|---:|
| Mean gain-corrected population-target MSE | 0.48315319393948364 | 0.33909463755161573 |
| Mean raw MSE | 0.5239750426980665 | 0.3767493017789699 |
| Mean reported squared error | 0.3557354206159417 | 0.34908495774354753 |
| Actual/reported RMS ratio | 1.1654103543938663 | 0.9855868340182047 |
| Stable recommendation over recorded ICA starts | 9/10 | 9/10 |
| Maximum native/export replay error | 0 / 0 | 0 / 0 |

Primary MSE decreases **29.82%**, raw MSE **28.10%**. There are **2 improved, 1 worse, 7 tied** participants. The mean paired difference is **−0.1440585564**, with the frozen 10,000-resample participant-bootstrap 95% interval **[−0.4516028726, +0.0211657481]**. It crosses zero: an improvement in this sample is observed, but a general advantage is not established. The frozen retention rule uses the observed mean and raw-error guardrail, not statistical significance.

| Participant | Original MSE | Variant MSE | Difference (variant − original) |
| S001 | 1.393875181 | 1.393875181 | +0.000000000 |
| S002 | 0.017795836 | 0.017795836 | +0.000000000 |
| S003 | 0.240220901 | 0.240220901 | +0.000000000 |
| S004 | 1.750526001 | 0.250978241 | -1.499547760 |
| S005 | 0.256920126 | 0.333267769 | +0.076347643 |
| S006 | 0.416617982 | 0.399232535 | -0.017385447 |
| S007 | 0.012026996 | 0.012026996 | +0.000000000 |
| S008 | 0.083435929 | 0.083435929 | +0.000000000 |
| S009 | 0.043610381 | 0.043610381 | +0.000000000 |
| S010 | 0.616502606 | 0.616502606 | +0.000000000 |

The gain is concentrated in **S004** (1.750526 → 0.250978). Excluding S004 as a sensitivity calculation, the other nine participants' mean difference is **+0.0065513551**, slightly worse. This sensitivity calculation does not replace the prespecified ten-person analysis. **S005 worsens** (0.256920 → 0.333268), while S006 improves modestly (0.416618 → 0.399233).

The aggregate RMS ratio moves closer to one, but ten development records do not establish uncertainty calibration; matching the aggregate ratio can conceal individual miscalibration. The known analytic-weight semantics issue is unchanged in this isolated experiment.

The two fixed controls reproduce their estimates within each pair. Their mean primary MSEs remain approximately 2.352045 (minimal) and 1.172762 (common ICA). The variant's lower point average does not by itself prove broad superiority over those methods.

## What changed and where the gain came from

Only the variance estimate inside `Measure.eyeBias` changes: pooled residual variance within blink/non-blink groups replaces variance around the overall mean. The blink/non-blink mean difference no longer contributes to the noise estimate used to discount that same difference. BiasZ=2, factorial shrinkage, and original weight/count semantics remain unchanged. The preceding analytic-weight correction was neutral on selected accuracy and is not combined with this variant.

The observed selected recipe changes for S002, S004, S005 and S006; S002's estimates and MSE nevertheless tie exactly. S004 switches from trial weighting to an EOG-regression recipe and loses the large original error on this constructed target. S005 switches filter settings in its EOG-regression recipe and worsens. S006 changes its low-pass setting and improves modestly. The single scoring change also affects the adaptive search, so this is a comparison of complete procedures, not a reranking of an identical fixed candidate list. These outcomes do not prove a causal neural/artifact separation mechanism.

There is a consequential stability tradeoff: S004 changes from stable recommendations [1,1,1] to [25,25,1] across the recorded ICA starts; S006 changes from unstable [6,25,2] to stable [2,2,2]. The overall count stays 9/10, but the participant driving the accuracy benefit is now unstable. The seed-check output does not provide ground-truth accuracy for every alternative start, so the 29.82% gain must not be presented as seed-robust.

## Verification and limits

All ten record jobs and the collector succeeded. Every worker logged `WITHIN_VARIANCE_CHECKS_PASSED meanNullPenalty=0.000655009416944`. The check covers an independently calculated within-group pooled variance, common-offset invariance, noiseless-group behavior, integer-frequency expansion, absent/all-blink behavior, 2,000 null simulations with the existing `<0.05/40` tolerance, and the unchanged existing factorial-shrinkage selection simulation.

The original `testWhatBlinksAddToTheMean` contains expectations for the old overall-variance formula and was not run for this experiment. The new closed-form and null checks validate the experimental formula; they are not a claim that the full legacy MATLAB suite passed. Production promotion would require legitimate updates to those formula-specific expectations and appropriate regression checks. They were not silently weakened or removed.

The cloud collector validates all ten paired records, matching sources/truth/options, matching EEGLAB, unchanged fixed-control estimates, baseline reproduction to 1e-7, and replay/export to 1e-6. Downloaded estimates/gains/truth were independently recomputed and match every reported ERP/control MSE and raw MSE. Native and executable-export differences are exactly zero for both arms in all ten records. The variant's file hash matches its frozen plan. Identical deterministic unit checks repeated on ten workers are one repeated check, not ten independent statistical validations.

Unresolved limits include analytic weights versus frequency counts, the pooled-variance assumptions, uncertainty with very few trials, and genuine neural differences correlated with blinking. The semisynthetic protocol includes a neural amplitude increase on blink trials, so a blink contrast cannot identify artifact causally. None of the existing short-fault detection or real-task validation limits is resolved here.

## Retained evidence and next decision

The exact tested variant, checker, frozen plan, workflow and ten paired records are retained on PR #51. Production stays original; no merge, release, or independent-confirmation opening is performed. This preserves the observed improvement without treating a one-participant-driven result as ready for default use.

The next focused question is whether S004's gain survives ICA restarts and independent short records, while explicitly tracking S005's regression. Freeze this candidate before that validation; do not choose pipelines per person using known truth. The present experiment is complete and no additional experiment is silently started by this report.

[Summary](../benchmarks/results/within_trial_20261104/summary.json) · [All paired records](../benchmarks/results/within_trial_20261104/records.json) · [Frozen plan](../benchmarks/results/within_trial_20261104/plan.json) · [Verification](../benchmarks/results/within_trial_20261104/verification.json) · [Exact retained candidate](../benchmarks/within_variant/Measure.m)

Executed-source identity: baseline `7151821edf8a5a54c3167d7f18b66c530911b568`; trial execution `f6265967ed5538629795ad308ee0baf3660dbf20`; the variant worker copies `benchmarks/within_variant/Measure.m` into the production path before starting a fresh MATLAB session. The variant record's algorithmCommit identifies this recipe, not an untouched checkout. Later evidence-only commits are not claimed as a newly executed full-platform validation.
