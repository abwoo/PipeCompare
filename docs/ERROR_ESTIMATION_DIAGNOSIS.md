# Error-estimation diagnosis: two hypotheses checked

[Cloud execution](https://github.com/abwoo/PipeCompare/actions/runs/41641303598), MATLAB R2026a Update 5, unchanged production `7151821edf8a5a54c3167d7f18b66c530911b568`. Execution recipe `a7d01ebaa0c3ff042159d3f550a94d210c52aec4`. No production change or merge.

## Design

The actual `Measure.eyeBias` implementation was called in six prespecified cases, 4,000 repetitions each, seed 520113. Forty independent unit-variance Gaussian scores have fixed labels (20 blink, 20 other). Group shifts are 0, 3 or 8. Analytic weights are either equal or have one 9x-weight trial in each group, normalized to mean one as production weights are. Known bias is weighted blink share times group shift. Exact SE is calculated from the fixed weights and known unit variance, independently of the implementation.

A within-group, analytic-weight-aware pooled-variance SE was evaluated as a diagnostic alternative only. Its weighted residual degrees of freedom and variance of weighted group means use squared weights. This assumes equal within-group variances and fixed, exogenous weights. Real data-derived weights, label error and signal/eye confounding are not covered.

## Results

SE columns report root-mean-square estimated standard error, not end-to-end ERP reconstruction error.

| Weights | Group shift | Exact SE | Empirical SE | Original RMS SE | Diagnostic alternative RMS SE |
|---|---:|---:|---:|---:|---:|
| Equal | 0 | 0.158114 | 0.156079 | 0.157919 | 0.157966 |
| Equal | 3 | 0.158114 | 0.156552 | 0.287540 | 0.157970 |
| Equal | 8 | 0.158114 | 0.155871 | 0.659274 | 0.158288 |
| Unequal, normalized | 0 | 0.252538 | 0.249770 | 0.155229 | 0.253249 |
| Unequal, normalized | 3 | 0.252538 | 0.254813 | 0.285001 | 0.252315 |
| Unequal, normalized | 8 | 0.252538 | 0.252459 | 0.659841 | 0.253091 |

**First hypothesis confirmed under the stated model:** overall variance includes between-group signal. With equal weights and shift 3, the original SE is 82% above the exact SE; shift 8 raises it over fourfold. This decreases the soft-thresholded bias penalty. The behavior is consistent with an intentionally conservative design, not by itself proof that removing the conservatism improves selection. Under the equal-weight null, mean squared penalty increases from 0.000186 to 0.000309 with the alternative. That tradeoff must be retained in an end-to-end evaluation.

**Second hypothesis confirmed for analytic weights:** sums of weights do not represent independent sample counts. Normalization to mean one does not fix this: under the unequal-weight null the original SE is about 39% too low. A diagnostic interval `bias +/- 1.96 SE` covers the known zero bias in 78.03% of replications, compared with 96.8% for the alternative. These intervals are diagnostic probes, not a claim about an existing product confidence-interval display, and the alternative is not proven exactly calibrated.

Multiplying every analytic weight by 0.1, 1 or 10 leaves the estimated bias at 1.5, yet changes the original squared penalty from 0 to 1.039616 to 1.821993. A fixed weighted mean's uncertainty should be invariant to such scaling.

**Frequency counts are different:** the original function exactly matches physically expanding integer-count observations. Bootstrap counts must keep that interpretation. Replacing every combined weight/count matrix with an effective sample size would change bootstrap semantics incorrectly. A repair should carry analytic weights and resampling multiplicities separately through point scores, threshold search, bootstrap and crossfit.

## Relation to the ten-participant pilot

The selected pipeline explicitly weights trials only for participant 4. Other participants still compare weighted candidates, so this does not exclude indirect selection effects. From the saved selected calibration rows, 19 of 20 condition readouts have zero soft-thresholded eye-bias penalty (after removing the separately reported shape component); only participant 4 condition 2 has a nonzero term. This establishes where the penalty is active, not whether its zeros are erroneous.

The existing artifacts contain aggregate calibration rows but not all trial-level score/weight/label arrays needed to recompute corrected scores. No improved ten-participant accuracy is claimed from these probes. A corrected implementation must be freshly fitted and compared with the original on all ten records before retention, with the same revert-on-worsening rule and no PR merge.

Real neural responses can covary with blinks, including in the frozen pilot generator. A blink/non-blink contrast cannot identify artifact bias alone. Fixing SE arithmetic does not resolve this scientific identification limit.

## Evidence and disposition

[Executable diagnostic](../benchmarks/diagnose_error_estimation.m) and [complete numeric results](../benchmarks/results/error_diagnosis_20261103/diagnosis.json). Both hypotheses were tested; production is retained unchanged. The most clearly specified next implementation target is separating analytic weights from frequency counts. The within-group variance alternative requires its own false-penalty and signal-preservation evaluation rather than automatic promotion.
