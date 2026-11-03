# Analytic-weight correction: paired result and decision

[Cloud run](https://github.com/abwoo/PipeCompare/actions/runs/41645985391), execution `7061918c44f23c1a0172503506ea4d280360d447`. All ten pairs completed after retrying participant 3's MATLAB installer crash. The crash occurred before algorithm execution; the missing-record collector failure was consequential. Successful original artifacts were preserved and the missing pair was supplied on retry.

## Decision

**Retain original production under the frozen improvement rule.** The correction passes weight-semantics and relevant existing statistical checks, but has no observed selected-accuracy benefit on these ten records. Its experimental code remains available for further review; this is not evidence that the diagnosed weight-semantics defect is imaginary or unimportant in other inputs.

| Metric | Original | Weight correction |
|---|---:|---:|
| Mean gain-corrected population-target MSE | 0.48315319393948364 | 0.48315319393948364 |
| Mean raw MSE | 0.5239750426980665 | 0.5239750426980665 |
| Mean reported squared error | 0.3557354206159417 | 0.3554235314423144 |
| Actual/reported RMS ratio | 1.1654103543938663 | 1.1659215742987483 |

All ten selected pipeline strings are identical between versions. Zero participants improved, zero worsened, ten tied on the primary outcome. Observed paired differences are all zero; the bootstrap interval [0,0] is a consequence of identical outputs on these ten inputs, not a universal equivalence guarantee.

The correction is active: participant 4's reported squared error changes from 0.3944216755847224 to 0.3913027838484498 while its actual MSE remains 1.7505260007935346. Thus arithmetic/reporting changed without changing the selected estimate. Calibration is not established by ten records and this result shows no calibration gain.

## Checks

The cloud collector validates all ten pairs, unchanged source/truth/options, fixed controls, exact algorithm recipe, previous-baseline MSE reproduction, native replay and executable export. Error values were additionally recomputed from the downloaded estimates, gains and known truth. The corrected function passed analytic-weight scaling invariance, mixed weight/frequency expansion equivalence and pure-frequency expansion equivalence; the existing blink-bias and factorial-selection statistical tests passed. No statistical assertion was weakened.

The isolated correction separates frequency counts from the combined weights, using sums of squared analytic weights in variance and group effective sizes. It retains overall variance and the original BiasZ=2 and factorial shrinkage. The prior ten-record data are reused for development; no independent confirmation is claimed.

[Summary](../benchmarks/results/weight_trial_20261104/summary.json), [full paired outcomes](../benchmarks/results/weight_trial_20261104/records.json), [frozen plan](../benchmarks/results/weight_trial_20261104/plan.json), [executed source](https://github.com/abwoo/PipeCompare/tree/7061918c44f23c1a0172503506ea4d280360d447/benchmarks/weight_variant).

The next isolated within-group variance experiment uses the original production baseline, because this trial did not satisfy the improvement gate. PR #51 remains unmerged and production is unchanged.
