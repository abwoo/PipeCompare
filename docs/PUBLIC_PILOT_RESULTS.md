# Results: ten short public EEG recordings

Executed on GitHub Actions, 2026-11-03, MATLAB R2026a Update 5 + EEGLAB 2026.0.0.

[Completed cloud run](https://github.com/abwoo/PipeCompare/actions/runs/41604029052) | [Frozen protocol](../benchmarks/public_pilot_protocol.json) | [Machine-readable evidence](../benchmarks/results/public_pilot_20261103/evidence.json)

## What this pilot supports

The unchanged main implementation completed all 10 participant excerpts and produced 10 assessable Standard recommendations. All 10 unaltered resting-alpha searches also produced a recommendation. Native replay and the generated executable script reproduced every selected estimate exactly in this run.

The known-target ERP experiment supports lower selected error than the minimal fixed pipeline on this pilot. An advantage over the common fixed ICA pipeline is not established by the paired interval. Calibration remains unestablished, and one participant changes recommendation under alternative ICA starts. These are separate findings; successful execution does not establish stable scientific performance.

The primary accuracy experiment uses a declared semisynthetic P3/blink injection into real resting EEG background. It is not a real P3 task, certified clean neural data, clinical validation, or a group-prior generalization study.

## Prespecified paired accuracy

MSE is in uV^2 against the known population-target P3 after the declared gain correction. Each participant has equal weight after its two condition errors are squared and averaged. Difference = Standard minus control; negative favors Standard. Intervals resample 10 participant pairs, not candidates or conditions.

| Comparison | Pairs | Standard mean MSE | Control mean MSE | Mean difference | 95% participant-bootstrap interval | Interpretation |
|---|---:|---:|---:|---:|---|---|
| Minimal fixed | 10 | 0.483153 | 2.352045 | -1.868891 | [-3.256516, -0.705841] | pilot-improvement |
| Common fixed ICA | 10 | 0.483153 | 1.172762 | -0.689609 | [-1.944103, 0.169792] | not-established |

Mean raw, uncorrected ERP MSE: Standard 0.523975; minimal 2.353385; common ICA 1.448850. Raw and gain-corrected outcomes are not interchangeable.

Both fixed controls produced a readout for all 10 inputs, but each passed the tool's separate signal/retention constraints on 7/10. Their prespecified primary comparisons retain all ten outputs; failures were not removed after seeing outcomes. This is a comparison with these fixed recipes, not every possible well-tuned conventional workflow.

## Error-report calibration

Actual/reported RMS ratio: **1.165410**, 95% interval **[0.670970, 1.596618]** over 10 selected readouts. Actual RMS = 0.695092 uV; reported RMS = 0.596436 uV.

The interval overlaps and extends beyond the prespecified [0.8,1.2] agreement band: **not-established**. Relative half-width is 39.71%, above the 10% precision target. The existing floor of 30 assessable independent records is also unmet. A point ratio of 1.17 is not a universal correction multiplier or proof of calibrated confidence intervals.

## Stability, faults and execution

- Full Standard fitting, adaptive cutoffs/EOG checks and two extra ICA starts were retained on every recording.
- ICA recommendation unchanged across the extra starts: **9/10**. Participant **6** changes recommendation; the original recommendation was not replaced post hoc.
- Native replay and executable-export maximum estimate difference: **0 uV on all 10**.
- Mean end-to-end data computation per participant: **103.3 seconds**, excluding MATLAB installation and native-backend setup.
- Native ICLabel and original MATLAB inference agreed on all checked winning classes; maximum probability difference was **1.848e-6**.

Strong injected Pz channel faults use 100 uV noise. Outside-injection marking is measured over all un-injected channel-samples; natural background is not known clean, so this is not a certified false-positive rate.

| Injected duration | Assessed inputs | Mean injected-sample recall | Mean outside-injection marking |
|---|---:|---:|---:|
| 1 s | 10/10 | 0.0% | 0.000% |
| 5 s | 10/10 | 0.0% | 0.000% |
| 10 s | 10/10 | 100.0% | 0.316% |
| 20 s | 10/10 | 100.0% | 0.000% |

The tested long-window detector misses the 1/5-second injections. Detection of these strong 10/20-second injections does not establish sensitivity for weak, adjacent-channel, or natural faults.

## Every planned participant

| Participant | Blink prevalence | Standard MSE | Minimal MSE | Common ICA MSE | ICA recommendation stable |
|---|---:|---:|---:|---:|---|
| 1 | 0.0 | 1.393875 | 1.171418 | 1.152529 | true |
| 2 | 0.1 | 0.017796 | 0.898050 | 0.175181 | true |
| 3 | 0.3 | 0.240221 | 0.768514 | 0.768514 | true |
| 4 | 0.5 | 1.750526 | 2.099136 | 0.971575 | true |
| 5 | 0.7 | 0.256920 | 5.638506 | 6.085790 | true |
| 6 | 0.0 | 0.416618 | 1.107339 | 0.171195 | false |
| 7 | 0.1 | 0.012027 | 0.831888 | 0.013549 | true |
| 8 | 0.3 | 0.083436 | 1.545504 | 1.545504 | true |
| 9 | 0.5 | 0.043610 | 6.298838 | 0.642550 | true |
| 10 | 0.7 | 0.616503 | 3.161254 | 0.201231 | true |

All ten planned sources were present, checksum-verified, execution-complete and assessable. There were no production execution errors and no dropped participants. All source, protocol, algorithm and execution identities agree. Artifact JSON was compared with its cloud log for every participant; point MSE and RMS summaries were independently rechecked.

## Reproducibility and limits

The production algorithm is 7151821edf8a5a54c3167d7f18b66c530911b568. The executed benchmark/checker revision is cca3b44e94ee5c4779f3f4d4f6ed4c3c6ec41b1e; later result and documentation commits do not imply a fresh compatibility run. No PR50 production changes were imported, no outcome-driven method tuning was performed, and its confirmation cohort remains unused.

Committed evidence: [per-record outcomes](../benchmarks/results/public_pilot_20261103/records.json), [compact table](../benchmarks/results/public_pilot_20261103/records.csv), [summary](../benchmarks/results/public_pilot_20261103/summary.json), [fault measurements](../benchmarks/results/public_pilot_20261103/faults.csv), and [provenance](../benchmarks/results/public_pilot_20261103/evidence.json). The run artifacts additionally retain candidate tables, calibration rows, executable exports and native logs for 90 days. Raw EEG and fitted ICA matrices are not committed.

This supports small-sample operational usability for the declared 16-channel, 160-Hz excerpts and these two measurement presets. It does not certify all channel counts, ERP components, acquisition systems, preprocessing starting points, operating systems or MATLAB releases. It does not justify promoting an experimental algorithm, merging PR50, or declaring a general accuracy/calibration guarantee.

## Follow-up scoring experiment

A paired replay of all ten records tested removing factorial smoothing from ranking. Mean MSE increased from 0.483153 to 0.487013 uV² (0 improved, 1 worsened, 9 unchanged), and the existing known-truth simulation regressed. The variant was rejected; the production algorithm remains unchanged. See [the complete trial results](RANKING_TRIAL_RESULTS.md).
