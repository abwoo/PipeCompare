# Held-out scoring comparison: result

Prospective held-out comparison of the retained within-group variance scoring
variant against unchanged production, on ten EEGMMIDB participants (S011–S020)
who were not used to develop the variant. Same 60-second excerpts and
known-target injection protocol as the original pilot. Same-source
generalization only; not a real P3 task, clinical, or independently executed
validation.

[Completed cloud run](https://github.com/abwoo/PipeCompare/actions/runs/42227299459),
execution `0f59671`. An earlier run (42225506767) stopped on subject 11 because
the newly added support recorder raised MATLAB `heterogeneousStrucAssignment`;
the recorder was fixed and the run was restarted from scratch. No partial
outcomes from the failed run are mixed into this result.

## Decision: inconclusive — do not promote

All ten pairs are exact ties. The variant changes nothing on these ten
participants: identical selected pipeline, identical estimates, identical MSE,
identical reported error, identical gain, in all ten cases.

| Metric | Original (baseline) | Variant |
|---|---:|---:|
| Mean gain-corrected population-target MSE | 1.182348 | 1.182348 |
| Mean raw MSE | 1.224095 | 1.224095 |
| Mean difference (variant − baseline) | 0.0 | 95% CI [0, 0] |
| Wins / losses / ties | 0 / 0 / 10 | |
| Gates passed | mean reduction ≥5%: no; paired CI entirely below zero: no; raw mean not worse: yes; wins exceed losses: no; no material individual regression: yes; variant uncertainty supported: no | |

The 29.8% mean improvement seen on the ten development participants
(0.483153 → 0.339095) does not replicate on these ten held-out participants.
That development gain was concentrated in a single participant (S004) and
does not generalize here. The variant is therefore **not promoted to
production default**. It remains a retained exploratory candidate only, and
this result is evidence against treating it as a general improvement.

## Every held-out participant

| Participant | Baseline MSE | Variant MSE | Difference |
|---|---:|---:|---:|
| 11 | 1.025268 | 1.025268 | 0 |
| 12 | 1.662072 | 1.662072 | 0 |
| 13 | 0.585963 | 0.585963 | 0 |
| 14 | 0.866039 | 0.866039 | 0 |
| 15 | 0.348138 | 0.348138 | 0 |
| 16 | 0.455505 | 0.455505 | 0 |
| 17 | 1.369654 | 1.369654 | 0 |
| 18 | 2.942403 | 2.942403 | 0 |
| 19 | 2.383795 | 2.383795 | 0 |
| 20 | 0.184644 | 0.184644 | 0 |

## Fixed comparators on the held-out ten

On these ten fresh participants the original production algorithm does not
establish an advantage over either fixed control:

| Comparison | Standard MSE | Control MSE | Mean difference | 95% CI | Verdict |
|---|---:|---:|---:|---|---|
| Minimal fixed | 1.182348 | 1.194369 | −0.012021 | [−1.118557, 0.937190] | not-established |
| Common fixed ICA | 1.182348 | 1.901495 | −0.719147 | [−2.078950, 0.475519] | not-established |

Note the contrast with the development ten, where Standard beat the minimal
control clearly (0.483 vs 2.352). The held-out ten are harder overall
(Standard mean MSE 1.182 vs 0.483 on the development ten), and the minimal
control is nearly tied here. Ten same-source records per cohort cannot settle
which method is generally better.

## Calibration on the held-out ten

Actual/reported RMS ratio: **1.537936**, 95% CI **[1.223575, 1.989480]**.
Actual RMS 1.087358 µV; reported RMS 0.707024 µV. The interval lies entirely
above the prespecified [0.8, 1.2] agreement band: **adverse**, not merely
unestablished. On these fresh participants the tool's reported error
understates the actual error by about 54%. This is worse than the
development ten (ratio 1.165410, CI [0.670970, 1.596618]). The variant's
calibration numbers are identical to baseline because every selected output
is identical. Ten records do not establish stable calibration; the direction
here is that reported uncertainty is too small.

## Verification and limits

- All ten planned pairs valid; 11/11 cloud jobs (10 participants + collector)
  completed successfully.
- Every paired MSE, pipeline string, estimate, reported error and gain was
  independently recomputed from the downloaded records; all match the
  reported values exactly.
- Production algorithm: `7151821edf8a5a54c3167d7f18b66c530911b568`. Executed
  benchmark/checker: `0f5967166851883eed0ac975e404ff85e3125752`. Frozen
  protocol SHA256 `c41e8456f598b5bc766f102bd3f4bbc0d23226aeebae53d839cb0f536b67d329`.
- Retention and support: every pair retains trial-level support rows in
  `records.json`; no pair was dropped, imputed, or treated as zero error.
- Limits: same-source semi-synthetic validation (real resting EEG background
  plus injected known P3/blinks), 60-second excerpts, 16 channels, 160 Hz.
  It does not certify real task ERP, other channel counts, acquisition
  systems, or clinical use. A zero difference on ten participants rules out
  a large general effect of this variant on this protocol, but it does not
  prove the two scorings are identical on all inputs.

## What this changes

- The within-group variance candidate stays exploratory and is **not**
  promoted. No further promotion case can be built from the development ten
  alone.
- The open question is no longer "does the variant help" but whether the
  original production algorithm itself is reliably better than simple fixed
  controls on fresh data, and why its reported error understates actual
  error by roughly half on fresh participants. Those need a larger fresh
  cohort and a calibration-focused design, not more tuning on the same
  twenty participants.
- PR #51 remains open and unmerged. Production is unchanged.

Committed evidence: [summary](../benchmarks/results/heldout_scoring_20261106/summary.json),
[paired records](../benchmarks/results/heldout_scoring_20261106/records.json),
[collector summary](../benchmarks/results/heldout_scoring_20261106/summary.md).
Frozen protocol: [HELDOUT_SCORING_PROTOCOL.md](HELDOUT_SCORING_PROTOCOL.md).
