# External HMM diagnostic results — 2026-11-05

## Decision

Keep the unchanged production HMM. This diagnostic does not justify disabling it or tuning its threshold to historical S005 trials. It shows detection of prescribed synthetic blink pulses and increased marking during real instructed-blink tasks, while natural-blink false-positive rates remain unresolved. No algorithm modification, P3 accuracy improvement, promotion, merge, or full-pipeline independent confirmation is claimed.

## Completion and missingness

[Cloud run 41957880741](https://github.com/abwoo/PipeCompare/actions/runs/41957880741) succeeded at execution commit `3140252b2e871bd8233ac3d2e285c576f35dbb9b`. MATLAB and Python checks passed. All ten participants have a 60-second blink-task excerpt; seven have a 60-second open-baseline excerpt. Subjects **4, 9 and 10** lack open_base annotations and were not replaced. Total usable input duration is **17 minutes**, not 20. Across one/two-channel configurations, 34 measured outputs and six explicit missing slots account for all 40 planned slots. All 34 external and 40 synthetic HMM fits are valid.

These are new PhysioMotion participants, not the EEGMMIDB participants previously called S004/S005. Numeric IDs must not be treated as the same people.

## Known-truth controls

Ten fixed synthetic backgrounds, each paired with no pulse or ten compact raised-cosine pulses. Pulse peaks are 4/8/16 relative to unit-SD AR noise, not empirically calibrated natural-blink severity levels. A hit means a padded state overlaps ±100 ms of the known pulse center.

| Pulse peak | Hits / events | Mean marked fraction outside exact pulse support |
|---|---:|---:|
| 4 | 98 / 100 | 0.117778% |
| 8 | 100 / 100 | 0.109630% |
| 16 | 100 / 100 | 0.575556% |

Both misses occur in seed 8, amplitude 4. The three amplitude conditions reuse the same backgrounds and event timings; they are not 300 independent trials. These are center-detection rates under the specified simulation, not duration recovery, precision, or general real-EEG accuracy.

In the no-pulse controls, the mean marked fraction is **0.118667%**, or **0.712 seconds across 600 seconds**. Three of ten null traces have nonzero marking; seven have none. Thus false marking is small in total duration under this noise model, but not absent. Outside-support fractions in pulse conditions include the fixed 50 ms dilation and should not be presented as event-level false-positive rates.

## Real-data descriptive check

Primary input is Fp1. The values below are percentages of time marked after the fixed dilation, not blink recall or false-positive rates.

| Available block | Participants | HMM | Robust single-point rule | HMM-only time beyond robust rule |
|---|---:|---:|---:|---:|
| Open baseline | 7 | 13.1286% | 6.5702% | 6.7983% |
| Instructed blink | 10 | 33.4933% | 18.5628% | 16.2293% |

The different denominators must not be treated as a paired comparison. On the **same seven complete participants** (1,2,3,5,6,7,8), HMM marked time rises from **13.1286% to 35.2976%**; each of the seven rises. This is consistent with responsiveness to the instructed-blink condition. More marking does not establish higher accuracy.

Prespecified two-channel sensitivity check: HMM marks 13.8167% of open blocks (n=7) and 35.4750% of blink blocks (n=10). Robust-rule marking is 10.2433% and 19.8872%, respectively. Adding Fp2 changes the masks; it is not demonstrated to improve accuracy.

## Interpretation and limits

The real data have human task/segment labels, not per-blink onset/offset truth. They contain no dedicated EOG channel; Fp1/Fp2 are EEG proxies. Open-baseline EEG may contain spontaneous blinks and other artifacts. Consequently the 13.13% baseline marking cannot be called a 13.13% false-positive rate, and 33.49% blink-task marking cannot be called 33.49% sensitivity. HMM segment counts likewise are not counts of correctly identified blinks.

This check supplies no evidence that historical S005 IDs 3/31/38 are false positives, and does not validate trial-23 rejection policy. It evaluates only the detector, not adaptive pipeline selection, artifact removal or P3 estimation. Synthetic truth is simplified; the true real-data balance between added sensitivity and unwanted marking remains unknown.

## Verification and provenance

- Cloud checker recomputed external fractions and segment counts from all exported masks. Synthetic count/range checks passed; synthetic raw masks were not exported, so independent sample-level reconstruction of synthetic hit decisions is not claimed.
- Downloaded ZIP SHA256 verified: `bbbff94182545a9e237ac87c7f12dfc0c2e51d55b0cd043b0444a0044a5cbc9c`.
- All 34 external result rows and 40 synthetic rows exactly match the cloud log. Local mask recount also passed. Complete-pair arithmetic was independently recomputed.
- Artifact `detector-validation-results`, ID **13383513198**, contains extracted inputs, per-subject source/range digests, annotations, masks, summaries and missingness records.
- The earlier run 41957688403 stopped before detector execution on absent subject-4 baseline annotation. Missingness handling was amended before any detector outcome, and the successful run retains all missing slots explicitly.
- Protocol and code: `benchmarks/detector_validation/`. Production/PR51 code head remains `e8a3c5322b3ebcdbaf56ef4441c8841c5b82b0ef`; PR51 remains unmerged.

## Next decision

No speculative detector repair is warranted from these results alone. The remaining decisive evidence is event-level independent eye-state truth (video or properly adjudicated EOG), including background intervals, before modifying HMM thresholds. If the aim returns to overall tool accuracy, compare candidates on a prespecified held-out P3 benchmark; detector occupancy cannot substitute for that endpoint. No additional experiment is automatically started here.
