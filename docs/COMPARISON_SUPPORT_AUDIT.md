# Artifact flag and comparison-support audit

## Outcome
S005's sparse comparison support has three observed contributors: a high realized injected-blink burden (33/41 trials), three additional continuous-EOG HMM flags, and removal of one remaining unflagged trial by the post-EOG candidate's amplitude-threshold rejection. This is a traced support limitation; no trial-index mismatch or unexplained data loss was found.

No production repair is justified as a proven implementation-bug fix by this audit alone. No accuracy improvement is claimed. Original defaults remain unchanged and PR51 is unmerged.

## Evidence and checks
- [Cloud run 41827407681](https://github.com/abwoo/PipeCompare/actions/runs/41827407681), execution `50990f6ec5d50a4cba56d67bb21a576b7b95f54d`, completed successfully.
- Re-fetched SHA256-checked original S004/S005 sixty-second inputs; reused the frozen simulation generator with observation-only truth recording and identical random draws.
- Original and instrumented detectors produce exactly identical flags. Recreated population/realized targets match saved targets within 1e-12; reference identities and flags match all 60 saved subject/start/arm/stage snapshots.
- Audited 1,270 successful candidate snapshots and 60 selected snapshots. Every selected candidate's absent trial IDs are exactly explained by its rejection records.
- Independently recomputed all 1,270 snapshots' retention, injection and unflagged counts from trial-ID CSVs using Python. Verified detector-source unions, selected rejection identities and identical S005 support transitions at all five ICA starts.
- Artifact `comparison-support-audit`, ID 13320092418; ZIP SHA256 `169c19dd5347a18ecc7966b1eb5e9b2105b2f74bbcc74ebe4def054a07653353`.
- No ICA fitting, new ranking, threshold change or ten-person accuracy rerun was performed.

## S005 support accounting
The frozen generator has nominal blink probability 0.7 for S005. The realized random draw injected blinks into 33/41 trials (80.49%). This describes simulated injections into the public EEG background, not a measured blink rate of the real participant.

| Stage / definition | Condition 11 (21 trials) | Condition 31 (20 trials) |
|---|---:|---:|
| No injected blink | 4 | 4 |
| Unflagged by the original artifact detector, before processing | 2 | 3 |
| Unflagged in the initial selected candidate | 2 | 3 |
| Unflagged after cutoff-search selection | 2 | 3 |
| Unflagged in the final single-EOG recommendation | 1 | 3 |
| Unflagged in the final expanded-EOG recommendation | 1 | 3 |

The last four rows describe candidates selected at different search stages, not sequential deletion on one fixed candidate. Actual deletion attribution uses each candidate's own saved rejection decisions. Both arms use the retained within-variance experimental formula; “single-EOG” is the incumbent of the earlier expansion comparison, not original main.

### Detector contribution
All 33 injected-blink trials are flagged. Three non-injected trials are also flagged:
- ID 3, condition 31;
- IDs 31 and 38, condition 11.

All three are flagged **only by the continuous EOG HMM branch** (including its existing 50 ms dilation and window/baseline overlap logic). They are not flagged by the trial-level robust blink rule, sample-level 5-SD rule, horizontal-eye rules or muscle rule. Overall flags: 36/41. The horizontal detector flags 33/41 and muscle 0/41, but neither adds these three extras.

This is not a validated false-positive count: the resting EEG background, and the EOG proxy constructed from frontal background, lack clean natural-artifact labels. Additional HMM flags can reflect background activity, natural artifacts, or detector sensitivity/temporal extension. The audit identifies the responsible branch but does not establish which explanation is correct.

### Rejection contribution
Before EOG candidate selection, the unflagged identities are:
- condition 11: **23, 33**;
- condition 31: **17, 26, 27**.

Both final EOG recipes remove ID **23** at their `reject_threshold` step (step 10), leaving only ID 33 in condition 11. The selected single-EOG recipe keeps 15/21 and 13/20 total trials; expanded EOG keeps 16/21 and 14/20, but they have the same 1-and-3 unflagged support. All five ICA starts show this same pattern. The rejection-index records exactly explain the missing identities.

ID 23 has no injected blink and no detector flag. That does not prove it should survive amplitude rejection: other high-amplitude background activity may justify rejection. This audit does not establish the per-channel waveform reason for that trial crossing the selected threshold, nor that retaining it improves MSE.

### Connection to the previous bootstrap failure
One available comparator in 21 reference trials is missed with probability (20/21)^21 = 35.89% by the current IID bootstrap. The observed union across both conditions is 820/1999 = 41.02% invalid draws. Restoring ID 23 alone would still leave only two comparison trials in condition 11; missing both has probability (19/21)^21 = 12.22%. Thus merely protecting that one trial does not establish adequate comparison support.

## S004 control
S004 has 16/41 injected blinks and 20/41 artifact flags, leaving 13 and 8 reference comparators in its two conditions. All injected trials are flagged; four extra flags occur on IDs 2, 11, 15, 16. S004 has substantially more reference support than S005. Its previously diagnosed rare invalid draws in other candidates remain distinct from S005's severe selected-candidate support loss.

## Interpretation and next step
The main limitation here is sparse reference information in a short, high-blink-burden development case, compounded by detector flags and adaptive rejection. It is not evidence that EOG regression has no accuracy benefit, and does not justify disabling HMM, deleting failed bootstrap draws, or using injection truth for production selection.

A focused next diagnostic should inspect the HMM state intervals around IDs 3, 31 and 38, plus amplitude/channel evidence for rejected ID 23. Known-truth null/weak-blink/high-prevalence simulations can then test whether a detector change is warranted; real-background non-injected trials must not be declared clean truth. Any proposed repair must undergo paired accuracy testing under frozen criteria. If the data genuinely lack comparison support, the tool should report that insufficiency; numerical changes cannot manufacture missing observations.

No repair was made, no worse version was activated, and no PR was merged.
