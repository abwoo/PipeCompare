# PR51: frozen scoring comparison on ten held-out participants

Status: protocol fixed before running any new participant. PR51 remains open and unmerged; production remains the original implementation.

Execution amendment: run 42225506767 stopped on subject 11 because the new support recorder assigned a populated MATLAB structure into an empty structure without matching fields. Both arms reached their selected point estimates, but support recording and replay were incomplete; subjects 12–20 were skipped. The correction initializes the output array from the first actual record and adds a two-condition weighted/unweighted recorder check. Source selection, simulation, algorithm files and all decision gates are unchanged. The failed run remains preserved; the new execution is collected separately, with no mixing of its partial records. Because partial subject-11 outcomes were exposed by the failure logs, this amendment is explicitly recorded rather than represented as occurring before any outcome exposure.

Question: does the retained within-group variance scoring change improve selected P3 accuracy on new-to-PR51 participants, without unacceptable individual deterioration or unsupported uncertainty estimates?

## Frozen inputs and scope

- Original production: `7151821edf8a5a54c3167d7f18b66c530911b568`; experiment based on PR51 head `e8a3c5322b3ebcdbaf56ef4441c8841c5b82b0ef`.
- Variant: `benchmarks/within_variant/Measure.m`, SHA256 `ebf3252a8ca582effe2af4c5f2ec2f2eee8f512c25d5a430115ec94655b33c11`. Only pooled within-flag-group residual variance changes. Keep factorial smoothing, BiasZ=2 and legacy weighting semantics. Do not add rejected EOG expansion or weight fixes.
- Public EEGMMIDB 1.0.0, DOI `10.13026/C28G6P`, ODC-By-1.0. Subjects **11–20**, run R01, first 60 s, 160 Hz, same 16 channels. This is ten people and ten source minutes, yielding twenty paired algorithm evaluations; conditions, candidates and bootstrap replicates are not additional people.
- `benchmarks/heldout_protocol.json` pins each EDF checksum against the upstream manifest and all simulation/search/statistical parameters. Inject known P3 and blink signals using seed `951000 + subject`; prevalence cycles 0/10/30/50/70% twice. These are new PR51 development participants, not a claim to have opened or audited the PR50 confirmation cohort.
- Same real background, injected truth, fixed controls, search seed and full Standard procedure for both arms. Fresh MATLAB sessions; variant activated only in the disposable Actions checkout. Fresh minimal/common fixed controls provide the usual standard-preprocessing comparisons.
- This tests same-source generalization of a semi-synthetic task. It cannot establish clean-neural truth in natural EEG, real P3-task validity, cross-dataset performance, clinical utility or independently conducted confirmation.

## Minimal run

1. Verify unchanged production paths, original/variant SHA256 and source checksums. Run the native dependency sanity check required by the cloud environment.
2. Run participant 11 as an execution smoke check, then 12–20 at at most three concurrent jobs. Whether participant 11 improves is **not** a scheduling condition. The only stop condition is an execution failure.
3. Per participant: full original Standard search, fixed controls, native/export replay; then the frozen variant with the same inputs in a fresh session. Preserve rankings and detailed trial scores/weights/reference flags.
4. Reconstruct the existing 1999 bootstrap score columns for observation only; record every nonfinite score without dropping draws. Save compact scoring state for later diagnosis without rerunning preprocessing.
5. Recompute MSEs independently in Python and verify the whole ten-person coverage, fixed-control invariance, execution identities and replay agreement.

Do not repeat historical Python/MATLAB engineering suites, subjects 1–10 development searches, earlier ICA/EOG experiments, alpha searches, channel-fault batteries, PhysioMotion/HMM or optical EyeLink validation. These are not missing evidence for this narrowly scoped comparison. New collector checks are limited to preventing incorrect conclusions or incomplete coverage.

## Prespecified decision

Primary endpoint: per-person mean squared gain-corrected population-target error over two conditions, then the mean of the ten paired differences (variant minus original), in µV². Resample whole people together across arms: 10,000 bootstrap draws, seed 510003, the existing order-statistic 95% interval. Also report raw MSE, conditional error, retention, pipeline choice, ICA stability, fixed-control differences and calibration with explicit denominators.

For **retain-heldout-supported**, require all of:

- 10/10 complete, valid pairs; no source/protocol/control/replay mismatch.
- Mean primary MSE falls at least **5%**, paired 95% interval lies entirely below zero, and wins exceed losses (ties within 1e-9 µV²).
- Mean raw MSE does not rise beyond numerical tolerance.
- Each person's primary MSE increase is no larger than `max(0.05 µV², 20% of original MSE)`.
- Variant selected estimates have finite reported errors and assessed finite artifact bias. Per condition require effective unflagged support >=5 if positively weighted flagged trials remain, or >=2 if none remain. Effective size is `(sum weights)^2 / sum(weights^2)`. All feasible-family bootstrap scores must be finite. Baseline failures remain visible but are not used to excuse variant failures.

The 5%, individual deterioration and effective-support thresholds are deliberately chosen practical screening rules for this pilot, **not** established physiological or clinical cutoffs. They are fixed before outcomes and never change the algorithm's choices.

- **Reject** if the complete cohort's mean primary or raw MSE worsens: retain original production; preserve the adverse experiment, do not erase evidence.
- **Inconclusive/no promotion** if point estimates improve but any support, interval or individual guardrail fails. Report accuracy and uncertainty separately.
- **Incomplete** if any pair is missing, failed or inconsistent. Never substitute another person, pool old development subjects, choose a favorable attempt, or report a reduced denominator as ten.
- No outcome automatically changes production or merges PR51. Even a passing result is limited to these ten same-source recordings; stable calibration retains the existing >=30-person floor.

## Checkpoints and recovery

- Save baseline results immediately and complete pair artifacts with subject and run-attempt IDs. Preserve failure logs and partial files.
- For infrastructure failures rerun only failed jobs at the **same commit**. Preserve earlier attempts; collector uses the highest attempt even if incomplete, never the most favorable result. Successful participants are not rerun for a better outcome.
- The collector rejects mixed execution commits/runs/EEGLAB versions, mismatched within-pair attempts, modified source/protocol digests, changed fixed controls or replay error >=1e-6 µV.
- Baseline-only checkpoints are diagnostic and never count as completed pairs. Any source/code/protocol change requires a clearly labeled new run rather than mixing records.
- Scientific rejection/inconclusiveness is a valid completed result, not a red Actions job. Execution/evidence corruption produces an incomplete report and a failing collector.
