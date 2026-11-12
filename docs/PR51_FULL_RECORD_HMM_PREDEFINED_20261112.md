# PR51 stage 2B — full-record source HMM support test, frozen development scope

Prepared before recalculating the full-record EEG HMM masks. **This is an exploratory use of previously observed participants 4 and 8**, not prospective independent inference. The original 300-second, 20-participant ERP CORE P3 protocol and 15/20 PARTIAL report remain immutable. No new held-out participant is used or requested.

## Why this experiment is justified

The event-only, SHA256-verified read in [run 43179038504](https://github.com/abwoo/PipeCompare/actions/runs/43179038504) found:
- Subject **4**: original 300 s has 18 valid targets / 95 valid standards, complete source **544 s** has 30 / 159.
- Subject **8**: original 300 s has 22 / 95, complete source **460 s** has 30 / 126.

These are correct/RT-valid stimuli, NOT clean EEG trials. Subject 4 had only **4 / 3** HMM-unflagged source trials at 300 s, subject 8 **5 / 20**. More recorded events may add usable support without weakening constraints, but only source HMM recomputation can test it.

## Predefined paired source-replay experiment

1. Fetch the *identical original source files* and verify existing pinned `.set`, `.fdt`, `events.tsv` SHA256; no new source cohort.
2. For each participant, prepare two data inputs with the same original 16 EEG channels, recorded VEOG, event relabeling, 256 Hz downsample, correct-response and RT 200–1000 ms rules: a **300-second baseline** and a **complete-record** excerpt (subject 4: 544 s; subject 8: 460 s). An opt-in `--duration-seconds` switch is isolated to this development branch; the original fetch default stays 300 s. All source and duration choices are documented and fixed before HMM results.
3. Build the same EEGLAB dataset as the original pilot. Run `pipecompare.eval.Measure.reference(EEG,contract)` with exactly the existing target/standard codes 11/31, epoch -200..800 ms, baseline -200..0 ms, Pz P300 window 300..600 ms; **no ICA, no ranking, no altered HMM thresholds and no EOG regression**.
4. First reproduce the original 300-second source HMM-unflagged counts **exactly**: subject 4 = 4/3, subject 8 = 5/20. Failure is an experiment/semantics discrepancy; do not assert an improvement from the longer window in that case.
5. For the complete source report source target and standard counts, HMM flags, *actual unflagged IDs*, unflagged N, and whether each condition reaches the unchanged floor of 5. Preserve all outcomes, including failure. HMM source-unflagged N is only an upper bound on downstream effective support; it is **not** candidate-level post-ICA effective N or an estimator accuracy result.

## Interpretation

If full-record source HMM support remains below five in any condition, the original source comparison shortage remains unsolved for that participant. If both conditions meet five, that is only evidence that there are enough HMM-unflagged source events to justify **another separately pinned full-pipeline exploratory test**, including weighting, trial rejection, P300 shape and calibration. No automatic promotion, no modification of the original 300-second frozen cohort, no past-run pooling, and no 14-person new frozen study.

Original production defaults and PR #51 remain unchanged and unmerged.
