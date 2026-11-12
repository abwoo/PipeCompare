# PR #51 ICA scientific readout — six-start development-only protocol

**Execution branch:** `codex/pr51-ica-science-recorder-20261112`. **Raw run to compare against:** [43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528) at `5c9388562974c8305da95cb59e7249982788d589`. The observer run was launched at [43184910782](https://github.com/abwoo/PipeCompare/actions/runs/43184910782). This is a development-repeat experiment, not a new held-out sample or the original prospective validation. Original 20-person ERP CORE P3 result remains **15/20, PARTIAL**, and PR #51 remains unmerged.

## What is unchanged and what is new

The source protocol manifest, original `Measure.m` and `Rank.m`, reference HMM labels, five-trial effective support floor, 300–600 ms Pz target-minus-standard estimand, actual filter/ICA/IC removal/rejection decisions and all seeds 0–5 are unchanged. Source recordings are downloaded in cloud by the existing checksum-pinned fetcher. All changes are isolated research recorder code in `+pipecompare/+research/IcaScience.m`, opt-in through `pop_pipecompare(...,'sciencecapture',true)`, plus read-only candidate snapshots after numerical evaluation and a source-identity checker.

Observation does not change a candidate's original status or ranking. The observer stores only derived science metadata and ERP arrays. If it fails, the original candidate still completes under the old scientific rules, but the observer separately labels its science record as failed, and GitHub Actions validation cannot claim complete scientific evidence.

## Readouts and provenance

- Each actual fitted ICA solution: component count, channel labels and the mixing/topographic matrix `pinv(icaweights*icasphere)`; this has arbitrary component sign and scale and requires matching before comparing.
- ICLabel probability matrix and class names, plus exactly which components the original `icremove` step discarded. These are recorded before the final EEG candidate object is forgotten.
- The actual evaluated two-condition Pz ERP on its complete timebase (not merely 300–600 ms). Full-wave condition averages are based on the *same retained urevent trial identities and positive weights* used by original `Measure.candidate`. The integral over P300 must reproduce the original uncorrected estimate; division by the original per-seed gain must reproduce the previously reported gain-corrected amplitude.
- Full candidate/reference trial IDs, retained IDs, original source blink flags, per-trial scores and weights, effective N. The physical–effective trial support inequality remains enforced.
- Read-only independent P300 difference bootstrap SE, 1000 with-replacement trial draws in each condition, fixed RNG `520004 + subject + 1000*seed + candidateId`. This is *IID trial resampling*, not a new noise-free truth metric, not a bias measure and not independent participant validation.
- Within a **fixed candidate recipe and participant**, match components across starts against seed 0 using absolute centered scalp-map correlations and Hungarian assignment (sign-invariant). Report matched component count, median/minimum matched correlation, unmatched component count and matched removed-component Jaccard. Do not compare raw component indices across ICA runs.
- Full gain-corrected Pz target-minus-standard difference-wave RMS change between each alternate start and seed 0. This is descriptive morphology stability, not accuracy against an unknown neural ground truth.

## Required consistency and stop gates

Compare observer re-run with the immutable prior cloud artifact for the **same subject**: algorithm/protocol/Measure identities, completion, per-seed recommended candidate IDs, original scores, scientific Rank feasibility/strata, gain-corrected P300, two-condition effective/physical support and seed alignment. Any difference outside recorded numerical tolerance must be reported as an observer parity failure; no original evaluation results may be overwritten or silently substituted. A science record is complete only when channel/mixing matrices, labels, ICA removal evidence, Pz time coordinates and actual kept-trial identities all align. Archive all positive/negative raw outputs even when the independent validator rejects them.

The evaluated subjects are previously observed ERP CORE development subjects **11, 15 and 21** only. Predeclared pipeline IDs are 11→12/25, 15→5/2/6, 21→25/4 (unchanged Standard, Minimax and stability-first choices). Selecting a new candidate based on the resulting topography/waveform/SME data is forbidden in this experiment. The untouched 14-subject group stays sealed until a distinct selection criterion and evaluation gate have been set in advance and scientifically qualified on development data.

No activation of a production ranking rule, PR #51 merge, successful 14-subject experiment or independent clinical conclusion follows from a successful observer run alone.
