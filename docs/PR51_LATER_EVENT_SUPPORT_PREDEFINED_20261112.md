# PR51 Stage 2 follow-up — later task-event availability, subjects 4 and 8

Prepared before downloading any new later-event windows. This is **development-only** and may not revise the original ERP CORE P3 frozen 300-second excerpt, 20-subject PARTIAL 15/20 verdict, original HMM eligibility, ICLabel thresholds, or paired precision metrics.

## Context and purpose

Subject 4 had just 4 unflagged targets and 3 unflagged standards under the frozen blink reference; subject 8 had 5 and 20 and a fixed comparator lacked paired support. No weighted ranking or eye regression can restore more clean **source trials** than exist under the original evidence definition. This stage asks a narrower prerequisite: do original, pinned ERP CORE recordings contain *additional correct/RT-valid task events* after 300 seconds?

## Frozen before checking extended windows

Fetch **only** original, SHA256-verified `.set` metadata and `events.tsv` for the same participants 4 and 8; do not download or process full EEG `.fdt` signal. Derive exact original recording duration from `EEG.pnts/EEG.srate`, check `.fdt` metadata size, and leave data selection unchanged.

Count event-filter-eligible targets and standards for 300, 360, 420, 480 seconds where source duration supports them, plus the maximum full integer-second window. Apply the original correct-response and 200–1000 ms RT rules and the exact epoch boundary rule; explicitly retain the frozen 300 s row as reference. Record additional events versus 300 s, source SHA identities and unavailable longer windows. Do not choose a single new duration based on outcome or call the extra events HMM-clean, epoch-retained or scientifically usable.

Even a higher raw correct-trial count only warrants **consideration** of a new prospectively pinned longer-window scientific feasibility test. That separate test would need full new EEG source and unchanged HMM/VEOG support assessments, measured-waveform preservation, its own versioned protocol, and a cohort/selection policy. It cannot be pooled with the original 20 subject results or used to change the original scientific eligibility gates.

This task must never launch the 14 completely unseen ERP CORE subjects, activate an algorithm, or merge PR #51. Code and run are isolated in `codex/pr51-support-extension-diagnostic-20261112`.
