# ERP CORE / PR51: frozen-run root cause and next validation (2026-11-12)

Status: **diagnosis only; no merge, no main change, no threshold relaxation, no scientific promotion**.

## Frozen identities

- Original real-P3 execution: [GitHub Actions 42809789847](https://github.com/abwoo/PipeCompare/actions/runs/42809789847), execution commit a5aab337d5f8bbe7f9dbdaa41b4588b80dc91f8b.
- Source cohort: 20 prespecified ERP CORE P3 participants, 300 s each, 16 EEG + real VEOG, 256 Hz.
- Frozen protocol SHA256: bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b.
- Production reference algorithm commit: 7151821edf8a5a54c3167d7f18b66c530911b568.
- Original collector summary: [2026-11-11](../benchmarks/results/erpcore_recollect_20261111/summary.md), 15/20 valid pairs, verdict PARTIAL.
- Independent reproduction surface: \`benchmarks/diagnose_erpcore_logs.py\` downloads original individual-job logs, validates subjects/identity, and writes per-subject reason counts, support, VEOG flags, and ICA seed choices. **It neither refits nor recomputes scientific outcomes.**

## Failure taxonomy from the preserved raw MATLAB job logs

| Subject(s) | Observed reason | Interpretation | Priority |
|---|---|---|---|
| 18, 23 | Each has 18 candidate-step errors in \`icremove(threshold=auto)\`: "Unable to perform assignment because the left and right sides have a different number of elements"; six other candidates fail the 16% filter-amplitude guardrail | Software exception + separate morphology constraint. Not justified to label as universally bad EEG | P0 |
| 17 | 24 candidates fail; all include 31 ms filter-only peak latency shift; some also have 17–18% amplitude change and 6–7% artificial deflection | A frozen shape-validity exclusion; could be real distortion **or** a peak estimator unstable on ambiguous waves | P1 |
| 4 | Standard feasible 19/25 but fixed control readout unavailable; fixed control effective unflagged support 1 target/1 standard; Standard effective unflagged support 4 target/3 standard | Primary support-collapse case, even where final trial counts look large | P1 |
| 8 | Standard effective unflagged 5 target/20 standard; both fixed controls have 4 target/18 standard | Comparator-specific support failure; cannot be solved by counting retained epochs alone | P1 |
| 11, 15, 21 | ICA three-seed candidate IDs respectively 12/25/25, 5/3/3, 25/4/4 | Different selected pipeline IDs; not yet evidence that P3 amplitude or its SE materially changes | P2 |

The 18/23 exception **must be localized to an exact stack frame** before altering production source. The exception is displayed at the \`icremove\` step in the current executor, so do not assert that \`IcThreshold.choose\` itself is already proven at fault. A plausible source is dimension/orientation handling in candidate-specific ICLabel paths, but that remains a hypothesis.

Additional evidence from original diagnostics:
- Subject 4 selected Standard: 18 target and 95 standard epochs kept, but only 4/3 effective unflagged support; the HMM/VEOG diagnostic reports 92 standard flags and median flagged VEOG 241.76 uV. The 3 remaining unflagged standard epochs have median max VEOG 219.75 uV; a binary blink flag is not proof of optical truth.
- Subject 8 fixed arms: 4 target epochs retained, effective unflagged target support 4 (< frozen floor 5). Standard: 16 retained target epochs, effective unflagged support 5.
- Subjects 18 and 23 have 0 target epochs in both fixed arms; the simple published-pipeline approximation also rejects all targets by its fixed ±100 uV criterion. This is independent of the \`icremove\` runtime exception, so a code repair cannot by itself guarantee a complete valid pair.
- Subject 17 fixed comparators return finite P3 estimates. The shape gate is therefore a scientific validity decision, not an inability to load EEG.
- The current Standard readout initializes \`passedConstraints=false\` and never fills it in the \`readout()\` helper. Do not use this field as an actual Standard constraint verdict; instead inspect \`r.ranking.table.status\`, \`result.feasible\`, and support. Any recorder repair must not alter the frozen scientific result.

## Statistical and reproducibility limits

The prespecified full-cohort stability gate is at least 18/20. The current score "13/16 same" is **conditional on the 16 available Standard fits**, not a passing 13/20. Because 4 fits lack valid pairs, **the count gate cannot currently pass regardless of the remaining seeds**. Do not change this 18/20 requirement post hoc; additionally report conditional 13/16 and coverage 16/20 as separate diagnostics. The original completeness gate is 20/20, separate from precision and calibration.

The 15 valid-pair precision results are positive but conditional on inclusion. The 5 unassessable/incomplete cases may not be missing at random; avoid upgrading population-level superiority. The reference-pipeline comparison is **only an approximation**, not bitwise reproduction of ERP CORE's published EEGLAB/ERPLAB workflow (original P9/P10 reference, subject-specific rejection, and ICA differ).

## Stage-wise, falsifiable follow-up protocol (no frozen protocol edits)

**Stage A — forensic replay without EEG**: run the read-only GitHub Actions log auditor on all original 20 jobs. Assert frozen commit/digest/cohort, exact reason counts, support denominators, and ICA identities; save all 20 source-log SHA256 values. Every omitted or mixed job fails closed. Run unit tests. This does not prove a code change works.

**Stage B — root-cause reproduce 18/23**: on an isolated branch preserving the original EEGLAB and MATLAB versions, add *diagnostic-only* extended stack capture around the failed candidate step and run only subjects 18 and 23. Determine stack line, affected array sizes and channel/IC counts. Write synthetic unit regression tests reproducing the exact failure. Only then propose a minimal code repair, compare fixed vs original on these two subjects, and ensure all original candidates on unaffected subjects replay unchanged. Never lower a signal guardrail during this experiment.

**Stage C — distinguish filter failure from noisy peak (17)**: export the frozen injected-known-signal waveform and filter-only difference curves (300–600 ms), and bootstrap the peak location with **trial** resamples. Measure width and multimodality of the peak-latency distribution, mean-window distortion, integrated area and waveform correlation. A new shape-validity metric may be proposed as a *prospective alternative*, but the original 10 ms limit and comparison results stay intact. A noisy peak is not automatic permission to pass.

**Stage D — data-quality and support audit (4, 8, 18, 23)**: fixed stage-by-stage epoch IDs, rejection reason, thresholds, VEOG peak/derivative, blink HMM interval, ICLabel decisions, target/standard support and effective-weight denominator. Distinguish EOG signals, saturation, baseline drift, true blink contamination, overly aggressive flags, and all-channel amplitude rejection. Evaluate candidate fixes with independent external EOG/optical blink references where available; avoid overwriting trial flags or force-retaining rejected trials.

**Stage E — ICA uncertainty (11, 15, 21)**: for each pre-specified ICA start, save selected pipeline ID, P3 difference amplitude, SME/bootstrap SE, trial set, ICs removed, ICLabel scores, and scalp-activation matching across runs. Report both label agreement and measurement-equivalence (prespecified µV margins chosen without inspecting hold-out outcomes). Add ICA split-half/component matching or a RELICA-like stability diagnostic as research only. The current identity-based \`same\` remains unaltered.

**Stage F — prospective confirmation**: after repairs pass existing regression tests, pre-register a new protocol/hash, run both original and repaired pipelines on the same held-out ERP CORE subjects and identical trial definitions. Retain 20/20 coverage accounting, non-finite readout failures, negative results, and individual regressions. A longer full-task window may be evaluated separately as a new design; do not silently replace the frozen 300 s experiment or its original conclusion.

## Literature basis (not proof for this repository)

1. Kappenman ES et al., [ERP CORE](https://pubmed.ncbi.nlm.nih.gov/33099010/), NeuroImage 2021. Real P3 protocols have trial exclusion and reference choices that matter; reproducibility requires reporting them.
2. Luck SJ et al., [Standardized measurement error](https://pubmed.ncbi.nlm.nih.gov/33782996/), Psychophysiology 2021. Trial count alone is not precision; report the SE for the actual P3 score.
3. Pontifex MB et al., [Variability of ICA decomposition may impact EEG signals](https://pubmed.ncbi.nlm.nih.gov/28026876/), Psychophysiology 2017. Repeated ICA can perturb reconstructed P3 amplitude.
4. Clayson PE et al., [How does noise affect amplitude and latency measurement of ERPs?](https://doi.org/10.1111/psyp.12001), Psychophysiology 2013. Mean amplitude is often less noise-sensitive than isolated peaks.
5. Groppe DM et al., [Identifying reliable independent components via split-half comparisons](https://pmc.ncbi.nlm.nih.gov/articles/PMC3062525/). A split-half design probes training-data uncertainty rather than only random initialization.
6. Pion-Tonachini L et al., [ICLabel](https://pmc.ncbi.nlm.nih.gov/articles/PMC6592775/), NeuroImage 2019. Artifact prevalence and separation quality limit automatic IC classification reliability.

The CMS, Medicare, FDA public-health sources do **not** supply individual EEG signals or ERP CORE trial-level ground truth; they cannot resolve this methodological question.
