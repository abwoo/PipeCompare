# S005 waveform and rejection audit

Completed 2026-11-05. Diagnostic only; production defaults and PR51 code remain unchanged. PR51 remains unmerged.

## Decision

Both requested diagnostics are complete. The extra HMM flags are not caused solely by temporal dilation; trial 23 exceeds each selected recipe's frozen final threshold even when EOG regression is skipped with other decisions fixed. Neither finding demonstrates an implementation defect warranting a production repair. Do not disable HMM, force-retain trial 23, or tune a threshold to these observed outcomes.

This is a seed-zero reconstruction of existing development results, not a new adaptive search, new participant sample, independent scientific confirmation, or an accuracy-improvement experiment. “Independent verification” below means recomputation of exported sample evidence in Python.

## 1. HMM waveforms: non-injected trial IDs 3, 31, 38

Times are relative to event onset; baseline is −200 to 0 ms and the P3 measurement window is 300 to 600 ms.

| Trial ID | Raw HMM active interval (ms) | Overlapping measurement region | Raw / padded overlapping samples | Maximum feature robust z |
|---|---:|---|---:|---:|
| 3 | −256.25 to −75 | Baseline | 21 / 29 | 6.165412 |
| 31 | 456.25 to 562.5 | P3 | 18 / 32 | 6.990485 |
| 38 | 456.25 to 600 | P3 | 24 / 32 | 7.653848 |

All three overlap before the 50 ms temporal dilation. The robust single-sample detector does not flag them. Removing all injected eye signal from the continuous input still produces HMM overlap in all three cases. Original and no-injected-eye EOG samples are exactly equal inside their baseline/P3 windows (maximum difference zero). These findings support background activity as the source of the overlap, rather than injected-blink spillover or dilation alone.

**Ground-truth limitation:** EOG1 in this benchmark is a synthetic proxy, 0.5 × frontal Fp1 background plus injected blink, not an independently recorded EOG channel. Non-injected does not mean artifact-free. We cannot distinguish natural eye activity, other frontal background fluctuations, and false-positive HMM labels from these data alone. The smoothed HMM feature is evaluated after downsampling by six from 160 Hz; its 3.5 peak criterion is a segment qualification criterion, not a pointwise blink classifier. A high feature z is not proof of a real blink.

## 2. Trial 23: exact amplitude rejection

The following are the two previously selected within-variance experimental recipes, not a comparison of original main against the variant.

| Selected recipe | Crossing channel | Frozen final threshold (µV) | Peak before EOG (µV) | Peak after EOG (µV) | Peak time (ms) | Excess after EOG (µV) |
|---|---|---:|---:|---:|---:|---:|
| Single-EOG, 0.1–20 Hz | F7 | 21.8 | 24.022066 | 23.399380 | 556.25 | 1.599380 |
| Expanded-EOG, 0.1–40 Hz | F8 | 26.0 | 26.090050 | 26.229204 | 56.25 | 0.229204 |

Each recipe has exactly one crossing channel. EOG regression reduces the F7 peak by 0.622686 µV and increases the F8 peak by 0.139154 µV. Both already cross their respective frozen final thresholds before EOG. The counterfactual that skips EOG while keeping other decisions fixed also rejects trial 23. Neither recipe removes ICA components; post-EOG and pre-rejection samples therefore agree.

The F7 peak falls inside the P3 measurement window. The F8 peak is outside both baseline and P3, at 56.25 ms. The rejection rule tests the whole −200 to 800 ms epoch across EEG channels (excluding EOG1), whereas the P3 objective is measured at Pz over 300–600 ms. Thus the expanded recipe exposes a design tradeoff between whole-epoch quality control and preserving trials for a narrower measurement. It is not by itself a coding error. A small 0.229204 µV excess is not justification for retrospective threshold adjustment.

**Causality limit:** these comparisons freeze the final threshold and other saved decisions. They do not rerun the entire adaptive procedure without EOG, and do not exclude EOG's influence on threshold selection or upstream search decisions.

## Validation completed

- Frozen inputs and saved targets reproduced.
- Both selected recipes exactly reproduce saved rejected trial identities and per-trial P3 estimates; maximum replay error = 0.
- MATLAB diagnostic marker: `WAVEFORM_AUDIT_CHECKS_PASSED hmmCases=3 selectedRecipes=2 seed=0`.
- Independent Python recomputation checks original/control EOG identity, HMM overlap counts, all channel/stage maxima, peak times and threshold crossings against exported CSV samples.
- Python marker: `WAVEFORM_SAMPLE_VERIFICATION_PASSED`; three HMM cases and two recipes passed.
- Both resulting plots were visually inspected. No algorithm modification or ten-person rerun was performed.
- Prior support audit found the same trial-23 rejection across five ICA starts; detailed waveforms here reconstruct seed zero only.

## Provenance and deliverables

| Item | Record |
|---|---|
| Production commit | `7151821edf8a5a54c3167d7f18b66c530911b568` |
| PR51 code head at completion | `e8a3c5322b3ebcdbaf56ef4441c8841c5b82b0ef` |
| MATLAB execution commit | `3628ad32b1d8d295712b481ece42e24d001eb7a4` |
| MATLAB run | [41830489631](https://github.com/abwoo/PipeCompare/actions/runs/41830489631) |
| MATLAB artifact | `waveform-audit-results`, ID 13323257137 |
| MATLAB ZIP SHA256 | `70ba43d9c206370d98fedcd469b4d1dcc6bc8cbb1bd3b45a8fb73b7b79f9ada0` |
| Python execution commit | `0378b78e0ec6f172b2ca176d9082bb1868864089` |
| Python run | [41951157537](https://github.com/abwoo/PipeCompare/actions/runs/41951157537) |
| Python artifact | `verified-waveform-plots`, ID 13380742991 |
| Python ZIP SHA256 | `fabffd27fdcf5438c440e3320e90ed69b23c4cce2498a407ec85c426d994ae28` |

Small audit JSONs and verified SVG plots are committed alongside this report under `benchmarks/results/waveform_audit_20261105/`. Complete continuous and stage-specific sample CSVs remain in the MATLAB run artifact. Source diagnostics are under `benchmarks/waveform_audit/` on the isolated `codex/pr51-waveform-audit-20261104` branch. Documentation is committed after execution with [skip ci]; no latest-documentation-head full-suite pass is claimed.

## What remains unresolved

No natural-blink truth is available for classifying the three additional HMM labels. Comparison support in S005 remains sparse; restoring only trial 23 does not establish sufficient support. If further work is authorized, first validate detector behavior with independent EOG/blink labels or a prespecified known-truth control. Any alternative rejection window needs a scientific rationale fixed before evaluating held-out data, including checks for retained artifact and P3 bias. Do not use these same trial identities to design a passing rule.
