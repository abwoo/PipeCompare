# Optical blink-reference diagnostic results

Completed 2026-11-06. [Cloud run 42146097502](https://github.com/abwoo/PipeCompare/actions/runs/42146097502) succeeded at commit `2d68c97db346a319d26f5e299684333082ca79d4`.

## Decision

Keep unchanged production HMM. The detector matches all 21 binocular EyeLink blink-reference events in this short example. Two additional HMM segments overlap optical saccade intervals, consistent with the implementation's stated scope of blinks AND eye movements. This is not proof that all extra marks are correct: overlap alone does not establish causation, and optical reference validity for background is unavailable. No threshold adjustment, production repair, P3 improvement claim or merge.

## Sample and alignment

One EYE-EEG natural-reading demonstration participant, not ten people. Analyzed 86.552734375 seconds between the first and last shared trigger (EEG time 12.447265625–98.998046875 s). Endpoint-only clock alignment used two of 20 shared triggers; the 18 interior triggers had maximum residual **1.939135 ms**, below the frozen 10 ms gate. No timing parameters were fitted using detected blinks or EEG amplitudes. There are 21 overlapping left/right optical blink events.

EEG was already filtered and average-referenced by the source authors. Primary inputs are fixed IO1−Fp1 / IO2−Fp2 bipolar proxies; Fp1-only is a prespecified sensitivity check. Default HMM parameters and 50 ms dilation are unchanged. One-to-one event matching permits ±100 ms around reference intervals. Optical events are independent of EEG amplitudes but are automated EyeLink classifications, not adjudicated video truth.

## Results

| Detector | Predicted segments | Matched optical blinks | Missed optical blinks | Unmatched segments | Marked time |
|---|---:|---:|---:|---:|---:|
| Raw HMM, before dilation | 23 | 21 | 0 | 2 | 16.4730% |
| HMM with 50 ms dilation | 23 | 21 | 0 | 2 | 19.1718% |
| Robust 5-MAD point rule, padded | 44 | 21 | 0 | 23 | 11.9192% |
| HMM + robust union | 32 | 21 | 0 | 11 | 20.5032% |
| Fp1-only HMM, padded | 25 | 21 | 0 | 4 | 10.1094% |

Reference-event recall is 21/21 for all tested configurations in this one record. HMM dilation is not required to obtain these 21 matches. More/fewer segments do not alone establish better performance: brief threshold crossings can fragment a waveform, while HMM can mark longer continuous intervals. No predicted segment exceeds two seconds.

The two unmatched padded HMM intervals are **47.630859375–47.88671875 s** and **75.755859375–76.36328125 s** in the original EEG time axis. Both overlap ESACC intervals. All 11 unmatched union intervals likewise overlap saccades. This exploratory breakdown was added separately from the original executed protocol and does not convert these intervals into certified correct detections.

## Critical negative-reference limitation and reporting correction

Inspection of the source ASC found **101,635 sample rows, with both pupil-area fields zero in every row**. Under the frozen rule requiring positive binocular pupil areas, zero non-blink background seconds qualify as valid optical observation. The 2.408203125 seconds of eligible time consist only of the explicitly added binocular blink-reference intervals.

Thus all unmatched segments remain **indeterminate**, not certified false positives. Specificity, true precision and physiological false-positive rate are not estimable. In particular, the initial cloud summary's zero unmatched-per-eligible-minute value is misleading when background eligibility is absent. The follow-up evaluator changes that rate to **null / unavailable** when eligible background seconds are zero. It also exports eligibleBackgroundSeconds explicitly. No detector output, matching identity, threshold, reference interval or synchronization parameter was changed.

The positive reference matches remain reportable as agreement with recorded optical blink events. They do not establish general real-blink sensitivity or complete annotation coverage. The historical S005 extra flags still cannot be classified using this different participant's record.

## Verification and reproducibility

- Native MATLAB HMM execution and cloud Python one-to-one matching completed.
- Source ZIP SHA256: `54dd72106c5d019f7bd95d18a3a613a18f5be72f929f9aa87da5036e90544651`.
- Downloaded run artifact SHA256 verified: `c817c44f863643ec83727bcbaf6ce60ddeb9cfbdd1592904d0a8e59e32d91db3`.
- Artifact ID **13469733380**, name `eye-reference-results`, includes original source ZIP, inputs, optical reference, detector masks, fitted-model metadata and initial summary.
- Local reference reconstruction matches all original reference fields and the exported eligibility mask exactly.
- Local evaluator reproduces every original match/count/mask metric; only the explicitly corrected unavailable rate and additional exploratory saccade/background fields differ.
- One-to-one assignment checks cover merged predictions and cases where greedy matching loses valid pairs.
- The follow-up Python analysis is local postprocessing of the cloud masks. No rerun of MATLAB or production algorithm is claimed for the documentation revision.
- PR51 code head remains `e8a3c5322b3ebcdbaf56ef4441c8841c5b82b0ef`; this work is isolated on `codex/pr51-eye-reference-20261106`.

Source documentation: https://www.eyetracking-eeg.org/testdata.html . Test dataset: https://www.eyetracking-eeg.org/testdata/reading.zip .

## Practical implication

The evidence supports retaining HMM while investigating the tool's other established limitations. It does not support removing HMM to rescue historical comparison trials, and it does not show that the HMM/robust union should be replaced. A larger independent dataset with valid background eye-state measurements is still required for a genuine false-positive assessment. This diagnostic does not measure final P3 error or repair sparse comparison groups.

