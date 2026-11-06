# Frozen optical-reference diagnostic

Small feasibility follow-up to PR51. No production changes or merge. One openly downloadable EYE-EEG natural-reading demonstration participant (~102 s total); do not split into pseudo-independent people. Dataset source https://www.eyetracking-eeg.org/testdata.html and https://www.eyetracking-eeg.org/testdata/reading.zip . Source ZIP SHA256 54dd72106c5d019f7bd95d18a3a613a18f5be72f929f9aa87da5036e90544651. Download limited to 14.4 MB; raw data retained in run artifact, not committed.

## Reference and synchronization

Use the full interval between first and last shared trigger 3. Pair all 20 ordered EEG triggers with EyeLink MYKEYWORD 3 messages. Fit clock slope/offset from endpoints only; require each of 18 interior residuals <=10 ms and no interior EEG boundary. Timing is not adjusted using blinks or EEG amplitudes. Metadata-only inspection before detector execution found max interior residual ~1.94 ms. No detection outcomes have yet been inspected.

Use overlapping left/right EBLINK intervals, union their bounds and merge overlapping resulting intervals. Retain only full events within synchronized bounds. These are EyeLink's optical classifications, independent of EEG amplitudes, but NOT human-video adjudicated truth. Unilateral labels are not negative labels. Exclude ±100 ms around all optical labels from background eligibility; add binocular reference intervals back. Otherwise require both pupil areas >0 and a recorded optical sample within 3 ms. Unknown/dropout time is not certified negative. The site documents prior 0.1–100 Hz filtering and average reference; this is not unprocessed EEG.

## Detector and metrics

Run unchanged production Hmm.events with default PeakZ 3.5 / MinSeconds .1. Primary input explicitly defined as two bipolar proxies IO1-Fp1 and IO2-Fp2; report Fp1-only HMM as sensitivity check, not optimized channel selection. Compare raw HMM, 50 ms padded HMM, robust 5-MAD point rule with same padding, and their union. No fitting to optical labels. Fit each detector on the synchronized continuous excerpt only.

Match predicted intervals one-to-one to optical intervals if they overlap after allowing ±100 ms on the reference interval, maximizing match count. Report matched/missed reference events and reference recall. For unmatched predictions, distinguish >=80% reference-eligible samples from indeterminate intervals. Report unmatched eligible predictions per eligible minute, all predicted intervals, eligible fractions and segments >2 seconds. Unmatched is reference disagreement, NOT certified physiological false positive; missing optical labels and EEG artifact other than blink remain possibilities. Do not report unqualified accuracy or precision from incomplete reference coverage.

## Gates and interpretation

Reject interpretation if source digest, sample count, trigger correspondence, clock residual, continuity, finite EEG values or HMM fit validity fails. No threshold relaxation to manufacture a pass. Export inputs, optical timestamps/intervals, all detector masks and match identities. Independently evaluate masks in Python and keep all unmatched events available for review. Sparse reference events or large disagreements are reported, not grounds for choosing a different time window. This is feasibility/reference-agreement evidence in one person, not ten-person confirmation. No HMM repair or P3 rerun is authorized by a favorable result alone.

The SET metadata references a differently named FDT; the archive's actual reading_eeg.fdt is read explicitly after exact 72 × 52225 float32 size verification. No manual file rewriting or EEG resampling is used.

Before detector output inspection, add a descriptive breakdown of unmatched eligible predictions by overlap with optical ESACC intervals (zero additional tolerance). Hmm.events is explicitly designed for blinks AND eye movements. Thus saccade-related predictions are not necessarily erroneous blink detections. The blink-reference recall remains the same frozen endpoint; saccade overlap is explanatory, not an alternate pass gate.

## Post-run reporting correction

The source ASC pupil-area fields were subsequently found to be zero for every sample. Therefore no non-blink background passes the frozen optical-availability rule. The initial cloud summary's zero unmatched-per-eligible-minute value is not interpretable as a false-positive estimate. The follow-up evaluator reports this rate as null whenever eligible non-blink background is zero. Detector masks, event matching, time alignment and all thresholds are unchanged; the original cloud artifact is preserved. Saccade overlap is an exploratory annotation added outside the original executed commit, not a new independently validated endpoint.
