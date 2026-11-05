# Frozen detector diagnostic, 2026-11-05

No production edits, no merge, no post-result threshold tuning. Execute native `pipecompare.eval.Hmm.events` (PeakZ 3.5, MinSeconds 0.1) from the previously audited branch. This is a detector diagnostic, not a full Standard/P3 accuracy comparison.

## External data

PhysioMotion Artifact, OpenNeuro ds006386, DOI 10.18112/openneuro.ds006386.v1.0.1, CC0. Source repository tree inspected: 647f61744019ad41801598384bb461b7fdc06560. Public S3 objects are fetched at execution; headers, downloaded byte ranges, annotation and event bytes are SHA256 recorded. These range digests identify downloaded excerpts, NOT a full-file checksum or immutable S3 version guarantee.

Subjects 1–10, run 01, no outcome-based substitution. Two 60-second slices each: first `open_base` annotation start + 2 seconds, and first event `trial_type=4` (instructed blinking) onset. 20 minutes total. No filtering/resampling beyond unchanged detector internals. Primary input Fp1; prespecified sensitivity input both Fp1/Fp2. EEG channels are used as proxies, not claimed independently recorded EOG. EDF physical-unit conversion is explicit. Download only covering EDF records, not the entire dataset. Missing inputs fail visibly; do not silently drop subjects.

Compare HMM alone, robust 5-MAD point rule and their union, with fixed 50 ms dilation. Report fraction marked, HMM-only fraction and segment counts per person/block. Main descriptive contrast: blink minus open block fraction. These quantities are NOT sensitivity/specificity, nor are HMM segments individual blinks.

Critical source limitation: the human annotations span tasks and may contain multiple artifacts; an entire blink task is not continuously blink-contaminated. Open baseline is not certified blink-free. No event-level false-positive/false-negative estimates can be made from these blocks. There are no dedicated EOG channels. Human annotations were made on processed bipolar EEG with task guidance, not independent video eye-state measurements. This external check cannot resolve whether historical S005 IDs 3/31/38 are false positives.

## Known-truth controls

10 fixed random seeds (1000–1009), 60 seconds at 250 Hz: stationary AR(1) coefficient .8 noise standardized to unit SD plus unit-amplitude 10 Hz sinusoid. Pair each background with compact raised-cosine pulses centered at 5,10,...,50 seconds, half-width .3 s and peak amplitudes 4,8,16. Null arm amplitude zero. Fit the unchanged HMM separately per trace; report valid fit status.

Event detection means any padded HMM state within ±.1 s of each known center. Report hit/total events (100 per amplitude pooled), occupancy outside the exact injected support, and null occupancy. Outside-support occupancy includes the fixed dilation; it is not event-level precision. These controls have synthetic truth but do not establish performance on real natural blinks. No changes are promoted from this diagnostic.

## Decision and recovery

Completion requires 10/10 subjects, 40 external outputs (2 blocks × 2 channel configurations), 40 synthetic outputs and finite metrics within [0,1]. Invalid HMM fits remain in outputs and must be discussed, not removed. Failures are infrastructure/data failures, not evidence of bad scientific performance. Low detection or high occupancy is recorded without weakening assertions to manufacture a scientific pass. Numerical execution success does not imply the detector is accurate.

Run independent Python count/metric checks. Each subject's inputs/provenance and each trace's masks are retained as artifacts; rerun only after identifying failures, verify identical source hashes before pooling across attempts. No adaptive threshold change, no repeated ten-person P3 search. An independently video/EOG-labeled event-level dataset remains needed for natural-blink error rates; Georgia Tech BLINK was considered but its published download currently requires a contact form, and no personal information was submitted.

Sources: https://www.nature.com/articles/s41597-026-07120-7 ; https://github.com/OpenNeuroDatasets/ds006386 ; https://gnan.ece.gatech.edu/eeg-eyeblinks/

## Data-availability amendment before any detector outcomes

Initial run 41957688403 stopped while parsing subject 4: no open_base annotation exists. No HMM outcomes were generated. Keep all ten participant slots; record missing block/configuration slots explicitly and run available blocks without substituting participants or inferring baseline from another task. Completion checks require 40 accounted-for slots, not falsely claim 40 measured outputs. Available-block averages have explicit denominators and are descriptive, not paired full-cohort contrasts. Missing comparisons remain unresolved.
