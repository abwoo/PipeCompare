# PR51 P1–P2: weak-component matching confidence and ERP association (diagnostic-only protocol)

Date: 2026-11-12. This protocol is committed **before running the new matching-confidence code**; however the existing 3 development participants and their initial low-correlation results were previously observed. This is exploratory, never a fresh validation.

## Frozen source and scope
- Use only archived opt-in six-start real-EEG ICA science records from GitHub run **43184910782**, which have already passed independent original-Run parity checks. Subjects **11, 15, 21**. Candidate roster frozen as 11→[12,25], 15→[5,2,6], 21→[25,4].
- Seeds 0–5 are exploratory. Seed-0 anchored cross-start comparisons do not provide statistically independent samples; do not pool candidate arms sharing an ICA fit as independent biological replications.
- No new EEG fitting, no production source mutation, no IC-label/removed-component changes, no original Rank re-ranking, no p-value, no 14-person confirmation cohort, no PR #51 merge.

## P1: explain whether weak maps track measured ERP changes, without causal attribution
Per predeclared candidate, for each seed 1–5 compare its actual gain-corrected Pz difference waveform to seed 0 using full-epoch RMS and 300–600 ms mean-amplitude difference. Record archived bootstrap SE and its change relative to seed 0, along with counts of weak matched scalp maps (absolute centered correlation below **0.90**, the existing unchanged screening value), broken removal agreement and ICLabel top-class changes. Report each start separately and positive/negative findings; do not claim a retained component's *causal contribution* to ERP. The saved matrices are mixing-only; **component activations are absent**. The sensor-loading ratio of a weak component, if reported, is geometry-only and cannot estimate a time-locked neural contribution.

## P2: determine whether the weak pairs could be assignment artifacts
For each seed-0-to-seed-S decomposition, align electrodes by exact channel label, demean each mixing-map column, normalize to unit L2 norm and form absolute correlation. Compare the original one-to-one global Hungarian assignment to **unconstrained nearest topographic neighbour** in both directions. A weak matched pair is:
- **intrinsic-weak** when the seed-zero component has no counterpart with |r| >= 0.90 anywhere among seed-S components;
- **assignment-conflict** if it has at least one strong counterpart, but one-to-one assignment chooses a weaker match;
- **mutual-nearest** if the assigned seed-zero and seed-S components are each other's strongest absolute-correlate partners (tie means ambiguous).

Also compute an **assignment score-margin diagnostic**: global matched-correlation sum minus the best global matching forbidding a given pair (not a probability or validated threshold). A zero/near-zero margin signals non-uniqueness/ambiguity. Compare the direct seed0→seed2 correspondence with the composed seed0→seed1→seed2 assignment (cycle consistency); a mismatch is evidence of assignment/path dependence, not a diagnosed biological failure. All metrics, including unavailable cases, must be retained; **do not relax the existing min-correlation 0.90 gate**.

## Stop conditions
- Missing source identity, decomposition, channel-label alignment, finite mixing maps, removed-component provenance, complete waveform or bootstrap estimate → fail closed and preserve prior negative decision.
- The prospective 14-subject roster must remain untouched. Even a more internally consistent match on these subjects does **not** justify selecting a new winning recipe or merging PR #51.
- If P1 cannot attribute components to P300 from existing snapshots, design a separate optional pre-removal activation/readout stage with parity tests before any new cloud EEG fit. Avoid reconstructing component contributions from mixing alone, or treating scalar amplitude changes as signal truth.

**Primary reporting:** Component-level weak-match explanation and observed waveform/precision changes, plus method limitations. An upstream minimum 0.90 failure remains a failure, irrespective of matching ambiguity.
