# PR51 Stage 3C — read-only ICA / ERP scientific recorder

**Research-only branch** `codex/pr51-science-recorder-20261112`. No production default change, no PR51 merge, no access to the untouched fourteen-subject holdout. Original preregistered ERP CORE verdict remains **15/20, PARTIAL**.

## Intent and observational constraints

Earlier analyses found the all-seed-feasible Minimax and stability-first candidates on subjects 11/15/21 without an improvement in unseen-seed P300 amplitude stability relative to original ranking. Candidate-ID stability cannot substitute for independent ICA component recovery or ERP preservation. The next step is to observe the complete scientific state *without feeding it back into ranking*.

The research switch `PIPECOMPARE_PR51_SCIENCE_RECORD=1` is disabled by default. When enabled, the existing `Executor.evaluateLeaf` computes its original Measure and Rank-relevant fields *before* a separately guarded `ScienceRecorder.capture` call. A recorder exception is reported in `cand.science`, **never converted into a failed scientific pipeline**. Each seed's observation is stored alongside the original candidate's seed-specific measurement; no new seed is chosen and no original gate is changed.

The capture uses the **existing evaluated candidate** and final epoched EEGLAB array, not a fresh ICA computation or a second signal processing route.

## Captured measurements and limitations

- **ICA:** per-path fitted unmixing (`icaweights × icasphere`), inverse mixing/topography, channel indices, component count and actual `icremove` step's removed IC indices. Only original fitted matrices are represented. Component ordering, sign and scale are not invariant and require aligned spatial comparisons before interpretation.
- **ERP:** full baseline-corrected **Pz** condition-mean waveforms over the complete -200…800 ms epoch and their target-minus-standard difference. They are explicitly marked **uncorrected for gain**; use per-seed `gain` separately and never equate them with the gain-corrected point estimate.
- **Trials:** original urevent-based identities, condition reference IDs, kept IDs, original per-trial scalar P300 window scores and trial weights. This supports a second independent alignment and support check without replacing retained trial decisions.
- **Bootstrap:** independently calculated by a Python collector from each condition's stored per-trial scalar scores and weights, with **1000 deterministic resamples** and fixed dedicated random state. The recorder deliberately does not invent a MATLAB bootstrap SE; it records that independent recomputation is needed.
- **Missing at this stage:** independent ICLabel per-component probability vectors, trial-specific entire waveform matrices, full spatiotemporal maps, and per-subject cross-site stability. Do not claim them as recorded.

## Verification gates

1. The native full Standard search, three fixed outputs and replay/export continue to run with the pinned source and MATLAB/EEGLAB identities.
2. All six ICA initialization seeds (0–5) and all measured candidates get a science snapshot, with original candidate ID and seed, valid ICA topography dimensions and removal indices.
3. The complete waveform time axis, condition count, original urevent ID alignment and trial-score/weight lengths are checked independently.
4. Bootstrap point estimates from the recorded original trial scores must equal the gain-corrected MATLAB scalar P300 to within 1e-5 μV; bootstrap SE must be finite. The historical trial floor remains unchanged and **technical recorder success does not authorize scientific promotion**.
5. Separate cross-run noninterference replay compares the new six-seed experiment to original [Actions 43169562528](https://github.com/abwoo/PipeCompare/actions/runs/43169562528), including original recipe IDs, Rank feasibility, objective scores, source fingerprints and P300 values.

The research workflow is [`.github/workflows/pr51-science-recorder.yml`](https://github.com/abwoo/PipeCompare/blob/codex/pr51-science-recorder-20261112/.github/workflows/pr51-science-recorder.yml). It runs the same already-studied **three** ERP CORE subjects [11,15,21] in separate real-MATLAB Actions jobs; CI is not a new frozen 14-person study.

**Decision gate:** If recording fails or changes candidate/ERP results, treat this as an instrumentation bug and fix the recorder before any new algorithm. If recording succeeds, use recorded spatial components, removed-component identities, waveform morphology and effective trial support to form a *predeclared* candidate that must show scientific benefit on genuinely independent evidence before a new 14-person frozen evaluation. No production adoption or PR51 merge is implied by the recorder's existence.
