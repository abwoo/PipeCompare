# PR51 P1b: passive pre-IC-removal component event projection — locked research protocol

Created 2026-11-12 before any new real-data recorder run. Scope: already explored real ERP CORE P3 development subjects **11, 15, 21**, ICA seeds **0–5** and predeclared candidate IDs (11: 12,25; 15: 5,2,6; 21: 25,4). No prospectively independent subjects and no production algorithm changes.

## Question and existing gap

The old six-start science run contains mixing maps, ICLabel probabilities, IC removal identities, retained-trial IDs and final Pz ERP waveforms. It **does not** contain the removed ICs' event-locked activations. Mixing maps alone cannot establish whether a weakly correlated, retained or removed component contributes to P300. A separate run must capture data *before* `pop_subcomp`, without affecting the actual EEG passed to the original pipeline.

## Recorder specification, no treatment effect inference

An opt-in observer, under the pre-existing `scienceCapture=true` **research-only** flag, runs on the in-memory pre-`icremove` EEG immediately before the original cleaning step. Its input EEG is immutable. Observer failure must be represented as a science failure; never alter the main candidate status, ICLabel decision, weights or original Rank computation.

Compute the original ICA linear operator `W=icaweights*icasphere` and mixing map `A=pinv(W)` after aligning EEG channel indices and Pz. On a private copy of the pre-removal data, select only pre-declared target `11` and standard `31` ERP epochs using the **existing contract** and original urevent time-lock IDs. For each eligible physical event, apply the original component time-series projection `A(Pz,k)*(W*EEG.data(ICA_channels))_k`. Subtract the same private-copy baseline average from each component, then average 300–600 ms, storing the per-event component × trial values in μV.

For each event, also compute the pre-removal sensor Pz baseline-corrected 300–600 ms mean amplitude and the numerical residual between that sensor value and the sum of all linear component projections. **Do not silently force a zero residual:** rank-reduced decompositions may leave data outside the ICA subspace, which is scientifically relevant. Record rank, dimensions, event urevent IDs, condition identities, the actual labels and all residuals. No ground-truth assumption.

The subsequent candidate's final retained trial IDs and original per-trial weights can later be joined to these pre-removal signatures. However any such weighted sum is a **pre-removal stage-specific quantity**, not the final post-filter/post-removal/rejection P300 estimate. No full component-to-final-ERP causal attribution without an explicit counterfactual/pipeline linearity verification.

## Scientific integrity and execution

- Only code paths gated by the existing `scienceCapture=true` are modified; when false, original pipeline behavior stays identical.
- The observer executes before the actual `icremove` step but its return values do not influence that step. On failure, a science-failure record is saved without rejecting a candidate.
- Original `Measure.m`, `Rank.m`, source protocol, effective trial floor 5, ICA initialization seeds, calibration, artifact decisions, candidate ranking and all originally frozen study gates **remain unchanged**.
- Pure MATLAB unit tests must show no input EEG mutation, exact projection for known synthetic mixing data, honest non-zero rank-deficient residual, and missing/misaligned data fail closed. The existing P300/science source-parity validator must confirm the real re-run preserves all prior scientific Rank/P300 outputs.
- This is an *additional real fit on previously used development subjects*, not a new frozen experiment or evidence of improved recommendation performance. The 14 untouched subjects remain sealed.

**No PR #51 merge or production promotion from readout alone.**
