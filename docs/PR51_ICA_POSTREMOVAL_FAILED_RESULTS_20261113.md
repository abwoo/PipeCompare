# PR51 P2 initial post-removal hypothesis rejected — real EEG forensic record

Date 2026-11-13. Original protocol `docs/PR51_ICA_POSTREMOVAL_PROTOCOL_20261113.md` was committed **before** new evidence. It predicted `actual Pz before−after pop_subcomp = sum(pre-observer removed components)` from a pre-observer using `pinv(icaweights*icasphere)` as its mixing matrix. This was not demonstrated on real EEG.

## Failed real scientific experiment, gates never adjusted

[Real six-start GitHub Actions 43405859824](https://github.com/abwoo/PipeCompare/actions/runs/43405859824): Python and MATLAB isolated science unit tests passed. In the actual development EEG results for subjects 15 and 21, original frozen ERP CORE P300/Rank science parity and original pre-removal recorder both **passed** independently, but the new post-removal subtraction validation **failed**:
- Subject 15, example original candidate 5 seed0: max event-wise Pz prediction discrepancy **1.695616 μV**.
- Subject 21, example original candidate 25 seed0: **4.807143 μV**.

The numeric absolute tolerance was **0.0001 μV**, unchanged. This is an unequivocal failed scientific model, not a slight floating-point residual. Subject 11 is not counted as validated unless its own test completed.

An independent **archived-data only** [Actions 43407447651](https://github.com/abwoo/PipeCompare/actions/runs/43407447651) compared actual post-step Pz, original ICA `pinv(W)`-based removed predictions, opposite signs and incorrectly keeping rather than deleting components. All 18 comparisons in subject15 and 12 comparisons in subject21 remained above the frozen tolerance; flipping sign and summing kept instead of removed made the error worse. The observed effect roughly tracked the predicted signal but contained systematic non-negligible residuals; this alone does **not** identify the physical origin.

## Source-level mechanistic hypothesis, to test prospectively

The official EEGLAB `pop_subcomp` source **reconstructs the kept components using the in-memory `EEG.icawinv` matrix** and then replaces the original ICA-channel EEG data with that reconstruction. It does not guarantee equality to `originalEEG.data - sourceRemovedProjection` as estimated using an independently recomputed pseudoinverse. In rank-reduced/partially captured data, the discarded out-of-ICA-subspace residual and any difference between the recorder's `pinv(W)` and EEGLAB's actual stored `icawinv` can matter.

Reference implementation: `https://github.com/sccn/eeglab/blob/develop/functions/popfunc/pop_subcomp.m`, in particular its `component_keep`, `compproj=EEG.icawinv(:,component_keep)*eeg_getdatact(...)`, `EEG.data(EEG.icachansind,:,:)=compproj` statements. The real experiment pins EEGLAB 2026.0.0; verifying precisely that version's runtime matrix is a separate requirement.

**Next specific falsifiable hypothesis:** record the original actual in-memory `EEG.icawinv` at the very step, together with `W*EEG.data` activations, the ICA-subspace reconstruction residual, and the resulting sensor waveform. For a nonempty removal, predict the **measured post-removal Pz** from the sum of retained components using the *actual recorded native* mixing coefficients, not `pinv(W)`. With zero removed components, predict identity/no-op. Keep the original failed `pinv(W)` hypothesis and its threshold unchanged and explicitly separate; do not retroactively recode it as a pass.

This is a **research measurement semantics** correction to be tested in another isolated branch, **not** a claim that an algorithm improves neural physiology. No holdout, recipe promotion, fallback threshold relaxation, or PR #51 merge. Original 15/20 PARTIAL remains.
