# PR51 P2d prospective analysis: exact EEGLAB native ICA removal, fixed-event P300 and final ERP

**Frozen before any real P2c native `EEG.icawinv` result has been examined** (2026-11-13). Use only historical-real-development subjects 11,15,21, six actual ICA starts, and the immutable original Standard candidate IDs (11:12;15:5;21:25), plus previously predeclared alternatives 11:25;15:2/6;21:4. No independent-subject confirmation or promotion can come from these repeatedly reviewed participants.

## Required source and provenance gate

Use only P2c [real native EEGLAB run 43408389159](https://github.com/abwoo/PipeCompare/actions/runs/43408389159), and only if its entire workflow is **SUCCESS**, with all participant-level original source/Rank/P300 identity checks and direct `nativeRealEEGLABProjectionVerified` tests passed. Reject any one-participant partial run, source SHA mismatch or independent trial/weight mismatch. This is a new physical equation; the earlier failed `pinv(W)` subtraction study stays failed.

## Fixed per-stage computations

For each predeclared candidate, seed0 vs each of seeds1..5, use matched original physical `urevent` IDs and conditions, and the original recorded final retained event/positive weight sets F0, Fs. Let `D(sensor,F)` denote baseline-corrected 300–600ms **Pz target-minus-standard at that physical stage**, not calibrated final P300.

- Pre-step fixed-trial signal `P=D(pre_s,F0)-D(pre_0,F0)`;
- Actual native post-`pop_subcomp` fixed-trial signal `R=D(post_s,F0)-D(post_0,F0)`;
- Immediate ICA-removal contribution `C=R-P`, verified from separately stored **actual** native retained component reconstruction/no-op step (not independent mathematical evidence, but a separate source parity check);
- Trial retention/weight term `Q=D(post_s,Fs)-D(post_s,F0)`, total post-step change `T=R+Q`;
- Distinct original *final gain-corrected* ERP P300 drift, full-epoch waveform RMS and original bootstrap precision from previously independently validated scientific readout. Never treat `T` or `C` as equal to the final ERP estimand; they are not automatically commensurate, because filters, re-reference, epoching, artifact rejection and gain calibration can occur after the stage.
- Native IC pseudoinverse-to-actual-mixing discrepancy and pre-step native subspace-residual maximum; disclose missing effect identities and whether changes in removed-component set correlate with signal changes. Observed correlation is not physiologic causal truth.

Compute original Standard vs each candidate's absolute cross-seed stage drift and the unchanged final ERP stability, but treat this as **exploratory**, with no prechosen winner and no newly lowered noninferiority bounds. Preserve all good and bad arms including seed0 and both target and standard. Subject-level repeated recipe arms are not biologically independent units.

## Fail-closed checks / disposition

Original Rank/Measure/experiment protocol and calibrated P300 parity, source unflagged and effective-trial constraints, native EEGLAB full source reconstruction, matched physical events and weights must pass. If any absent, no result; no substitution from earlier failed artifacts. No threshold relaxation of topology r>=0.90, positive weighted trial floor≥5 or previous multimetric v1 criterion. Never open unused 14-person cohort or merge PR #51 because a source equation works.

Success is an **inspectable method-level physics/signal attribution advance** and a means of predeclaring a better stability metric later on new subjects, *not* a validated superior ICA selector.
