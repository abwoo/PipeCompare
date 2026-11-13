# PR51 P2c: actual EEGLAB 2026.0.0 native IC-subspace reconstruction — prospective correction protocol

**Precommitted 2026-11-13 before implementing or running this second hypothesis.** This is a method-measurement correction on the already repeatedly inspected development subjects 11,15,21; not a clinical/algorithm efficacy claim or independent experiment. Retain P2's earlier failed subtraction protocol and its real negative results; **do not retroactively loosen that protocol's 0.0001 μV gate**.

## Why the preceding test failed and why this hypothesis differs

The previous pre-deletion recorder derived an ICA mixing matrix `pinv(icaweights*icasphere)`, then predicted after-removal raw Pz as `original Pz - sum(removed original component projections)`. This assumption failed at real ERP CORE P3 by 1.6–5.2 μV in the 15/21 development archives, despite original Rank/P300 parity and zero source-ERP drift. Sign flipping and replacing removed by kept components did not fix it.

**Crucial source-verified EEGLAB semantics:** the *pinned actual EEGLAB release* `sccn/eeglab@2026.0.0`, `functions/popfunc/pop_subcomp.m` lines ~174–176 and 208, computes
`kept=setdiff(1:nICA,removed);`
`compproj = EEG.icawinv(:,kept)*eeg_getdatact(EEG,'component',kept,...);`
`EEG.data(EEG.icachansind,:,:)=compproj`.
Thus its output is the **native stored-mixing kept-component reconstruction**, not necessarily original sensor minus projections from a newly calculated pseudoinverse. This can also omit residual data not represented by the particular native ICA subspace. It is an issue of **scientific source provenance and data units**, not a license to alter original thresholds.

## New passive native provenance (same real ICA, no new decision)

Only when existing research `scienceCapture=true`:
- At the **actual pre-`icremove` in-memory EEG**, capture `EEG.icawinv` and the original `W=EEG.icaweights*EEG.icasphere` before any code changes the component-removal choice. Require native mixing dimension = [nICAchannels,nICs], finite values, same Pz channel in the ICA montage, complete original urevent IDs. No guessed or posthoc computed inverse substitution.
- Retain both older `pinv(W)` component Pz projections and, separately labeled, `EEG.icawinv(Pz,:).*W*EEG.data` native component projections on the **same original event window** with exact baseline correction. Record the native sensor reconstruction residual = physical original Pz minus sum(all native projected components), with per-event residual and maximum norm.
- Optionally retain a compact native-vs-pseudoinverse mixing matrix max-absolute difference and ICA matrix rank. **No changes** to ICA fit, `pop_subcomp`, removed components, filtering, channel reference, original per-trial weights, Rank or gain calibration.
- The existing after-removal observer reads the actual EEG produced by `pop_subcomp` on the same physical events, before next pipeline operation.

## Frozen falsifiable predictions, choose by actual original removed IDs

1. **If at least one IC is removed**: for every original event, `measured post-removal Pz == sum(native stored-mixing projections of all retained ICs)`; independently check `measured pre-removal Pz - measured post-removal Pz == sum(native removed IC projections) + original native subspace residual`. These are equivalent but must both be explicitly checked, not treated as two independent pieces of evidence.
2. **If no IC is removed**: `Steps.icRemove` does not call `pop_subcomp`, so `measured post==measured pre` (no-op), even with reduced-rank ICA.
3. Physical original `urevent` IDs and condition assignments exactly align pre/post and final positive-weight retained trial IDs are drawn only from those sources. No imputation of a missing event.
4. Both equations must pass a frozen `1e-4 μV` absolute and `1e-6` relative tolerance. **If they fail, preserve failure and diagnose whether `EEG.icawinv` changed during `icRemove` (e.g. ICLabel) or whether `eeg_getdatact` uses a different actual activation than `W*data`**. Do not increase tolerance or use original holdout.

## Execution and scientific gates

MATLAB and Python isolated tests must cover actual projection, a deliberately nonzero native-vs-`pinv(W)` mismatch, rank-deficient subspace residual, no-op removal, missing native mixing, exact time-locked event IDs, reconstruction failure and original EEG input immutability. New *separate* scientific validator uses the correct source-verified EEGLAB equation. The first P2 failed validator remains in GitHub history and its failed verdict is not rewritten.

Cloud real ERP CORE P3: three *previously inspected* participants (11,15,21), all six actual starts, original six-start frozen source control run `43169562528`, exact source identity, raw event/physical source/support, untouched `Measure.m` and `Rank.m`, full true ERP waveform and scientific bootstrap parity. No fresh subject, no more than original starts, no new production candidate, no main/PR #51 change.

Success constitutes a **methodological advance in directly validating real EEGLAB ICA projection semantics and correcting scientific attribution**. It does **not** prove a superior ICA recommendation or authorize PR #51 merge. The original ERP CORE 20-subject outcome remains **15/20 PARTIAL**, multimetric-v1 remains ABSTAIN, all 14 untouched subjects remain sealed.
