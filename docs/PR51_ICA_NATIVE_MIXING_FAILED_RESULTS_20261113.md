# PR51 P2c native mixing hypothesis also rejected by real ERP CORE EEG

Date 2026-11-13. [Native `EEG.icawinv` real-EEG workflow 43408389159](https://github.com/abwoo/PipeCompare/actions/runs/43408389159) independently passed original Rank/P300/ICA source parity and archived the missing native-mixing fields but failed the source-verified post-subcomp reconstruction hypothesis on subject15: maximum actual-vs-predicted Pz event amplitude error **1.695618 μV**, not within predeclared 0.0001 μV. No numeric gate was relaxed.

[Read-only GitHub archive-forensic 43409512876](https://github.com/abwoo/PipeCompare/actions/runs/43409512876) inspected **all 18** source candidate × seed pairs for subject15, without new EEG fits:
- `max(abs(EEG.icawinv - pinv(W))) == 0` in all, and component-level Pz amplitude difference exactly **0** for native vs pinv, so the matrix identity **does not explain the discrepancy**.
- Native all-component Pz residual before deletion stays about **4.6×10⁻⁶ μV**; this is far smaller than the actual source discrepancy, which ranges about **1.62–1.70 μV** in subject15.
- The source effect is therefore **not** explained by a different mixing matrix, a large pre-removal source residual, a flipped sign, or simply reusing kept instead of removed IC indices.

**Next falsifiable candidate:** EEGLAB's real `pop_subcomp` uses `eeg_getdatact(EEG,'component',kept,...) `. The pinned `eeg_getdatact.m` (2026.0.0) may use `EEG.icaact` cached component activations (or on-disk cached ICA activities) rather than recomputing `W*EEG.data`. The earlier observer used `W*EEG.data` on a private epoch copy and could not prove those activations matched the ones consumed at deletion time. An independent, passive **in-step** probe just before and after actual `pop_subcomp` is needed; it must record whether the cached projection matches actual native output and compare cached vs fresh activation without changing the EEG or removed components.

This is a reason to strengthen scientific measurement provenance, not a basis to change the code's ICA recommendation behavior, relax r≥0.90 topology stability, touch the unused 14-person holdout, or merge PR #51. Original ERP CORE **15/20 PARTIAL** unchanged, v1 multi-metric selector still ABSTAIN.
