# PR51 independent ICA initialization stress-test — frozen *development* rules

Prepared before accessing any new ICA fits with seed indices 3, 4, 5. This document freezes the **analysis comparison within the original three development subjects**, NOT a new held-out-subject validation and NOT a new release. All changes are isolated on `codex/pr51-ica-unseen-seeds-20261112`.

## Roster, source and experiment

- **Participants:** ERP CORE P3 `[11,15,21]` already used to diagnose unstable candidate IDs. No additional participants are used.
- **Source:** the exact original 300-second P3 subset, the same per-file SHA256-pinned 20-subject source manifest in `benchmarks/erpcore_protocol.json`. EEG/ICLabel setup and event conventions remain those of the existing frozen study.
- **Independent ICA initialization:** original baseline seed 0 plus five other predefined starts 1..5.
- **Training ICA starts:** only 0, 1, 2. Candidate ranking uses the maximum original numerical objective over those three starts, restricted to candidates whose `Rank.status` is `feasible` at all three starts, same baseline-selected stratum, finite gain-corrected P300, and at least 5 HMM-unflagged effective trials per condition at all three starts. No ties guessed and no missing results imputed; method abstains.
- **Unseen ICA starts:** 3, 4, 5, **never used** to choose the candidate. For the prespecified training-selected candidate, report per-unseen-seed original `Rank.status`, objective, gain-corrected P300 target minus standard amplitude, and both conditions' HMM-unflagged effective support; evaluate whether the candidate survives all three original scientific gates. Report the original baseline recipe on exactly the same starts.
- **No threshold optimization on validation starts.** The 5-trial floor and original scientific constraints are unchanged; raw P300 and gain-corrected values must not be compared as if equal. Any failure is a negative result and must be kept.

## Interpretation constraints

This is **initialization-level generalization within observed development participants**; it cannot establish participant-level or site-level generalization. Stable candidate ID by construction in the training starts is NOT a scientific result. An improvement here only motivates formal design and data source pinning for the 14 fully disjoint ERP CORE holdout participants: `[24,25,26,27,28,29,31,32,33,34,36,37,38,39]`. None of those should be loaded before the prospective code and protocol are frozen.

The original 20-person record remains 15/20 valid paired measurements, verdict PARTIAL, independently of all newer development runs. PR #51 remains unmerged.
