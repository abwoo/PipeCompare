# PR51 14-subject ERP CORE P3 holdout: source provenance locked, scientific experiment NOT yet frozen

**Date:** 2026-11-12.

## Fully verified source-identity work

- The original ERP CORE P3 metadata-only OSF crawl `https://github.com/abwoo/PipeCompare/actions/runs/43163072602` completed successfully, with 5/5 tests. **No heldout EEG signal bytes were downloaded.**
- All 14 predeclared unseen subjects **24,25,26,27,28,29,31,32,33,34,36,37,38,39** have three files each: `_eeg.set`, `_eeg.fdt`, `_events.tsv`, total **42**.
- GitHub branch `codex/pr51-holdout14-design-20261112` now contains the exact source-identity snapshot at `benchmarks/results/holdout14_osf_candidate_manifest_20261112.json`. The JSON bytes SHA-256 is **`7cabbd348c0267a1b42e140be7e627f1c0ba7873ac6b87997544d6035f61aeb6`**. Each of 42 entries records public OSF URL, individual OSF metadata SHA-256, file size, OSF file ID, and subject identity. The fully validated cloud provenance commit workflow succeeded: https://github.com/abwoo/PipeCompare/actions/runs/43164198449.
- The 14-person cohort selection and original 6 P3 exclusions have previously passed 6/6 regression tests and an independent nonoverlap audit: https://github.com/abwoo/PipeCompare/actions/runs/43162284848.
- OSF-provided hashes are **source metadata claims**. At subsequent freeze-run time, the actual downloaded file bytes MUST be independently SHA-256 checked against this manifest before the corresponding participant is admitted to analysis.

## Exact distinction between provenance and confirmatory freeze

This is an **inspectable candidate source manifest**, not yet a complete scientific freeze. `benchmarks/erpcore_holdout14_preanalysis_20261112.json` deliberately retains `status=PRE-FREEZE-DESIGN-NOT_EXECUTED`, source `status=UNPINNED` in the scientific protocol and algorithm `status=UNSELECTED` because measurement-preserving ICA recommendation and two independent seed triplets have **not** yet been demonstrated/approved. This metadata document does not silently toggle those statuses.

Before running any of the 14 *unseen* EEG cases, scientifically choose and pin:
- Single exact algorithm implementation Git commit; explicitly compare original pipeline ranking, common-feasible minimax, and estimator-family-restricted alternative as different methods on development data only.
- Two independent ICA initialization groups, and a ranking/recipe ID comparability protocol that prevents same-seed mechanical determinism from masquerading as independent recommendation stability.
- Immutable 14/14 participant-by-participant missingness/abstention reporting, separate original and fixed controls, a same-scale P300 gain/SME/shape audit, original HMM-unflagged effective trial floor >=5 each arm/condition.
- Immutable acceptance, bootstrap and calibration gates and independent evidence evaluator; no outcome-based tuning of the heldout 14.

**No 14-subject EEG result has been produced.** The original 20-person PR51 remains 15/20 valid pairs and **PARTIAL**; PR51 remains unmerged.
