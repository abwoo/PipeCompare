# PR51 prospective held-out ERP CORE validation — design only, NOT frozen (2026-11-12)

**No held-out EEG has been loaded, scored, inspected or used to select algorithms under this plan.** This file is a planning specification; it is **not** a preregistration or frozen implementation. The experiment must not be started before the independent ICA robustness method and scientific gates are fully specified and the engineering suites pass. PR #51 remains open and unmerged. Original 20-subject ERP CORE verdict remains PARTIAL.

## Disjoint participant cohort

The existing frozen ERP CORE P3 protocol defines the first twenty of 34 eligible participants as development/evaluation data:

`[1,2,3,4,5,7,8,11,12,13,14,15,16,17,18,19,20,21,22,23]`.

Its published P3 summary exclusion set is `[6,9,10,30,35,40]` (artifact rejection, target count or task accuracy), which is fixed independently of this algorithm.

**Proposed remaining 14 held-out P3-usable subjects:**

`[24,25,26,27,28,29,31,32,33,34,36,37,38,39]`.

This roster is the set difference `{1..40} \ (exclusions ∪ prior-20)`; neither sorted outcomes nor future model performance are used to choose it. This is **within ERP CORE**, not a claim of external-dataset or independent-site generalizability. Source data file identity and SHA256 for each remaining participant are **not yet pinned**. Do not run a freeze until every `.set`, `.fdt` and event file is independently resolved, checksummed and preserved in a separately versioned manifest.

## Experimental sequence and required gates

1. Finish development research on the 20 **non-held-out** participants, including 4/8 insufficient-support and 11/15/21 unstable-ID results. A numerically stable minimax rank that fails any seed's pre-existing `Rank.status=feasible` is **not valid**. Keep original `minimum_effective_unflagged=5` on each arm/condition, do not lower the ocular-artifact or signal-preservation criteria.
2. Lock the new robust-selection algorithm (code hash, exact tie-breaks, strata, failed-seed abstention, IC random seeds, candidate set, model complexity controls), all auxiliary dependencies and P300 measurement scale. Full waveform preservation and removed-component audit are required before promotion. The frozen method must be able to abstain if no all-seed eligible candidate exists.
3. Freeze a new **held-out-specific decision protocol BEFORE running any of these 14 participants**. Original 20-subject thresholds should be copied as scientifically comparable as possible; cohort-sized integer counts need a prospective, explicit justification and cannot be silently substituted for original 18/20 or 20/20. Prespecify analysis of missing pairs, precision/bootstrap CI, gain calibration, P300 signal-shape constraints, per-subject regressions, overall ICA selection stability, runtime and technical failures. Do not modify any threshold after seeing held-out outcomes.
4. Execute Standard original baseline, revised robust candidate, Minimal-fixed and Common-fixed ICA in matched fixed source windows, and retain the original HMM/VEOG flags as observed proxy information, not ground truth. Independent tests verify each `source.json` checksum, matrix settings, EEG-LAB/Matlab pin, code SHA, seed schedule, 300 s excerpt, event IDs, target/standard trial levels and all-gate eligibility.
5. Evaluate **all 14** proposed hold-out subjects, **including invalid/abstained subjects**. Coverage is not quietly changed to the subset producing finite results. Report paired sample availability and confidence intervals, all exclusion reasons, source-to-trial matching, method-specific quality and technical cost. A successful Actions job merely means the process ran; it does not establish scientific performance.
6. New-algorithm versus the unchanged original baseline is the primary comparison; full per-seed recommendation identity alone is insufficient if the new robust decision always selects the same recipe by construction. Evaluate recommendation stability across independent fit starts on held-out data. No performance claims from previous 20 development subjects count as confirmatory evidence.

## Stopping rules / no-go

- No qualified new method that passes the **existing** source, signal-shape and trial-support checks → **do not freeze or execute held-out**.
- No complete SHA256-pinned source manifest for the 14 → **do not freeze or execute held-out**.
- No all-seed-feasible candidate for an individual → abstain and report incomplete/invalid, never impute another pipeline.
- Held-out results reveal a defect → report the failed frozen experiment separately, develop a versioned next protocol, do not silently rerun as if confirmatory.
- PR #51 may be merged only after explicit independent approval; this plan provides none.

The prospective holdout is reserved for independent verification after method development, not a way to rescue the original incomplete 20-subject cohort.
