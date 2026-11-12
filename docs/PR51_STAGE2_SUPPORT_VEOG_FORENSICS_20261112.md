# PR51 — stage 2: frozen trial-support and VEOG forensics

Date: 2026-11-12. This is exploratory forensics, **not** a new confirmatory experiment. PR51 remains OPEN and unmerged. The original 20-person P3 validation remains PARTIAL (15/20 valid pairs).

## Pinned sources and evidence

- Frozen original GitHub run: https://github.com/abwoo/PipeCompare/actions/runs/42809789847.
- Original algorithm SHA: `7151821edf8a5a54c3167d7f18b66c530911b568`.
- Original protocol SHA256: `bb65c9a522e9f45210e5ee941a669630b75587a5255664aa0f1402b572e89c6b`.
- 20-person / 60-arm derived-marker support audit (6/6 tests and artifacts): https://github.com/abwoo/PipeCompare/actions/runs/43149415431.
- Per-event source-SHA-verified VEOG probe (6/6 tests and 5/5 source-subject jobs): https://github.com/abwoo/PipeCompare/actions/runs/43150159789.
- Original selected Standard pipeline full per-event event IDs, HMM flags, weights and P3 scores for subjects 4/8 are preserved in `benchmarks/results/erpcore_support_trial_20261112/selected_004.json` and `selected_008.json`. These are original MATLAB-result-marker-derived arrays, **not raw EEG or manually adjudicated blink truth**.

## Independent VEOG per condition

The fixed 50 µV cut refers to baseline-corrected peak absolute VEOG during the -0.2 to +0.8 s event epoch. It is a diagnostic proxy, not an artifact ground truth label.

| Subject | Target VEOG >50 µV | Standard VEOG >50 µV | Target median µV | Standard median µV |
|---|---:|---:|---:|---:|
| 4 | 17/18 | 94/95 | 235.21 | 239.95 |
| 8 | 21/22 | 83/95 | 233.22 | 219.37 |
| 17 | 29/32 | 106/126 | 107.66 | 99.12 |
| 18 | 26/26 | 114/121 | 153.32 | 135.55 |
| 23 | 31/31 | 127/135 | 165.98 | 162.83 |

## Cross-tabulation using the original HMM flags for selected Standard pipeline

Subject 4: target HMM flagged 14/18 (14/14 VEOG above 50 µV); unflagged 4 (3/4 above 50). Standard HMM flagged 92/95 (92/92 above 50); unflagged 3 (2/3 above 50). Standard effective-unflagged target/standard support 4/3, fixed comparators 1/1.

Subject 8: target HMM flagged 17/22 (17/17 VEOG above 50); unflagged 5 (4/5 above 50). Standard HMM flagged 75/95 (75/75 above 50); unflagged 20 (8/20 above 50). Standard effective-unflagged support 5/20, fixed comparators 4/18. Frozen trial-support floor >=5 per condition and comparator remains unchanged.

Subjects 18 and 23: original Standard has no recommended pipeline so its saved selected-support arrays are **missing, not zero**. Both fixed comparator arms retained zero targets, and independently the target VEOG proxy exceeds 50 µV in all target events (26/26 and 31/31). Engineering exception fixes now complete the search but still return zero feasible Standard candidates. No valid pair has been restored.

Subject 17: Standard signal-shape failure was 31 ms P3 peak-latency deviation (frozen shape gate), fixed comparator support target4. Independent target VEOG >50 µV in 29/32. Do not relax shape or effective-trial gates to manufacture pass outcomes.

## Critical original MATLAB secondary diagnostic bug

The frozen `benchmarks/erpcore_pilot.m` local `hmmVeogDiagnostic` checks for the literal condition string `target`, but `selectedSupport` uses actual labels `11` and `31`. Thus the legacy aggregate *mis-grouped target events with standard events*; old `hmmVeog.nFlagged`, medians and `flaggedWithVeog` should NOT be interpreted as a correct per-condition breakdown. It did not change the original P300 estimate, pipeline ranking, candidate eligibility, raw EEG or selectedSupport arrays; it affected the secondary summary only.

A new, **versioned and non-mutating** `benchmarks/p3_hmm_veog_diagnostic.m` accepts numeric aliases `11` and `31`, preserves target and standard metrics separately, and validates indices. It does **not** rewrite the frozen original source recorder or previous results. It is covered by `tests/test_p3_hmm_veog_diagnostic.m` and `.github/workflows/pr51-hmm-veog-matlab-v2.yml`.

## Scientific decision: data observability, not post-hoc rescue

Original blink reference flags are measured on original EEG and reused for artifact-bias identifiability. In `Measure.eyeBias`, where all available target trials are flagged, there are zero artifact-free comparator events and bias is unassessable (`NaN`); a visually cleaned result alone cannot retrospectively supply measured clean reference trials. The absence of nonblink support is a *methodological missingness / observability limitation*. Every raw subject and failed arm remains in the record.

Do not change the pre-existing >=5 effective unflagged support floor, frozen signal-shape gate, or HMM threshold after viewing these results. A prospective eye-correction method would require separately assessed signal preservation, independent/held-out clean-reference data, explicitly versioned eligibility rules and a new freeze before an independent evaluation. Original 20-subject development sample cannot be reused to claim independent confirmation.

## Ordered next gate

1. Finish the wide bugfix `test_engine` MATLAB suite and the new versioned diagnostic unit tests.
2. Stage 2 evidence remains **diagnostic**, with no assertion of recovered clean trials. Preserve originals and missingness.
3. Proceed to ICA *measurement-level* robustness: subject 11/15/21 alternate seed final P300 estimates, SEs, IC removal, shape and trial support, not only candidate ID.
4. Freeze updated code+protocol+population before truly held-out ERP CORE or independent external P3 confirmation. Never merge PR #51 without explicit instruction.
