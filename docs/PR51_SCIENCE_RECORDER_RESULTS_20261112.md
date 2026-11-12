# PR #51 Stage 3C — six-seed real-EEG scientific recorder: completed development audit

**2026-11-12; status: DEVELOPMENT ONLY, scientific method NOT promoted.** PR #51 stays OPEN / UNMERGED. Original preregistered 20-person ERP CORE result is **15/20 valid pairs, PARTIAL**, unchanged. No untouched 14-person ERP CORE participants were accessed.

## Executions and instrument repair

- Research recorder and native six-start MATLAB EEG jobs: [Actions 43189687480](https://github.com/abwoo/PipeCompare/actions/runs/43189687480), source `1b7a7a896bdc157d499df8a4e4de179fa38ceb03`. All 3 subjects saved complete MATLAB result artifacts, but subject 21's post-hoc Python recorder validator failed because MATLAB `jsonencode` writes unavailable numerical trials as JSON `null`. This run is **NOT a fully successful workflow**; do not misstate it.
- Fix isolated on `codex/pr51-science-recorder-json-fix-20261112`: explicitly convert numerical scores, weights and waveforms to float arrays before NumPy `isfinite`. This preserves unavailable trial scores as NaN and never counts them as valid observations. Regression [43190760439](https://github.com/abwoo/PipeCompare/actions/runs/43190760439) **PASS**. No MATLAB, data, source events, signal processing or scientific threshold changes for the Python fix.
- No-refit, full saved-artifact recovery: [Actions 43191798394](https://github.com/abwoo/PipeCompare/actions/runs/43191798394) **SUCCESS**, with archived `pr51-science-reanalysis-43191798394` containing all scientific records, parity comparisons and topology/ERP analysis. This reuses the ORIGINAL all-three MATLAB artifacts, not recomputed EEG.

## Record identity and scientific validity

| Subject | Candidate × ICA start snapshots | Valid audited | Failed |
| --- | ---: | ---: | ---: |
| 11 | 120 | 120 | 0 |
| 15 | 114 | 114 | 0 |
| 21 | 114 | 114 | 0 |
| **Total** | **348** | **348** | **0** |

The scientific recorder captures each evaluated ICA decomposition's weights/sphere, per-channel mixing/scalp topographies, IC-removal indices, complete Pz -200..800 ms condition and difference ERP averages, source urevent identities, original P300 per-trial scalar scores and weights. Independent Python computation resamples these aligned scores **1000 times** for the gain-corrected P300 bootstrap SME, checking its mean point against the original MATLAB value. Artifact validity does not establish signal correctness or benefit.

Original IDs, seed recommendations, rank statuses and original score objectives are unchanged. Comparing the recorder-enabled and original six-seed runs, subjects 11 and 15 have **exactly identical P300 numbers** in all six starts. Subject 21 has numerical differences up to **5.610759e-7 µV**, below the original protocol's **1e-6 µV replay tolerance**. They are **not bit-identical**; both the measured differences and the borrowed, previously fixed numerical tolerance are disclosed. This is an engineering numerical-consistency check, **not a new or relaxed scientific performance threshold**.

## Full-waveform, IC-removal and independent bootstrap findings

Reported numbers are descriptive means across **unseen ICA starts 3, 4, 5**, each compared to the corresponding FIXED recipe's seed-0 ICA fit. The full-epoch RMS difference refers to the **raw (not gain-corrected) Pz target-minus-standard mean ERP waveform**. The SME is gain-corrected bootstrap SE for the P300 mean-amplitude difference in µV. Removal Jaccard compares removed components *after* permutation/sign-invariant scalp-map matching; 1.0 means identical matched removed sets, not proven ground-truth artifact removal. Recipe ID is fixed within each row; different recipes are NOT being reselected on held-out seeds.

| Subject | Fixed recipe | ID | Mean unseen bootstrap SME (µV) | Mean waveform RMS vs seed 0 (µV) | Lowest removed-IC Jaccard |
| --- | --- | ---: | ---: | ---: | ---: |
| 11 | Original | 12 | 1.0538 | 0.00000 | 1.000 |
| 11 | Minimax / Stability-first | 25 | 1.0261 | 0.00000 | 1.000 |
| 15 | Original | 5 | 1.2009 | 0.01854 | 0.667 |
| 15 | Minimax | 2 | 1.0983 | 0.01271 | 1.000 |
| 15 | Stability-first | 6 | 1.1094 | 0.02073 | 0.667 |
| 21 | Original | 25 | 1.5532 | 0.01248 | 1.000 |
| 21 | Minimax / Stability-first | 4 | 1.5755 | 0.01248 | 1.000 |

Across those ICA fit pairs, the median absolute matched IC scalp-map correlation ranges by participant approximately: subject 11 **0.94163–0.95982**; 15 **0.99546–0.99997**; 21 **0.99974–0.99993**. Correlation is descriptive and affected by component count, scalp map alignment, and choice of spatial channels; high correlation does not certify good brain signal preservation.

### Methodological inference

- Subject **15**: Minimax has the most consistent matched removed-IC set (1.0 versus 0.667), the lowest average whole-ERP-waveform RMS change and bootstrap SME of these three compared recipes. But its **P300 amplitude range** from the previous audit, **0.014652 µV**, is worse than the original pipeline's **0.006436 µV**. Thus the scalar P300-range criterion disagrees with two distinct evidence dimensions; neither automatically establishes unbiased ERP or true ocular artifact removal.
- Subject **11**: Both original and robust recipe produce unchanged Pz difference waveform across the compared starts even though scalp ICA maps differ. Robust recipe's bootstrap SME decreases mildly; no external truth or superiority claim.
- Subject **21**: Original and robust recipe have the same within-recipe waveform RMS and stable matched removal sets; the Minimax/bootstrap SME is slightly larger than the original. This defeats any claim of across-subject Pareto improvement.

**Scientific verdict: continue development, do not promote** Minimax, stability-first or any new composite selector. The recorders are now sufficiently complete for *descriptive* spatial/temporal hypothesis generation; they do NOT establish a general causal improvement with three already-observed participants. Do not overfit a weighted objective to these exploratory outcomes, re-label the old holdout, loosen support rules or merge PR51.

## Next controlled sequence

1. Inspect individual IC-removal identities and actual matched topographies, particularly subject 15's seed4/5 difference, and confirm whether components removed correspond to independent VEOG/HEOG artifact evidence without treating proxy thresholds as ground truth.
2. Improve reporting: explicitly flag identical EEG waveforms despite IC-map drift, and separate uncertainty about the P300 mean estimate from cross-initialization stability and sample support.
3. Prospectively specify a selection rule and *decision gate* that includes bootstrap precision, full-waveform preservation, original scientific feasibility and abstention. Define source/version hashes and exact comparison/equivalence handling **before any unseen participant data**.
4. Only a candidate with meaningful, defensible development evidence and a protocol frozen beforehand may enter a **new** 14-person study; all unavailable/non-feasible outcomes count against it. Production remains the original algorithm until successful independent evaluation.

This report uses only existing cloud GitHub artifacts. There is no new 14-participant experiment or production change.
