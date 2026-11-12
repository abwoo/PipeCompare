# PR51 Stage 3C — pre-specified ICA rank-lock experiment on subject 11

**Date 2026-11-12. Prepared before any ICA re-fits under the rank-locked implementation.** This is a targeted *experimental sensitivity harness correction*, not a modification to the production Standard selector, the frozen historical ERP CORE data/protocol, or PR #51. Retain prior negative outcomes.

## Demonstrated confound from authentic source

On original development subject 11, the two previously analyzed ICA-containing recipes 12 and 25 have 17 channels × **15** components for seed 0, but 17 × **16** components for each seed 1–5. All six original fingerprints exist; these differences are observed real data, not just JSON artifacts. Zero removed ICs made the scalar P300 point estimate exactly unchanged but did not establish decomposition stability.

Both the repository's optional MATLAB `runIca` implementation and upstream EEGLAB 2026.0.0 `pop_runica` estimate PCA rank from the **first 3,000 samples** of the fitting copy. Current `ctx.icaSeed>0` alters sample order **before** this estimation; the rank estimator therefore sees a different subset, so PCA dimensionality and ICA optimization both vary. This is a confound in the existing ICA initialization sensitivity test, not evidence that default (seed-zero) production rank estimation is broken.

Upstream source: [EEGLAB pop_runica.m, version 2026.0.0](https://github.com/sccn/eeglab/blob/2026.0.0/functions/popfunc/pop_runica.m), `getrank(double(tmpdata(:,1:min(3000,size(tmpdata,2)))))` and `pca_opt` branches.

## Single-variable intervention and test

- Preserve original subject 11, original 300-second pinned recording, MATLAB R2026a, EEGLAB 2026.0.0, Standard search, P300 estimand, trial floor 5 and all ranking rules.
- Preserve **seed 0** exactly: no PCA rank-lock option added; run unmodified Standard path and verify previously observed source candidate ID12 and its original P300 estimate within 1e-8 μV (else mark failed).
- For seed indices 1..5 only, compute the **same** 3,000-column rank-estimator rule on the unshuffled prepared ICA fitting copy **before** sample permutation, pass the resulting dimension via native `'pca',r` option, then shuffle in the originally pinned seed order. This isolates sample permutation from the PCA dimension decision. The fallback path must not append a second contradictory PCA option.
- Record all six actual ICA fingerprints and statuses; verify dimension equality on the originally problematic subject11 source recipes 12 and 25 (or report failure). Validate original Rank scientific gates and effective HMM-unflagged trial support; do not equate matching rank with topographic or P300 stability.
- Retain original unmended run [43178065954](https://github.com/abwoo/PipeCompare/actions/runs/43178065954) as comparison, rather than overwriting or relabeling results.
- Run pure MATLAB numeric-rank regression for first-3000-sample subset under permutation; run original Python ICA leakage tests and latest source fingerprint tests before real EEG. Keep this on a separate development branch.

## Promotion boundary

A 15-component fit in all starts would remove one *confound* and enable scientifically meaningful cross-start topology comparisons. It would not establish a superior recommendation or prove P300 quality/accuracy. Candidate identity, ERP waveform, trial support and calibration remain separate endpoints; negative effects remain reported.

No thresholds or original 20-subject report changed (15/20, PARTIAL); no 14-person new freeze, production activation or PR #51 merge.
