# PR #51 P1b — component-resolved stage P300 vs downstream ERP association, locked analysis protocol

**2026-11-13.** This is a *retrospective, read-only secondary analysis* of the already explored ERP CORE P3 development subjects 11, 15, 21. A complete original Rank/P300-parity and pre-removal provenance capture is a prerequisite; **do not silently analyze incomplete or modified evidence**. The preceding hard-stop method `docs/PR51_ICA_MULTIMETRIC_PROTOCOL_20261112.md` remains unchanged and rejected.

## Source and analysis boundary

Only use the real cloud-only development rerecording [Actions 43397978997](https://github.com/abwoo/PipeCompare/actions/runs/43397978997), which must have archived `result.json`, `ica_science_validated.json` and `ica_contribution_validation.json` for each of subjects 11, 15, 21 with all original candidate/seed identities preserved. The immutable original six-start reference remains run 43169562528. Predeclared candidates are 11:[12,25], 15:[5,2,6], 21:[25,4]. No new real subjects, holds or ICA starts beyond seeds 0–5. Do not use any other branch/rerun or accidentally combine artifacts from attempts.

**No intervention:** only analyze the actual saved `componentUv` signature from the *pre-IC-removal stage*, plus the original full final ERP Pz waveform and source trial provenance. The sign-invariant, channel-aligned Hungarian component matching algorithm from the successful P1/P2 audit remains exactly fixed; keep **0.90** for weak scalp maps. ICA component sign changes do not alter the Pz projected contribution, provided the corresponding inverse weights/mixing are calculated together.

## Per-candidate, per-start comparisons

For each seed `s=1..5` match the real seed0 ICA mixing columns to those at seed s and compare stage-specific target-minus-standard Pz component contributions using the saved frozen **event identities and final positive trial weights**. This *stage projection* is not the final ERP counterfactual.

Store (and verify for every seed):

1. For each matched component: old/new IC index, absolute scalp-map correlation, removed/retained status in each actual pipeline, ICLabel top classes when present, and signed stage Pz contribution change in μV. Components with `r<0.90` remain explicitly labeled weak; do not use a median correlation to hide them.
2. Sum of signed matched weak and strong contribution changes, contributions from unmatched components when the ICA ranks differ, and the measured pre-removal Pz sensor/residual changes. Their sum must equal the observed stage pre-removal target-minus-standard Pz change within a fixed numerical tolerance. Do not force a zero residual in reduced-rank ICA.
3. The **absolute** contribution-change budget for retained weak components vs others, without interpreting an absolute budget as a physiological effect-size, source identity or causal attribution.
4. The final measured P300 gain-corrected difference change and full Pz waveform RMS vs seed0, extracted from original unchanged science readouts. **Do not subtract stage pre-removal μV from gain-corrected final ERP amplitudes** as though they share an estimand.
5. Whether original final **retained event identity AND original weights** are exactly comparable between starts. If either changes, label the across-start component difference as contaminated by selection/weighting differences. Report unflagged physical/effective trial counts per condition separately. No absent trial or value is treated as a zero.

No participant-level significance claims, no candidate ranking, no hypothesis testing, no reweighting after seeing answers. An exploratory case-wise association cannot establish that an individual ICA component causes final P300 change; after IC removal nonlinear subsequent operations can differ.

## Gates and disposition

The independent stage validator must confirm `sourcePz = sum(original ICA component projections)+reconstruction residual` for each actual pre-removal urevent (including non-zero rank-deficient residual) and original original-Run Rank/P300 score/strata/source trial support parity. A missing seed0 candidate is **FAIL**, not an innocent skipped observer. If the source run fails, postpone analysis and fix code without changing the scientific gate.

This analysis can identify whether low-correlation retained components have an appreciable *stage signature* and whether downstream processing/weights plausibly account for final ERP drift. It **cannot** promote Minimax, stability-first, multimetric-v1, change the >=5 effective-unflagged floor, loosen r>=0.90, run frozen 14-subject confirmation, or merge PR #51. Original 20-subject verdict **15/20 PARTIAL** is unchanged.
