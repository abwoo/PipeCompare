# PR51 multi-metric ICA v1 — frozen development result and scientific stop (2026-11-12)

**Evidence:** [GitHub Actions 43193549225](https://github.com/abwoo/PipeCompare/actions/runs/43193549225), all regression/unit jobs and archived-record science evaluation **SUCCESS** at SHA `e99cc01b57bdd6e8e08be9540732d774114fbd52`; complete result artifact `pr51-ica-multimetric-development-43193549225`. First attempt [43193303017](https://github.com/abwoo/PipeCompare/actions/runs/43193303017) failed solely because the new workflow had not installed NumPy for the *existing* ICA science tests; dependencies were pinned and the same tests then passed. No EEG refit.

## Strict protocol

The selection rule was written in [PR51_ICA_MULTIMETRIC_PROTOCOL_20261112.md](PR51_ICA_MULTIMETRIC_PROTOCOL_20261112.md) **before execution**, and has not been loosened after inspecting results. Original Rank scientific feasibility, baseline-matched strata, effective physical trial n >=5, ICA component matching and classification evidence, complete waveform, P300 amplitude range and bootstrap SE must pass. Candidate selection uses only starts 0–2, with starts 3–5 for descriptive assessment. The three subjects and aggregate values had been seen in previous development work; this is **not** independent participant confirmation.

## Finding: every baseline fails the predeclared topology screen

| ERP CORE subject | Original candidate ID | Lowest matched scalp-map correlation over selection starts 1–2 | Hard screen | Decision |
|---|---:|---:|---|---|
| 11 | 12 | 0.661017 | >=0.90 | Abstain |
| 15 | 5 | 0.579231 | >=0.90 | Abstain |
| 21 | 25 | 0.799610 | >=0.90 | Abstain |

The alternative recipes in the previously declared comparison sets also fail this per-component minimum correlation screen. The low *minimum* correlations coexist with much higher *median* correlations previously reported across all five alternate starts; those are different statistics and no contradiction. A high median must not erase the scientific risk of a badly matched component.

**This is a negative method-development result, not a bug to tune away.** The unit tests enforce that the original baseline's failure causes explicit abstention; no alternative is secretly promoted, no trial floor is reduced and no unreported winner is chosen.

## Next defensible diagnostic, not a post-hoc gate change

Characterize whether the low-correlation matched components are retained or removed in one/both ICA fits; the associated ICLabel classifications and component topographies may expose a scientifically plausible failure mode. The optimal matching algorithm must be treated as a descriptive approximation, not proof of identical neural sources. Only a new, separately documented protocol can evaluate a new ranking proposal after this investigation. The frozen v1 rejection remains part of the record.

**Disposition:** Production Standard unchanged. Neither Minimax, stability-first nor multi-metric v1 is authorized as a replacement. The original 20-subject ERP CORE result remains **15/20 PARTIAL**; PR #51 remains OPEN and UNMERGED; the untouched 14 subjects remain sealed.
