# PR51 stage 3 — pre-specified ICA-start sensitivity and robustness research plan

**Written before reviewing the new seed-specific P300 measurements.** This is an exploratory method-development design, not a revision to the original preregistered ERP CORE run. Do not merge PR #51 or apply any retrospective change to its original scientific gates.

## Measurement-level diagnosis

Initial real-P3 unstable subjects **11, 15, 21** have frozen original candidate recommendation IDs `[12,25,25]`, `[5,3,3]`, and `[25,4,4]` across ICA seed index `0,1,2`. Distinct IDs may or may not cause materially different P300 estimates.

The new per-seed audit stores:

- The selected candidate index for the current ICA start, selected on the same samples by the unchanged ranking.
- Target and standard P300 mean-amplitude estimates at Pz 300–600 ms, and their difference; per-condition SME/error when finite.
- Reference HMM-unflagged source trial counts, weighted effective HMM-unflagged n per condition; absence is explicit NaN / unavailable, **not zero**.
- Existing candidate objective values, absolute and relative P300 amplitude differences between seed0 and each alternate seed. No new arbitrary scientific pass/fail threshold will be introduced after inspection.

The initial run is: https://github.com/abwoo/PipeCompare/actions/runs/43156055411. Also require four existing MATLAB scientific test suites to pass on the isolated readout change: https://github.com/abwoo/PipeCompare/actions/runs/43156452318.

**Caution:** these scalar estimates alone do not characterize full waveform shape, removed independent component identity, or out-of-sample effects. Seed-specific recommendation IDs identify *pipeline recipes*, not necessarily identical independently fitted ICA components. Claims must be bounded accordingly.

## Prospective algorithm directions to test AFTER the diagnostic

**Control A:** unchanged Standard ranking and original single fitted ICA instance, plus existing stability audit. This is not a new algorithm.

**Candidate B — robust common-feasible path:** treat seed variation as a robustness stress test and consider only pipeline IDs feasible under every planned ICA initialization. For candidate `i`, compute a robust maximum across the seed-specific original *comparable, strictly finite* objective values:

```
J_robust(i) = max_{s in {0,1,2}} J_original(i,s)
eligible(i) = all_s original_scientific_constraints_pass(i,s)
i* = argmin_{eligible(i)} J_robust(i)
```

Ties must use the existing deterministic rank/complexity tie-breaking policy. If no common eligible candidate exists, **abstain**, not weaken constraints. This is a **research proposal**, not an implemented or accepted rule; assess whether cross-seed worst-case ranking is comparable across fits and whether the required extra computing time is acceptable before executing it. Seeds and all rule choices must be pinned before new independent outcomes.

**Candidate C — stability-aware uncertainty reporting:** retain the original chosen recipe but report the seed-resolved P300 estimate distribution and an explicit uncertainty/abstention state if the recommendation or measurement cannot be reproduced. This avoids pretending a fixed seed establishes stability; it does **not** claim to satisfy the original 18/20 recommendation consistency gate.

Do not optimize parameters to make subjects 11/15/21 satisfy the old 18/20 criterion; their outcomes are used as development evidence. Full waveform, altered ICs, induced artifact, support and independent-holdout P300 bootstrap SE must be checked for any promoted variant.

## Boundaries for stage 4

Only after a candidate method and all scientific thresholds are **frozen** should the new algorithm be run on genuinely held-out ERP CORE subjects or an external source. Independent subjects are not the same as an independent acquisition site. The original 20 subjects are development data and must never be re-labeled as a new confirmatory study. Publication claims require full missingness/abstention accounting and parallel fixed comparator processing.

## Do not confuse

- GitHub Actions success with scientific feasibility.
- Seed-ID stability with P300 measurement stability.
- Small test-retest variability with unbiased P300 estimation.
- HMM or VEOG eye proxies with ground-truth blink labels.
- The exploratory Stage2 EOG-regression failure with a general proof that no ocular correction can work.
