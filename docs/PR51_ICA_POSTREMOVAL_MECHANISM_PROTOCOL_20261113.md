# PR51 P2b — fixed-event actual IC-removal mechanism audit: preregistered exploration

**Recorded 2026-11-13 before reviewing any actual post-removal 3-subject results**. Development cohort 11,15,21 and seed0–5 has been examined many times; this is not independent confirmation. The known original selection/weight decomposition showed that pre-removal Pz was unchanged on frozen seed0 trials but final ERP sometimes changed; this experiment tests the exact intermediate physical signal.

## Hard prerequisite and immutable sources
Use only the actual verified GitHub new-EEG run [43405859824](https://github.com/abwoo/PipeCompare/actions/runs/43405859824) *if and only if all jobs succeed* and the artifact for each of the three subjects includes original Rank/ERP parity, validated pre-removal ICA projection `ica_contribution_validation.json`, validated actual post-removal Pz `ica_postremoval_validated.json`. Use original six-seed, subject-locked `result.json` from that same run, not separate existing runs or software simulation.

No new subjects, ICA fits in this analysis, source mutation or threshold recalibration. Candidate roster 11:[12,25], 15:[5,2,6], 21:[25,4]; compare each of six starts against original seed0. Standard choice and original source trial support remain fixed.

## Event-level arithmetic in physically matched stage units

Use original independent urevent identities and actual **sensor Pz** baseline-corrected 300–600ms mean at (1) immediately before real `icremove`, (2) immediately after real `icremove`, both **uncalibrated μV**, and original final selection/positive weights `F_s`. Define the weighted target-minus-standard difference operator `D(Y,F)`.

For seed s, compute and separately report:
- **Pre removal frozen-trial signal drift:** `P_s = D(pre_s,F_0) - D(pre_0,F_0)`.
- **Post removal frozen-trial signal drift:** `R_s = D(post_s,F_0) - D(post_0,F_0)`.
- **Direct change in IC-deletion effect:** `C_s = R_s - P_s`, with independent event-level prediction from original removed-component projections: `C_s = -D(removedComponents_s,F_0)+D(removedComponents_0,F_0)`, numeric parity at 1e-4 μV.
- **Post-removal selection/weight term:** `Q_s = D(post_s,F_s) - D(post_s,F_0)`, and total stage change `T_s=D(post_s,F_s)-D(post_0,F_0)`. Require exact identity `T_s = R_s + Q_s` at 1e-6 μV.
- Final original **gain-corrected** P300 difference from the original scientific measurement and full ERP waveform RMS vs seed0. Final endpoint is a different estimand; report it alongside stage values, not subtract stage from final as though calibrated identically.

For every start and candidate retain actual removed IC indices, counts, whether *physical* event identities, retained final event identity and trial weights are identical, original support count, original Rank feasibility. Group descriptive examples by unchanged-vs-changing event/weight choices. Do not infer counterfactual neural truth from differences in cleaned signals. Repeated candidate recipes and seeds are not biologically independent observations.

## Stop conditions
Any source, observed post-removal, original Rank/P300 parity, actual IC-removal component projection, original retained urevent identity or matched time window failure terminates the scientific comparison; no imputation. The existing per-matched-ICA minimum absolute topographic correlation r>=0.90 gate remains unmodified and **failure cannot be waived**. No alteration to original >=5 effective unflagged trials or original frozen ERP CORE baseline. No promotion from these previously investigated samples. 14-person holdout remains sealed, PR #51 unmerged.

## Success criterion
A useful **methodological research advance** is an evidence-based decomposition identifying whether cross-start differences in the *same* original physical event arise **at the IC-removal step**, rather than from later selective trial/weight changes. This is not automatically a superior recommender. Only after finding a reproducible, specific locus should we predeclare an updated method on a new isolated development branch and validate using genuinely independent subjects.
