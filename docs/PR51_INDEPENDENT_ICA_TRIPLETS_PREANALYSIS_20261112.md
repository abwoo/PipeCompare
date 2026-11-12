# PR #51: pre-specified independent ICA initialization-triplet replication

**Written in the new isolated GitHub branch before any 6-start participant-level results have been read.** This is research-methodology validation on original ERP CORE P3 **development** subjects 11, 15, 21, **not** a new held-out EEG experiment.

## Research question

Does the proposed original-estimator-family-locked, original-science-feasibility-gated minimax ICA selection recommend the same pipeline from **two non-overlapping groups of ICA starts**, rather than trivially choosing the same pipeline when presented twice with the same starts?

## Exact frozen-for-this-development-audit grouping

- Group A: original ICA start **0**, alternate independent starts **1 and 2**.
- Group B: alternate independent starts **3, 4, 5**.
- Seed-group overlap: **zero**. Each group independently uses its own first seed's original `Rank.run` recommended ID as the baseline P300 estimation family (weighted vs unweighted).
- The same ERP CORE acquisition record is used by both groups. **ICA starts are independent; participants, source EEG and statistical scores are not independent**. Agreement is specifically initialization-recommendation reproducibility, not cross-site clinical generalization.
- A candidate is eligible only if under **every** seed of that group: the original scientific `Rank.run` status is `feasible`, effective HMM-unflagged target and standard support are each **at least 5**, objective values are finite and nonnegative, scientific ranking strata agree, and the weighted/unweighted P300 estimation family is unchanged from that group's own first-seed baseline.
- Within the eligible set select minimum worst objective across that group's three starts; tied minima, incomparable strata and missing data **must abstain**. Do not change any original thresholds after viewing results.

## Recorded outcomes

Exactly one of:

- `same-across-disjoint-triplets`: both groups found a scientific-feasible candidate **with the same pipeline ID**.
- `different-across-disjoint-triplets`: both groups selected science-feasible, but **different** candidate IDs.
- `at-least-one-group-abstained`: either seed group cannot legitimately select a candidate.

The first category is a necessary but not sufficient condition for a robust ICA method. It does **not** establish same fitted independent components, full waveform equivalence, trustworthy P300 magnitude, or a statistical group-level improvement. Record all corrected P300 estimates, SME/effective HMM support and implementation identity. All three original development participants must be reported, even if negative; no participant dropping, no replacing failures with NaN=0, and no promoting group-to-group metric with zero independent test data.

## Implementation and reproducibility

The ordinary simple-mode `icaSeeds=2` remains unchanged. Only GitHub Actions exploratory branch `codex/pr51-independent-ica-triplets-20261112` opts into six-start audit with environment `PR51_DEV_SIX_SEED_AUDIT=1`, which raises the two additional ICA reruns from 2 to 5 **on that development branch only**. The original scientific ranking and fixed Minimal/Common pipelines are not changed.

- The independently validated algorithms are `+pipecompare/+eval/IcaRobust.m`, `IcaRobustFamily.m`, `IcaIndependentGroups.m`.
- 28/28 original and new MATLAB science/missingness tests PASS: https://github.com/abwoo/PipeCompare/actions/runs/43167448188
- 15/15 independent Python auditor tests PASS: https://github.com/abwoo/PipeCompare/actions/runs/43167764200
- Real-data six-ICA-start development workflow: https://github.com/abwoo/PipeCompare/actions/runs/43167949961. **No scientific success is declared merely because the workflow has started or GitHub Actions finishes.**

## Genuine new frozen study remains separate

ERP CORE's 14 predeclared unseen subjects [24,25,26,27,28,29,31,32,33,34,36,37,38,39] remain unexamined. All 42 public OSF source file metadata hashes are pinned at `benchmarks/results/holdout14_osf_candidate_manifest_20261112.json` on distinct `codex/pr51-holdout14-design-20261112`; metadata hash is not yet an actual file-byte SHA verification. Independent 14-person scientific confirmation **must not start** until the precise method version, all loss/abstention rules and the independent ICA initialization comparison criterion are fully frozen and source file checks pass. This workflow is not that frozen experiment.

**PR #51 remains OPEN and unmerged, original 20-person result 15/20 valid matched pairs, PARTIAL.**
