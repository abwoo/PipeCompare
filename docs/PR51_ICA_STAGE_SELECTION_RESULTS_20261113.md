# PR #51 P1b/P1c scientific results: pre-removal ICA stage versus downstream P300

Date: **2026-11-13**. All development-only, real ERP CORE P3 EEG, no new subject or ICA fit for these analyses. This report is a **negative-to-mechanistic observational result**, not a new working algorithm, clinical confirmation, or evidence authorizing PR #51 merge.

## Verified provenance

- Independent full three-subject, six-start scientific observer rerun [43397978997](https://github.com/abwoo/PipeCompare/actions/runs/43397978997), all jobs **SUCCESS** including actual MATLAB 11/15/21, source identity & frozen original Rank/P300 parity, reconstruction residual and retained-trial joins.
- Prior P1b matched-component / final ERP readout [43403766941](https://github.com/abwoo/PipeCompare/actions/runs/43403766941), unit and actual archived scientific analysis **SUCCESS**; full reproducible JSON artifact `pr51-ica-stage-impact-development-43403766941`.
- P1c fixed-reference stage-signal vs trial selection/weight algebra [43404342763](https://github.com/abwoo/PipeCompare/actions/runs/43404342763), separate protocol committed **before** running, all tests and original real EEG source audit **SUCCESS**; detailed JSON artifact `pr51-ica-stage-selection-decomposition-43404342763`.
- Three **previously inspected** subjects only, 7 candidate recipe arms (11:[12,25], 15:[5,2,6], 21:[25,4]), seeds 0–5, **35** seed0-versus-seed(1–5) comparisons. Reusing the *same* subject and one underlying fitted ICA across recipe arms is **not independent biological replication**.

## Finding: the stage target-minus-standard drift is not a fixed-trial signal shift

Let `D_s` be the **actual pre-IC-removal** baseline-corrected Pz 300–600ms target-minus-standard mean over the candidate's final selected event IDs/weights. For each subsequent initialization s, we computed the same stage measurement on seed0's **frozen exact retained event IDs and weights** and proved the arithmetic identity `total stage drift = fixed-trial signal difference + selection/weight arithmetic difference`, separately checking both conditions, all source sensor event identities and original stage-sensor parity.

**In all 35 cases, the fixed-seed0-trial/weight *target-minus-standard* 300–600 ms Pz difference change was numerically zero.** This is not a claim that each ICA component, every channel or the entire sensor waveform is invariant. The nonzero pre-removal stage difference when using each seed's own final trial/weight choices was explained arithmetically by those choices.

| Development case | Pre-removal target−standard Pz difference change (μV), seeds1 / 2 | Frozen-seed0 trial/weight Pz difference change (μV), seeds1 / 2 | Trial/weight arithmetic term (μV), seeds1 / 2 | Final **different-estimand** gain-corrected P300 drift (μV), seeds1 / 2 |
|---|---|---|---|---|
| Subject 11 original candidate12 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| Subject 15 original candidate5 | +0.000593 / −0.003550 | 0 / 0 | +0.000593 / −0.003550 | −0.005836 / −0.007143 |
| Subject 15 alternative candidate2 | 0 / 0 | 0 / 0 | 0 / 0 | +0.010678 / +0.015696 |
| Subject 21 original candidate25 | **+1.495006 / +1.415300** | **0 / 0** | **+1.495006 / +1.415300** | +0.013672 / +0.101204 |
| Subject 21 alternative candidate4 | +0.088072 / +0.061678 | 0 / 0 | +0.088072 / +0.061678 | +0.021635 / −0.013567 |

All values stage-only except final endpoint in the right-hand column. Do not compare stage and final deltas as if they were the same estimand; they are different points in the processing chain, and the final is gain-calibrated. All 35 algebraic residual/conservation and original event joins passed. The decomposition uses seed0's own fixed reference; it is **not an identified causal treatment effect** of trial selection.

## Where the remaining problem must be investigated

- The original v1 all-component **minimum 0.90 topographic correlation gate still fails**. The earlier component audit found these weak maps were retained in both compared ICA fits; even when Pz stage difference cancels, their contribution and post-removal signal quality are not established safe.
- For subject 15 alternative candidate2, the *final retained event identities and weights are identical across seeds*, the pre-removal Pz difference is numerically unchanged, **yet final gain-corrected P300 changes**. This directly identifies a **post-IC-removal or downstream processing/measurement** locus for further explanation (ICA projection changes, IC removal set, filter/reference/epoch path and Rank calibration all remain possibilities). It does not isolate which.
- For subject 21 original candidate25, seeds1/2, the large pre-removal stage difference is arithmetically explained by altered final event/weight selection, while final gain-corrected differences are much smaller; do not claim 1.5 μV of true evoked neural signal error.
- Both mechanisms can be present in one participant, even if the original raw sensor events were the same. These are subject-level development observations, not grounds for replacing Standard by candidate4 or changing trial support thresholds.

## Next experiment: isolate post-removal divergence before any new recommender

Implement a **second read-only observer immediately after actual IC removal and at final epoch/weight extraction**, with source-event IDs and fixed-reference trial/weight comparisons. Check step-by-step sensor consistency (pre-removal, post-removal, final baseline, final calibrated ROI) without inferring an IC's causal effect. In parallel, compare per-fit IC removal decisions using the existing sign/permutation-invariant component matches; block results with missing/misaligned source provenance. Report whether changed ICA removal, trial rejection, weights or later filters align with final P300 differences. Preserve unchanged Standard baseline.

Do not select a new threshold or candidate based on repeatedly inspected subjects. A separately frozen prospective experimental design with a genuinely independent participant sample would be required before scientific promotion.

**Disposition:** Original ERP CORE 20-subject result **15/20 PARTIAL** unchanged; multi-metric v1 **ABSTAIN**; PR #51 **OPEN and UNMERGED**; untouched 14-subject group sealed. No production change.
