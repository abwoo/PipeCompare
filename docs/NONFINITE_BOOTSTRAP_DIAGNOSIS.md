# Nonfinite bootstrap diagnosis — completed

## Finding
All 101,496 invalid candidate/condition/resample rows are caused by zero weighted mass in the **unflagged-by-blink-detector comparison group**. This is an identification/support failure in those resamples. No invalid SME rows were found. The original overall-variance formula is also invalid on every one of these rows; this observed failure is not introduced by the within-group-variance denominator.

"Unflagged" does not certify clean EEG. This diagnostic does not establish whether the blink detector overflags trials or whether rejection removes useful comparison trials. Those upstream causes remain to be checked.

## Reproduction and verification
- MATLAB diagnostic execution: `3816a28cbd6922f4beccc10b830a3f7656673d7f`.
- [MATLAB run 41825137589](https://github.com/abwoo/PipeCompare/actions/runs/41825137589): successful; all ten score matrices reconstructed with identical finite masks and finite-value error below 1e-10; previous eight fallback cases reproduced exactly.
- 36 feasible candidates and 1,999 paired draws per subject/start. Saved EOG experiment outputs reused; no ICA refitting.
- [Separate Python verification 41825720716](https://github.com/abwoo/PipeCompare/actions/runs/41825720716), execution `898cf373aeb6f082d23de1272d73c45a1fb29ac9`: successful. Verified exact candidate/condition/draw/cause correspondence between JSON and all 101,496 CSV rows, distinct score-cell counts and invalid draw-column counts.
- Raw artifact: `nonfinite-diagnosis-results`, ID 13319351381, ZIP SHA256 `58971f4a6c9495353a6bf1323b30a7709a9f75cb721762f8af9089671a1b3b6d`. Includes diagnosis.json and invalid_draws.csv with counts and analytic weights.
- Verified summary artifact: `verified-nonfinite-summary`.
- Counts are repeated computational evaluations, not independent observations. The same paired-bootstrap seed/counts are reused across ICA starts.

## Results
| Subject / ICA seed | Invalid draw columns among any candidate / 1999 | Incumbent invalid draws | Challenger invalid draws | Invalid score cells |
|---|---:|---:|---:|---:|
| S004 / 0 | 0 | 0 | 0 | 0 |
| S004 / 1 | 1 | 0 | 0 | 6 |
| S004 / 2 | 0 | 0 | 0 | 0 |
| S004 / 3 | 1 | 0 | 0 | 4 |
| S004 / 4 | 1 | 0 | 0 | 6 |
| S005 / 0 | 820 | 820 | 820 | 19600 |
| S005 / 1 | 820 | 820 | 820 | 19600 |
| S005 / 2 | 820 | 820 | 820 | 19600 |
| S005 / 3 | 820 | 820 | 820 | 19600 |
| S005 / 4 | 820 | 820 | 820 | 19600 |

### S004
Starts 1, 3 and 4 each contain only one invalid draw column, affecting respectively 6, 4 and 6 other candidates. Incumbent and challenger have zero invalid scores in every start. Affected candidates retain six unflagged comparison trials in condition 2 (20 reference trials); missing all six under IID resampling has probability (14/20)^20 = 0.000797923. The blanket all-candidates finite-evidence rule turns this rare failure in other candidates into a complete fallback.

This does not show that replacement should be allowed: finite starts 0 and 2 already failed the simultaneous bound, and the expanded recommendation's actual MSE is worse at starts 0, 1 and 3. Invalid-draw handling alone is not an established fix for the S004 tradeoff or its high-error start 2.

### S005
Both incumbent and challenger retain only **one unflagged trial** in condition 1 and **three** in condition 2, at all five ICA starts.
- Condition 1: 21 reference trials, 767/1999 invalid draws. Incumbent retains 15 trials (14 flagged, 1 unflagged); challenger retains 16 (15 flagged, 1 unflagged). Probability of missing the one usable comparator: (20/21)^21 = 35.894%.
- Condition 2: 20 reference trials, 84/1999 invalid draws. Incumbent retains 13 (10 flagged, 3 unflagged); challenger 14 (11 flagged, 3 unflagged). Comparator absence probability: (17/20)^20 = 3.876%.
- Union across conditions: **820/1999 = 41.0205%** invalid draws for both recommended candidates. The two condition failure sets overlap in 31 draws.
- Across all 36 candidates: 19,600 invalid score cells and 20,296 invalid condition rows per ICA start.

The expanded output's lower known-truth MSE (0.333268 -> 0.176831) remains an observed development result, but the current resampling evidence cannot certify it. Repeating ICA starts does not create additional unflagged observations.

## Decision
Diagnosis complete. Replacement gate remains rejected; expanded-search candidate remains unpromoted. Original production defaults and PR51 code head remain unchanged and PR51 remains unmerged. No invalid draws dropped, no alpha relaxation, no threshold retuning, no new-start or ten-person follow-up launched.

## Next justified work (not executed in this diagnostic)
1. Audit reference blink flags and stage-wise retention in S005. Count flagged/unflagged support before preprocessing, after adaptive rejection, and in both selected EOG recipes. Use the injection schedule only for diagnostic validation of flags, never as ranking input.
2. Distinguish excessive flagging from loss of genuine comparison trials and from genuinely high blink burden; change only a demonstrated implementation defect.
3. A computational fallback must report insufficient comparison support. It must not impute zero eye bias, silently omit failed draws, or present forced stratification around one comparator as newly gained information.
4. If support is genuinely insufficient, evaluate a prespecified additional short window with more comparison trials; freeze selection and criteria before observing error. Any resampling/gate revision needs known-truth coverage checks before replaying these development cases.
5. Continue to reject regressions under the previously frozen paired accuracy rule; independent confirmation remains unopened.
