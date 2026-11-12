# PR #51 stage 2 — real-data trial-support and EOG remediation verdict (2026-11-12)

**Scope:** exploratory research and diagnostic engineering only. PR #51 is OPEN and **must not be merged**. Preserve the original 20-person ERP CORE P3 verdict: 15/20 valid matched pairs, **PARTIAL**. No original scientific eligibility rule, effective-trial floor, or ranking criterion has been retrospectively relaxed.

## Verifiable GitHub Actions evidence

- Technical MATLAB exception fix and broad scientific regression (4/4 suites PASS): https://github.com/abwoo/PipeCompare/actions/runs/43148371209
- Complete stage 2 trial and eye-proxy reconciliation (14/14 Python tests and 4/4 cloud jobs PASS): https://github.com/abwoo/PipeCompare/actions/runs/43154764226
- Exploratory **fixed Minimal versus fixed Minimal + EOG regression** comparison (subject 4 and 8 jobs and collector PASS): https://github.com/abwoo/PipeCompare/actions/runs/43155328846
- Relevant source: `benchmarks/audit_erpcore_support.py`, `benchmarks/validate_veog_support.py`, `benchmarks/erpcore_pilot.m`, and separately versioned `benchmarks/pr51_eog_trial_support.m`.

All the runs above completed technically. **Success in Actions is not scientific eligibility.**

## Identified event-condition bug (diagnostic only)

The previous `hmmVeogDiagnostic` compared reference labels with literal `target`, but the Standard recorder stores `11` / `31`. This caused target trials to be matched to the standard condition or silently skipped by source-event ID. The diagnostic reporter has been repaired using explicit `11` / `31` and `target` / `standard` mapping; unknown labels now fail. No EEG data, source epochs, preprocessing choices or original P300 gate were changed by this fix.

## Subject 4

- Source reference trials: target 18, standard 95.
- HMM-marked eye artifact: target 14, standard 92. Thus only **4 and 3 unflagged source trials**, respectively, before any further rejection or weighting.
- Exploratory VEOG baseline-corrected absolute peak >50µV: **17/18** target and **94/95** standard. The 50µV proxy is *not ground truth* and must not replace the preregistered HMM criterion.
- Frozen Standard effective unflagged support: 4 / 3. Fixed Minimal and Common: 1 / 1.
- Existing EOG regression experiment: fixed Minimal 1 / 1; fixed Minimal+EOG-regression 1 / 1. Both rejected; neither yields a scientifically valid two-condition P300.
- **Decision:** insufficient support under the original HMM source flags. No threshold-selection or weighting algorithm can yield effective n≥5 from 4 and 3 unflagged original reference trials without altering the underlying label/evidence definitions. Do not rescue by loosening gates.

## Subject 8

- Source reference trials: target 22, standard 95.
- HMM-marked artifact: target 17, standard 75. Thus 5 / 20 unflagged reference trials.
- Independent baseline-corrected VEOG proxy >50µV: **21/22** target and **83/95** standard; 4 of the 5 target trials not HMM-marked still exceed 50µV.
- Frozen Standard effective unflagged support: 5 / 20; fixed Minimal and Common: 4 / 18, so the pair fails.
- Exploratory fixed Minimal: target/standard support **4 / 18**, P3 difference point estimate **2.1101598µV**; scientific status rejected.
- Fixed Minimal + pre-existing `eogregress` before rereferencing: support **0 / 7**, target P3 estimate undefined; scientific status rejected.
- **Decision:** the tested EOG-regression insertion worsens support and does not pass the fixed paired comparison. Do not promote this particular treatment. A different prospective EOG correction would require new design, a P300 shape-preservation audit and independent confirmation.

## Causal / methodological interpretation

The support shortfall here is **not** fixed by the stage-1 numerical exception repair. Retained epoch count, positive-weight epoch count, HMM-unflagged effective support and VEOG amplitude are separate quantities. VEOG is an independent but imperfect proxy for ocular artifacts. **No** evidence establishes that a >50µV event alone is a blink or that HMM flags are true labels.

Stage-2 disposition: keep insufficient-data abstention, preserve all negative cases and trial-level artifacts, do not combine the exploratory runs with the original preregistered 20-subject evaluation. Proceed with stage-3 measurement-level ICA stability, then consider any prospective new method before held-out frozen confirmation.

## Stage-3 and stage-4 boundary

The old ICA-ID consistency results for 11,15,21 were `[12,25,25]`, `[5,3,3]`, `[25,4,4]`. Stage 3 must compare actual P300 amplitude, SME and support **from each seed's newly recommended pipeline**, not only ID. Stage 4 must be a new frozen protocol on truly unseen cases, not a retrospective relabeling of the first 20.
