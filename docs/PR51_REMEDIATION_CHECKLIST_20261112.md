# PR #51 — staged remediation and acceptance criteria (2026-11-12)

**Status:** exploratory isolated branch only. Original PR #51 remains unmerged. The original 20-participant ERP CORE P3 outcome (15 valid paired analyses; decision PARTIAL) remains immutable. Changes in this branch are not a replacement for the original preregistered analysis.

## Order and stop/go decisions

1. **Execution correctness**: eliminate verified MATLAB exceptions without changing scientific eligibility requirements.
2. **Effective trial support**: determine why condition-specific clean/effective support collapses, then prospectively evaluate repairs without changing the original cutoffs after looking at results.
3. **ICA stability**: quantify robustness of the *measurements* and recommended pipelines across seeds, not just the numerical candidate IDs.
4. **New frozen validation**: freeze an updated algorithm/protocol and run a genuinely independent confirmation. A retrospective rerun on subjects used to develop the fixes remains exploratory.

### Stage 1 — execution and numerical correctness

- ERP CORE 18 and 23 original failure: `IcThreshold.choose:241` / `scoreSet` returned empty threshold and objective if the count gate passed but all selection objectives were non-finite.
- `Threshold.curve` silently turned NaN trial variance into zero via MATLAB `max(NaN,0)`. The revised computation clips negative numerical roundoff but preserves missing and non-finite evidence.
- `Threshold.choose` used to try to bootstrap a non-existent best candidate. It must instead report the non-finite objective, and must never silently count it as feasible.
- Once search returned no recommendation, `Presets.icaText` compared a nonempty candidate-index column with `[]`; this should not crash an otherwise completed search.
- Unit acceptance: non-finite-score and auto-rejection cases remain non-finite; normal finite case remains finite; absent recommendation is supported. Both paths must be covered by MATLAB GitHub Actions.
- Integration acceptance: both 18 and 23 full real-data jobs complete with `completed=true`, `errors=[]`, `notRun=0`; every remaining rejection is a *scientific* constraint, not a MATLAB execution failure. At least one previously successful control subject must be checked for regressions.
- The exploratory smoke must record the implementation SHA/diff and original immutable protocol + dataset SHA. Remove temporary exception instrumentation before final clean acceptance.
- Broad suite acceptance: `test_methodology`, `test_statistics`, `test_engine`, `test_simple` run against the revised source. Report any test failures honestly.

### Stage 2 — trial-support accounting (not a threshold-relaxation exercise)

- Subjects 4 and 8: trace *reference*, *retained*, *blink-flagged*, and *effective unflagged* target/standard trial counts separately for Standard, Minimal and Common.
- Subjects 18 and 23: all 24 original Standard paths lacked scientific feasibility; the fixed comparators retained **zero target trials**. Restoring executable code cannot by itself create a valid pair.
- For each trial retain stable source urevent index, condition, baseline, Pz/ROI metrics, VEOG/HEOG proxy measurements when available, rejection reason, selected threshold, blink mask, clean weight and effective-support contribution. Reconcile aggregates to original frozen raw results and candidate tables.
- Distinguish data loss, conservatively flagged blinks, hard-amplitude rejection, and undefined bias assessment; report sensitivity against an **external** eye-signal proxy (not a gold-standard blink label). Do not tune eye thresholds or effective-trial floor to rescue failures after observing outcomes.
- Any new preprocessing candidate must be designed before a fresh validation, be evaluated with identical fixed comparator treatment, and be audited for preserved P300 amplitude/shape and unchanged real-data source identity.

### Stage 3 — ICA recommendation stability

- Study subjects 11, 15, 21 with frozen original per-seed recommended candidate IDs `[12,25,25]`, `[5,3,3]`, `[25,4,4]`.
- Recompute per-seed final P300 target-standard mean-amplitude estimate, SME/bootstrapped SE, kept/effective support, removed ICs, and signal-loss guard; compare those measures rather than assuming that distinct candidate IDs imply scientifically important differences.
- Evaluate seed sensitivity under a fixed prespecified equivalence criterion and held-out trials. A single fixed random seed is a reproducibility control, **not** proof of algorithm robustness.
- Preserve original >=18/20 participant recommendation-consistency gate unless a wholly new protocol is approved before independent evaluation.

### Stage 4 — independent frozen P300 experiment

- Freeze code SHA, source-data hashes, toolchain versions, candidate recipes, bootstrap seeds, primary endpoint, equivalence margin, artifact/support floor, pass/fail gates, and planned missingness handling *before* outcomes are seen.
- Keep original 20 subjects and their original scientific decision as exploratory development evidence. Prefer a held-out set of ERP CORE participants (subject to availability and nonoverlap) or an independent external cohort for confirmation.
- Run every planned subject and all fixed comparators. Publish both intent-to-evaluate counts and valid-pair counts. No result is silently omitted.
- Carry forward the originally specified >=5 effective unflagged trials per condition and fixed comparator arm, >=18/20 ICA stability criterion, and 20/20 complete valid pairs **for the original protocol only**; an independently preregistered protocol must state its own immutable thresholds and sample size before collection.

## Reproducibility / claims policy

- A successful GitHub Actions *job* means computation executed. It does not mean an EEG participant met scientific gates.
- The 18/23 exploratory smoke jobs are evidence about **technical correctness**, not additional successful subject pairs.
- No edits or merges to PR #51 until the user explicitly authorizes them; maintain separate branches for new methodology variants.
- Keep failing and excluded subjects, error traces, trial-support insufficiency, and raw result identities in all analyses.
