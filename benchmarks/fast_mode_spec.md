# Round 10 fast mode — specification (diagnostic only, never evidence)

Baseline: head `6035b50ec32a500236fa45e8c5b6d743b3c6e01d`, branch
`codex/methodology-audit-20261023`. This package is local only:
no push, no merge, no main, no MATLAB run locally. CI is the first
MATLAB verification of every MATLAB change in this package.

## The loop the user ordered (2026-10-28 23:42 and 23:45)

- Cancel a shadow that is destined to time out; do not wait it out.
  Change immediately, push immediately, then run fast, screen
  immediately, fix immediately, and loop again. There is no day/night
  distinction in this loop; a finished fast run is screened at once,
  not held for the morning.
- Every fast run must fit in 1 hour. A fast run that needs more
  than 1 hour is not a fast run; it is a sliced run that must stop
  at its slice, upload, and be screened as partial.
- Fast tests and micro searches use the existing 4000-series sample
  data (design seeds 40101-40104) to FIND problems. A design seed
  is never held-out evidence: a result on a design seed must never
  be reported as a held-out verdict, never merged into a paired
  mean, and never used to self-certify a fix. The final conclusion
  still requires held-out evidence (20031+, frozen criteria).

## Three layers, all diagnostic, none is full evidence

### Layer 1 — seconds: fake-candidate tests (matlab-tests, < 1 hour)

No search, no EEG. Fake candidates, fake fragments, fake
calibration records. This layer proves logic on known inputs:
contractMatch on struct arrays, table column names from
`makeValidName`, guard notes, bestSet flags, crossfit finite
guards, eyeBias representation, checker guards (filename seed,
duplicate seed-arm, planned vs present, CSV pick/selectedKey),
fast-config validation, resume-guard rejections, variant B gate
closed by default. A Layer 1 pass is only a ticket to Layer 2.

### Layer 2 — minutes: offline replay (pure Python stdlib, < 10 minutes)

Inputs are already-fitted families and already-completed fragments
only. Outputs are per-seed old_pick / new_pick / changed / reason /
flag, with flags restricted to `bestSetUndefined /
guardUnassessed / eyeBiasNotApplicable / calibrationCoverage /
timeoutSliced / partialSearch / none`. Replayable: ranking-rule
changes with the family unchanged (guard, bestSet visibility,
eyeBias representation, calibration reporting, report columns).
NOT replayable, and a "no change" replay counts for nothing:
family changed, ICA changed, variance definition changed, data
changed, guard pass-through changed (variant B must refit;
replay on a fitted family is variant A only).
Calibration output always carries ratio AND coverage together;
a ratio without coverage is incomplete, never a result.
`benchmarks/fast_replay.py` classifies every fragment and
`require_paired_eligible()` refuses any non-held-out fragment
for a paired computation. The MATLAB merge and the Python
checker carry the same refusal independently.

### Layer 3 — 1 hour: micro search, slice-diagnostic (one CI job < 60 min)

Frozen recipe (the config is frozen; the execution gap below is
part of the recipe, not a footnote to hide):
- 1 seed, design seed 40101 by default (40102-40104 for later
  fast iterations only), 1 stratum (one bench_synth record),
  1 arm (Balanced, the slow arm, by default), production family.
- Checkpoint non-empty, finite maxSeconds <= 2400, sliceSeconds
  == 2400 (40 minutes). Empty checkpoint + Inf time is forbidden.
- Job budget: setup 10-15 min + search slice 30-40 min +
  upload/collect 5 min = < 60 min. At 40 minutes the job slices,
  uploads `manifest.mat`, `root.mat`, `leaf_*.mat`, `ica_cache/*`,
  and writes a fragment marked `partialSearch` with the reason.
  A sliced fragment is a partial diagnostic, not a completion.
- Micro picks are NEVER compared by MSE with full-family picks.
  A micro run that reproduces a bug is reported as "reproduced
  under the micro recipe", with the recipe difference stated.
- Three tables are mandatory, missing one fails the run:
  per-step timing (median/max/count per step type + slowest-5
  candidates; wall clock only from tic/toc, never summed across
  candidates), failure-visibility (per arm assessable / reason /
  flag), coverage (planned vs present vs assessable vs sliced;
  planned != present must be shown, never dressed as planned).

### Recipe deviation, stated at the top (do not push past this)

The frozen micro recipe ALSO names a reduced family, a reduced
trial subset, `icaSeeds=1`, `refineCutoffs` off, `eogCheck` off.
That part is NOT wired in this package, and this package must
not be pushed as if it were: at head 6035b50, pop_pipecompare's
documented option list does not include icaSeeds, refineCutoffs
or eogCheck, and Presets.searchOptions() hardcodes
refineCutoffs=true, icaSeeds=2, eogCheck=true. Threading those
overrides through pop_pipecompare is a production change that
needs its own design, pre-registration of the reduction rule,
and MATLAB proof. Until it lands, the Layer 3 job in this
package is a slice-diagnostic run on the production family,
labelled `slice-diagnostic`, and `validate_fast_config`
returns `executionReady=false` for the reduced recipe.
Anyone who describes a slice-diagnostic fragment as a
reduced-family micro result is mislabelling it.

## Timeout no longer means evidence lost (slicing)

- Every slice uploads its checkpoint artifacts and writes a
  readable fragment (even a partial one) with `partialSearch` +
  completed arms + unfinished reason. Waiting for a whole seed
  to finish before writing anything is forbidden in fast mode.
- Resume guard (strict, stricter than Executor.resume itself):
  versions pinned (manifest MATLAB/EEGLAB must equal the current
  run exactly; Executor.resume only warns, this guard errors),
  artifactVariance recorded explicitly and asserted (never
  defaulted from an old checkpoint), adaptive-leaf re-derivation
  logged. Even a passing guard is NOT an exact resume:
  Executor.resume re-derives adaptive leaves beyond the manifest.
- Sliced results are labelled `sliced validation runs`. They
  are never merged into a previous paired mean, never counted
  in a previous coverage, and never called completed full runs.
  A completed sliced validation run is a new, separately
  labelled run, judged on its own.

## Variant B: validation first, arm second

Variant B is NOT an arm in this package. The gate
(`benchmarks/variant_b_gate.m`) is closed by default and only
opens on a recorded 1-seed guard on/off validation via the
bypass-pop_pipecompare-optimize path (M3) or the production
pass-through path (M1), with plan-construction reuse verified
in MATLAB and an explanation when on/off picks agree.
Three runs of pop_pipecompare, and a variant A replay
relabelled as variant B, are explicitly rejected by the gate:
pop_pipecompare has no protectReferenceCleanup option and the
first Rank.run defaults the guard to true, so three runs are
identical and their guarded-minus-unguarded difference is
identically zero — false evidence. Production Standard/Balanced
definitions are never swapped, renamed, or redefined to make
variant B look good.

## Visibility items, in priority order (report-only unless marked)

1. bestSet incomplete: REPORT-ONLY in this package. Rank
   already flags `bestSetUndefined` and notes it in the table
   (Round 8); this package adds an assertion test that the
   flag is set and the pick is unchanged. The strict
   keep-only-best variant is a prototype, NOT wired into
   Rank.run: wiring it WOULD change picks and therefore needs
   its own held-out revalidation, separately labelled.
2. Calibration coverage: REPORT-ONLY. `Calibration.summarize`
   already returns ratio + coverage + plannedRecords (Round 8);
   this package adds an assertion test that a ratio is never
   reported without coverage, and that a stratum collapse
   cannot hide inside a pooled ratio (summarizeByStratum).
   Changing the dropping rule itself would change results and
   is NOT in this package.
3. Hard-data failure visibility: REPORT-ONLY in this package.
   Every fast fragment must carry assessable / reason / flag;
   a silent normal-looking pick on incomplete data is a
   failure of the run, not a result. Changing what is
   assessable (the gate, the dropping rule) would change
   picks and is NOT in this package.
4. Diffusion / calibration-gate report wiring into shadow JSON:
   NOT in this package (helpers exist, wiring is unfinished).
   Report columns are diagnostic only until wired and proven.
5. Successive halving / surrogate screening / CI cache /
   matlab-tests suite split: NOT in this package. Halving and
   surrogate need frozen rules + a named arm + held-out
   pairing; cache and split were measured as not worthwhile
   (setup ~89s vs searches of hours; split setup ~693s total).

## What counts as "problem found", and what must return to full held-out

Problem found = at least 2 of 3 agree: (1) fake test reproduces
and passes after the fix; (2) offline replay changes a pick with
the reason landing on a named flag; (3) micro search reproduces
with timing/visibility tables locating the step and arm.
One item alone is a lead, not a finding.
Must return to full held-out (fast loop can never substitute):
final proof, a new arm, family/ICA/variance/data changed,
calibration-gate upgrade, variant B arm creation. Fast results
are screening and diagnosis; held-out paired evidence is proof.

## Continuous loop (no waiting for a user prompt)

After a fast run finishes (complete or sliced), screen it
immediately: assessable/reason, fingerprint, pick/selectedKey,
finite MSE, provenance commit == current head, timing present.
If a problem is found, fix it in the next local package at
once and push it as soon as the previous fast job is final —
a push cancels a pull_request shadow, so a fast push must
never be made while a fast job on the same branch is still
in_progress, unless the user explicitly orders that
cancellation. The loop does not pause for morning, for a
user message, or for a "no change" summary: no change still
gets screened, and the next fix is prepared in parallel.
