# PR #50: scope, current behavior and evidence

This is an implementation/evidence index. It does not establish publication
readiness or independent whole-procedure confirmation.

## Scope

Production methodology corrections and reproducible validation infrastructure
remain in PR #50. The consolidation changes no Rank, Measure, HMM, spectral or
empirical-Bayes algorithm. Its baseline is
`6b26931619b8b11f630326638cdf2a4cb478e0fa`.

Eight disconnected experimental helpers, their three unit tests and four completed
one-off audit workflows are preserved
on [archive/pr50-prototypes-20261030](https://github.com/abwoo/PipeCompare/tree/archive/pr50-prototypes-20261030).
That branch points to the complete original baseline, so no source, test or
historical result is lost. They are absent from the delivery branch.
This includes halving, surrogate screening, group-relative scores, diffusion and
calibration report prototypes, the strict best-set prototype and variant-B design.
The archived one-off runs stay accessible through their original run links.
The variant-B admission check remains a diagnostic helper, not a validated arm.

## Current recommendation rule

Introduced at `d908041ff5aa5f8c62cf111e5b75eae397d2f950`:

| Evidence under the fitted score | Decision |
|---|---|
| Baseline significantly better: baseline - proposal + critical < 0 | Use baseline |
| Indistinguishable scores | Keep proposal; maintain |
| Missing usable comparison | Keep proposal; unassessed |
| Proposal already is baseline | Keep baseline |
| Reference repair comparator ineligible or absent | No repair substitution |

The automatic average-reference comparator and explicit upgrade baseline share
this rule. It controls a conditional fitted score, not unknown true MSE.
An explicit baseline takes precedence. This consolidation does not reverse the
rule again. The earlier b27b435 guard instead retained the baseline when a proposed
upgrade was uncertain; evidence for it must not be attributed to the current rule.

## Evidence ledger

| Algorithm revision | Evidence | Interpretation |
|---|---|---|
| 4116ea11f3d5d39e4d3d282894dcf487172bd50d | [40282756637](https://github.com/abwoo/PipeCompare/actions/runs/40282756637) | Four development pairs; Standard MSE 0.487085, Balanced 0.749338. Failed the original engineering admission tolerance. |
| b27b4351846a8f53b7472e80529b2f86cf8a60cd | [40458110149](https://github.com/abwoo/PipeCompare/actions/runs/40458110149) | Earlier guard, cached-family replay only: 40102 worsened; 40103 1.352610 to 0.218319. Does not validate current d908041 rule. |
| 87f8e2680f1af817cfecad1456c4c6881c9540dc | [40460170075](https://github.com/abwoo/PipeCompare/actions/runs/40460170075) | Seven fixed-operator source cases; no ICA fitting. Explains residual blink/movement and cancellation, not current algorithm accuracy. |
| 6b26931619b8b11f630326638cdf2a4cb478e0fa | [40805638630](https://github.com/abwoo/PipeCompare/actions/runs/40805638630) | Numerical tests: 180 passed, 0 failed, 2 skipped. |
| 6b26931619b8b11f630326638cdf2a4cb478e0fa | [40805638787](https://github.com/abwoo/PipeCompare/actions/runs/40805638787) | Native/GUI compatibility: 177 passed, 0 failed, 2 skipped. |
| 6b26931619b8b11f630326638cdf2a4cb478e0fa | [40805638604](https://github.com/abwoo/PipeCompare/actions/runs/40805638604) | Eight Standard fragments, one Balanced fragment; one complete pair. Six Balanced jobs timed out; one failed without an available job log. The old checker crashed on singleton seeds. |

The arithmetic blend 0.495976 included three old replay values plus one
carried-forward error. It was never a complete current-version benchmark.
One-pair bootstrap intervals are degenerate: a criterion label from that
calculation is not evidence of stable calibration or population generalization.
A summary workflow's green check may validate historical artifacts only.

## Execution and recovery

- Full `shadow-selection.yml` is manual only and never cancels a running job.
- `pr50-consolidation.yml` checks the affected contracts and re-aggregates the
  nine artifacts from 40805638604. Original algorithm SHA and source run stay
  attached; processing SHA is recorded separately. No EEG is fitted for aggregation.
- `shadow-checkpoint.yml` runs a single explicit seed/arm with a finite budget,
  uploads the checkpoint even after interruption, and accepts an artifact from
  the exact same code/toolchain/seed/arm for continuation.
- The passport, root/manifest identity, completed leaves and ICA random state
  are checked before reuse. A completed retry reads its saved result.
- Checkpointed fragments are a separate execution, rejected by the old unsliced
  aggregator. Adaptive decisions can be re-derived on resume; exact full-process
  parity and independent confirmation are not certified by infrastructure tests.
- Old failed jobs that never wrote checkpoints cannot be resumed from final JSON.

## Remaining work, in order

1. Consolidation checks and saved-output verification are complete; see below.
2. Use the tested checkpoint path for missing work in an explicitly labelled
   execution; keep completed old fragments and failed coverage visible. Reuse old
   controls only when input, contract and backend provenance match; never silently
   combine different execution protocols.
3. Before changing a recommendation rule, freeze its objective, comparator and
   whole-procedure validation protocol. No new default follows from a prototype.
4. Independent cohort calibration and full-process fit/validation separation
   remain unresolved scientific tasks. Confirmation remains unopened.

## Completed consolidation verification

- [41052090671](https://github.com/abwoo/PipeCompare/actions/runs/41052090671):
  101 affected cases passed. One new native fixture failed and was corrected;
  that case is superseded by the proof below. The overall run is red, and is
  not represented as a successful run. Already passing cases were not repeated.
- [41053098369 native proof step](https://github.com/abwoo/PipeCompare/actions/runs/41053098369/job/124053760559):
  the remaining native checkpoint test passed and the original-fragment MATLAB
  merge completed. The overall run later failed at Python empty-pick comparison;
  this is superseded by the saved-output check below.
- [41053921374](https://github.com/abwoo/PipeCompare/actions/runs/41053921374):
  13 Python cases passed. Independent verification reports 8 seeds and 10
  assessable arm records; all paired values, intervals, verdicts and CSVs match.
  It reused completed MATLAB outputs and started no MATLAB or ICA.
- Across these scoped executions, 102 MATLAB cases have passed. The two
  remaining native guard tests passed in the first set; the interrupted/reused
  search proof passed in the later set. This is not a new full compatibility run.
- All 42 production implementation blobs match the 6b26931 baseline. The
  checkpoint helper and native proof test match their tested 298744a blobs.
- [Machine-readable provenance](../benchmarks/results/pr50_consolidation.json)
  records source, merge and checking revisions separately. The final integration
  contains documentation, workflow registration and archiving changes and uses
  `[skip ci]` to avoid restarting completed validation.

The native proof checks interrupted/uninterrupted recommendation and score
agreement, byte-preserved completed leaves and ICA files, cached completion
retry, and rejection of mismatched identities. The manual cross-worker workflow
is prepared but a full scientific run through it has not been completed.
No missing full-record experiment or independent confirmation was launched
during consolidation.
