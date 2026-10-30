# Diagnostic and scientific execution boundaries

The maintained status and run ledger is [PR50_STATUS](../docs/PR50_STATUS.md).
Original Round 8–10 plans, prototypes and tests are preserved on
[archive/pr50-prototypes-20261030](https://github.com/abwoo/PipeCompare/tree/archive/pr50-prototypes-20261030).

- Fake-candidate tests and saved-output replay diagnose implementation behavior.
- Design seeds 40101–40104 are diagnostic only.
- The micro workflow uses the production family with a finite checkpointed
  slice; reduced-family/trial overrides are not implemented.
- A completed diagnostic is not independent scientific evidence.
- Full shadow searches are explicit; normal PR updates do not restart them.
- One seed/arm may use the separate slice/resume path. Its passport binds code,
  seed, arm, variance and toolchain. A checkpoint is not a completed result.
- Cooperative budgets are checked between operations. An individual ICA fit
  can exceed the budget, so a hard runner timeout is still possible.
- Checkpointed executions do not supplement an old unsliced paired mean.
- A changed family, variance, preprocessing or ranking rule requires validation
  at its own frozen revision; replay of old decisions cannot substitute.
