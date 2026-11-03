# Ten short public EEG records: frozen exploratory pilot

This PR evaluates main algorithm 7151821edf8a5a54c3167d7f18b66c530911b568.
PR #50 remains closed and none of its production changes or Balanced recipe is imported.

## Frozen source and computation

- PhysioNet EEG Motor Movement/Imagery Dataset 1.0.0, participants 1–10, R01, first 60 seconds each.
- Ten distinct participants, 10 minutes total, 160 Hz; the same 16 declared scalp channels for every source.
- Each complete source EDF is SHA-256 verified against the frozen official manifest before extraction.
- MATLAB R2026a and EEGLAB 2026.0.0 run on GitHub Actions. Native ICLabel is verified against its original MATLAB inference; the compiler compatibility changes do not alter model weights or inference mathematics.
- An unchanged resting alpha filter search checks execution on unmodified real EEG.
- The whole main Standard procedure is freshly fitted per source: bad-channel handling, ICA, automatic IC/rejection choices, searched filters, trial weights, cutoff/EOG checks, and two extra ICA starts. Candidate/condition rows never count as independent records.
- The source channels, source digest, source excerpt, frozen protocol digest, production revision and execution revision are retained separately.
- No raw EEG, ICA matrices or derived per-trial subject data are committed.

## Known-target comparison

Resting EEG does not supply P3 ground truth. A separate, declared semisynthetic copy uses a known P3 field, randomized event schedule/amplitude/latency, and injected blink field with declared prevalence. The natural background is centered, randomly signed/circularly shifted and scaled by a fixed rule. It is not certified clean neural activity. Genuine neural amplitude covaries with blinking in the generator and remains part of the target.

The population target includes amplitude and latency variation; the realized clean-trial target is reported separately. A single P3 mean at Pz, 300–600 ms, is scored. The canonical reference is the average of the 16 scalp sites, excluding EOG. Raw error and gain-corrected error are both retained.

Standard is compared with two prespecified controls on the exact same recording:
1. Minimal: 0.1–30 Hz filtering, average reference, epoch/baseline, 100 uV rejection.
2. Common: the same with fixed bad-channel handling, ICA and Eye/Muscle ICLabel threshold 0.8.

No candidate or comparator is chosen by true error. Truth is used only after fitting to calculate outcome error and a descriptive oracle regret within the fitted feasible family. That regret is not an independent selection evaluation.

Native decision replay and the generated executable script are both checked against the stored selected estimate at 1e-6 uV. Known channel-fault durations 1/5/10/20 s are a separate detector diagnostic. Marking outside an injected interval on natural EEG is not an independently certified false-positive rate.

## Prespecified interpretation

All ten planned inputs must be present and execution-complete, with at least eight assessable Standard recommendations, for the pilot usability screen. Missing recommendations and execution failures are visible for all planned participants. Green computation is not a calibration or accuracy conclusion.

The primary comparison is the participant-mean squared gain-corrected population-target error. Paired differences are Standard minus each fixed control, with 10,000 participant-bootstrap draws, seed 510003. A 95% upper bound below zero supports improvement on this pilot; an interval containing zero is not-established.

Calibration is the RMS actual error divided by RMS reported error, after averaging condition squares within each participant. The [0.8,1.2] interval is a prespecified preliminary agreement band. Interval overlap is not-established; an interval entirely outside is adverse. The existing floor of 30 assessable independent records and relative half-width <=10% remains unmet by construction. Ten records cannot establish stable calibration, clinical performance, universal EEG improvement, or unseen-participant group-prior transfer.

The Python collector independently recomputes errors and MSEs, validates provenance and unchanged Standard options, rejects mixed/duplicate records, and retains all bootstrap mass, including undefined/zero-denominator calibration draws.

## Reproduction

The workflow public-pilot.yml runs the new branch on pushes affecting this benchmark, or can be dispatched manually. It verifies the new parser/checker contracts, runs participant 1 before the remaining nine, retains a separate artifact per participant, and collects coverage even after failures. No broad historical compatibility matrix is needed for this benchmark-only change.

Commands for a single source after installing the declared dependencies:
- python benchmarks/fetch_public_pilot.py --subject 1 --out input
- In MATLAB with EEGLAB initialized: addpath('benchmarks'); public_pilot('input','result',1)
- python benchmarks/summarize_public_pilot.py fragments pilot_summary

Source: Schalk, G. (2009), EEG Motor Movement/Imagery Dataset (1.0.0), PhysioNet, https://doi.org/10.13026/C28G6P, Open Data Commons Attribution License 1.0. Original acquisition: Schalk et al. (2004), BCI2000, IEEE Transactions on Biomedical Engineering 51(6), 1034–1043.

The frozen input/metric manifest is benchmarks/public_pilot_protocol.json. Executed outcomes, cloud run links, adverse findings and participant-level records are in [the completed results](PUBLIC_PILOT_RESULTS.md). No confirmation cohort is opened and no merge or release is performed.

## Cloud execution

Frozen benchmark execution: [Actions run 41604029052](https://github.com/abwoo/PipeCompare/actions/runs/41604029052), code revision cca3b44e94ee5c4779f3f4d4f6ed4c3c6ec41b1e. The run completed successfully; [results and limitations](PUBLIC_PILOT_RESULTS.md) separate execution from scientific conclusions. Later result/documentation revisions must keep this execution revision distinct from the production algorithm revision.
