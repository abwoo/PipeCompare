# Frozen held-out scoring comparison

Prospective held-out participant comparison for the PR51 scoring variant, using real resting EEG plus known injected P3/blinks. Same-source generalization only; not real ERP-task, clinical, or independently executed validation.

Decision: **inconclusive-no-promotion**. Valid pairs: 10/10.

Mean comparison and prespecified gates:
```json
{
  "meanComparison": {
    "baselineMSE": 1.1823481131639453,
    "variantMSE": 1.1823481131639453,
    "difference": 0.0,
    "ci95": [
      0.0,
      0.0
    ],
    "relativeReduction": 0.0,
    "baselineRawMSE": 1.2240953395504908,
    "variantRawMSE": 1.2240953395504908,
    "rawDifferenceCI95": [
      0.0,
      0.0
    ],
    "wins": 0,
    "losses": 0,
    "ties": 10
  },
  "gates": {
    "meanReductionAtLeastFivePercent": false,
    "pairedCIEntirelyBelowZero": false,
    "rawMeanNotWorse": true,
    "winsExceedLosses": false,
    "noMaterialIndividualRegression": true,
    "variantUncertaintySupported": false
  }
}
```

| Participant | Baseline MSE | Variant MSE | Difference | Individual guardrail | Variant uncertainty |
|---|---:|---:|---:|---|---|
| 11 | 1.02527 | 1.02527 | +0 | True | True |
| 12 | 1.66207 | 1.66207 | +0 | True | True |
| 13 | 0.585963 | 0.585963 | +0 | True | False |
| 14 | 0.866039 | 0.866039 | +0 | True | False |
| 15 | 0.348138 | 0.348138 | +0 | True | False |
| 16 | 0.455505 | 0.455505 | +0 | True | False |
| 17 | 1.36965 | 1.36965 | +0 | True | True |
| 18 | 2.9424 | 2.9424 | +0 | True | False |
| 19 | 2.3838 | 2.3838 | +0 | True | True |
| 20 | 0.184644 | 0.184644 | +0 | True | False |

All ten planned subjects, including missing/invalid pairs:
```json
[
  {
    "subject": 11,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-11-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 12,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-12-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 13,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-13-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 14,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-14-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 15,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-15-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 16,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-16-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 17,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-17-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 18,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-18-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 19,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-19-attempt-1"
    ],
    "error": null
  },
  {
    "subject": 20,
    "present": true,
    "validPair": true,
    "attempts": [
      "heldout-pair-20-attempt-1"
    ],
    "error": null
  }
]
```

Point accuracy, uncertainty support and engineering completion are separate outcomes. Ten same-source recordings cannot establish stable calibration or real-task superiority. No production change or merge. See summary.json for retention, support, calibration denominators, fixed comparators and reasons.
