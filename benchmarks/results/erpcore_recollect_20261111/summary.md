# ERP CORE real-P3 validation — collector summary

Decision: **partial**

Records: 15 valid pairs of 20 planned.
Missing: none; problems: 0.
Excluded pairs (support < 5): [8].
Mixed execution identity: False.

Incomplete pairs (Amendment A5, listed not dropped): 4 of 20 planned.
  subject 4: minimal arm unavailable: fixed readout not finite.
  subject 17: Standard arm unavailable: no feasible pipeline (0 of 24 candidates).
  subject 18: Standard arm unavailable: no feasible pipeline (0 of 24 candidates).
  subject 23: Standard arm unavailable: no feasible pipeline (0 of 24 candidates).

Precision Standard vs minimal: mean diff -0.597 uV, 95% CI [-0.941, -0.342], wins 14, losses 1, individual regressions 0.
Precision Standard vs common: mean diff -0.392 uV, 95% CI [-0.583, -0.219], wins 14, losses 1, individual regressions 0.
Calibration (Standard, reported vs bootstrap SE): ratio 0.973, 95% CI [0.956, 0.988], band [0.8, 1.2].
ICA-seed stability: 13/16 same recommendation.
Reference agreement (Standard vs published pipeline): mean |diff| 3.892 uV, correlation 0.760.

Gates: {"precision_vs_minimal": true, "precision_vs_common": true, "calibration": true, "stability": false, "individual_regressions": true, "complete_pairs": false}
