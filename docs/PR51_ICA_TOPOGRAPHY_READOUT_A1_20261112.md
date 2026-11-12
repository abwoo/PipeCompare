# ICA topology audit A1 — singleton component indices in MATLAB JSON

Written while raw six-start ICA topology workflow [43178065954](https://github.com/abwoo/PipeCompare/actions/runs/43178065954) is still running and before reading newly fitted component maps.

The source generator legitimately stores a scalar removed component index when one IC is removed. MATLAB `jsonencode` serializes that 1×1 numeric as a JSON number, not an array. The first Python topology comparison accepted only list-shaped removal identities and would mark this scientifically valid singleton as missing. The A1 parser accepts a finite positive integral JSON number as a one-item set, while preserving original fail-closed checks for invalid indices, duplication, absent evidence and component identity. It does not alter ICA fits, topology matching, scientific metrics, thresholds, selected candidates, or historic results.

Regression: `tests/test_ica_topography.py::test_matlab_singleton_removed_indices_are_not_missing`. The completed historical raw result will be re-audited via source SHA and workflow run identity. **No EEG refit is required for this JSON correction.** Neither a source data success nor CI pass changes the original PR51 PARTIAL 15/20 or authorizes the 14-subject new frozen study.
