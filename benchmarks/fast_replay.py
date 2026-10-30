#!/usr/bin/env python3
"""Fast-loop replay guards (Round 10, diagnostic only, stdlib only).

This file does NOT run MATLAB, does NOT run a search, and does NOT
change any pick. It only classifies fragments so that a fast,
diagnostic fragment can never be mistaken for held-out evidence.

Evidence classes (exactly one per fragment):
  held-out           fastMode absent or mode 'full', labels fullyHeldOut
                     true, designSeedDiagnostic false.
  design-diagnostic  design seed 40101-40104, fastMode.mode 'micro' or
                     'slice-diagnostic', designSeedDiagnostic true.
                     Fast problem-finding only (user instruction
                     2026-10-28 23:45); never a held-out verdict.
  slice-diagnostic   checkpoint-sliced single-arm run, 40-minute slice,
                     fragment marked partialSearch when the slice
                     interrupted the search. Never merged into a
                     previous paired mean.
  sliced-validation  a completed run after a validated resume guard.
                     Still labelled sliced, never merged into a
                     previous paired mean, and never called an exact
                     resume (adaptive leaves are re-derived).
  partial            any other incomplete fragment (assessable false
                     with a reason, or partialSearch true).

require_paired_eligible() raises ValueError on any non-held-out
fragment passed to a paired computation. The held-out checker
(check_shadow_selection.py) has the same guard in MATLAB merge;
this Python guard is defense in depth, not a replacement.

calibration_report() always returns ratio AND coverage together.
A ratio without coverage is returned as incomplete, never as a
complete calibration result.
"""
from __future__ import annotations

DESIGN_SEEDS = frozenset((40101, 40102, 40103, 40104))


def fast_mode_of(frag: dict) -> dict:
    fm = frag.get("fastMode") or {}
    if not isinstance(fm, dict):
        return {}
    return fm


def evidence_class(frag: dict) -> str:
    fm = fast_mode_of(frag)
    mode = fm.get("mode") or "full"
    if fm.get("slicedValidation"):
        return "sliced-validation"
    if fm.get("partialSearch"):
        return "partial"
    if fm.get("designSeedDiagnostic") or mode in ("micro", "slice-diagnostic"):
        return "design-diagnostic" if frag.get("seed") in DESIGN_SEEDS else "slice-diagnostic"
    labels = frag.get("labels") or {}
    if labels.get("fullyHeldOut") is True:
        return "held-out"
    return "partial"


def require_paired_eligible(frags: list[dict]) -> None:
    """Raise on any fragment that is not held-out evidence."""
    for frag in frags:
        cls = evidence_class(frag)
        if cls != "held-out":
            raise ValueError(
                f"fragment seed {frag.get('seed')} is {cls}, not held-out "
                "evidence; it must never enter a paired mean or verdict."
            )


def calibration_report(records: list[dict]) -> dict:
    """records: [{assessable: bool, actual_sq: float, reported_sq: float}].

    Returns ratio AND coverage together, always. complete is True
    only when every planned record is assessable and the ratio is
    finite; otherwise the result is incomplete and the caller must
    report coverage alongside the ratio, never the ratio alone.
    """
    planned = len(records)
    ok = [r for r in records if r.get("assessable")]
    coverage = (len(ok) / planned) if planned else 0.0
    if not ok:
        return {"planned": planned, "assessable": 0, "coverage": coverage,
                "ratio": None, "complete": False}
    import math
    actual = math.sqrt(sum(r["actual_sq"] for r in ok) / len(ok))
    reported = math.sqrt(sum(r["reported_sq"] for r in ok) / len(ok))
    ratio = (actual / reported) if reported else None
    finite = ratio is not None and ratio == ratio and ratio not in (float("inf"), float("-inf"))
    return {"planned": planned, "assessable": len(ok), "coverage": coverage,
            "ratio": ratio, "complete": (len(ok) == planned and finite)}
