#!/usr/bin/env python3
"""Fail-closed audit of the scientific feasibility of robust ICA candidates.

NOT confirmatory. The original recommendation remains unchanged, and the
effective *unflagged* per-condition trial threshold is fixed at 5.
"""
import argparse
import json
import math
from pathlib import Path

FROZEN_SUPPORT_FLOOR = 5.0


def finite(value):
    return isinstance(value, (float, int)) and not isinstance(value, bool) and math.isfinite(value)


def validate(r):
    if not r.get("completed") or r.get("errors") or r.get("notRun") != 0:
        raise ValueError("incomplete real-P3 run, cannot report ICA results")
    s = r.get("icaSeedCheck")
    if not isinstance(s, dict):
        raise ValueError("missing ICA check")
    ids, obj = s["ids"], s["objective"]
    statuses, mask, strata = s["statusBySeed"], s["supportPassBySeed"], s["stratumBySeed"]
    robust = s["robust"]
    if not (len(ids) == len(obj) == len(statuses) == len(mask) == len(strata)):
        raise ValueError("ICA candidate axes differ")
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate candidate IDs")
    nseed = len(s["recommended"])
    if nseed != 3 or len(s["measurements"]) != nseed:
        raise ValueError("expected original plus two ICA initialization seeds")
    all_eligible = []
    for idx, objective, status, support, group in zip(ids, obj, statuses, mask, strata):
        if any(len(x) != nseed for x in (objective, status, support, group)):
            raise ValueError(f"misaligned seed axes for candidate {idx}")
        is_good = (all(x == "feasible" for x in status)
                   and all(x is True for x in support)
                   and all(finite(v) and v >= 0 for v in objective)
                   and len(set(group)) == 1)
        if is_good:
            all_eligible.append((idx, max(objective), group[0]))
    groups = {v[2] for v in all_eligible}
    assert robust["originalRecommendationUnchanged"] is True
    if robust["status"] == "selected":
        if not all_eligible:
            raise ValueError("selected candidate fails at least one scientific gate")
        if len(groups) != 1:
            raise ValueError("selected candidate across incomparable scientific strata")
        if robust.get("allSeedsScientificallyFeasible") is not True:
            raise ValueError("selector falsely claimed all-seed scientific feasibility")
        if robust.get("allSeedsSupportFloorPassed") is not True:
            raise ValueError("selector falsely claimed effective-unflagged support")
        candidate = robust.get("candidateId")
        selected = [v for v in all_eligible if v[0] == candidate]
        if len(selected) != 1:
            raise ValueError("selected candidate fails at least one scientific gate")
        best = min(v[1] for v in all_eligible)
        winners = [v for v in all_eligible if abs(v[1]-best) <= 1e-12*max(1, abs(best))]
        if len(winners) != 1 or candidate != winners[0][0]:
            raise ValueError("candidate not unique minimax among feasible candidates")
        if not math.isclose(robust["worstObjective"], best, rel_tol=0, abs_tol=1e-9):
            raise ValueError("reported worst objective differs from candidate matrix")
        measurements = robust.get("seedMeasurements", [])
        if len(measurements) != nseed:
            raise ValueError("robust candidate measurements absent")
        for seed, m in enumerate(measurements):
            if m["seed"] != seed or m["candidateId"] != candidate:
                raise ValueError("chosen recipe not re-evaluated at each ICA seed")
            eff = m["effectiveUnflagged"]
            if len(eff) != 2 or not all(finite(v) and v >= FROZEN_SUPPORT_FLOOR for v in eff):
                raise ValueError("scientific trial support below fixed floor")
            raw = m["differenceEstimate"]
            corrected = m["correctedDifferenceEstimate"]
            gain = m["gain"]
            if not (finite(raw) and finite(corrected) and finite(gain) and gain > 0
                    and math.isclose(raw/gain,corrected,rel_tol=0,abs_tol=1e-8)):
                raise ValueError("invalid gain-corrected per-seed P300")
    elif robust["status"] == "abstain":
        if robust.get("candidateId") is not None:
            raise ValueError("abstention cannot expose a selected candidate")
        if robust.get("seedMeasurements"):
            raise ValueError("abstention cannot invent scientific measurements")
        if not all_eligible:
            expected = "no-common-feasible-candidate"
        elif len(groups) != 1:
            expected = "multiple-feasible-strata-cannot-be-ranked-together"
        else:
            best=min(v[1] for v in all_eligible)
            ties=[v for v in all_eligible if abs(v[1]-best) <= 1e-12*max(1,abs(best))]
            expected=("indistinguishable-minimax-tie-requires-predeclared-tiebreak"
                      if len(ties)>1 else None)
        if robust["reason"] != expected:
            raise ValueError("abstained despite unique scientifically feasible candidate")
    else:
        raise ValueError("unrecognized scientific robust selector state")
    return {
        "subject": r["subject"], "selectedOriginalIds": s["recommended"],
        "originalSame": s["same"], "scientificFeasibleAcrossSeeds": len(all_eligible),
        "robustState": robust["status"], "robustId": robust["candidateId"],
        "robustReason": robust["reason"],
        "robustWorstObjective": robust["worstObjective"],
        "robustMeasurements": robust.get("seedMeasurements", []),
        "scientificGateChanged": False, "confirmation": False,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    x=validate(json.loads(args.input.read_text()))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x, indent=2))
    print("PR51_ICA_ROBUST_FEASIBILITY",json.dumps(x))


if __name__=="__main__":
    main()
