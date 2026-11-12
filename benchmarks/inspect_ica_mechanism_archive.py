#!/usr/bin/env python3
"""Read-only inspection of the archived ERP CORE scoring audit schema.

This job never refits ICA or chooses a method. Do not make scientific claims
from missing MAT fields; inspect availability before designing downstream
readouts.
"""
import json
from pathlib import Path
import scipy.io

def items(x):
    if x is None:
        return []
    import numpy as np
    return list(x.flat) if isinstance(x, np.ndarray) else [x]

for mat in sorted(Path("evidence").rglob("scoring_audit.mat")):
    r=Path(str(mat).replace("scoring_audit.mat","result.json"))
    record=json.loads(r.read_text())
    a=scipy.io.loadmat(str(mat), squeeze_me=True, struct_as_record=False)["audit"]
    cands=items(a.cands)
    print("PR51_MAT_SCHEMA",json.dumps({
        "subject":record["subject"],"nFeasibleBaseline":len(cands),
        "firstCandFields":list(cands[0]._fieldnames) if cands else [],
        "selectedIds":[record["icaSeedCheck"]["recommended"][0]],
        "stepsFirstType":str(type(cands[0].steps)) if cands else None,
        "firstSteps":[{"type":str(getattr(x,"type","")),
                       "icsRemoved":str(getattr(x,"icsRemoved","")),
                       "comps":str(getattr(x,"comps",""))[:50]}
                      for x in items(cands[0].steps)][:12] if cands else [],
        "firstCandidateId":float(cands[0].id) if cands else None,
    }))
