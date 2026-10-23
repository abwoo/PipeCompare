import csv
import json
import re
from pathlib import Path
root=Path('timing_source')
result=[]
with (root/'recordings.csv').open(encoding='utf-8-sig',newline='') as stream:
    rows=list(csv.DictReader(stream))
for row in rows:
    subject=int(float(row['subject']))
    transcript=(root/f'S{subject:03d}R01.mat.log').read_text(errors='replace')
    values=[float(v) for v in re.findall(r'Finished in ([0-9.]+) s',transcript)]
    reported=float(row['seconds'])
    lower=max(values) if values else None
    result.append(dict(subject=subject,reportedFitSeconds=reported,engineElapsedReports=values,
        explicitEngineLowerBoundSeconds=lower,
        implicitTimerContradiction=lower is not None and lower>reported+1))
print('PIPECOMPARE_TIMING_JSON '+json.dumps(result))
