"""Remove an unused GPU include lookup from the pinned CPU-only MEX build.
This edits the compiler driver, never ICLabel features, weights or inference.
"""
import hashlib
from pathlib import Path
import subprocess
import sys

root=Path(sys.argv[1]).resolve()
expected='092ecb8c04031ddcb1d9a7ba6ebe76a250c4672c'
actual=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
if actual!=expected:
    raise RuntimeError('Unexpected MatConvNet source: '+actual)
path=root/'matlab'/'vl_compilenn.m'
text=path.read_text()
line="  sprintf('-I\"%s\"',fullfile(toolboxdir('distcomp'),'gpu','extern','include')), ...\n"
if text.count(line)!=1:
    raise RuntimeError('Pinned compiler CPU compatibility anchor changed')
before=hashlib.sha256(path.read_bytes()).hexdigest()
path.write_text(text.replace(line,''))
print('CPU-only compiler compatibility patch; source='+actual+'; originalSHA256='+before)
