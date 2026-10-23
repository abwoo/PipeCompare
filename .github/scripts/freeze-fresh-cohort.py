import hashlib
import json
import re
import urllib.request
urls=[
'https://physionet.org/files/eegmmidb/1.0.0/SHA256SUMS.txt',
'https://physionet.org/files/eegmmidb/1.0.0/SHA256SUMS',
'https://physionet-open.s3.amazonaws.com/eegmmidb/1.0.0/SHA256SUMS.txt']
entries={};used=None;raw=None
for url in urls:
    try:
        with urllib.request.urlopen(url,timeout=60) as response:raw=response.read()
        for line in raw.decode().splitlines():
            parts=line.split()
            if len(parts)==2 and re.fullmatch('[0-9a-f]{64}',parts[0]):
                entries[parts[1].lstrip('*/.')]=parts[0]
        if 'S021/S021R01.edf' in entries:
            used=url;break
    except Exception as ex:print(type(ex).__name__+': '+str(ex))
assert used and entries['S021/S021R01.edf']=='983781fb39c3a747051a04f23484164de485b97ba43a7fd63fd8fe3d7f6f6f3c'
files=[]
for subject in range(51,81):
    path=f'S{subject:03d}/S{subject:03d}R01.edf'
    files.append(dict(subject=subject,path=path,sha256=entries[path]))
print('PIPECOMPARE_FRESH_MANIFEST_JSON '+json.dumps(dict(manifestURL=used,manifestSHA256=hashlib.sha256(raw).hexdigest(),files=files)))
