#!/usr/bin/env python3
"""Instrument a disposable checkout of Executor.m to reveal caught icremove stack traces.

Do NOT commit patched production source. This adds stderr logging only; scientific
outputs from this temporary instrumented run are diagnostic, never confirmation.
"""
from pathlib import Path
import argparse
import subprocess

ANCHOR = "pipecompare.utils.log('FAILED %s: %s (%d candidate(s) affected)', in.label, ME.message, numel(under));"
ADDITION = """
if strcmp(in.type, 'icremove')
    fprintf(2, 'PR51_ICREMOVE_STACK subject=%s\\n%s\\n', getenv('PILOT_SUBJECT'), ...
        getReport(ME, 'extended', 'hyperlinks', 'off'));
end
"""

def patch(src):
    if src.count(ANCHOR) != 1:
        raise ValueError(f"expected 1 frozen Executor.failUnder anchor, got {src.count(ANCHOR)}")
    return src.replace(ANCHOR, ANCHOR + ADDITION, 1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    args = parser.parse_args()
    root = args.repo
    target = root / "+pipecompare/+run/Executor.m"
    original = target.read_text()
    changed = patch(original)
    assert len(changed) > len(original)
    target.write_text(changed)
    source = subprocess.run(["git","diff","--","+pipecompare/+run/Executor.m"],
                            cwd=root,capture_output=True,text=True,check=True).stdout
    if source.count("+    fprintf(2, 'PR51_ICREMOVE_STACK") != 1 or source.count("@@") == 0:
        raise ValueError("unexpected diagnostic patch diff; refuse to run")
    (root.parent/"result").mkdir(exist_ok=True)
    (root.parent/"result"/"diagnostic_instrumentation.diff").write_text(source)
    print("PR51 temporary execution logging instrumentation applied; original branch source unchanged.")

if __name__=="__main__":
    main()
