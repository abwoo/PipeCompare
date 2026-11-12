#!/usr/bin/env python3
"""Patch a disposable cloud checkout solely to capture MATLAB exception origins.

The temporary changed Executor.m is never committed. No scientific parameters,
rankings, production source or frozen protocol are changed in the Git branch.
MATLAB evalc captures normal transcript, so stack frames are also attached to
candidate failure messages in diagnostic result.json (not scientific outputs).
"""
from pathlib import Path
import argparse
import subprocess

ANCHOR = "c.status = failStatus(ME); c.message = sprintf('%s: %s', in.label, ME.message);"
INJECTION = """
    if strcmp(in.type, 'icremove') && ~isempty(ME.stack)
        frames = arrayfun(@(s) sprintf('%s:%d', s.name, s.line), ...
            ME.stack, 'UniformOutput', false);
        c.message = sprintf('%s [PR51_STACK %s]', c.message, strjoin(frames, ' > '));
    end
"""

def patch(src):
    if src.count(ANCHOR) != 1:
        raise ValueError(f"Expected one Executor.failUnder message anchor, got {src.count(ANCHOR)}")
    return src.replace(ANCHOR, ANCHOR + INJECTION, 1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("repo", type=Path)
    args = parser.parse_args()
    target = args.repo / "+pipecompare/+run/Executor.m"
    original = target.read_text()
    target.write_text(patch(original))
    diff = subprocess.run(["git", "diff", "--", "+pipecompare/+run/Executor.m"],
                          cwd=args.repo, capture_output=True, text=True, check=True).stdout
    if diff.count("+    if strcmp(in.type, 'icremove')") != 1 or "PR51_STACK" not in diff:
        raise ValueError("Unexpected patch; abort before EEG processing")
    (args.repo.parent / "result").mkdir(exist_ok=True)
    (args.repo.parent / "result" / "diagnostic_instrumentation.diff").write_text(diff)
    print("PR51 temporary diagnostic message enrichment applied; no scientific parameter changes.")

if __name__ == "__main__":
    main()
