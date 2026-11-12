#!/usr/bin/env python3
"""Log reject_threshold MATLAB stack frames only in a disposable cloud checkout.

Do not commit the patched Executor.m; this is forensic instrumentation, not
the algorithm under validation. The recorder audits the production diff
before patching and archives this temporary instrumentation diff.
"""
from pathlib import Path
import argparse
import subprocess

ANCHOR = "c.status = failStatus(ME); c.message = sprintf('%s: %s', in.label, ME.message);"
EXTRA = """
    if strcmp(in.type, 'reject_threshold') && ~isempty(ME.stack)
        frames = arrayfun(@(s) sprintf('%s:%d', s.name, s.line), ...
            ME.stack, 'UniformOutput', false);
        c.message = sprintf('%s [PR51_REJECT_STACK %s]', c.message, strjoin(frames, ' > '));
    end
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("repository", type=Path)
    root = parser.parse_args().repository
    target = root / "+pipecompare/+run/Executor.m"
    content = target.read_text()
    assert content.count(ANCHOR) == 1, "Cannot uniquely identify frozen failure marker"
    target.write_text(content.replace(ANCHOR, ANCHOR + EXTRA, 1))
    diff = subprocess.run(["git", "diff", "--", "+pipecompare/+run/Executor.m"],
                          cwd=root, capture_output=True, text=True, check=True).stdout
    if "PR51_REJECT_STACK" not in diff or diff.count("@@") != 2:
        raise ValueError("unexpected disposable diagnostic patch")
    (root.parent / "result").mkdir(exist_ok=True)
    (root.parent / "result" / "reject_threshold_trace.diff").write_text(diff)
    print("PR51 temporary reject_threshold failure stack instrumentation applied")

if __name__ == "__main__":
    main()
