"""LLM-assisted review prototype, NOT a search component (Round 8).

Status: PROTOTYPE / AUXILIARY ONLY. No LLM is called by this file,
no external API, no user data sent anywhere, no model weights are
embedded. It is a local, template-based summariser for CI logs and
diffs: it extracts failed test names, timeout lines, and changed
files, so a human (or a future LLM review step, with the user's
explicit data approval) can review faster. It never decides a pick,
never changes a ranking, never writes to the shadow evidence, and
must not be wired into Executor, Rank, or the shadow script.
If a future step does call an external LLM, that step needs its own
approval, data scope, and a rule that LLM output is never evidence:
only CI + held-out paired results are evidence.
"""
import re
import sys


def summarise_log(text):
    failed = sorted(set(re.findall(r"FAILED[:\s]+([A-Za-z0-9_\.]+)", text)))
    timeouts = text.count("timed out")
    return {"failed_tests": failed, "timeout_lines": timeouts, "chars": len(text)}


if __name__ == "__main__":
    data = sys.stdin.read()
    print(summarise_log(data))
