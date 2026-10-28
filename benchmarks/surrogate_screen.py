"""Surrogate screen prototype, NOT production code (Round 8 bundle).

Status: PROTOTYPE / DESIGN. No MATLAB implementation, no shadow-arm
wiring, no LLM, no external API, no user data sent anywhere.
What it is: a tiny stdlib-only linear predictor that maps cheap
candidate features (interpolatedFraction, icsRemoved, filter cutoffs,
etc.) to a predicted MSE, for SCREENING ONLY. Survivors must still be
ranked by full-fidelity Rank.run; a low-fidelity / predicted score is
never the final score.
Rules (frozen before any held-out use):
  - Train only on historical seeds (< 20031). Test only on 20031+.
  - Tuning the feature set or threshold on held-out results is double
    use and invalidates the arm.
  - Report, per held-out seed: rank correlation, survival rate of the
    full-fidelity best, paired MSE of the screened pick vs the
    exhaustive pick, wall time, and pick stability. A screen that
    drops the true best is a failure, not a speedup.
  - This file is not called by Executor, Rank, or the shadow script
    in Round 8. Wiring it as a new arm requires the gates above.
"""
import math


def featurize(candidate):
    """candidate: dict with cheap fields. Returns a feature vector."""
    return [
        1.0,
        float(candidate.get("interpolatedFraction", 0.0) or 0.0),
        float(candidate.get("icsRemoved", 0) or 0),
        float(candidate.get("lowpass", 0) or 0),
        float(candidate.get("highpass", 0) or 0),
    ]


def _solve(A, b):
    # Gaussian elimination with partial pivoting, stdlib only.
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-12:
            continue
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                for c in range(col, n + 1):
                    M[r][c] -= f * M[col][c]
    return [M[i][n] / M[i][i] if abs(M[i][i]) > 1e-12 else 0.0 for i in range(n)]


def fit(candidates, mses, ridge=1e-6):
    """Fit a linear predictor. candidates/mses: training rows only."""
    X = [featurize(c) for c in candidates]
    d = len(X[0])
    AtA = [[sum(row[i] * row[j] for row in X) for j in range(d)] for i in range(d)]
    Atb = [sum(row[i] * y for row, y in zip(X, mses)) for i in range(d)]
    for i in range(d):
        AtA[i][i] += ridge
    return _solve(AtA, Atb)


def predict(weights, candidate):
    feats = featurize(candidate)
    return sum(w * f for w, f in zip(weights, feats))


def screen(weights, candidates, keep_fraction=0.25):
    """Return the kept candidate indices (lowest predicted MSE first)."""
    scored = sorted(range(len(candidates)), key=lambda i: predict(weights, candidates[i]))
    n_keep = max(1, int(math.floor(len(candidates) * keep_fraction)))
    return scored[:n_keep]
