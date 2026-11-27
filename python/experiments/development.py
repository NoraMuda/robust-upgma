"""Development study (Section 3.3): the data on which the steps of Robust UPGMA were
chosen. It uses its own seed stream, SeedSequence([99001, condition, replicate]),
never the seeds of the confirmatory simulation (20260523), and none of the FCPS,
real or null data. 171 conditions x 10 replicates = 1,710 data sets.

Round 1 (exploratory): variants of the triangle repair (Step 1) with trimming (Step 3).
Round 2 (Supplementary Table S5): six pipelines combining the four steps; the best,
pipeline C (Steps 1-4), was frozen as Robust UPGMA.

    python development.py --round 2
Output: results/python/development_round{1,2}.csv.gz
"""
from __future__ import annotations

import argparse
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform

from common import (RESULTS, STRUCTURES, conditions, contaminate_rows, corrupt_distances, cut_tree, evaluate,
                    generate, run_methods)
from robust_upgma import robust_upgma, row_outliers, winsorise_distances

ROUND2 = {
    "A: Steps 1+3": dict(winsorise=False, size_cut=False),
    "B: Steps 1+3+4": dict(winsorise=False),
    "C: Steps 1+2+3+4 (chosen)": dict(),
    "D: Steps 1+2+3": dict(size_cut=False),
    "E: Steps 3+4": dict(triangle=False, winsorise=False),
    "F: Steps 1+3+4, c = 2.5": dict(winsorise=False, c=2.5),
}
ROUND1 = {
    "1a: Step 1 (medians) + 3": dict(ql=0.5, qu=0.5),
    "1b: Step 1 (quantiles 0.75/0.25) + 3": dict(ql=0.75, qu=0.25),
    "1c: Step 1 (quantiles 0.6/0.4) + 3": dict(ql=0.6, qu=0.4),
    "1d: Step 1 (upper bounds only) + 3": dict(sides="upper"),
    "1e: Step 1 applied twice + 3": dict(passes=2),
    "1f: Step 2 + 3": dict(passes=0, winsorise=True),
}


# ---- round-1 helpers: triangle bounds with quantiles other than the median ----
def _bounds_q(D, s, rng, ql, qu):
    n = D.shape[0]
    R = np.arange(n) if n <= s else np.sort(rng.choice(n, s, replace=False))
    A = D[:, R]
    L, U = np.zeros((n, n)), np.zeros((n, n))
    for i in range(n - 1):
        B = A[i + 1:]
        L[i, i + 1:] = np.quantile(np.abs(A[i] - B), ql, axis=1)
        U[i, i + 1:] = np.quantile(A[i] + B, qu, axis=1)
    return L + L.T, U + U.T


def round1_variant(d, K, ql=0.5, qu=0.5, sides="both", passes=1, winsorise=False, s=64, c=3.0):
    D = squareform(d)
    rng = np.random.default_rng(0)
    for _ in range(passes):
        L, U = _bounds_q(D, s, rng, ql, qu)
        if sides == "upper":
            L = np.zeros_like(L)
        D = np.minimum(np.maximum(D, L), U)
    if winsorise:
        D = squareform(winsorise_distances(squareform(D, checks=False), 2.5)[0])
    n = D.shape[0]
    out = row_outliers(D, c)                     # Step 3 with a plain cut at K
    inl = np.flatnonzero(~out)
    if len(inl) < K + 1:
        inl, out = np.arange(n), np.zeros(n, bool)
    lab_in = cut_tree(linkage(squareform(D[np.ix_(inl, inl)], checks=False), "average"), K)
    lab = np.empty(n, int)
    lab[inl] = lab_in
    for i in np.flatnonzero(out):
        lab[i] = int(np.argmin([D[i, inl[lab_in == g]].mean() for g in range(K)]))
    return lab


def one(task):
    rnd, cid, r = task
    struct, sep, corr, eps = conditions()[cid]
    rng = np.random.default_rng(np.random.SeedSequence([99001, cid, r]))
    X, y = generate(struct, sep, rng)
    K, mask = STRUCTURES[struct]["K"], None
    if corr == "none":
        d = pdist(X)
    elif corr == "row":
        Xc, out = contaminate_rows(X, eps, rng)
        d, mask = pdist(Xc), ~out
    else:
        d = corrupt_distances(pdist(X), corr, eps, rng, y=y, n=len(y))
    labs = run_methods(d, K, ["UPGMA", "W-UPGMA", "PAM"])
    for name, kw in (ROUND1 if rnd == 1 else ROUND2).items():
        labs[name] = round1_variant(d, K, **kw) if rnd == 1 else robust_upgma(d, K, **kw).labels
    return [dict(structure=struct, separation=sep, corruption=corr, level=eps, rep=r, method=m,
                 ari=evaluate(y, lab, mask=mask)["ari"]) for m, lab in labs.items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--round", type=int, choices=[1, 2], default=2)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    tasks = [(a.round, cid, r) for cid in range(len(conditions())) for r in range(a.reps)]
    t0 = time.time()
    with Pool(a.workers) as pool:
        res = pool.map(one, tasks, chunksize=10)
    out = RESULTS / f"development_round{a.round}.csv.gz"
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row for rows in res for row in rows]).to_csv(out, index=False)
    print(f"{len(tasks)} data sets in {time.time() - t0:.0f} s -> {out}")


if __name__ == "__main__":
    main()
