"""Confirmatory Monte Carlo simulation (Sections 3.2 and 4).

171 conditions (3 structures x 3 separations x [6 corruption models x 3 levels + clean])
x 100 replicates = 17,100 data sets. In every replicate all methods receive the same
corrupted dissimilarities; the ablation variants are run on the first 30 replicates.
Seeds: numpy SeedSequence([20260523, condition id, replicate]).

    python simulation.py                 # 100 replicates, all cores
    python simulation.py --reps 5        # quick check
Output: results/python/sim_trials.csv.gz (one row per data set x method)
"""
from __future__ import annotations

import argparse
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from common import (ABLATIONS, METHODS, RESULTS, STRUCTURES, conditions, contaminate_rows,
                    corrupt_distances, evaluate, generate, run_methods)

N_ABLATION_REPS = 30


def make_dataset(cid, r):
    """The (corrupted) data set of condition cid, replicate r."""
    struct, sep, corr, eps = conditions()[cid]
    rng = np.random.default_rng(np.random.SeedSequence([20260523, cid, r]))
    X, y = generate(struct, sep, rng)
    mask = None
    if corr == "none":
        d = pdist(X)
    elif corr == "row":
        Xc, out = contaminate_rows(X, eps, rng)
        d, mask = pdist(Xc), ~out                # scored on the uncontaminated observations
    else:
        d = corrupt_distances(pdist(X), corr, eps, rng, y=y, n=len(y))
    return X, y, d, mask, STRUCTURES[struct]["K"]


def one_trial(task):
    cid, r = task
    struct, sep, corr, eps = conditions()[cid]
    X, y, d, mask, K = make_dataset(cid, r)
    methods = METHODS + (list(ABLATIONS) if r < N_ABLATION_REPS else [])
    D_clean = squareform(pdist(X))
    return [dict(cond=cid, structure=struct, separation=sep, corruption=corr, level=eps, rep=r, method=m,
                 **evaluate(y, lab, D_eval=D_clean, mask=mask)) for m, lab in run_methods(d, K, methods).items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    tasks = [(cid, r) for cid in range(len(conditions())) for r in range(a.reps)]
    t0 = time.time()
    with Pool(a.workers) as pool:
        res = pool.map(one_trial, tasks, chunksize=10)
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row for rows in res for row in rows]).to_csv(RESULTS / "sim_trials.csv.gz", index=False)
    print(f"{len(tasks)} data sets in {time.time() - t0:.0f} s -> {RESULTS / 'sim_trials.csv.gz'}")


if __name__ == "__main__":
    main()
