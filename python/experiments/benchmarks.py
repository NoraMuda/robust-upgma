"""Benchmarks on FCPS and labelled real data (Sections 3.4 and 5).

Settings: clean, 10% salt-and-pepper, 10% mixed and 5% observation outliers; 30
realisations of each corrupted setting (and of the clean setting for EngyTime and
MAGIC, which are class-stratified subsampled to 2,000 observations in every
realisation). All methods are run on the same realisation.
Seeds: numpy SeedSequence([20260524, data-set index, setting index, realisation]).

    python benchmarks.py                 # 30 realisations, all cores
Output: results/python/bench_trials.csv.gz
"""
from __future__ import annotations

import argparse
import os
import pickle
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from common import (FCPS, METHODS, REAL, RESULTS, contaminate_rows, corrupt_distances, evaluate,
                    load_dataset, run_methods)

DATASETS = FCPS + list(REAL)
SETTINGS = [("clean", None, 0.0), ("salt_pepper10", "salt_pepper", 0.10), ("mixed10", "mixed", 0.10),
            ("row05", "row", 0.05)]
MAX_N = 2000


def one(task):
    di, si, r = task
    name, (sname, kind, eps) = DATASETS[di], SETTINGS[si]
    X, y, src = load_dataset(name)
    rng = np.random.default_rng(np.random.SeedSequence([20260524, di, si, r]))
    if len(y) > MAX_N:                           # class-stratified subsample
        idx = np.concatenate([rng.choice(np.flatnonzero(y == g), int(round(MAX_N * np.mean(y == g))),
                                         replace=False) for g in np.unique(y)])
        X, y = X[idx], y[idx]
    n, K = len(y), len(np.unique(y))
    mask = None
    if kind is None:
        d = pdist(X)
    elif kind == "row":
        Xc, out = contaminate_rows(X, eps, rng)
        d, mask = pdist(Xc), ~out
    else:
        d = corrupt_distances(pdist(X), kind, eps, rng, y=y, n=n)
    D_clean = squareform(pdist(X))
    return [dict(dataset=name, source=src, n=n, K=K, setting=sname, rep=r, method=m,
                 **evaluate(y, lab, D_eval=D_clean, mask=mask)) for m, lab in run_methods(d, K, METHODS).items()]


def _run(task):
    return task, one(task)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    tasks = []
    for di, name in enumerate(DATASETS):
        n = len(load_dataset(name)[1])
        for si, (_, kind, _) in enumerate(SETTINGS):
            reps = a.reps if (kind is not None or n > MAX_N) else 1
            tasks += [(di, si, r) for r in range(reps)]
    tasks.sort(key=lambda t: -min(len(load_dataset(DATASETS[t[0]])[1]), MAX_N))   # big jobs first
    RESULTS.mkdir(parents=True, exist_ok=True)
    cache = RESULTS / "bench_partial.pkl"        # lets an interrupted run resume
    done = pickle.loads(cache.read_bytes()) if cache.exists() else {}
    todo = [t for t in tasks if t not in done]
    t0 = time.time()
    with Pool(a.workers) as pool:
        for i, (task, rows) in enumerate(pool.imap_unordered(_run, todo), 1):
            done[task] = rows
            if i % 25 == 0 or i == len(todo):
                cache.write_bytes(pickle.dumps(done))
                print(f"{len(done)}/{len(tasks)} runs, {time.time() - t0:.0f} s", flush=True)
    res = [done[t] for t in sorted(done)]
    pd.DataFrame([row for rows in res for row in rows]).to_csv(RESULTS / "bench_trials.csv.gz", index=False)
    cache.unlink()
    print(f"{len(tasks)} runs in {time.time() - t0:.0f} s -> {RESULTS / 'bench_trials.csv.gz'}")


if __name__ == "__main__":
    main()
