"""Export test cases for R/verify_against_python.R: for 47 cases (43 simulated
conditions and 4 benchmark sets) the exact dissimilarities, the reference set that
Step 1 uses, a checksum after Steps 1-2, and the partition of every method.

    python export_r_fixtures.py
Output: R/fixtures/
"""
from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from common import ABLATIONS, METHODS, REPO, corrupt_distances, load_dataset, run_methods
from simulation import make_dataset
from robust_upgma import triangle_repair, winsorise_distances

OUT = REPO / "R" / "fixtures"


def cases():
    from common import conditions
    for cid in range(0, 171, 4):
        struct, sep, corr, eps = conditions()[cid]
        _, y, d, _, K = make_dataset(cid, 0)
        yield f"sim{cid:03d}_{struct}_{sep}_{corr}_{eps}", d, y, K
    for name, kind, eps in [("Hepta", None, 0), ("Iris", "salt_pepper", .1), ("Penguins", None, 0),
                            ("Lsun3D", "mixed", .1)]:
        X, y, _ = load_dataset(name)
        from scipy.spatial.distance import pdist
        d = pdist(X)
        if kind:
            d = corrupt_distances(d, kind, eps, np.random.default_rng(7), y=y, n=len(y))
        yield f"bench_{name}_{kind or 'clean'}", d, y, len(np.unique(y))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta, labs = [], []
    for name, d, y, K in cases():
        n = len(y)
        ref = np.arange(n) if n <= 64 else np.sort(np.random.default_rng(0).choice(n, 64, replace=False))
        d12, _ = winsorise_distances(triangle_repair(d)[0], 2.5)
        with gzip.open(OUT / f"{name}_d.csv.gz", "wt") as f:
            f.write("\n".join(repr(float(v)) for v in d))
        meta.append(dict(case=name, n=n, K=K, ref=" ".join(str(i + 1) for i in ref),
                         y=" ".join(str(int(v) + 1) for v in y), checksum_step12=float(d12.sum())))
        for m, lab in run_methods(d, K, METHODS + list(ABLATIONS)).items():
            labs.append(dict(case=name, method=m, labels=" ".join(str(int(v)) for v in lab)))
    pd.DataFrame(meta).to_csv(OUT / "cases.csv", index=False)
    pd.DataFrame(labs).to_csv(OUT / "python_labels.csv.gz", index=False)
    print(len(meta), "cases ->", OUT)


if __name__ == "__main__":
    main()
