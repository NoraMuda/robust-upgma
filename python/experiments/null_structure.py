"""Data without cluster structure (Sections 3.5 and 4.4).

Null data: NULL-Gauss (150 x 5 from N(0, I)), NULL-Unif (150 x 5 from U[0,1]^5),
GolfBall (FCPS); positive controls: B3 well separated and Hepta. GolfBall and Hepta
enter as random subsamples of 150 points, so every dip test uses the same number of
distances. Each set: 100 replicates, clean and with 10% salt-and-pepper corruption.
Recorded: dip-test p-values on the raw and on the winsorised distances, and the
partition of Robust UPGMA asked for K = 3 (group shares, clean-geometry silhouette).
Seeds: numpy SeedSequence([20260525, replicate]).

    python null_structure.py
Output: results/python/null_trials.csv.gz
"""
from __future__ import annotations

import diptest
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score

from common import RESULTS, corrupt_distances, generate, load_dataset
from robust_upgma import robust_upgma, winsorise_distances

R = 100


def record(name, setting, X, d, rep):
    d_win, _ = winsorise_distances(d, 2.5)
    lab = robust_upgma(d, 3).labels
    sizes = np.sort(np.bincount(lab, minlength=3))[::-1] / len(lab)
    return dict(dataset=name, setting=setting, rep=rep, dip_p=diptest.diptest(d)[1],
                dip_p_winsorised=diptest.diptest(d_win)[1], largest=sizes[0], smallest=sizes[-1],
                sil=silhouette_score(squareform(pdist(X)), lab, metric="precomputed")
                if len(np.unique(lab)) > 1 else np.nan)


def main():
    rows = []
    golf, hepta = load_dataset("GolfBall")[0], load_dataset("Hepta")[0]
    for r in range(R):
        rng = np.random.default_rng(np.random.SeedSequence([20260525, r]))
        sets = {"NULL-Gauss": rng.standard_normal((150, 5)),
                "NULL-Unif": rng.random((150, 5)),
                "B3-well (control)": generate("B3", "well", rng)[0],
                "GolfBall": golf[rng.choice(len(golf), 150, replace=False)],
                "Hepta (control)": hepta[rng.choice(len(hepta), 150, replace=False)]}
        for name, X in sets.items():
            d = pdist(X)
            rows.append(record(name, "clean", X, d, r))
            rows.append(record(name, "salt_pepper10", X, corrupt_distances(d, "salt_pepper", 0.10, rng), r))
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(RESULTS / "null_trials.csv.gz", index=False)
    print("->", RESULTS / "null_trials.csv.gz")


if __name__ == "__main__":
    main()
