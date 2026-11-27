"""Basic tests:  cd python && pytest"""
import gzip
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score

from robust_upgma import (madn, min_cluster_size, robust_upgma, triangle_repair, upgma, w_upgma,
                          winsorise_distances)

FIXTURES = Path(__file__).resolve().parents[2] / "R" / "fixtures"


def blobs(n_per=40, seed=0):
    rng = np.random.default_rng(seed)
    centres = np.array([[0, 0], [5, 0], [2.5, 4.3]])
    y = np.repeat([0, 1, 2], n_per)
    return centres[y] + rng.normal(size=(3 * n_per, 2)), y


def test_madn_normal():
    v = np.random.default_rng(0).normal(scale=2.0, size=200_000)
    assert madn(v) == pytest.approx(2.0, rel=0.02)


def test_step1_leaves_clean_metric_data_unchanged():
    X, _ = blobs()
    D = squareform(pdist(X))
    D1, n_changed = triangle_repair(D)
    assert n_changed == 0
    np.testing.assert_array_equal(D1, D)


def test_winsorised_values_within_bounds():
    X, _ = blobs()
    d, info = winsorise_distances(pdist(X), k=2.5)
    assert d.min() >= info["lower"] - 1e-12 and d.max() <= info["upper"] + 1e-12


def test_square_and_condensed_input_agree():
    X, _ = blobs()
    d = pdist(X)
    np.testing.assert_array_equal(robust_upgma(d, 3).labels, robust_upgma(squareform(d), 3).labels)


def test_recovers_clean_clusters():
    X, y = blobs()
    res = robust_upgma(pdist(X), 3)
    assert adjusted_rand_score(y, res.labels) > 0.95
    assert len(np.unique(res.labels)) == 3


def test_robust_to_corrupted_entries():
    X, y = blobs(seed=1)
    rng = np.random.default_rng(2)
    d = pdist(X)
    bad = rng.random(d.size) < 0.10
    d[bad] = rng.choice([0.0, 10 * d.max()], size=bad.sum())
    ari_r = adjusted_rand_score(y, robust_upgma(d, 3).labels)
    ari_u = adjusted_rand_score(y, upgma(d, 3))
    assert ari_r > 0.8 and ari_r > ari_u


def test_min_cluster_size():
    assert min_cluster_size(300, 3, 0.2) == 20
    assert min_cluster_size(5, 3, 0.2) == 2


def test_all_steps_off_equals_upgma():
    X, _ = blobs()
    d = pdist(X)
    res = robust_upgma(d, 3, triangle=False, winsorise=False, trim=False, size_cut=False)
    np.testing.assert_array_equal(res.labels, upgma(d, 3))


def test_only_step2_equals_w_upgma():
    X, _ = blobs()
    d = pdist(X)
    res = robust_upgma(d, 3, triangle=False, trim=False, size_cut=False)
    np.testing.assert_array_equal(res.labels, w_upgma(d, 3))


@pytest.mark.skipif(not (FIXTURES / "cases.csv").exists(), reason="R/fixtures not found")
def test_reference_partitions():
    """Partitions stored in R/fixtures (used to check the R port) are reproduced exactly."""
    import pandas as pd
    cases = pd.read_csv(FIXTURES / "cases.csv")
    labs = pd.read_csv(FIXTURES / "python_labels.csv.gz")
    labs = labs[labs.method == "R-UPGMA"].set_index("case").labels
    for row in cases.itertuples():
        with gzip.open(FIXTURES / f"{row.case}_d.csv.gz", "rt") as fh:
            d = np.array([float(v) for v in fh.read().split()])
        got = robust_upgma(d, int(row.K)).labels
        exp = np.array([int(v) for v in labs[row.case].split()])
        np.testing.assert_array_equal(got, exp, err_msg=row.case)
