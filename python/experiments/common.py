"""Shared pieces of all experiments: simulation design, corruption models
(Table 2), competing methods (Section 3.2.5), metrics (Section 3.6) and data.

Method names are those of the paper:
  R-UPGMA       Robust UPGMA (Algorithm 1)
  W-UPGMA       Step 2 alone (MADn winsorisation + UPGMA)
  UPGMA, Single, Complete        cut at K
  UPGMA+cut, Complete+cut        the size-constrained cut of Step 4 on raw dissimilarities
  PAM, HDBSCAN
  w/o Step 1 (triangle), w/o Step 2 (MADn bounds), w/o Step 3 (trimming), w/o Step 4 (cut)
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.metrics.cluster import contingency_matrix

PYTHON_DIR = Path(__file__).resolve().parents[1]
REPO = PYTHON_DIR.parent
DATA = REPO / "data"
RESULTS = REPO / "results" / "python"
sys.path.insert(0, str(PYTHON_DIR))
from robust_upgma import (cut_tree, min_cluster_size, robust_upgma, size_constrained_cut,  # noqa: E402
                          w_upgma)

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------------- simulation design
STRUCTURES = {
    "B3": dict(K=3, p=5, sizes=(50, 50, 50), sd=(1, 1, 1)),
    "U4": dict(K=4, p=5, sizes=(100, 50, 30, 20), sd=(1, 0.5, 1.5, 0.75)),
    "H3": dict(K=3, p=50, sizes=(50, 50, 50), sd=(1, 1, 1)),
}
SEPARATIONS = {"well": 6.0, "moderate": 4.0, "overlap": 2.5}
CORRUPTIONS = ["salt_pepper", "amplification", "mixed", "gaussian", "row", "structured"]
LEVELS = [0.05, 0.10, 0.20]


def conditions():
    """The 171 conditions of the simulation, in a fixed order (condition id = index)."""
    out = []
    for s in STRUCTURES:
        for sep in SEPARATIONS:
            out.append((s, sep, "none", 0.0))
            out += [(s, sep, c, e) for c in CORRUPTIONS for e in LEVELS]
    return out


def simplex_centres(K, p, delta):
    """K centres with all pairwise distances equal to delta (needs p >= K)."""
    C = np.zeros((K, p))
    C[np.arange(K), np.arange(K)] = delta / np.sqrt(2)
    return C


def generate(struct, sep, rng):
    cfg = STRUCTURES[struct]
    C = simplex_centres(cfg["K"], cfg["p"], SEPARATIONS[sep])
    X = np.vstack([C[g] + cfg["sd"][g] * rng.standard_normal((nk, cfg["p"]))
                   for g, nk in enumerate(cfg["sizes"])])
    return X, np.repeat(np.arange(cfg["K"]), cfg["sizes"])


# ----------------------------------------------------------------------------- corruption (Table 2)
def corrupt_distances(d, kind, eps, rng, y=None, n=None):
    """Corrupted copy of the condensed distance vector d (models a-e of Table 2)."""
    d = d.copy()
    m = d.size
    if kind == "structured":                     # (e) amplification of pairs touching cluster 0
        i, j = np.triu_indices(n, k=1)
        elig = np.flatnonzero((y[i] == 0) | (y[j] == 0))
        n_sel = min(int(round(eps * m)), elig.size)
        S = rng.choice(elig, size=n_sel, replace=False)
        d[S] = d[S] * rng.uniform(2, 5, n_sel)
        return d
    n_sel = int(round(eps * m))
    S = rng.choice(m, size=n_sel, replace=False)
    if kind == "salt_pepper":                    # (a)
        q01, q99 = np.quantile(d, [0.01, 0.99])
        salt = rng.random(n_sel) < 0.5
        d[S] = np.where(salt, q99 * rng.uniform(1.5, 3.0, n_sel), q01 * rng.uniform(0.0, 1.0, n_sel))
    elif kind == "amplification":                # (b)
        d[S] = d[S] * rng.uniform(2, 5, n_sel)
    elif kind == "mixed":                        # (c)
        F = rng.uniform(2, 5, n_sel)
        up = rng.random(n_sel) < 0.5
        d[S] = np.where(up, d[S] * F, d[S] / F)
    elif kind == "gaussian":                     # (d)
        d[S] = np.abs(d[S] + rng.normal(0, 0.5 * d.std(), n_sel))
    else:
        raise ValueError(kind)
    return d


def contaminate_rows(X, eps, rng):
    """(f) Observation outliers: returns (X_contaminated, outlier_mask)."""
    X = X.copy()
    n, p = X.shape
    n_out = int(round(eps * n))
    idx = rng.choice(n, size=n_out, replace=False)
    u = rng.normal(size=(n_out, p))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    delta = rng.uniform(5, 10, size=(n_out, 1)) * X.std() * np.sqrt(p)
    X[idx] = X[idx] + delta * u
    mask = np.zeros(n, bool)
    mask[idx] = True
    return X, mask


# ----------------------------------------------------------------------------- methods (Section 3.2.5)
METHODS = ["R-UPGMA", "W-UPGMA", "UPGMA", "UPGMA+cut", "Single", "Complete", "Complete+cut", "PAM", "HDBSCAN"]
ABLATIONS = {"w/o Step 1 (triangle)": dict(triangle=False), "w/o Step 2 (MADn bounds)": dict(winsorise=False),
             "w/o Step 3 (trimming)": dict(trim=False), "w/o Step 4 (cut)": dict(size_cut=False)}


def run_methods(d, K, methods=METHODS):
    """Partitions of every method on the same condensed dissimilarity vector d."""
    n = squareform(d).shape[0]
    D = squareform(d)
    out = {}
    for m in methods:
        if m == "R-UPGMA":
            out[m] = robust_upgma(d, K).labels
        elif m in ABLATIONS:
            out[m] = robust_upgma(d, K, **ABLATIONS[m]).labels
        elif m == "W-UPGMA":
            out[m] = w_upgma(d, K)
        elif m in ("UPGMA", "Single", "Complete"):
            out[m] = cut_tree(linkage(d, {"UPGMA": "average"}.get(m, m.lower())), K)
        elif m in ("UPGMA+cut", "Complete+cut"):
            Z = linkage(d, "average" if m == "UPGMA+cut" else "complete")
            out[m] = size_constrained_cut(Z, D, K, min_cluster_size(n, K))
        elif m == "PAM":
            import kmedoids
            out[m] = np.asarray(kmedoids.pam(D, K, random_state=0).labels)
        elif m == "HDBSCAN":                     # K not used; noise points (-1) form one group
            from sklearn.cluster import HDBSCAN
            out[m] = HDBSCAN(min_cluster_size=10, metric="precomputed", copy=True).fit_predict(D)
        else:
            raise ValueError(m)
    return out


# ----------------------------------------------------------------------------- metrics (Section 3.6)
def matched_accuracy(y, lab):
    """Accuracy after optimal one-to-one matching (Hungarian algorithm)."""
    C = contingency_matrix(y, lab)
    r, c = linear_sum_assignment(-C)
    return C[r, c].sum() / len(y)


def evaluate(y, lab, D_eval=None, mask=None):
    """ARI, matched accuracy, number of clusters, share of the largest cluster and the
    silhouette in the clean Euclidean geometry D_eval (optionally on mask only)."""
    if mask is not None:
        y, lab = y[mask], lab[mask]
    res = dict(ari=adjusted_rand_score(y, lab), acc=matched_accuracy(y, lab),
               n_clusters=len(np.unique(lab)), largest=np.bincount(lab - lab.min()).max() / len(lab))
    if D_eval is not None:
        Dm = D_eval if mask is None else D_eval[np.ix_(mask, mask)]
        nl = len(np.unique(lab))
        res["sil"] = silhouette_score(Dm, lab, metric="precomputed") if 1 < nl < len(lab) else np.nan
    return res


# ----------------------------------------------------------------------------- data (Section 3.4)
FCPS = ["Hepta", "Tetra", "Lsun3D", "Chainlink", "Atom", "EngyTime", "TwoDiamonds", "WingNut", "Target"]
REAL = {"Iris": "iris", "Wine": "wine_recognition", "Breast Cancer": "breast_cancer_wisconsin_diagnostic",
        "E. coli": "ecoli", "Penguins": "penguins", "MAGIC": "magic"}


def load_dataset(name):
    """(X, y, source). FCPS data as published; real data z-scored, constant features
    removed, Penguins restricted to the four morphometric variables (complete cases)."""
    if name in FCPS or name == "GolfBall":
        import rdata
        r = rdata.read_rda(str(DATA / "fcps" / f"{name}.rda"))[name]
        X = np.asarray(r["Data"], float)
        y = np.unique(np.asarray(r["Cls"]).astype(int).ravel(), return_inverse=True)[1]
        return X, y, "FCPS"
    df = pd.read_csv(DATA / "pmlb" / f"{REAL[name]}.tsv.gz", sep="\t")
    if name == "Penguins":
        df = df[["bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "target"]].dropna()
    y = np.unique(df.pop("target").to_numpy(), return_inverse=True)[1]
    X = df.to_numpy(float)
    X = X[:, X.std(0) > 0]
    return (X - X.mean(0)) / X.std(0), y, "Real"
