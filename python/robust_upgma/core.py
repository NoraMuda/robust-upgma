"""Robust UPGMA (Algorithm 1 of the paper) and its building blocks.

Notation
--------
D        n x n symmetric dissimilarity matrix, d_ii = 0
d        condensed vector of the m = n(n-1)/2 entries d_ij, i < j
MADn(v)  1.4826 * med|v - med(v)|

Algorithm 1 (defaults s = 64, k = 2.5, c = 3, f = 0.2)
  Step 1  triangle repair: with a reference set R of min(s, n) observations,
          l_ij = med_r |d_ir - d_jr|,  u_ij = med_r (d_ir + d_jr),
          d~_ij = min(max(d_ij, l_ij), u_ij).  Leaves uncorrupted metric data unchanged.
  Step 2  MADn winsorisation: L = max(0, med - k MADn), U = med + k MADn of all d~_ij;
          d'_ij = min(max(d~_ij, L), U).
  Step 3  MADn trimming: rho_i = med_{j != i} d'_ij; observation i is set aside when
          rho_i > med(rho) + c MADn(rho).
  Step 4  UPGMA on the remaining observations and a size-constrained cut: the smallest
          K* >= K whose cut contains K groups with >= ceil(f n_inliers / K) members;
          all other observations join the cluster with the smallest mean dissimilarity.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

MADN_CONSTANT = 1.4826  # 1 / Phi^{-1}(3/4)


# ----------------------------------------------------------------------------- helpers
def madn(v) -> float:
    """Normalised median absolute deviation, 1.4826 * med|v - med(v)|."""
    v = np.asarray(v, dtype=float).ravel()
    return float(MADN_CONSTANT * np.median(np.abs(v - np.median(v))))


def as_square(D) -> np.ndarray:
    """Square dissimilarity matrix from a square matrix or a condensed vector."""
    D = np.asarray(D, dtype=float)
    if D.ndim == 1:
        return squareform(D)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square matrix or a condensed distance vector")
    if not np.allclose(D, D.T):
        raise ValueError("D must be symmetric")
    return D.copy()


def as_condensed(D) -> np.ndarray:
    D = np.asarray(D, dtype=float)
    return D.copy() if D.ndim == 1 else squareform(as_square(D), checks=False)


def cut_tree(Z: np.ndarray, n_clusters: int) -> np.ndarray:
    """Cut a SciPy linkage matrix into (at most) n_clusters groups, labels 0..K-1."""
    return fcluster(Z, t=n_clusters, criterion="maxclust") - 1


# ----------------------------------------------------------------------------- Step 1
def triangle_bounds(D: np.ndarray, s: int = 64, seed: int = 0):
    """Matrices (L, U) of the triangle-median bounds over a reference set of
    min(s, n) observations drawn with numpy.random.default_rng(seed)."""
    D = np.asarray(D, float)
    n = D.shape[0]
    R = np.arange(n) if n <= s else np.sort(np.random.default_rng(seed).choice(n, s, replace=False))
    A = D[:, R]
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    sR, h = A.shape[1], A.shape[1] // 2
    for i in range(n - 1):                       # pairs j > i; mirrored below
        B = A[i + 1:]
        lo = np.abs(A[i] - B)
        up = A[i] + B
        if sR % 2:
            L[i, i + 1:] = np.partition(lo, h, axis=1)[:, h]
            U[i, i + 1:] = np.partition(up, h, axis=1)[:, h]
        else:                                    # median of an even number of values
            pl = np.partition(lo, [h - 1, h], axis=1)
            pu = np.partition(up, [h - 1, h], axis=1)
            L[i, i + 1:] = 0.5 * (pl[:, h - 1] + pl[:, h])
            U[i, i + 1:] = 0.5 * (pu[:, h - 1] + pu[:, h])
    return L + L.T, U + U.T


def triangle_repair(D, s: int = 64, seed: int = 0):
    """Step 1. Returns the repaired square matrix and the number of entries changed."""
    D = as_square(D)
    L, U = triangle_bounds(D, s, seed)
    n_changed = int(np.triu((D < L) | (D > U), 1).sum())
    return np.minimum(np.maximum(D, L), U), n_changed


# ----------------------------------------------------------------------------- Step 2
def winsorise_distances(D, k: float = 2.5):
    """Step 2 (alone: W-UPGMA). Clip every dissimilarity to [L, U].
    Returns the condensed clipped vector and diagnostics."""
    d = as_condensed(D)
    med = float(np.median(d))
    s = madn(d)
    L, U = max(0.0, med - k * s), med + k * s
    info = dict(k=k, median=med, madn=s, lower=L, upper=U,
                n_low=int((d < L).sum()), n_high=int((d > U).sum()))
    return np.clip(d, L, U), info


# ----------------------------------------------------------------------------- Step 3
def row_outliers(D: np.ndarray, c: float = 3.0) -> np.ndarray:
    """Step 3. Boolean mask of observations with rho_i > med(rho) + c MADn(rho)."""
    n = D.shape[0]
    rho = np.median(D[~np.eye(n, dtype=bool)].reshape(n, n - 1), axis=1)
    return rho > np.median(rho) + c * madn(rho)


# ----------------------------------------------------------------------------- Step 4
def size_constrained_cut(Z: np.ndarray, D: np.ndarray, K: int, min_size: int) -> np.ndarray:
    """Step 4 cut: the smallest K* >= K whose cut has K groups of >= min_size members;
    the K largest are kept and every other observation joins the kept group with the
    smallest mean dissimilarity. Falls back to the plain cut at K."""
    n = D.shape[0]
    for Kp in range(K, min(n, 20 * K) + 1):
        lab = cut_tree(Z, Kp)
        sizes = np.bincount(lab)
        big = np.flatnonzero(sizes >= min_size)
        if len(big) >= K:
            keep = big[np.argsort(-sizes[big])][:K]
            new = -np.ones(n, int)
            for g, b in enumerate(keep):
                new[lab == b] = g
            rest = np.flatnonzero(new < 0)
            if len(rest):
                M = np.stack([D[np.ix_(rest, np.flatnonzero(new == g))].mean(1) for g in range(K)], 1)
                new[rest] = M.argmin(1)
            return new
    return cut_tree(Z, K)


def min_cluster_size(n: int, K: int, f: float = 0.2) -> int:
    return max(2, int(np.ceil(f * n / K)))


# ----------------------------------------------------------------------------- Algorithm 1
@dataclass
class RobustUPGMAResult:
    labels: np.ndarray        # partition of all n observations, labels 0..K-1
    linkage: np.ndarray       # SciPy linkage matrix (dendrogram) of the inliers
    inliers: np.ndarray       # indices of the observations used for the dendrogram
    outliers: np.ndarray      # indices set aside in Step 3
    info: dict = field(default_factory=dict)

    def __repr__(self):
        return (f"RobustUPGMAResult(K={len(np.unique(self.labels))}, n={len(self.labels)}, "
                f"outliers={len(self.outliers)}, step1_changed={self.info.get('step1_n_changed')}, "
                f"step2_changed={self.info.get('step2_n_changed')})")


def robust_upgma(D, K: int, s: int = 64, k: float = 2.5, c: float = 3.0, f: float = 0.2,
                 triangle: bool = True, winsorise: bool = True, trim: bool = True,
                 size_cut: bool = True, seed: int = 0) -> RobustUPGMAResult:
    """Robust UPGMA (Algorithm 1).

    D : square dissimilarity matrix or condensed vector; K : number of clusters.
    Set triangle=False for non-metric dissimilarities. The switches triangle,
    winsorise, trim and size_cut turn the four steps on or off (ablation).
    """
    D = as_square(D)
    n = D.shape[0]
    info = {"step1_n_changed": 0, "step2_n_changed": 0}
    if triangle:
        D, info["step1_n_changed"] = triangle_repair(D, s, seed)
    if winsorise:
        d, w = winsorise_distances(D, k)
        D = squareform(d)
        info.update(lower=w["lower"], upper=w["upper"], step2_n_changed=w["n_low"] + w["n_high"])
    out = row_outliers(D, c) if trim else np.zeros(n, bool)
    inl = np.flatnonzero(~out)
    if len(inl) < K + 1:
        inl, out = np.arange(n), np.zeros(n, bool)
    Din = D[np.ix_(inl, inl)]
    Z = linkage(squareform(Din, checks=False), "average")
    m = min_cluster_size(len(inl), K, f)
    lab_in = size_constrained_cut(Z, Din, K, m) if size_cut else cut_tree(Z, K)
    lab = np.empty(n, int)
    lab[inl] = lab_in
    for i in np.flatnonzero(out):
        lab[i] = int(np.argmin([D[i, inl[lab_in == g]].mean() for g in range(K)]))
    info["min_cluster_size"] = m
    return RobustUPGMAResult(lab, Z, inl, np.flatnonzero(out), info)


def w_upgma(D, K: int, k: float = 2.5) -> np.ndarray:
    """W-UPGMA: Step 2 alone followed by UPGMA cut at K (first-submission method)."""
    d, _ = winsorise_distances(D, k)
    return cut_tree(linkage(d, "average"), K)


def upgma(D, K: int) -> np.ndarray:
    """Ordinary UPGMA (average linkage) cut at K clusters."""
    return cut_tree(linkage(as_condensed(D), "average"), K)
