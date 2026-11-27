"""Robust UPGMA: average-linkage hierarchical clustering for corrupted dissimilarities.

>>> from robust_upgma import robust_upgma
>>> res = robust_upgma(D, K=3)       # D: square dissimilarity matrix or condensed vector
>>> res.labels, res.linkage, res.outliers, res.info
"""
from .core import (MADN_CONSTANT, RobustUPGMAResult, cut_tree, madn, min_cluster_size,
                   robust_upgma, row_outliers, size_constrained_cut, triangle_bounds,
                   triangle_repair, upgma, w_upgma, winsorise_distances)

__version__ = "2.0.0"
__all__ = ["robust_upgma", "w_upgma", "upgma", "RobustUPGMAResult", "triangle_bounds", "triangle_repair",
           "winsorise_distances", "row_outliers", "size_constrained_cut", "min_cluster_size", "cut_tree",
           "madn", "MADN_CONSTANT"]
