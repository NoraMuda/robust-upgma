"""Minimal example: three groups, 10% of the dissimilarities corrupted.

    python examples/example.py
"""
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score

from robust_upgma import robust_upgma, upgma

rng = np.random.default_rng(1)
centres = np.array([[0, 0], [4, 0], [2, 3.5]])
y = np.repeat([0, 1, 2], 50)
X = centres[y] + rng.normal(size=(150, 2))

d = pdist(X)                                   # condensed dissimilarities, 11,175 entries
bad = rng.random(d.size) < 0.10                # corrupt 10% of them
d[bad] = rng.choice([0.0, 10 * d.max()], size=bad.sum())

res = robust_upgma(d, K=3)
print(res)
print("ARI, Robust UPGMA:", round(adjusted_rand_score(y, res.labels), 3))
print("ARI, UPGMA       :", round(adjusted_rand_score(y, upgma(d, 3)), 3))
print("observations set aside in Step 3:", res.outliers)

# Dendrogram of the inliers (optional):
# from scipy.cluster.hierarchy import dendrogram; import matplotlib.pyplot as plt
# dendrogram(res.linkage, no_labels=True); plt.show()
