# Robust UPGMA

Average-linkage hierarchical clustering (UPGMA) for **corrupted dissimilarity matrices**, with
implementations in **Python** and **R** and the code that reproduces every table and figure of

> Muda, N. Robust unweighted pair group method with arithmetic mean clustering for corrupted
> dissimilarities using triangle-inequality repair and median absolute deviation bounds.
> *Scientific Reports* (under review).

## The method (Algorithm 1)

Given an n x n dissimilarity matrix D and the number of clusters K:

1. **Triangle repair.** With a reference set R of min(s, n) observations, each d_ij is clipped to
   [med_r |d_ir - d_jr|, med_r (d_ir + d_jr)]. For uncorrupted metric data nothing changes.
2. **MADn winsorisation.** Dissimilarities are clipped to med ± k·MADn.
3. **MADn trimming.** Observations whose median dissimilarity exceeds med + c·MADn of these
   medians are set aside.
4. **UPGMA and a size-constrained cut.** The dendrogram of the remaining observations is cut at
   the smallest level that gives K groups of at least ceil(f·n/K) members; all other observations
   join the nearest cluster (smallest mean dissimilarity).

Defaults: s = 64, k = 2.5, c = 3, f = 0.2. Use `triangle = False` for non-metric dissimilarities.

## Quick start

**Python** (3.10 or later)

```bash
cd python
pip install .                 # the method only needs numpy and scipy
python examples/example.py
```

```python
from scipy.spatial.distance import pdist
from robust_upgma import robust_upgma

res = robust_upgma(pdist(X), K=3)   # square matrix or condensed vector
res.labels      # partition of all observations
res.linkage     # SciPy linkage matrix of the inliers (scipy.cluster.hierarchy.dendrogram)
res.outliers    # observations set aside in Step 3
res.info        # bounds and number of entries changed in Steps 1 and 2
```

**R** (4.1 or later)

```r
setwd("path/to/robust-upgma/R")
source("robust_upgma.R")
res <- robust_upgma(as.matrix(dist(X)), K = 3)
res$labels; plot(res$tree); res$outliers
```

The R code gives the same partitions as the Python code for the same input (see below).
Details: [`R/README_R.md`](R/README_R.md).

## Reproducing the paper

Python is the reference implementation: the numbers in the paper come from it.

```bash
cd python
pip install -r requirements.txt && pip install .
cd experiments
python simulation.py      # 17,100 simulated data sets (Tables 4-6, Figs 1-2)
python development.py --round 1 && python development.py --round 2   # development study (Table S5)
python benchmarks.py      # FCPS and real data (Table 8, Figs 3-4)
python null_structure.py  # data without cluster structure (Table 7)
python analysis.py        # all tables and figures -> results/python/tables, results/python/figures
```

or simply `bash run_all.sh`. The full simulation takes a few hours on 4 cores
(`--reps 10` gives a quick check); all scripts accept `--workers`.

In R, open `R/robust_upgma_all_in_one.R`, choose what to run at the top, and click *Source*,
or run the separate scripts listed in `R/README_R.md`. `R/verify_against_python.R` checks that
R and Python give identical partitions on 47 test cases.

### Seeds

| Study | Seed (NumPy `SeedSequence`) | Size |
|---|---|---|
| Development study | `[99001, condition, replicate]` | 171 conditions x 10 = 1,710 data sets |
| Confirmatory simulation | `[20260523, condition, replicate]` | 171 conditions x 100 = 17,100 data sets |
| Benchmarks | `[20260524, data set, setting, realisation]` | 15 data sets x 4 settings, 30 realisations |
| Null data | `[20260525, replicate]` | 100 replicates |

The development and confirmatory simulations use different seeds, so the method was evaluated on
data that were not used to choose it. R uses its own generator (`set.seed`), so R results for
simulated and corrupted data agree with Python within Monte Carlo error, not digit for digit;
clean-data results are identical.

## Repository layout

```
python/robust_upgma/     the method (core.py)
python/experiments/      simulation, development, benchmarks, null data, analysis
python/tests/            pytest tests  (cd python && pytest)
python/examples/         minimal example
R/                       R implementation, experiments, single-file script, fixtures for the R-Python check
data/                    FCPS and real data sets (see data/README.md)
results/python/          trial-level results, tables and figures used in the paper
results/R/               the same from the R implementation
```

## Citation

See [`CITATION.cff`](CITATION.cff). Please cite the paper if you use this code.

## Licence

Code: MIT (see [`LICENSE`](LICENSE)). The data in `data/` keep their original licences
(see [`data/README.md`](data/README.md)).
