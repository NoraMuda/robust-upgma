# Robust UPGMA in R

R version of the method and of every experiment in the manuscript. Run all scripts **from this `R/` folder**
(in RStudio: *Session → Set Working Directory → To Source File Location*, or `setwd(".../robust-upgma/R")`).
If you prefer a single script, open `robust_upgma_all_in_one.R`, choose what to run at the top, and click *Source*.

## 1. Install (once)

```r
install.packages(c("cluster", "dbscan", "diptest", "clue", "matrixStats"))
```
`parallel` and `stats` come with R. Tested with R 4.3.3, cluster 2.1.6, dbscan 1.1-12, diptest 0.77-0,
clue 0.3-65, matrixStats 1.2.0.

## 2. Use the method on your own data

```r
source("robust_upgma.R")
D   <- as.matrix(dist(X))            # or any metric dissimilarity matrix
res <- robust_upgma(D, K = 3)        # defaults: s = 64, k = 2.5, c = 3, f = 0.2
res$labels                           # partition of all observations
plot(res$tree)                       # average-linkage dendrogram of the inliers (res$inliers)
res$outliers                         # observations flagged in Step 3
res$info                             # MADn bounds, number of entries changed in Steps 1-2
```
Use `triangle = FALSE` for non-metric dissimilarities. `w_upgma(D, K)` gives Step 2 alone (W-UPGMA).
See `example.R`.


## 3. Files

| File | Contents |
|---|---|
| `robust_upgma.R` | Algorithm 1 (Steps 1–4), W-UPGMA |
| `scipy_compat.R` | SciPy-compatible average/complete linkage and `fcluster`; NumPy-compatible means |
| `experiment_utils.R` | simulation design, corruption models (Table 2), competitors, ARI, matched accuracy, silhouette |
| `data_loader.R` | FCPS and real data (from `../data/`), preprocessing |
| `run_simulation.R`, `run_benchmarks.R`, `run_null.R` | experiments |
| `analysis.R` | tables, figures, comparison with the manuscript |
| `verify_against_python.R`, `fixtures/` | R-versus-Python check |
| `example.R` | minimal example |
| `robust_upgma_all_in_one.R` | everything above in one file (generated from the other files): choose what to run, click Source |
