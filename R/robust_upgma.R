# =============================================================================
# Robust UPGMA in R  (Algorithm 1 of the manuscript)
#
#   Step 1  triangle repair of individual dissimilarities
#   Step 2  MADn winsorisation of all dissimilarities          (= W-UPGMA alone)
#   Step 3  MADn trimming of outlying observations
#   Step 4  UPGMA (average linkage) + size-constrained cut
#
# This is a line-by-line port of the Python reference implementation
# (python/robust_upgma/core.py).  Given the same dissimilarities and the same reference
# set R of Step 1, it returns the same partition (see verify_against_python.R).
# Requires: stats (base R); matrixStats is used for speed if installed.
# Hierarchical clustering uses the SciPy-compatible port in scipy_compat.R, so
# that ties among winsorised dissimilarities are broken exactly as in Python.
# =============================================================================
if (!exists("scipy_linkage")) source("scipy_compat.R")   # run scripts from the R/ folder

MADN_CONSTANT <- 1.4826

madn <- function(v) {
  v <- as.numeric(v)
  MADN_CONSTANT * median(abs(v - median(v)))
}

as_dist_matrix <- function(D) {
  if (inherits(D, "dist")) return(as.matrix(D))
  D <- as.matrix(D)
  if (nrow(D) != ncol(D)) stop("D must be a square dissimilarity matrix or a 'dist' object")
  if (max(abs(D - t(D))) > 1e-10) stop("D must be symmetric")
  D
}

.row_medians <- function(M) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowMedians(M)
  else apply(M, 1, median)
}

# ---- Step 2: MADn winsorisation (also the whole of W-UPGMA) -----------------
winsorise_distances <- function(D, k = 2.5) {
  D <- as_dist_matrix(D)
  d <- D[upper.tri(D)]
  med <- median(d); s <- madn(d)
  L <- max(0, med - k * s); U <- med + k * s
  Dw <- pmin(pmax(D, L), U)
  diag(Dw) <- 0
  list(D = Dw, lower = L, upper = U, median = med, madn = s,
       n_low = sum(d < L), n_high = sum(d > U))
}

# ---- Step 1: triangle repair -------------------------------------------------
# ref: indices (1-based) of the reference set R.  If NULL, min(s, n) indices are
# drawn with set.seed(seed).  (Python draws them with numpy; to reproduce the
# Python partitions exactly, pass the Python indices + 1 as `ref`.)
reference_set <- function(n, s = 64, seed = 0) {
  if (n <= s) return(seq_len(n))
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)
  R <- sort(sample.int(n, s))
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  R
}

triangle_bounds <- function(D, ref) {
  n <- nrow(D)
  A <- D[, ref, drop = FALSE]                    # n x s
  L <- matrix(0, n, n); U <- matrix(0, n, n)
  for (i in seq_len(n - 1)) {
    j <- (i + 1):n
    B <- A[j, , drop = FALSE]
    ai <- matrix(A[i, ], nrow = length(j), ncol = ncol(A), byrow = TRUE)
    L[i, j] <- .row_medians(abs(B - ai))         # med_r |d_ir - d_jr|
    U[i, j] <- .row_medians(B + ai)              # med_r (d_ir + d_jr)
  }
  list(L = L + t(L), U = U + t(U))
}

triangle_repair <- function(D, s = 64, seed = 0, ref = NULL) {
  D <- as_dist_matrix(D)
  if (is.null(ref)) ref <- reference_set(nrow(D), s, seed)
  b <- triangle_bounds(D, ref)
  up <- upper.tri(D)
  n_low <- sum((D < b$L)[up]); n_high <- sum((D > b$U)[up])
  D2 <- pmin(pmax(D, b$L), b$U)
  list(D = D2, n_changed = n_low + n_high, ref = ref)
}

# ---- Step 3: MADn trimming of observations -----------------------------------
row_outliers <- function(D, c = 3) {
  n <- nrow(D)
  off <- D; diag(off) <- NA
  r <- apply(off, 1, median, na.rm = TRUE)       # rho_i = med_{j != i} d_ij
  r > median(r) + c * madn(r)
}

# ---- Step 4: UPGMA + size-constrained cut ------------------------------------
upgma_tree <- function(D) scipy_linkage(D, "average")
cut_k <- function(tree, K) scipy_fcluster_maxclust(tree, K)

size_constrained_cut <- function(tree, D, K, min_size) {
  n <- nrow(D)
  for (Kp in K:min(n, 20 * K)) {
    lab <- cut_k(tree, Kp)
    sizes <- tabulate(lab, nbins = Kp)
    big <- which(sizes >= min_size)
    if (length(big) >= K) {
      keep <- big[order(-sizes[big])][seq_len(K)]
      new <- match(lab, keep)                    # 1..K or NA
      rest <- which(is.na(new))
      if (length(rest)) {
        M <- sapply(seq_len(K), function(g)
          np_rowmeans(D[rest, which(new == g), drop = FALSE]))
        M <- matrix(M, nrow = length(rest))
        new[rest] <- max.col(-M, ties.method = "first")
      }
      return(new)
    }
  }
  cut_k(tree, K)
}

# ---- Algorithm 1 --------------------------------------------------------------
# Returns a list with: labels (1..K), tree (hclust object of the inliers, for
# plot()), inliers,
# outliers, and diagnostics.  triangle = FALSE for non-metric dissimilarities.
robust_upgma <- function(D, K, s = 64, k = 2.5, c = 3, f = 0.2, triangle = TRUE,
                         seed = 0, ref = NULL, winsorise = TRUE, trim = TRUE,
                         size_cut = TRUE) {
  D <- as_dist_matrix(D)
  n <- nrow(D)
  info <- list(step1_n_changed = 0L)
  if (triangle) {
    tr <- triangle_repair(D, s, seed, ref)
    D <- tr$D; info$step1_n_changed <- tr$n_changed
  }
  if (winsorise) {
    w <- winsorise_distances(D, k)
    D <- w$D
    info$lower <- w$lower; info$upper <- w$upper
    info$step2_n_changed <- w$n_low + w$n_high
  }
  out <- if (trim) row_outliers(D, c) else rep(FALSE, n)
  inl <- which(!out)
  if (length(inl) < K + 1) { inl <- seq_len(n); out <- rep(FALSE, n) }
  Din <- D[inl, inl, drop = FALSE]
  tree <- upgma_tree(Din)
  m <- max(2, ceiling(f * length(inl) / K))
  lab_in <- if (size_cut) size_constrained_cut(tree, Din, K, m) else cut_k(tree, K)
  lab <- integer(n)
  lab[inl] <- lab_in
  for (i in which(out)) {
    lab[i] <- which.min(sapply(seq_len(K), function(g) np_mean(D[i, inl[lab_in == g]])))
  }
  info$min_cluster_size <- m
  list(labels = lab, tree = scipy_to_hclust(tree), inliers = inl, outliers = which(out), info = info)
}

# W-UPGMA: Step 2 only (the method of the first submission)
w_upgma <- function(D, K, k = 2.5) {
  w <- winsorise_distances(D, k)
  cut_k(upgma_tree(w$D), K)
}
