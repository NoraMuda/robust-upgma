# =============================================================================
#  ROBUST UPGMA  -  everything in one R file
#  (generated from the R files in this folder; edit those, not this file)
#
#    Part 0  settings (EDIT THIS PART)
#    Part 1  the method (Algorithm 1) and SciPy/NumPy-compatible helpers
#    Part 2  simulation design, corruption models, competitors, metrics, data
#    Part 3  experiments: verification, example, simulation, benchmarks,
#            null data, tables/figures + comparison with the manuscript
#    Part 4  run
#
#  How to use (RStudio): open this file from the R/ folder of the repository,
#  choose what to run in Part 0, then click "Source". The repository folder is
#  found automatically. Results are written to <repository>/results/R/
# =============================================================================

# ============================ Part 0: SETTINGS ==============================
PROJECT_DIR <- ""   # leave empty: found automatically (or type the repository path, using "/")

RUN_VERIFY     <- TRUE    # check that R = Python on 47 test cases      (~30 s)
RUN_EXAMPLE    <- TRUE    # small example with a dendrogram plot        (~2 s)
RUN_SIMULATION <- FALSE   # Tables 4-6, Figs 1-2                        (~30-60 min)
RUN_BENCHMARKS <- FALSE   # Table 8, Figs 3-4                           (~15 min with 10 reps)
RUN_NULL       <- FALSE   # Table 7                                     (~5 min)
RUN_ANALYSIS   <- FALSE   # tables, figures, comparison with manuscript (~1 min)

SIM_REPS   <- 100         # replicates per condition (manuscript: 100)
BENCH_REPS <- 30          # realisations per corrupted setting (manuscript: 30)
NULL_REPS  <- 100         # replicates of the null study (manuscript: 100)
CORES      <- max(1L, parallel::detectCores() - 1L)

# ---- packages (installed automatically if missing) --------------------------
for (pkg in c("cluster", "dbscan", "diptest", "clue", "matrixStats")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}


# ============================ Part 1-2: method, design, utilities ===========
# =============================================================================
# SciPy-compatible hierarchical clustering in R
#
# After MADn winsorisation many dissimilarities are exactly equal (all clipped
# values equal U), so the order in which tied pairs are merged matters.
# stats::hclust and SciPy break such ties differently.  To reproduce the
# Python results exactly, this file ports
#   * scipy.cluster.hierarchy.linkage(method = "average" / "complete"), which
#     uses the nearest-neighbour-chain algorithm (scipy/_hierarchy.pyx: nn_chain),
#   * scipy.cluster.hierarchy.fcluster(Z, K, criterion = "maxclust"), including
#     its cluster numbering (cluster_monocrit) and its behaviour with tied
#     heights (it may return fewer than K clusters).
# For single linkage stats::hclust is used: the clusters at a given height are
# unique for single linkage, so ties do not matter there.
# =============================================================================

# Returns list(Z = (n-1) x 2 matrix of SciPy node ids (0-based, label()-ed),
#              h = merge heights sorted ascending, n = n)
scipy_linkage <- function(D, method = c("average", "complete")) {
  method <- match.arg(method)
  n <- nrow(D)
  Dm <- D
  size <- rep(1L, n)
  chain <- integer(n); clen <- 0L
  mx <- integer(n - 1); my <- integer(n - 1); mh <- numeric(n - 1)
  for (k in seq_len(n - 1)) {
    if (clen == 0L) { clen <- 1L; chain[1] <- which(size > 0L)[1] }
    repeat {
      x <- chain[clen]
      cand <- which(size > 0L); cand <- cand[cand != x]
      dv <- Dm[x, cand]
      if (clen > 1L) { y <- chain[clen - 1L]; cur <- Dm[x, y] } else { y <- NA_integer_; cur <- Inf }
      mn <- min(dv)
      if (mn < cur) { y <- cand[which(dv == mn)[1]]; cur <- mn }   # first index wins ties
      if (clen > 1L && y == chain[clen - 1L]) break
      clen <- clen + 1L; chain[clen] <- y
    }
    clen <- clen - 2L
    if (x > y) { tmp <- x; x <- y; y <- tmp }
    nx <- size[x]; ny <- size[y]
    mx[k] <- x; my[k] <- y; mh[k] <- cur
    size[x] <- 0L; size[y] <- nx + ny
    others <- which(size > 0L); others <- others[others != y]
    if (length(others)) {
      newd <- if (method == "average") (nx * Dm[others, x] + ny * Dm[others, y]) / (nx + ny)
              else pmax(Dm[others, x], Dm[others, y])
      Dm[others, y] <- newd; Dm[y, others] <- newd
    }
  }
  o <- order(mh, method = "radix")                 # stable, as numpy mergesort
  mx <- mx[o]; my <- my[o]; mh <- mh[o]
  # scipy label(): replace slot indices by union-find cluster ids
  parent <- seq_len(2 * n - 1)                     # 1-based ids; SciPy id = R id - 1
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  Z <- matrix(0L, n - 1, 2)
  for (i in seq_len(n - 1)) {
    rx <- find(mx[i]); ry <- find(my[i])
    Z[i, ] <- c(min(rx, ry), max(rx, ry)) - 1L     # 0-based like SciPy
    parent[rx] <- n + i; parent[ry] <- n + i
  }
  list(Z = Z, h = mh, n = n)
}

# Port of fcluster(Z, K, "maxclust"): smallest threshold t with at most K
# clusters; numbering as SciPy's cluster_monocrit depth-first traversal.
scipy_fcluster_maxclust <- function(L, K) {
  Z <- L$Z; h <- L$h; n <- L$n
  if (K >= n) return(seq_len(n))
  nclust <- n - findInterval(h, h)                 # clusters if cut at each height
  cutoff <- h[which(nclust <= K)[1]]
  # max height in each subtree (monotone linkage: equals the node height)
  T <- integer(n); n_cluster <- 0L; leader <- -1L
  visited <- logical(2 * n - 1)
  stack <- integer(n); k <- 1L; stack[1] <- 2L * n - 2L   # root id (0-based)
  while (k >= 1L) {
    root <- stack[k] - n + 1L                       # 1-based row of Z
    lc <- Z[root, 1]; rc <- Z[root, 2]
    if (leader == -1L && h[root] <= cutoff) { leader <- root; n_cluster <- n_cluster + 1L }
    if (lc >= n && !visited[lc + 1L]) { visited[lc + 1L] <- TRUE; k <- k + 1L; stack[k] <- lc; next }
    if (rc >= n && !visited[rc + 1L]) { visited[rc + 1L] <- TRUE; k <- k + 1L; stack[k] <- rc; next }
    if (lc < n) { if (leader == -1L) n_cluster <- n_cluster + 1L; T[lc + 1L] <- n_cluster }
    if (rc < n) { if (leader == -1L) n_cluster <- n_cluster + 1L; T[rc + 1L] <- n_cluster }
    if (leader == root) leader <- -1L
    k <- k - 1L
  }
  T
}

# Single linkage: components at the smallest height giving at most K clusters
single_linkage_cut <- function(D, K) {
  hc <- stats::hclust(stats::as.dist(D), method = "single")
  n <- nrow(D)
  nclust <- n - findInterval(hc$height, hc$height)
  cutoff <- hc$height[which(nclust <= K)[1]]
  stats::cutree(hc, h = cutoff + 0)                # same partition as SciPy (unique for single)
}

# Convert to an hclust object (for plotting the dendrogram)
scipy_to_hclust <- function(L) {
  n <- L$n
  merge <- matrix(0L, n - 1, 2)
  conv <- function(id) if (id < n) -(id + 1L) else id - n + 1L
  for (i in seq_len(n - 1)) merge[i, ] <- c(conv(L$Z[i, 1]), conv(L$Z[i, 2]))
  hc <- list(merge = merge, height = L$h, order = seq_len(n), labels = NULL,
             method = "average", call = match.call(), dist.method = NULL)
  class(hc) <- "hclust"
  hc$order <- stats::order.dendrogram(stats::as.dendrogram(hc))
  hc
}

# NumPy-compatible means.  numpy.mean uses pairwise summation (blocks of 128,
# eight accumulators), whereas R's mean()/rowMeans() accumulate in long double.
# The results differ in the last bit, which decides exact ties when an
# observation is equally far from two clusters (common after winsorisation).
np_rowsums <- function(M) {          # pairwise sum of each row, as numpy
  n <- ncol(M)
  if (n < 8) { s <- numeric(nrow(M)); for (j in seq_len(n)) s <- s + M[, j]; return(s) }
  if (n <= 128) {
    r <- M[, 1:8, drop = FALSE]
    m8 <- n - n %% 8
    if (m8 > 8) for (i in seq(9, m8, by = 8)) r <- r + M[, i:(i + 7), drop = FALSE]
    res <- ((r[, 1] + r[, 2]) + (r[, 3] + r[, 4])) + ((r[, 5] + r[, 6]) + (r[, 7] + r[, 8]))
    if (m8 < n) for (j in (m8 + 1):n) res <- res + M[, j]
    return(res)
  }
  n2 <- n %/% 2; n2 <- n2 - n2 %% 8
  np_rowsums(M[, 1:n2, drop = FALSE]) + np_rowsums(M[, (n2 + 1):n, drop = FALSE])
}
np_rowmeans <- function(M) np_rowsums(M) / ncol(M)
np_mean <- function(x) np_rowsums(matrix(x, nrow = 1)) / length(x)

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

# =============================================================================
# Shared pieces for the R experiments: data generation, corruption models,
# competing methods and evaluation metrics.  Mirrors experiments/common.py and
# experiments/simulation.py of the Python code.
# =============================================================================

# ---- simulation design -------------------------------------------------------
STRUCTURES <- list(
  B3 = list(K = 3, p = 5,  sizes = c(50, 50, 50),       sd = c(1, 1, 1)),
  U4 = list(K = 4, p = 5,  sizes = c(100, 50, 30, 20),  sd = c(1, 0.5, 1.5, 0.75)),
  H3 = list(K = 3, p = 50, sizes = c(50, 50, 50),       sd = c(1, 1, 1))
)
SEPARATIONS <- c(well = 6, moderate = 4, overlap = 2.5)
CORRUPTIONS <- c("salt_pepper", "amplification", "mixed", "gaussian", "row", "structured")
LEVELS <- c(0.05, 0.10, 0.20)

conditions <- function() {
  out <- list()
  for (s in names(STRUCTURES)) for (sep in names(SEPARATIONS)) {
    out[[length(out) + 1]] <- list(structure = s, separation = sep, corruption = "none", level = 0)
    for (cr in CORRUPTIONS) for (e in LEVELS)
      out[[length(out) + 1]] <- list(structure = s, separation = sep, corruption = cr, level = e)
  }
  out   # 171 conditions, same order as Python
}

simplex_centres <- function(K, p, delta) {
  C <- matrix(0, K, p)
  for (g in seq_len(K)) C[g, g] <- delta / sqrt(2)
  C
}

generate <- function(struct, sep) {
  cfg <- STRUCTURES[[struct]]
  C <- simplex_centres(cfg$K, cfg$p, SEPARATIONS[[sep]])
  X <- do.call(rbind, lapply(seq_len(cfg$K), function(g)
    sweep(cfg$sd[g] * matrix(rnorm(cfg$sizes[g] * cfg$p), cfg$sizes[g], cfg$p), 2, C[g, ], "+")))
  list(X = X, y = rep(seq_len(cfg$K), cfg$sizes))
}

# ---- corruption models (Table 2) --------------------------------------------
# D: full distance matrix; returns a corrupted, symmetric matrix
corrupt_distances <- function(D, kind, eps, y = NULL) {
  n <- nrow(D)
  up <- which(upper.tri(D))                       # the m pairs i < j
  d <- D[up]
  m <- length(d)
  if (kind == "structured") {
    ij <- arrayInd(up, dim(D))
    elig <- which(y[ij[, 1]] == 1 | y[ij[, 2]] == 1)
    n_sel <- min(round(eps * m), length(elig))
    S <- elig[sample.int(length(elig), n_sel)]
    d[S] <- d[S] * runif(n_sel, 2, 5)
  } else {
    n_sel <- round(eps * m)
    S <- sample.int(m, n_sel)
    if (kind == "salt_pepper") {
      q <- quantile(d, c(0.01, 0.99), type = 7, names = FALSE)
      salt <- runif(n_sel) < 0.5
      d[S] <- ifelse(salt, q[2] * runif(n_sel, 1.5, 3), q[1] * runif(n_sel, 0, 1))
    } else if (kind == "amplification") {
      d[S] <- d[S] * runif(n_sel, 2, 5)
    } else if (kind == "mixed") {
      F <- runif(n_sel, 2, 5); upw <- runif(n_sel) < 0.5
      d[S] <- ifelse(upw, d[S] * F, d[S] / F)
    } else if (kind == "gaussian") {
      d[S] <- abs(d[S] + rnorm(n_sel, 0, 0.5 * sd_pop(d)))
    } else stop(kind)
  }
  Dc <- matrix(0, n, n); Dc[up] <- d; Dc + t(Dc)
}

sd_pop <- function(x) sqrt(mean((x - mean(x))^2))   # numpy's std (ddof = 0)

contaminate_rows <- function(X, eps) {
  n <- nrow(X); p <- ncol(X)
  n_out <- round(eps * n)
  idx <- sample.int(n, n_out)
  u <- matrix(rnorm(n_out * p), n_out, p); u <- u / sqrt(rowSums(u^2))
  delta <- runif(n_out, 5, 10) * sd_pop(X) * sqrt(p)
  X[idx, ] <- X[idx, ] + delta * u
  mask <- rep(FALSE, n); mask[idx] <- TRUE
  list(X = X, outlier = mask)
}

# ---- competing methods --------------------------------------------------------
METHODS <- c("R-UPGMA", "W-UPGMA", "UPGMA", "UPGMA+cut", "Single", "Complete",
             "Complete+cut", "PAM", "HDBSCAN")
ABLATIONS <- c("w/o Step 1 (triangle)", "w/o Step 2 (MADn bounds)",
               "w/o Step 3 (trimming)", "w/o Step 4 (cut)")

run_methods <- function(D, K, methods = METHODS, ref = NULL, hdb_min_size = 10) {
  n <- nrow(D)
  out <- list()
  cutm <- function(method) cut_k(scipy_linkage(D, method), K)
  cut_size <- function(method) size_constrained_cut(scipy_linkage(D, method), D, K, max(2, ceiling(0.2 * n / K)))
  for (m in methods) {
    out[[m]] <- switch(m,
      "R-UPGMA"      = robust_upgma(D, K, ref = ref)$labels,
      "W-UPGMA"      = w_upgma(D, K),
      "UPGMA"        = cutm("average"),
      "UPGMA+cut"    = cut_size("average"),
      "Single"       = single_linkage_cut(D, K),
      "Complete"     = cutm("complete"),
      "Complete+cut" = cut_size("complete"),
      "PAM"          = cluster::pam(stats::as.dist(D), K, diss = TRUE)$clustering,
      "HDBSCAN"      = dbscan::hdbscan(stats::as.dist(D), minPts = hdb_min_size)$cluster,  # 0 = noise group
      "w/o Step 1 (triangle)"    = robust_upgma(D, K, triangle = FALSE)$labels,
      "w/o Step 2 (MADn bounds)" = robust_upgma(D, K, winsorise = FALSE, ref = ref)$labels,
      "w/o Step 3 (trimming)"    = robust_upgma(D, K, trim = FALSE, ref = ref)$labels,
      "w/o Step 4 (cut)"         = robust_upgma(D, K, size_cut = FALSE, ref = ref)$labels,
      stop("unknown method ", m))
  }
  out
}

# ---- metrics ------------------------------------------------------------------
adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  comb2 <- function(v) sum(v * (v - 1) / 2)
  a <- comb2(tab); b <- comb2(rowSums(tab)); c <- comb2(colSums(tab))
  n2 <- length(x) * (length(x) - 1) / 2
  expected <- b * c / n2; maximum <- (b + c) / 2
  if (maximum == expected) return(1)     # both partitions trivial (as scikit-learn)
  (a - expected) / (maximum - expected)
}

matched_accuracy <- function(y, lab) {
  tab <- unclass(table(y, lab))
  k <- max(dim(tab))
  sq <- matrix(0, k, k); sq[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  a <- clue::solve_LSAP(sq, maximum = TRUE)
  sum(sq[cbind(seq_len(k), as.integer(a))]) / length(y)
}

mean_silhouette <- function(lab, D) {
  if (length(unique(lab)) < 2 || length(unique(lab)) >= length(lab)) return(NA_real_)
  s <- cluster::silhouette(as.integer(factor(lab)), dmatrix = D)
  mean(s[, "sil_width"])
}

evaluate <- function(y, lab, D_eval = NULL, mask = NULL) {
  if (!is.null(mask)) {
    y <- y[mask]; lab <- lab[mask]
    if (!is.null(D_eval)) D_eval <- D_eval[mask, mask]
  }
  c(ari = adjusted_rand_index(y, lab), acc = matched_accuracy(y, lab),
    n_clusters = length(unique(lab)), largest = max(table(lab)) / length(lab),
    sil = if (is.null(D_eval)) NA_real_ else mean_silhouette(lab, D_eval))
}

# Data sets and settings of the benchmark study (shared by run_benchmarks.R and run_null.R)
# Paths are relative to the R/ folder: data files live in ../data/.

FCPS <- c("Hepta", "Tetra", "Lsun3D", "Chainlink", "Atom", "EngyTime", "TwoDiamonds", "WingNut", "Target")
REAL <- c("Iris" = "iris", "Wine" = "wine_recognition", "Breast Cancer" = "breast_cancer_wisconsin_diagnostic",
          "E. coli" = "ecoli", "Penguins" = "penguins", "MAGIC" = "magic")
SETTINGS <- list(clean = NULL, salt_pepper10 = c("salt_pepper", 0.10), mixed10 = c("mixed", 0.10),
                 row05 = c("row", 0.05))
MAX_N <- 2000

load_data <- function(name) {
  if (name %in% c(FCPS, "GolfBall")) {
    e <- new.env(); load(file.path("..", "data", "fcps", paste0(name, ".rda")), envir = e)
    obj <- get(name, envir = e)
    return(list(X = unname(as.matrix(obj$Data)), y = as.integer(factor(obj$Cls)), source = "FCPS"))
  }
  df <- read.delim(gzfile(file.path("..", "data", "pmlb", paste0(REAL[[name]], ".tsv.gz"))))
  if (name == "Penguins") {
    df <- df[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g", "target")]
    df <- df[complete.cases(df), ]
  }
  y <- as.integer(factor(df$target)); df$target <- NULL
  X <- as.matrix(df)
  X <- X[, apply(X, 2, sd) > 0, drop = FALSE]                     # drop constant features
  X <- scale(X, center = TRUE, scale = apply(X, 2, sd_pop))       # z-score (ddof = 0, as Python)
  list(X = unname(X), y = y, source = "Real")
}

# ============================ Part 3: experiments ===========================
run_verification <- function() {

  cases <- read.csv("fixtures/cases.csv", stringsAsFactors = FALSE)
  py <- read.csv(gzfile("fixtures/python_labels.csv.gz"), stringsAsFactors = FALSE)
  parse_int <- function(s) as.integer(strsplit(s, " ")[[1]])

  res <- list()
  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    d <- scan(gzfile(sprintf("fixtures/%s_d.csv.gz", cs$case)), quiet = TRUE)
    D <- as.matrix(structure(d, Size = cs$n, class = "dist"))   # same pair order as Python
    ref <- parse_int(cs$ref)
    # Steps 1-2 numerically
    D12 <- winsorise_distances(triangle_repair(D, ref = ref)$D)$D
    chk <- abs(sum(D12[upper.tri(D12)]) - cs$checksum_step12) / cs$checksum_step12
    labs <- run_methods(D, cs$K, methods = c(METHODS, ABLATIONS), ref = ref)
    for (m in names(labs)) {
      p <- parse_int(py$labels[py$case == cs$case & py$method == m])
      res[[length(res) + 1]] <- data.frame(case = cs$case, method = m,
                                            agreement_ARI = adjusted_rand_index(p, labs[[m]]),
                                            step12_rel_diff = chk)
    }
  }
  res <- do.call(rbind, res)
  summ <- aggregate(agreement_ARI ~ method, res, function(a) c(identical = sum(a > 1 - 1e-12), cases = length(a),
                                                                min_ARI = round(min(a), 3)))
  cat("\nSteps 1-2 numerical agreement: max relative difference of the repaired",
      "distances =", signif(max(res$step12_rel_diff), 3), "\n\n")
  print(do.call(data.frame, summ), row.names = FALSE)
  write.csv(res, "../results/R/verification_results.csv", row.names = FALSE)
  cat("\nDetails written to ../results/R/verification_results.csv\n")
  invisible(NULL)
}

run_example <- function() {
  set.seed(1)
  X <- rbind(matrix(rnorm(50 * 5, 0), 50), matrix(rnorm(50 * 5, 4), 50), matrix(rnorm(50 * 5, 8), 50))
  D <- as.matrix(dist(X))

  # corrupt 10% of the distances: half far too large, half far too small
  up <- which(upper.tri(D)); S <- sample(up, round(0.1 * length(up)))
  D[S[1:(length(S) %/% 2)]] <- D[S[1:(length(S) %/% 2)]] * 4
  D[S[-(1:(length(S) %/% 2))]] <- D[S[-(1:(length(S) %/% 2))]] * 0.05
  D[lower.tri(D)] <- t(D)[lower.tri(D)]

  res <- robust_upgma(D, K = 3)
  print(table(res$labels, rep(1:3, each = 50)))            # recovered vs true groups
  print(unlist(res$info))
  upgma <- cut_k(upgma_tree(D), 3)
  print(table(upgma, rep(1:3, each = 50)))                 # ordinary UPGMA for comparison
  plot(res$tree, labels = FALSE, main = "Robust UPGMA dendrogram (inliers)")
  invisible(NULL)
}

run_simulation <- function(reps = 100L, cores = 2L) {
  dir.create("../results/R", recursive = TRUE, showWarnings = FALSE)

  one_trial <<- function(task) {
    cid <- task[1]; r <- task[2]
    cond <- conditions()[[cid]]
    set.seed(20260523L + 1000L * cid + r)          # R's own seed stream (see README_R.md)
    g <- generate(cond$structure, cond$separation)
    X <- g$X; y <- g$y; K <- STRUCTURES[[cond$structure]]$K
    D_clean <- as.matrix(dist(X))
    mask <- NULL
    if (cond$corruption == "none") {
      D <- D_clean
    } else if (cond$corruption == "row") {
      cr <- contaminate_rows(X, cond$level)
      D <- as.matrix(dist(cr$X)); mask <- !cr$outlier
    } else {
      D <- corrupt_distances(D_clean, cond$corruption, cond$level, y = y)
    }
    methods <- if (r < 30) c(METHODS, ABLATIONS) else METHODS
    labs <- run_methods(D, K, methods = methods)
    do.call(rbind, lapply(names(labs), function(m) {
      ev <- evaluate(y, labs[[m]], D_eval = D_clean, mask = mask)
      data.frame(cond = cid, structure = cond$structure, separation = cond$separation,
                 corruption = cond$corruption, level = cond$level, rep = r, method = m,
                 t(ev), stringsAsFactors = FALSE)
    }))
  }

  tasks <- as.matrix(expand.grid(r = 0:(reps - 1), cid = seq_along(conditions()))[, c("cid", "r")])
  tasks <- lapply(seq_len(nrow(tasks)), function(i) tasks[i, ])
  cat(length(tasks), "datasets on", cores, "cores ...\n")
  t0 <- Sys.time()
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed"), envir = globalenv())
  res <- parallel::parLapplyLB(cl, tasks, one_trial)
  parallel::stopCluster(cl)
  res <- do.call(rbind, res)
  write.csv(res, gzfile("../results/R/sim_trials.csv.gz"), row.names = FALSE)
  cat("done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  invisible(NULL)
}

run_benchmarks <- function(reps = 30L, cores = 2L) {
  dir.create("../results/R", recursive = TRUE, showWarnings = FALSE)


  one_run <<- function(task) {
    name <- task$name; st <- task$setting; r <- task$rep
    dat <- load_data(name)
    set.seed(20260524L + 100000L * task$di + 1000L * task$si + r)
    X <- dat$X; y <- dat$y
    if (length(y) > MAX_N) {
      idx <- unlist(lapply(sort(unique(y)), function(g) {
        w <- which(y == g); w[sample.int(length(w), round(MAX_N * mean(y == g)))] }))
      X <- X[idx, , drop = FALSE]; y <- y[idx]
    }
    K <- length(unique(y)); n <- length(y)
    D_clean <- as.matrix(dist(X)); mask <- NULL
    s <- SETTINGS[[st]]
    if (is.null(s)) D <- D_clean
    else if (s[1] == "row") { cr <- contaminate_rows(X, as.numeric(s[2])); D <- as.matrix(dist(cr$X)); mask <- !cr$outlier }
    else D <- corrupt_distances(D_clean, s[1], as.numeric(s[2]), y = y)
    labs <- run_methods(D, K)
    do.call(rbind, lapply(names(labs), function(m) {
      ev <- evaluate(y, labs[[m]], D_eval = D_clean, mask = mask)
      data.frame(dataset = name, source = dat$source, n = n, K = K, setting = st, rep = r, method = m,
                 t(ev), stringsAsFactors = FALSE)
    }))
  }

  names_all <- c(FCPS, names(REAL))
  tasks <- list()
  for (di in seq_along(names_all)) for (si in seq_along(SETTINGS)) {
    nm <- names_all[di]
    big <- nm %in% c("EngyTime", "MAGIC")
    nrep <- if (!is.null(SETTINGS[[si]]) || big) reps else 1L
    for (r in 0:(nrep - 1)) tasks[[length(tasks) + 1]] <- list(name = nm, setting = names(SETTINGS)[si],
                                                                 rep = r, di = di, si = si)
  }
  cat(length(tasks), "runs on", cores, "cores ...\n")
  t0 <- Sys.time()
  cl <- parallel::makeCluster(cores)
  parallel::clusterExport(cl, setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed"), envir = globalenv())
  res <- do.call(rbind, parallel::parLapplyLB(cl, tasks, one_run))
  parallel::stopCluster(cl)
  write.csv(res, gzfile("../results/R/bench_trials.csv.gz"), row.names = FALSE)
  cat("done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  invisible(NULL)
}

run_null <- function(reps = 100L) {
  dir.create("../results/R", recursive = TRUE, showWarnings = FALSE)

  golf <- load_data("GolfBall")$X; hepta <- load_data("Hepta")$X
  record <- function(name, setting, X, D, r) {
    d <- D[upper.tri(D)]
    p_raw <- diptest::dip.test(d)$p.value
    W <- winsorise_distances(D)$D
    p_win <- diptest::dip.test(W[upper.tri(W)])$p.value
    lab <- robust_upgma(D, 3)$labels
    sz <- sort(tabulate(lab, 3), decreasing = TRUE) / length(lab)
    data.frame(dataset = name, setting = setting, rep = r, dip_p = p_raw, dip_p_winsorised = p_win,
               largest = sz[1], smallest = sz[3], sil = mean_silhouette(lab, as.matrix(dist(X))))
  }
  rows <- list()
  for (r in 0:(reps - 1)) {
    set.seed(20260525L + r)
    sets <- list("NULL-Gauss" = matrix(rnorm(150 * 5), 150), "NULL-Unif" = matrix(runif(150 * 5), 150),
                 "B3-well (control)" = generate("B3", "well")$X,
                 "GolfBall" = golf[sample.int(nrow(golf), 150), ],
                 "Hepta (control)" = hepta[sample.int(nrow(hepta), 150), ])
    for (nm in names(sets)) {
      X <- sets[[nm]]; D <- as.matrix(dist(X))
      rows[[length(rows) + 1]] <- record(nm, "clean", X, D, r)
      rows[[length(rows) + 1]] <- record(nm, "salt_pepper10", X, corrupt_distances(D, "salt_pepper", 0.10), r)
    }
  }
  res <- do.call(rbind, rows)
  write.csv(res, gzfile("../results/R/null_trials.csv.gz"), row.names = FALSE)
  res$reject_raw <- res$dip_p < 0.05; res$reject_win <- res$dip_p_winsorised < 0.05
  print(cbind(aggregate(cbind(dip_p, largest, smallest, sil) ~ dataset + setting, res, median)[, 1:6],
              aggregate(cbind(reject_raw, reject_win) ~ dataset + setting, res, mean)[, 3:4]))
  invisible(NULL)
}

run_analysis <- function() {
  dir.create("../results/R/tables", recursive = TRUE, showWarnings = FALSE)
  dir.create("../results/R/figures", recursive = TRUE, showWarnings = FALSE)
  TAB <- "../results/R/tables"; FIG <- "../results/R/figures"
  MAIN <- c("R-UPGMA", "W-UPGMA", "UPGMA", "UPGMA+cut", "Complete", "Complete+cut", "PAM", "HDBSCAN")
  FIGM <- c("R-UPGMA", "UPGMA", "UPGMA+cut", "Complete+cut", "PAM", "HDBSCAN")
  COLS <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300")
  SEP <- c("well", "moderate", "overlap"); CORR <- c("none", CORRUPTIONS[c(1, 2, 3, 4, 6, 5)])
  rd <- function(f) read.csv(gzfile(f), stringsAsFactors = FALSE)

  # Wilcoxon signed-rank test, zero differences split between ranks, normal
  # approximation with tie correction (scipy.stats.wilcoxon, zero_method="zsplit")
  wilcoxon_zsplit <- function(d) {
    n <- length(d); r <- rank(abs(d))
    rp <- sum(r[d > 0]) + sum(r[d == 0]) / 2; rm <- sum(r[d < 0]) + sum(r[d == 0]) / 2
    T <- min(rp, rm); mn <- n * (n + 1) / 4
    tt <- table(r); se <- sqrt(n * (n + 1) * (2 * n + 1) / 24 - sum(tt^3 - tt) / 48)
    if (se == 0) return(1)
    2 * pnorm(-abs((T - mn) / se))
  }

  violin <- function(data, cols, ylim = c(-0.1, 1.05), labels = NULL, main = "") {
    plot(NA, xlim = c(0.5, length(data) + 0.5), ylim = ylim, xaxt = "n", xlab = "", ylab = "", main = main,
         cex.main = 0.8, las = 1, cex.axis = 0.7)
    abline(h = seq(0, 1, 0.25), col = "#e6e6e6")
    for (i in seq_along(data)) {
      x <- data[[i]]
      if (length(unique(x)) > 1) {
        de <- density(x, bw = 0.15 * sd(x) + 1e-6, from = min(x), to = max(x))
        w <- 0.4 * de$y / max(de$y)
        polygon(c(i - w, rev(i + w)), c(de$x, rev(de$x)), col = cols[i], border = "white")
      } else segments(i - 0.35, x[1], i + 0.35, x[1], col = cols[i], lwd = 3)
      points(i, median(x), pch = 19, cex = 0.5)
    }
    if (!is.null(labels)) axis(1, at = seq_along(data), labels = labels, las = 2, cex.axis = 0.6)
  }

  out <- c()
  say <- function(...) { s <- paste0(...); cat(s, "\n"); out <<- c(out, s) }

  # ------------------------------------------------------------------ simulation
  if (file.exists("../results/R/sim_trials.csv.gz")) {
    sim <- rd("../results/R/sim_trials.csv.gz")
    main <- sim[sim$method %in% c(MAIN, "Single"), ]
    ov <- do.call(rbind, lapply(split(main$ari, main$method), function(a)
      c(mean = mean(a), median = median(a), fail = mean(a < 0.5), n = length(a))))
    write.csv(round(ov, 4), file.path(TAB, "sim_overall.csv"))
    sub <- main[main$level %in% c(0, 0.10), ]
    t4 <- aggregate(ari ~ separation + corruption + method, sub, median)
    t4w <- reshape(t4, idvar = c("separation", "corruption"), timevar = "method", direction = "wide")
    names(t4w) <- sub("^ari\\.", "", names(t4w))
    t4w <- t4w[order(match(t4w$separation, SEP), match(t4w$corruption, CORR)), c("separation", "corruption", MAIN)]
    write.csv(t4w, file.path(TAB, "T4_median_ari_eps10.csv"), row.names = FALSE)

    key <- c("structure", "separation", "corruption", "level", "rep")
    wide <- reshape(main[, c(key, "method", "ari")], idvar = key, timevar = "method", direction = "wide")
    names(wide) <- sub("^ari\\.", "", names(wide))
    t5 <- list()
    for (comp in c("PAM", "UPGMA+cut")) {
      for (sp in SEP) for (cr in CORR) {
        g <- wide[wide$separation == sp & wide$corruption == cr, ]
        d <- g[["R-UPGMA"]] - g[[comp]]
        t5[[length(t5) + 1]] <- data.frame(comparator = comp, separation = sp, corruption = cr, n = length(d),
          mean = mean(d), better = mean(d > 0.05), worse = mean(d < -0.05), p = wilcoxon_zsplit(d))
      }
    }
    t5 <- do.call(rbind, t5)
    t5$p_holm <- ave(t5$p, t5$comparator, FUN = function(p) p.adjust(p, "holm"))
    write.csv(t5, file.path(TAB, "T5_paired_vs_PAM_and_UPGMAcut.csv"), row.names = FALSE)

    ab <- sim[sim$rep < 30 & sim$method %in% c("R-UPGMA", ABLATIONS), ]
    t6 <- tapply(ab$ari, list(ab$corruption, ab$method), mean)[CORR, c("R-UPGMA", ABLATIONS)]
    t6 <- rbind(t6, "all (mean)" = tapply(ab$ari, ab$method, mean)[c("R-UPGMA", ABLATIONS)],
                "all (failure rate)" = tapply(ab$ari < 0.5, ab$method, mean)[c("R-UPGMA", ABLATIONS)])
    write.csv(round(t6, 4), file.path(TAB, "T6_ablation.csv"))

    png(file.path(FIG, "Fig1_sim_ari_distributions.png"), width = 7.2, height = 9.6, units = "in", res = 300)
    par(mfrow = c(7, 3), mar = c(0.6, 2.4, 1.2, 0.3), oma = c(4.5, 1.5, 0, 0))
    for (cr in CORR) for (sp in SEP) {
      g <- sub[sub$corruption == cr & sub$separation == sp, ]
      violin(lapply(FIGM, function(m) g$ari[g$method == m]), COLS,
             labels = if (cr == "row") FIGM else NULL,
             main = if (cr == "none") c(well = "Well separated", moderate = "Moderate", overlap = "Overlapping")[sp] else "")
      if (sp == "well") mtext(cr, side = 2, line = 1.6, cex = 0.5)
    }
    dev.off()
    png(file.path(FIG, "Fig2_sim_paired_vs_PAM_ecdf.png"), width = 7.2, height = 2.9, units = "in", res = 300)
    par(mfrow = c(1, 3), mar = c(4, 3.5, 1.5, 0.5))
    cc <- c(none = "#999999", salt_pepper = "#0d366b", amplification = "#256abf", mixed = "#5598e7",
            structured = "#eb6834", row = "#444444")
    for (sp in SEP) {
      plot(NA, xlim = c(-0.8, 1), ylim = c(0, 1), xlab = "ARI(R-UPGMA) - ARI(PAM)", ylab = "Cumulative proportion",
           main = sp, cex.main = 0.9)
      abline(v = 0, col = "grey")
      for (cr in names(cc)) { g <- wide[wide$separation == sp & wide$corruption == cr, ]
        lines(ecdf(g[["R-UPGMA"]] - g$PAM), col = cc[cr], do.points = FALSE, verticals = TRUE) }
    }
    legend("bottomright", legend = names(cc), col = cc, lwd = 1.5, cex = 0.6, bty = "n")
    dev.off()

    say("=== Simulation: overall ARI over all held-out simulated data sets ===")
    py <- read.csv("../results/python/tables/sim_overall.csv", row.names = 1)
    for (m in intersect(c(MAIN, "Single"), rownames(ov))) {
      se <- sd(main$ari[main$method == m]) / sqrt(ov[m, "n"])
      say(sprintf("%-13s R mean %.3f (MC s.e. %.3f)   manuscript %.3f   difference %+.3f   failure R %.1f%% / manuscript %.1f%%",
                  m, ov[m, "mean"], se, py[m, "mean"], ov[m, "mean"] - py[m, "mean"],
                  100 * ov[m, "fail"], 100 * py[m, "fail"]))
    }
  }

  # ------------------------------------------------------------------ benchmarks
  if (file.exists("../results/R/bench_trials.csv.gz")) {
    b <- rd("../results/R/bench_trials.csv.gz")
    med <- aggregate(ari ~ dataset + setting + method, b, median)
    t8 <- reshape(med, idvar = c("dataset", "setting"), timevar = "method", direction = "wide")
    names(t8) <- sub("^ari\\.", "", names(t8))
    write.csv(t8, file.path(TAB, "T8_benchmarks_median_ari.csv"), row.names = FALSE)
    m8 <- as.matrix(t8[, MAIN]); top <- apply(m8, 1, max)
    best <- colSums(m8 >= top - 0.01)
    say("\n=== Benchmarks: number of the 60 data set x setting combinations in which each method was best or tied ===")
    say(paste(sprintf("%s %d", MAIN, best), collapse = ", "))
    say("(Python: R-UPGMA 39, UPGMA+cut 19, HDBSCAN 15, PAM 13, UPGMA 8, W-UPGMA 8, Complete+cut 4, Complete 3)")
    pyb <- read.csv("../results/python/tables/T8_benchmarks_median_ari.csv", check.names = FALSE)
    cl <- merge(t8[t8$setting == "clean", c("dataset", MAIN)], pyb[pyb$setting == "clean", c("dataset", MAIN)],
                by = "dataset", suffixes = c(".R", ".py"))
    det <- !(cl$dataset %in% c("EngyTime", "MAGIC"))       # deterministic clean analyses
    for (m in MAIN) {
      dd <- abs(cl[det, paste0(m, ".R")] - cl[det, paste0(m, ".py")])
      say(sprintf("clean data, %-13s: %d of %d deterministic data sets identical to manuscript (max |diff| %.3f)",
                  m, sum(dd < 0.005), sum(det), max(dd)))
    }
    png(file.path(FIG, "Fig3_bench_saltpepper10.png"), width = 7.2, height = 6.2, units = "in", res = 300)
    par(mfrow = c(3, 5), mar = c(4.2, 2.2, 1.6, 0.3))
    for (ds in unique(b$dataset)) { g <- b[b$dataset == ds & b$setting == "salt_pepper10", ]
      violin(lapply(FIGM, function(m) g$ari[g$method == m]), COLS, labels = FIGM, main = ds) }
    dev.off()
    png(file.path(FIG, "Fig4_bench_mixed10.png"), width = 7.2, height = 6.2, units = "in", res = 300)
    par(mfrow = c(3, 5), mar = c(4.2, 2.2, 1.6, 0.3))
    for (ds in unique(b$dataset)) { g <- b[b$dataset == ds & b$setting == "mixed10", ]
      violin(lapply(FIGM, function(m) g$ari[g$method == m]), COLS, labels = FIGM, main = ds) }
    dev.off()
  }

  # ------------------------------------------------------------------ null data
  if (file.exists("../results/R/null_trials.csv.gz")) {
    nl <- rd("../results/R/null_trials.csv.gz")
    t7 <- cbind(aggregate(cbind(dip_p, largest, smallest, sil) ~ dataset + setting, nl, median),
                aggregate(cbind(reject_raw = dip_p < .05, reject_win = dip_p_winsorised < .05) ~ dataset + setting,
                          nl, mean)[, 3:4])
    write.csv(t7, file.path(TAB, "T7_null.csv"), row.names = FALSE)
    say("\n=== Null data (Table 7) ===")
    say(paste(capture.output(print(t7, digits = 2, row.names = FALSE)), collapse = "\n"))
  }

  # ------------------------------------------------------------------ timing
  tm <- do.call(rbind, lapply(c(150, 1000, 2000), function(n) {
    set.seed(1); D <- as.matrix(dist(matrix(rnorm(n * 10), n)))
    data.frame(n = n, upgma_s = system.time(cut_k(scipy_linkage(D, "average"), 3))[["elapsed"]],
               rupgma_s = system.time(robust_upgma(D, 3))[["elapsed"]])
  }))
  write.csv(tm, file.path(TAB, "T9_timing_R.csv"), row.names = FALSE)
  writeLines(out, "../results/R/comparison_with_manuscript.txt")
  cat("\nTables in", TAB, "; figures in", FIG, "; comparison in ../results/R/comparison_with_manuscript.txt\n")
  invisible(NULL)
}

# ============================ Part 4: RUN ===================================
# ---- find the package folder (the one that contains R/, data/ and results/) ----
.is_project <- function(p) length(p) == 1 && !is.na(p) && nzchar(p) &&
  dir.exists(file.path(p, "R")) && dir.exists(file.path(p, "data"))
.find_project <- function() {
  cand <- character(0)
  # 1. the folder of this script, when run with the "Source" button
  f <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
  if (!is.na(f)) cand <- c(cand, dirname(dirname(f)), dirname(f))
  # 2. the file open in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    f <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(f)) cand <- c(cand, dirname(dirname(normalizePath(f))), dirname(normalizePath(f)))
  }
  # 3. the working directory and the usual places
  home <- normalizePath("~", winslash = "/")
  docs <- c(home, file.path(home, "Documents"), file.path(dirname(home), "Documents"),
            file.path(Sys.getenv("USERPROFILE"), "Documents"), file.path(Sys.getenv("USERPROFILE"), "Downloads"),
            file.path(Sys.getenv("USERPROFILE"), "Desktop"))
  cand <- c(cand, getwd(), dirname(getwd()),
            file.path(docs, "robust-upgma"), file.path(docs, "robust-upgma-main"),
            file.path(docs, "robust-upgma-main", "robust-upgma-main"))
  for (p in unique(cand)) if (.is_project(p)) return(normalizePath(p, winslash = "/"))
  # 4. ask
  cat("\nPlease select the folder 'robust-upgma' (it contains the folders R, data and results).\n")
  p <- if (.Platform$OS.type == "windows") utils::choose.dir(caption = "Select the robust-upgma folder")
       else readline("Path to the robust-upgma folder: ")
  if (.is_project(p)) return(normalizePath(p, winslash = "/"))
  NA_character_
}
if (!.is_project(PROJECT_DIR)) PROJECT_DIR <- .find_project()
if (is.na(PROJECT_DIR)) stop("Could not find the robust-upgma folder. Unzip the package first ",
                             "(right-click the zip > Extract All) and run this script from its R folder.")
cat("Using the package folder:", PROJECT_DIR, "\n")

main <- function() {
  old_wd <- setwd(file.path(PROJECT_DIR, "R"))      # all paths below are relative to R/
  on.exit(setwd(old_wd), add = TRUE)
  if (RUN_VERIFY)     run_verification()
  if (RUN_EXAMPLE)    run_example()
  if (RUN_SIMULATION) run_simulation(SIM_REPS, CORES)
  if (RUN_BENCHMARKS) run_benchmarks(BENCH_REPS, CORES)
  if (RUN_NULL)       run_null(NULL_REPS)
  if (RUN_ANALYSIS)   run_analysis()
}
main()
cat("\nFinished. Results are in", file.path(PROJECT_DIR, "results", "R"), "\n")
