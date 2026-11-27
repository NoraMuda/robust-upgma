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
