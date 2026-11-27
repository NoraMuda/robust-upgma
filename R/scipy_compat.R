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
