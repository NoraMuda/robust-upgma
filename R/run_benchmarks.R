# =============================================================================
# Benchmarks on FCPS and labelled real data (Section 5), R version.
#   Settings: clean, 10% salt-and-pepper, 10% mixed, 5% observation outliers;
#   30 realisations of each corrupted setting.  Data sets with n > 2000
#   (EngyTime, MAGIC) are class-stratified subsampled to 2000 rows per run.
# Usage:  Rscript run_benchmarks.R [reps=30] [cores]
# Output: ../results/R/bench_trials.csv.gz
# =============================================================================
source("robust_upgma.R"); source("experiment_utils.R")
args <- commandArgs(trailingOnly = TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 30L
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)
dir.create("../results/R", recursive = TRUE, showWarnings = FALSE)

source("data_loader.R")

one_run <- function(task) {
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
invisible(parallel::clusterEvalQ(cl, { source("robust_upgma.R"); source("experiment_utils.R"); source("data_loader.R") }))
res <- do.call(rbind, parallel::parLapplyLB(cl, tasks, one_run))
parallel::stopCluster(cl)
write.csv(res, gzfile("../results/R/bench_trials.csv.gz"), row.names = FALSE)
cat("done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
