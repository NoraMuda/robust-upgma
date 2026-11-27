# =============================================================================
# Monte Carlo simulation (Sections 3.2 and 4 of the manuscript), R version.
#   171 conditions x `reps` replicates; every method on the same corrupted
#   dissimilarities; ablation variants on the first 30 replicates.
# Usage (from the R/ folder):
#   Rscript run_simulation.R            # 100 replicates, all cores but one
#   Rscript run_simulation.R 10 4       # 10 replicates, 4 cores (quick check)
# Output: ../results/R/sim_trials.csv.gz
# =============================================================================
source("robust_upgma.R"); source("experiment_utils.R")
args <- commandArgs(trailingOnly = TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 100L
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)
dir.create("../results/R", recursive = TRUE, showWarnings = FALSE)

one_trial <- function(task) {
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
invisible(parallel::clusterEvalQ(cl, { source("robust_upgma.R"); source("experiment_utils.R") }))
res <- parallel::parLapplyLB(cl, tasks, one_trial)
parallel::stopCluster(cl)
res <- do.call(rbind, res)
write.csv(res, gzfile("../results/R/sim_trials.csv.gz"), row.names = FALSE)
cat("done in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
