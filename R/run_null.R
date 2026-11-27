# =============================================================================
# Data without cluster structure (Sections 3.5 and 4.4), R version.
#   NULL-Gauss, NULL-Unif, GolfBall (subsample of 150) and two positive controls
#   (B3 well separated, Hepta subsample of 150); clean and 10% salt-and-pepper.
#   Dip test on raw and on winsorised distances; Robust UPGMA asked for K = 3.
# Usage:  Rscript run_null.R [reps=100]
# Output: ../results/R/null_trials.csv.gz
# =============================================================================
source("robust_upgma.R"); source("experiment_utils.R")
source("data_loader.R")
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[1]) else 100L
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
