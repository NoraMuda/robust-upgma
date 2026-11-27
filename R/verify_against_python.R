# =============================================================================
# Checks that the R port gives the same partitions as the Python reference code
# on identical inputs.  The fixtures (R/fixtures/) hold, for 47 test cases
# (43 simulated conditions + 4 benchmark sets), the exact dissimilarities,
# the reference set R that Python used in Step 1, and Python's partitions.
#
# Run from the R/ folder:   Rscript verify_against_python.R
# Agreement is measured by the ARI between the R and Python partitions
# (1 = identical partition, up to relabelling).
# =============================================================================
source("robust_upgma.R")
source("experiment_utils.R")

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
