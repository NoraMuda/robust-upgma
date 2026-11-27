# =============================================================================
# Tables and figures of the manuscript from the R results, plus a side-by-side
# comparison with the published (Python) numbers.
# Usage:  Rscript analysis.R          (after run_simulation.R, run_benchmarks.R,
#                                      run_null.R)
# Output: ../results/R/tables/*.csv, ../results/R/figures/*.png,
#         ../results/R/comparison_with_manuscript.txt
# =============================================================================
source("robust_upgma.R"); source("experiment_utils.R")
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
