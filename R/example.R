# Minimal example of Robust UPGMA (run from the R/ folder)
source("robust_upgma.R")
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
