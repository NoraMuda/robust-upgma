# Data sets and settings of the benchmark study (shared by run_benchmarks.R and run_null.R)
# Paths are relative to the R/ folder: data files live in ../data/.
if (!exists("sd_pop")) source("experiment_utils.R")

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

