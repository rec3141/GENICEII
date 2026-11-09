# Reference run of the R SPRING package for validate_spring.py (validation only; the network uses the Python port).
# SPRING itself is not on CRAN: its two source files are sourced from a checkout of https://github.com/GraceYoon/SPRING
# (R/SPRING.R, R/helpers.R). latentcor is installed in a scratch library (R_LIBS_USER); huge and pulsar are system packages.
# Usage: Rscript spring_ref.R <SPRING checkout> <data.tsv> <out prefix> [rep.num]
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; inp <- args[2]; out <- args[3]; reps <- if (length(args) > 3) as.integer(args[4]) else 20L
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({ library(latentcor); library(huge); library(pulsar) })
source(file.path(src, "R", "helpers.R")); source(file.path(src, "R", "SPRING.R"))
X <- as.matrix(read.delim(inp, row.names = 1, check.names = FALSE))
q <- mclr(X)
write.table(q, paste0(out, "_mclr.tsv"), sep = "\t", quote = FALSE, col.names = NA)
for (m in c("approx", "original")) {
  lc <- latentcor(q, types = "tru", method = m, use.nearPD = TRUE, nu = 0.001, ratio = 0.9, tol = 1e-6)
  write.table(lc$R, paste0(out, "_R_", m, ".tsv"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(lc$Rpointwise, paste0(out, "_Rpw_", m, ".tsv"), sep = "\t", quote = FALSE, col.names = NA)
  write.table(lc$K, paste0(out, "_K.tsv"), sep = "\t", quote = FALSE, col.names = NA)
}
lam <- exp(seq(log(0.6), log(0.6 * 1e-2), length.out = 20))
S <- read.delim(paste0(out, "_R_original.tsv"), row.names = 1, check.names = FALSE)
hm <- huge.mb(as.matrix(S), lambda = lam, sym = "or", verbose = FALSE)
path <- sapply(hm$path, function(a) sum(a[upper.tri(a)]))
saveRDS(hm$path, paste0(out, "_path.rds"))
for (i in seq_along(hm$path)) {
  a <- as.matrix(hm$path[[i]]); w <- which(upper.tri(a) & a != 0, arr.ind = TRUE)
  write.table(w, paste0(out, "_path_", i, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
fit <- SPRING(X, quantitative = FALSE, lambdaseq = lam, nlambda = 20, rep.num = reps, seed = 10010, ncores = 1, verbose = FALSE, Rmethod = "approx")
opt <- fit$output$stars$opt.index
a <- as.matrix(fit$fit$refit$stars); w <- which(upper.tri(a) & a != 0, arr.ind = TRUE)
write.table(w, paste0(out, "_refit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(lambda = lam, summary = fit$output$stars$summary, path_edges = path), paste0(out, "_stars.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("R SPRING: opt.index", opt, "edges", nrow(w), "\n")
