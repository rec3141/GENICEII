# SPIEC-EASI (Kurtz et al. 2015) on each layer, as a comparison method: spiec.easi(method = "mb") on the node counts
# (SpiecEasi applies clr with a pseudocount itself), lambda path of 20 values down to 0.01 x lambda.max, StARS with 20
# subsamples (SpiecEasi default threshold 0.05). Edge weight and sign = symmetrised MB coefficient (symBeta, mode "ave").
# Usage: Rscript run_spieceasi.R <network dir>   (reads cache/{dna,rna}_X.tsv, writes out/edges_spieceasi_{layer}.tsv)
args <- commandArgs(trailingOnly = TRUE); net <- if (length(args)) args[1] else "."
suppressPackageStartupMessages({ library(SpiecEasi); library(Matrix) })
for (layer in c("dna", "rna")) {
  X <- as.matrix(read.delim(file.path(net, "cache", paste0(layer, "_X.tsv")), row.names = 1, check.names = FALSE))
  set.seed(10010)
  se <- spiec.easi(X, method = "mb", lambda.min.ratio = 1e-2, nlambda = 20,
                   pulsar.params = list(rep.num = 20, seed = 10010, ncores = 8))
  B <- as.matrix(symBeta(getOptBeta(se), mode = "ave"))
  A <- as.matrix(getRefit(se))
  w <- which(upper.tri(A) & A != 0, arr.ind = TRUE)
  ids <- colnames(X)
  E <- data.frame(a = ids[w[, 1]], b = ids[w[, 2]], weight = round(B[w], 4))
  E$sign <- sign(E$weight)
  mrg <- se$select$stars$merge[[getOptInd(se)]]
  E$stability <- round(mrg[w], 3)
  sw <- E$a > E$b; tmp <- E$a[sw]; E$a[sw] <- E$b[sw]; E$b[sw] <- tmp
  write.table(E, file.path(net, "out", paste0("edges_spieceasi_", layer, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(layer, ncol(X), "nodes: StARS index", getOptInd(se), "lambda", signif(getOptLambda(se), 3), ",", nrow(E), "edges,", sum(E$sign > 0), "positive\n")
}
