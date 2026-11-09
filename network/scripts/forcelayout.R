# Force-directed layouts for forcelayout.py: igraph's Fruchterman-Reingold (layout_with_fr) on each edge set, edge weights
# multiplying the attraction, started from the given coordinates (the t-SNE map), 2000 iterations, no grid, fixed seed.
# Usage: Rscript forcelayout.R <edges.tsv> <init.tsv> <out.tsv>
#   edges.tsv: layout, i, j (1-based node rows of init.tsv), w;  init.tsv: x, y per node;  out.tsv: layout, node, x, y
# Only nodes with at least one edge in a layout are placed here; forcelayout.py puts the others on a ring.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(igraph))
E <- read.delim(args[1]); P <- as.matrix(read.delim(args[2])); n <- nrow(P)
out <- list()
for (k in unique(E$layout)) {
  e <- E[E$layout == k, ]
  v <- sort(unique(c(e$i, e$j)))
  g <- graph_from_edgelist(cbind(match(e$i, v), match(e$j, v)), directed = FALSE)
  set.seed(20261108)
  init <- P[v, , drop = FALSE] * sqrt(length(v))
  xy <- layout_with_fr(g, coords = init, weights = e$w, niter = 2000, grid = "nogrid")
  out[[k]] <- data.frame(layout = k, node = v, x = round(xy[, 1], 4), y = round(xy[, 2], 4))
  cat(k, length(v), "nodes,", nrow(e), "edges,", components(g)$no, "components\n")
}
write.table(do.call(rbind, out), args[3], sep = "\t", quote = FALSE, row.names = FALSE)
