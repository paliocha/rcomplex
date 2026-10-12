# Forman-Ricci curvature of the clique-lifted hypergraph per clique and hog.
# Design note: hypergraph-community-literature.md, section 9.2.
# Inputs, read from the directory in env var ENC_DIR:
#   E_<n>.rds   find_coexpressologs() edge table
#   g_<n>.rds   gene_cliques object;  cls_<n>.rds   its classification
# Not package code; no lint, no tests.
suppressPackageStartupMessages(library(collapse))
g <- readRDS(file.path(Sys.getenv("ENC_DIR"), "g_12000.rds")); cls <- readRDS(file.path(Sys.getenv("ENC_DIR"), "cls_12000.rds"))
m <- g$members; cl <- g$cliques
cat("cliques", nrow(cl), "member rows", nrow(m), "\n")
# node = (species, gene); d_k = number of cliques holding the node
node <- paste(m$species, m$gene, sep = "|")
nid <- match(node, unique(node)); d <- tabulate(nid)
# Forman-Ricci of a clique-as-hyperedge: F(e) = sum_{k in e} (2 - d_k)
frc <- fsum(2 - d[nid], m$clique_id)
cl$frc <- frc[as.character(cl$clique_id)]
cl$n <- tabulate(m$clique_id, nbins = nrow(cl))
cl$tier <- cls$classification[match(cl$clique_id, cls$clique_id)]
cat("\nF(e) > 0 means every member sits in few cliques (resolved); << 0 = combinatorial variant\n")
print(quantile(cl$frc, c(0, .01, .1, .25, .5, .75, .9, .99, 1)))
cat("\nshare of cliques with F > 0:", mean(cl$frc > 0), "\n")
cat("cliques with F > 0 per HOG: hogs", length(unique(cl$hog[cl$frc > 0])), "of", length(unique(cl$hog)), "\n")
cat("\nby tier: n, median F, share F > 0\n")
print(cl |> fgroup_by(tier) |> fsummarise(n = fnobs(frc), med_F = fmedian(frc), pos = fmean(frc > 0)))
# per HOG: how many cliques, and the max F — one resolved combination?
h <- cl |> fgroup_by(hog) |> fsummarise(n_cliques = fnobs(frc), max_F = fmax(frc), n_pos = fsum(frc > 0))
cat("\nHOGs by clique count (quantiles):\n"); print(quantile(h$n_cliques, c(.5, .9, .99, 1)))
cat("HOGs with >= 1 resolved clique (F > 0):", sum(h$n_pos > 0), "of", nrow(h), "\n")
cat("HOGs with > 100 cliques and 0 resolved:", sum(h$n_cliques > 100 & h$n_pos == 0), "\n")
# per copy: HCP-FRC mean over its cliques; spread within (hog, species) tells copies apart
m$frc <- cl$frc[match(m$clique_id, cl$clique_id)]
m$hog <- cl$hog[match(m$clique_id, cl$clique_id)]
cp <- m |> fgroup_by(hog, species, gene) |> fsummarise(n_cl = fnobs(frc), mean_F = fmean(frc))
ms <- cp |> fgroup_by(hog, species) |> fsummarise(copies = fnobs(gene), range_F = fmax(mean_F) - fmin(mean_F), max_F = fmax(mean_F))
cat("\n(hog, species) cells:", nrow(ms), " multi-copy:", sum(ms$copies > 1), "\n")
cat("multi-copy cells where copies differ in mean F (range > 0.5):", sum(ms$copies > 1 & ms$range_F > 0.5), "\n")
print(quantile(ms$range_F[ms$copies > 1], c(.5, .75, .9, .99)))
cat("\nexample multi-copy cells, largest range:\n")
ex <- ms[ms$copies > 1, ]; ex <- ex[order(-ex$range_F), ][1:3, ]
for (i in seq_len(nrow(ex))) print(cp[cp$hog == ex$hog[i] & cp$species == ex$species[i], ])
