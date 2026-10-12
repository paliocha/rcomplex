# False-twin classes among same-species paralogs; clique reduction.
# Design note: hypergraph-community-literature.md, section 9.3.
# Inputs, read from the directory in env var ENC_DIR:
#   E_<n>.rds   find_coexpressologs() edge table
#   g_<n>.rds   gene_cliques object;  cls_<n>.rds   its classification
# Not package code; no lint, no tests.
suppressPackageStartupMessages(library(collapse))
e <- readRDS(file.path(Sys.getenv("ENC_DIR"), "E_12000.rds")); g <- readRDS(file.path(Sys.getenv("ENC_DIR"), "g_12000.rds"))
e <- e[!is.na(e$q_value) & e$q_value < 0.1, ]       # gene-graph edges at alpha_graph
a <- paste(e$hog, e$species1, e$gene1, sep = "|"); b <- paste(e$hog, e$species2, e$gene2, sep = "|")
nodes <- unique(c(a, b)); ai <- match(a, nodes); bi <- match(b, nodes)
nb <- split(c(bi, ai), c(ai, bi))                   # neighbour sets
sig <- vapply(nb, function(x) rlang::hash(sort(unique(x))), "")
twin <- match(sig, unique(sig))[match(seq_along(nodes), as.integer(names(nb)))]
hog <- sub("\\|.*", "", nodes); sp <- sub("^[^|]*\\|([^|]*)\\|.*", "\\1", nodes)
cells <- data.frame(hog, sp, twin) |> fgroup_by(hog, sp) |>
  fsummarise(copies = fnobs(twin), classes = fndistinct(twin))
mc <- cells[cells$copies > 1, ]
cat("gene-graph nodes", length(nodes), " false-twin classes", length(unique(twin)), "\n")
cat("multi-copy (hog, species) cells", nrow(mc), " copies", sum(mc$copies), " -> classes", sum(mc$classes),
    "\n  cells where all copies are twins:", sum(mc$classes == 1),
    "  cells with some twin pair:", sum(mc$classes < mc$copies), "\n")
m <- g$members; key <- paste(g$cliques$hog[match(m$clique_id, g$cliques$clique_id)], m$species, m$gene, sep = "|")
tw_m <- twin[match(key, nodes)]; stopifnot(!anyNA(tw_m))
o <- order(m$clique_id, tw_m)
csig <- vapply(split(tw_m[o], m$clique_id[o]), function(x) rlang::hash(x), "")
cat("cliques", nrow(g$cliques), " -> distinct up to twin substitution:", length(unique(csig)), "\n")
ch <- g$cliques$hog[match(as.integer(names(csig)), g$cliques$clique_id)]
top <- names(sort(table(ch), decreasing = TRUE))[1:3]
for (h in top) cat(" ", h, ": cliques", sum(ch == h), "-> quotient", length(unique(csig[ch == h])),
  " copies/species:", paste(cells$copies[cells$hog == h], collapse = ","),
  " classes:", paste(cells$classes[cells$hog == h], collapse = ","), "\n")
