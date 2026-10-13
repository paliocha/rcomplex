# Is the 10-copy cap needed? Per-hog clique budget without it.
#
# For every hog in an edge table: genes with a called edge (q <
# alpha_graph) per species; twin classes = distinct partner sets; the
# product of class counts over species bounds the number of maximal
# cliques that span every species, and the sum over species subsets
# bounds all of them. Reports the hogs the current cap touches (> 10
# genes in a species), what their budget would be uncapped, after twin
# collapse, and under a product budget; and the exact maximal clique
# count (igraph) for hogs whose bound stays under a limit.
#
# Usage: Rscript copy-cap.R <edges.(tsv|rds)> [alpha_graph=0.1] [label]
args <- commandArgs(TRUE)
f <- args[1]
alpha <- if (length(args) > 1) as.numeric(args[2]) else 0.1
label <- if (length(args) > 2) args[3] else basename(f)
e <- if (grepl("\\.rds$", f)) readRDS(f) else data.table::fread(f, data.table = FALSE)
e <- e[
  !is.na(e$q_value) & e$q_value < alpha,
  c("hog", "gene1", "gene2", "species1", "species2")
]
cat(label, ": ", nrow(e), " called edges, ", length(unique(e$hog)), " hogs\n", sep = "")

per_hog <- function(d) {
  g <- c(d$gene1, d$gene2)
  s <- c(d$species1, d$species2)
  p <- c(d$gene2, d$gene1)
  nb <- split(p, g) # partner set per gene
  key <- vapply(nb, function(v) paste(sort(unique(v)), collapse = "|"), "")
  sp <- s[match(names(nb), g)]
  genes_sp <- table(sp)
  classes_sp <- tapply(key, sp, function(k) length(unique(k)))
  deg <- lengths(lapply(nb, unique))
  # the current cap: 10 most connected per species; ties by position
  tie <- vapply(names(genes_sp), function(x) {
    i <- which(sp == x)
    if (length(i) <= 10) {
      return(FALSE)
    }
    d10 <- sort(deg[i], decreasing = TRUE)[10]
    sum(deg[i] == d10) > 1 && sum(deg[i] >= d10) > 10
  }, NA)
  # budget: cliques spanning all species <= prod(classes); all cliques <=
  # prod(classes + 1) - 1 (every species in or out)
  data.frame(
    hog = d$hog[1], n_species = length(genes_sp), genes = sum(genes_sp),
    max_genes_sp = max(genes_sp), max_classes_sp = max(classes_sp),
    capped_now = any(genes_sp > 10), tie_at_cap = any(tie),
    capped_after_twins = any(classes_sp > 10),
    prod_genes = prod(as.numeric(genes_sp)),
    prod_classes = prod(as.numeric(classes_sp)),
    bound_all = prod(as.numeric(classes_sp) + 1) - 1
  )
}
t0 <- proc.time()[["elapsed"]]
H <- do.call(rbind, lapply(split(e, e$hog), per_hog))
cat("per-hog table in", round(proc.time()[["elapsed"]] - t0), "s\n")
H <- H[order(-H$bound_all), ]
cat(
  "\nhogs capped now (> 10 genes in a species):", sum(H$capped_now),
  " with a degree tie at the cut:", sum(H$tie_at_cap),
  "\nhogs that would still exceed 10 classes after twin collapse:",
  sum(H$capped_after_twins), "\n"
)
cat("\nclique bound (all species subsets) after twin collapse, no cap:\n")
print(quantile(H$bound_all, c(.5, .9, .99, .999, 1)))
cat(
  "hogs with bound > 1e4:", sum(H$bound_all > 1e4),
  " > 1e5:", sum(H$bound_all > 1e5), " > 1e6:", sum(H$bound_all > 1e6), "\n"
)
cat(
  "sum of bounds:", format(sum(H$bound_all), big.mark = ","),
  " sum excluding top 3:", format(sum(H$bound_all[-(1:3)]), big.mark = ","), "\n"
)
cat("\ntop 12 hogs by bound:\n")
print(H[1:12, c(
  "hog", "n_species", "genes", "max_genes_sp", "max_classes_sp",
  "capped_now", "prod_classes", "bound_all"
)], row.names = FALSE, digits = 3)

# exact maximal cliques for the capped hogs whose bound is manageable
cap_h <- H$hog[H$capped_now]
ex <- do.call(rbind, lapply(cap_h, function(h) {
  d <- e[e$hog == h, ]
  b <- H$bound_all[H$hog == h]
  if (b > 2e6) {
    return(data.frame(hog = h, bound = b, cliques_uncapped = NA, sec = NA))
  }
  g <- c(d$gene1, d$gene2)
  p <- c(d$gene2, d$gene1)
  nb <- split(p, g)
  key <- vapply(nb, function(v) paste(sort(unique(v)), collapse = "|"), "")
  s <- c(d$species1, d$species2)[match(names(nb), g)]
  rep_of <- names(nb)[match(paste(s, key), paste(s, key))]
  m <- unique(data.frame(
    a = rep_of[match(d$gene1, names(nb))],
    b = rep_of[match(d$gene2, names(nb))]
  ))
  gr <- igraph::graph_from_data_frame(m, directed = FALSE)
  t1 <- proc.time()[["elapsed"]]
  n <- igraph::count_max_cliques(gr, min = 3)
  data.frame(
    hog = h, bound = b, cliques_uncapped = n,
    sec = round(proc.time()[["elapsed"]] - t1, 2)
  )
}))
cat("\nexact twin-collapsed maximal cliques (>= 3) for hogs the cap touches now:\n")
ex <- ex[order(-ex$bound), ]
print(head(ex, 15), row.names = FALSE)
cat(
  "total uncapped cliques over capped hogs:", format(sum(ex$cliques_uncapped, na.rm = TRUE), big.mark = ","),
  " unresolved (bound > 2e6):", sum(is.na(ex$cliques_uncapped)), "\n"
)
out <- sub("\\.(tsv|rds)$", "", f)
saveRDS(list(H = H, exact = ex), file.path(dirname(f), paste0("copycap_", label, ".rds")))
