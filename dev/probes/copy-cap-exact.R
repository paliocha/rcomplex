# Exact twin-collapsed maximal cliques for the families copy-cap.R could not
# bound: uncapped against the 10 most connected classes per species.
# Usage: Rscript copy-cap-exact.R <edges.(tsv|rds)> <label>  (after copy-cap.R)
args <- commandArgs(TRUE)
f <- args[1]
label <- args[2]
e <- if (grepl("\\.rds$", f)) readRDS(f) else data.table::fread(f, data.table = FALSE)
e <- e[!is.na(e$q_value) & e$q_value < 0.1, c("hog", "gene1", "gene2", "species1", "species2")]
cc <- readRDS(file.path(dirname(f), paste0("copycap_", label, ".rds")))
big <- cc$exact$hog[is.na(cc$exact$cliques_uncapped)]
collapse_graph <- function(d, cap = NULL) {
  g <- c(d$gene1, d$gene2); p <- c(d$gene2, d$gene1); s <- c(d$species1, d$species2)
  nb <- split(p, g)
  key <- vapply(nb, function(v) paste(sort(unique(v)), collapse = "|"), "")
  sp <- s[match(names(nb), g)]
  cls <- paste(sp, key)
  rep_of <- names(nb)[match(cls, cls)]
  m <- unique(data.frame(a = rep_of[match(d$gene1, names(nb))], b = rep_of[match(d$gene2, names(nb))]))
  if (!is.null(cap)) {
    reps <- unique(rep_of); rsp <- sp[match(reps, names(nb))]
    deg <- table(c(m$a, m$b))[reps]; deg[is.na(deg)] <- 0
    keep <- unlist(lapply(split(reps, rsp), function(r) r[order(-deg[r], r)][seq_len(min(cap, length(r)))]))
    m <- m[m$a %in% keep & m$b %in% keep, ]
  }
  igraph::graph_from_data_frame(m, directed = FALSE)
}
count <- function(gr, limit = 180) {
  t1 <- proc.time()[["elapsed"]]
  n <- tryCatch({ setTimeLimit(elapsed = limit, transient = TRUE); igraph::count_max_cliques(gr, min = 3) },
    error = function(err) NA)
  setTimeLimit()
  c(n = n, sec = round(proc.time()[["elapsed"]] - t1, 1))
}
out <- do.call(rbind, lapply(big, function(h) {
  d <- e[e$hog == h, ]
  g0 <- collapse_graph(d); g1 <- collapse_graph(d, cap = 10)
  u <- count(g0); c10 <- count(g1)
  r <- data.frame(hog = h, classes = igraph::vcount(g0), edges = igraph::ecount(g0),
    uncapped = u[["n"]], sec = u[["sec"]], capped10 = c10[["n"]], sec10 = c10[["sec"]])
  print(r, row.names = FALSE); r
}))
cat("\n", label, "total uncapped:", format(sum(out$uncapped, na.rm = TRUE), big.mark = ","),
    " capped10:", format(sum(out$capped10, na.rm = TRUE), big.mark = ","),
    " timed out:", sum(is.na(out$uncapped)), "\n")
