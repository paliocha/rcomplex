# Probe 2: Ollivier-Ricci curvature of hog hyperedges and OT copy map.
# Design note section 9.6. Run from the repository root with ENC_DIR set.
# Needs leaf_se.rds, E_12000.rds, g_12000.rds, cls_B_12000.rds in ENC_DIR;
# nets12k.rds is written there on first run; orc_out.rds is the output.
# Not package code.
t0 <- Sys.time()
source("dev/probes/encodings-common.R")
n <- nrow(genes)
hogs_all <- sort(unique(genes$hog)) # no NA hogs in the 12k set
stopifnot(!anyNA(genes$hog))

# nu_g(h): mu_g uniform over within-species neighbours, projected on hogs.
# One genes x hogs matrix per species, stacked in genes order.
Nu <- do.call(rbind, lapply(sp_names, function(s) {
  a <- adj[[s]]
  dg <- rowSums(a)
  hix <- match(rowData(leaf[[s]])$hog, hogs_all)
  M <- sparseMatrix(
    i = seq_along(hix), j = hix, x = 1,
    dims = c(length(hix), length(hogs_all))
  )
  Diagonal(x = ifelse(dg > 0, 1 / dg, 0)) %*% a %*% M
}))
mass <- rowSums(Nu)
cat("genes", n, " zero-degree genes", sum(mass == 0), "\n")

# W1 with ground metric d(x, y) = 0 if hog(x) == hog(y), else 1.
# Any coupling pays 1 per unit of mass moved across hog buckets and 0 inside
# a bucket. Inside bucket h at most min(nu_g(h), nu_g'(h)) can stay at cost
# 0, and that bound is reached by matching inside every bucket at once, so
#   W1 = 1 - sum_h min(nu_g(h), nu_g'(h)).
# sum_h min(a, b) = (|a|_1 + |b|_1 - |a - b|_1) / 2, so W1 comes from one
# manhattan distance matrix per hog. ("none" bucket unused: no NA hogs.)
w1_mat <- function(idx) {
  X <- Nu[idx, , drop = FALSE]
  X <- X[, colSums(X) > 0, drop = FALSE]
  l1 <- as.matrix(dist(as.matrix(X), method = "manhattan"))
  m <- mass[idx]
  W <- 1 - (outer(m, m, "+") - l1) / 2
  dimnames(W) <- list(idx, idx)
  W
}

# hog hyperedges with members in >= 2 species
by_hog <- split(seq_len(n), genes$hog)
by_hog <- by_hog[vapply(by_hog, function(i) length(unique(genes$species[i])) >= 2, NA)]
cat(
  "hog hyperedges (>= 2 species)", length(by_hog), " max size",
  max(lengths(by_hog)), "\n"
)
t1 <- Sys.time()
Ws <- parallel::mclapply(by_hog, w1_mat, mc.cores = 6)
stopifnot(!any(vapply(Ws, inherits, NA, "try-error")))
cat("W1 matrices:", format(Sys.time() - t1), "\n")

# curvature per hog (cross-species pairs) and per copy
kap <- do.call(rbind, lapply(names(Ws), function(h) {
  W <- Ws[[h]]
  s <- genes$species[as.integer(rownames(W))]
  x <- W[outer(s, s, "!=") & upper.tri(W)]
  data.frame(
    hog = h, size = nrow(W), kappa_mean = 1 - mean(x),
    kappa_min = 1 - max(x)
  )
}))
copy_w1 <- unlist(lapply(Ws, function(W) {
  s <- genes$species[as.integer(rownames(W))]
  W[outer(s, s, "==")] <- NA
  setNames(rowMeans(W, na.rm = TRUE), rownames(W))
}), use.names = FALSE)
copy_ix <- as.integer(unlist(lapply(Ws, rownames), use.names = FALSE))
genes$copy_w1 <- NA_real_
genes$copy_w1[copy_ix] <- copy_w1
print(summary(kap[, 3:4]))

# Q4: kappa vs hog_class
kap$tier <- hog_tier[kap$hog]
kt <- kap[!is.na(kap$tier), ]
cat("\nQ4 kappa by hog_class (median)\n")
print(kt |> fgroup_by(tier) |> fsummarise(
  n = fnobs(kappa_mean),
  kappa_mean = fmedian(kappa_mean), kappa_min = fmedian(kappa_min),
  size = fmedian(size)
))
for (v in c("kappa_mean", "kappa_min")) {
  cc <- kt[kt$tier %in% c("complete_conserved", "unclassified"), ]
  cat(
    "AUROC", v, "complete_conserved vs unclassified:",
    round(auroc(cc[[v]], cc$tier == "complete_conserved"), 3),
    " any tier vs unclassified:",
    round(auroc(kt[[v]], kt$tier != "unclassified"), 3), "\n"
  )
}
cat("hogs with a tier but no kappa:", sum(!names(hog_tier) %in% kap$hog), "\n")

# Q2-analogue: copy with lowest mean W1 vs copy with most called edges
mc <- genes[genes$copies > 1 & !is.na(genes$copy_w1), ]
mc <- mc[ave(mc$n_called, mc$cell, FUN = max) > 0, ]
o <- order(mc$cell, -mc$n_called, mc$best_q)
top_call <- mc$node[o][!duplicated(mc$cell[o])]
o <- order(mc$cell, mc$copy_w1)
top_w <- mc$node[o][!duplicated(mc$cell[o])]
cat(
  "\ncopy with min mean W1 = most-called copy:", round(mean(top_w == top_call), 3),
  " baseline", round(mean(1 / table(mc$cell)), 3), " cells", length(top_call), "\n"
)

# Q5: twins, W1 of the twins' own projected measures (same species, hog)
w1_pair <- function(i, j) 1 - rowSums(pmin(Nu[i, , drop = FALSE], Nu[j, , drop = FALSE]))
mu_pair <- function(i, j) { # unprojected: exact W1 on genes, 0/1 metric
  s <- genes$species[i]
  out <- numeric(length(i))
  for (sp in unique(s)) {
    k <- which(s == sp)
    a <- adj[[sp]]
    dg <- rowSums(a)
    o <- match(sp, sp_names)
    base <- sum(sapply(adj, nrow)[seq_len(o - 1)])
    ii <- i[k] - base
    jj <- j[k] - base
    D <- Diagonal(x = ifelse(dg > 0, 1 / dg, 0)) %*% a
    out[k] <- 1 - rowSums(pmin(D[ii, , drop = FALSE], D[jj, , drop = FALSE]))
  }
  out
}
ps <- pair_sets()
q5 <- t(sapply(ps[c("twin", "nontwin", "random")], function(p) {
  set.seed(2)
  p <- p[sample(nrow(p), min(nrow(p), 20000)), , drop = FALSE]
  c(
    n = nrow(p),
    w1_hog = qtab(w1_pair(p[, 1], p[, 2]))[-1],
    w1_gene = qtab(mu_pair(p[, 1], p[, 2]))[-1],
    zero = mean(w1_pair(p[, 1], p[, 2]) < 0.01)
  )
}))
cat("\nQ5 within-species W1 (hog-projected and gene-level)\n")
print(q5)

# Copy map per (hog, A, B) with a multi-copy side
orth <- prepare_orthologs(genes[, c("species", "gene", "hog")])
stopifnot(all(c("gene1", "gene2", "hog") %in% names(orth)))
res_map <- do.call(rbind, lapply(combn(sp_names, 2, simplify = FALSE), function(pr) {
  m <- resolve_ortholog_map(orth, rownames(adj[[pr[1]]]), rownames(adj[[pr[2]]]),
    pr[1], pr[2],
    edges = E
  )
  m <- m[m$source == "coexpressolog", ]
  if (nrow(m)) data.frame(A = pr[1], B = pr[2], m) else NULL
}))
nid <- function(s, g) match(paste(s, g, sep = "|"), genes$node)
res_key <- paste(nid(res_map$A, res_map$gene1), nid(res_map$B, res_map$gene2))
Ec$i <- nid(Ec$species1, Ec$gene1)
Ec$j <- nid(Ec$species2, Ec$gene2)
Ec <- Ec[order(Ec$q_value, -Ec$effect_size), ]
best_q <- Ec[!duplicated(paste(
  Ec$hog, pmin(Ec$species1, Ec$species2),
  pmax(Ec$species1, Ec$species2)
)), ]
bq_key <- c(paste(best_q$i, best_q$j), paste(best_q$j, best_q$i))

t2 <- Sys.time()
cm <- do.call(rbind, parallel::mclapply(names(Ws), function(h) {
  W <- Ws[[h]]
  ix <- as.integer(rownames(W))
  s <- genes$species[ix]
  sp <- intersect(sp_names, s)
  out <- list()
  for (pr in if (length(sp) > 1) combn(sp, 2, simplify = FALSE) else list()) {
    a <- ix[s == pr[1]]
    b <- ix[s == pr[2]]
    kA <- length(a)
    kB <- length(b)
    if (kA == 1 && kB == 1) next
    C <- W[as.character(a), as.character(b), drop = FALSE]
    # OT with uniform marginals, scaled to integers: rows kB, cols kA
    plan <- lpSolve::lp.transport(
      C, "min", rep("=", kA), rep(kB, kA),
      rep("=", kB), rep(kA, kB)
    )$solution / (kA * kB)
    # heaviest cell, ties to the cheapest
    hc <- order(-round(plan, 12), C)[1]
    hk <- paste(a[row(C)[hc]], b[col(C)[hc]])
    mn <- order(C)[1]
    mk <- paste(a[row(C)[mn]], b[col(C)[mn]])
    rb <- which(C == apply(C, 1, min)[row(C)] & C == apply(C, 2, min)[col(C)])
    keys <- outer(a, b, paste)
    out[[length(out) + 1]] <- data.frame(
      hog = h, kA = kA, kB = kB,
      near_perm = mean(apply(plan, 1, max) >= 0.8 / kA),
      support = sum(plan > 1e-12), n_rbp = length(rb),
      c_range = diff(range(C)),
      ot_gain = mean(C) - sum(plan * C), # uniform coupling minus OT cost
      has_q = any(keys %in% bq_key),
      heavy_q = hk %in% bq_key, argmin_q = mk %in% bq_key,
      q_in_support = any(keys[plan > 1e-12] %in% bq_key),
      q_rbp = any(keys[rb] %in% bq_key),
      n_res = sum(keys %in% res_key),
      res_in_support = sum(keys[plan > 1e-12] %in% res_key),
      res_rbp = sum(keys[rb] %in% res_key)
    )
  }
  do.call(rbind, out)
}, mc.cores = 6))
cat("\ncopy map:", format(Sys.time() - t2), " triples", nrow(cm), "\n")
saveRDS(list(kap = kap, cm = cm, genes = genes), file.path(Sys.getenv("ENC_DIR"), "orc_out.rds"))

cm$shape <- ifelse(cm$kA > 1 & cm$kB > 1, "both multi", "one side 1")
cat("\nCopy map summary by shape\n")
print(cm |> fgroup_by(shape) |> fsummarise(
  n = fnobs(kA), near_perm_all = fmean(near_perm == 1),
  med_c_range = fmedian(c_range), med_ot_gain = fmedian(ot_gain),
  flat_C = fmean(c_range < 0.01)
))
q <- cm[cm$has_q, ]
cat("\nbest-q pair vs plan, triples with a called edge\n")
print(q |> fgroup_by(shape) |> fsummarise(
  n = fnobs(kA), heavy = fmean(heavy_q), argminC = fmean(argmin_q),
  base_cell = fmean(1 / (kA * kB)),
  in_support = fmean(q_in_support), base_support = fmean(support / (kA * kB)),
  rbp = fmean(q_rbp), base_rbp = fmean(n_rbp / (kA * kB))
))
r <- cm[cm$n_res > 0, ]
cat("\nresolve_ortholog_map coexpressolog pairs vs plan\n")
print(r |> fgroup_by(shape) |> fsummarise(
  triples = fnobs(kA), pairs = fsum(n_res),
  in_support = fsum(res_in_support) / fsum(n_res),
  base_support = fmean(support / (kA * kB)),
  rbp = fsum(res_rbp) / fsum(n_res), base_rbp = fmean(n_rbp / (kA * kB))
))
cat("resolved coexpressolog pairs total", nrow(res_map), "\n")
cat("total", format(Sys.time() - t0), "\n")
