# Probe 1: equal-edge random walk on the mixed hypergraph (12k leaf set).
# Design note section 9.6. Run from the repository root with ENC_DIR set.
# Needs leaf_se.rds, E_12000.rds, g_12000.rds, cls_B_12000.rds in ENC_DIR;
# nets12k.rds is written there on first run; rw_genes.rds is the output.
# Not package code.
t0 <- Sys.time()
source("dev/probes/encodings-common.R")
n <- nrow(genes)
off <- c(0, cumsum(sapply(adj, nrow)))[seq_along(sp_names)]
names(off) <- sp_names

# co-expression hyperedges (size 2): block-diagonal adjacency
A_c <- as(bdiag(adj[sp_names]), "generalMatrix")
# hog hyperedges: hogs with >= 2 member genes, members = all copies
hs <- table(genes$hog)
keep_h <- names(hs)[hs >= 2]
in_h <- genes$hog %in% keep_h
H <- sparseMatrix(
  i = which(in_h), j = match(genes$hog[in_h], keep_h),
  x = 1, dims = c(n, length(keep_h))
)
w <- 1 / (as.vector(hs[keep_h]) - 1) # W = diag(1 / (|e| - 1))
A_h <- H %*% Diagonal(x = w) %*% t(H)
diag(A_h) <- 0
A_h <- drop0(A_h)
A <- as(A_c + A_h, "generalMatrix") # H W H^T minus diagonal
d <- rowSums(A_c) + in_h # deg_H: incident hyperedges
dinv <- ifelse(d > 0, 1 / d, 0)
# P = Dv^-1 A, rows sum to 1 (sum_e (|e|-1)/(|e|-1) / deg); P_c = Dv^-1 A_c
stopifnot(abs(range((dinv * rowSums(A))[d > 0]) - 1) < 1e-9)
cat(
  "nodes", n, " hog hyperedges", length(keep_h), " coexp hyperedges",
  nnzero(A_c) / 2, " P nnz", nnzero(A), " max hog size", max(hs), "\n"
)

# Return probabilities without forming P^k. A symmetric, P = D^-1 A.
# S = A D^-1 A (symmetric), P^2 = D^-1 S, so
#   p2(i) = S_ii / d_i
#   p3(i) = sum_j (P^2)_ij P_ji = sum_j S_ij A_ij / (d_i d_j)
#   p4(i) = sum_j (P^2)_ij (P^2)_ji = sum_j S_ij^2 / (d_i d_j)
# Same with A_c for the co-expression-only walk (same D, sub-stochastic).
DA <- Diagonal(x = dinv) %*% A
DAc <- Diagonal(x = dinv) %*% A_c
ret <- function(blk, A, DA) {
  R1 <- t(A[, blk, drop = FALSE])
  S <- R1 %*% DA
  db <- dinv[blk]
  cbind(
    p2 = S[cbind(seq_along(blk), blk)] * db,
    p3 = as.vector((S * R1) %*% dinv) * db,
    p4 = as.vector((S * S) %*% dinv) * db
  )
}
blocks <- split(seq_len(n), ceiling(seq_len(n) / 1000))
t1 <- Sys.time()
res <- parallel::mclapply(blocks, function(b) {
  cbind(ret(b, A, DA), c = ret(b, A_c, DAc))
}, mc.cores = 6)
stopifnot(!any(vapply(res, inherits, NA, "try-error")))
R <- do.call(rbind, res)
colnames(R) <- c("p2", "p3", "p4", "c2", "c3", "c4")
genes <- cbind(genes, R, d = d)
genes$h4 <- pmax(genes$p4 - genes$c4, 0) # hog-mediated 4-step return
genes$h3 <- pmax(genes$p3 - genes$c3, 0)
cat("return probabilities:", format(Sys.time() - t1), "\n")
# check against explicit P^k rows for 3 multi-copy genes
P <- Diagonal(x = dinv) %*% A
for (i in which(genes$copies > 1)[c(1, 500, 5000)]) {
  x <- P[i, , drop = FALSE]
  pk <- numeric(3)
  for (k in 2:4) {
    x <- x %*% P
    pk[k - 1] <- x[1, i]
  }
  stopifnot(isTRUE(all.equal(pk, unname(R[i, 1:3]), tolerance = 1e-8)))
}
rm(P)
saveRDS(genes, file.path(Sys.getenv("ENC_DIR"), "rw_genes.rds"))

# nodes of interest: all genes in multi-copy cells + 5000 random single-copy
set.seed(1)
interest <- c(which(genes$copies > 1), sample(which(genes$copies == 1), 5000))
cat(
  "nodes of interest", length(interest), " multi-copy genes",
  sum(genes$copies > 1), "\n"
)
print(round(sapply(
  genes[interest, c("p2", "p3", "p4", "h3", "h4")],
  quantile, c(.1, .5, .9)
), 6))

# Q1: false twins. log floor 1e-12 (h4 = 0 when no hog hyperedge)
L <- function(x) log(pmax(x, 1e-12))
ps <- pair_sets()
q1 <- t(sapply(ps, function(p) {
  c(
    qtab(abs(L(genes$p4[p[, 1]]) - L(genes$p4[p[, 2]]))),
    h = qtab(abs(L(genes$h4[p[, 1]]) - L(genes$h4[p[, 2]])))[-1],
    share_h = signif(median(genes$h4[p[, 1]] / genes$p4[p[, 1]]), 3),
    ddeg = median(abs(log(genes$d[p[, 1]]) - log(genes$d[p[, 2]])))
  )
}))
cat("\nQ1 |log p4(a) - log p4(b)|, h = hog-mediated part\n")
print(q1)

# Q2: copy with max h4 vs copy with most called edges (ties: best q)
mc <- genes[genes$copies > 1, ]
mc <- mc[ave(mc$n_called, mc$cell, FUN = max) > 0, ]
o_call <- order(mc$cell, -mc$n_called, mc$best_q)
top_call <- mc$node[o_call][!duplicated(mc$cell[o_call])]
top_by <- function(v) {
  o <- order(mc$cell, -v)
  mc$node[o][!duplicated(mc$cell[o])]
}
cells <- unique(mc$cell[o_call])
agree <- function(v) mean(top_by(v) == top_call)
q2 <- c(
  cells = length(cells),
  h4 = agree(mc$h4), p4 = agree(mc$p4), low_deg = agree(-mc$d),
  high_deg = agree(mc$d),
  baseline = mean(1 / table(mc$cell))
)
cat("\nQ2 agreement with most-called copy\n")
print(round(q2, 3))

# Q3: hog mean h4 vs hog_class
hm <- fmean(genes$h4, genes$hog)
hm <- data.frame(hog = names(hm), h4 = hm, tier = hog_tier[names(hm)])
hm <- hm[!is.na(hm$tier), ]
cat("\nQ3 hog mean h4 by hog_class\n")
print(hm |> fgroup_by(tier) |> fsummarise(n = fnobs(h4), median = fmedian(h4)))
cc <- hm[hm$tier %in% c("complete_conserved", "unclassified"), ]
cat(
  "AUROC complete_conserved vs unclassified:",
  round(auroc(cc$h4, cc$tier == "complete_conserved"), 3), "\n"
)
cat(
  "AUROC any non-unclassified tier vs unclassified:",
  round(auroc(hm$h4, hm$tier != "unclassified"), 3), "\n"
)
cat("total", format(Sys.time() - t0), "\n")
