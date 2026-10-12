# Shared setup for encodings-rw.R and encodings-orc.R: 12k leaf set,
# networks, twins, tiers. Design note section 9.6. Run from the repository root.
# Inputs, read from the directory in env var ENC_DIR: leaf_se.rds,
# E_12000.rds, g_12000.rds, cls_B_12000.rds.
# Writes nets12k.rds to ENC_DIR on first run. Not package code.
# Shared setup for rw.R and orc.R: 12k leaf set, networks, twins, tiers.
suppressPackageStartupMessages({
  library(Matrix)
  library(collapse)
  library(SummarizedExperiment)
  library(rcomplex)
})
leaf <- lapply(readRDS(file.path(Sys.getenv("ENC_DIR"), "leaf_se.rds")), function(s) s[seq_len(12000), ])
sp_names <- names(leaf)
if (!file.exists(file.path(Sys.getenv("ENC_DIR"), "nets12k.rds"))) {
  nets <- lapply(leaf, function(s) compute_network(assay(s, "vst"), n_cores = 8))
  saveRDS(nets, file.path(Sys.getenv("ENC_DIR"), "nets12k.rds"))
}
nets <- readRDS(file.path(Sys.getenv("ENC_DIR"), "nets12k.rds"))
# binary adjacency: neighbor_lists.h keeps value >= thr, row != col
adj <- lapply(nets, function(n) {
  a <- n$network >= n$threshold
  diag(a) <- FALSE
  as(drop0(a) * 1, "generalMatrix")
})
genes <- do.call(rbind, lapply(sp_names, function(s) {
  data.frame(
    species = s, gene = rownames(leaf[[s]]), hog = rowData(leaf[[s]])$hog
  )
}))
genes$node <- paste(genes$species, genes$gene, sep = "|")
stopifnot(all(sapply(sp_names, function(s) {
  identical(rownames(nets[[s]]$network), rownames(leaf[[s]]))
})))

# called edges and twin classes (as enc/twins2.R)
E <- readRDS(file.path(Sys.getenv("ENC_DIR"), "E_12000.rds"))
Ec <- E[!is.na(E$q_value) & E$q_value < 0.1, ]
a <- paste(Ec$hog, Ec$species1, Ec$gene1, sep = "|")
b <- paste(Ec$hog, Ec$species2, Ec$gene2, sep = "|")
gnodes <- unique(c(a, b))
ai <- match(a, gnodes)
bi <- match(b, gnodes)
nb <- split(c(bi, ai), c(ai, bi))
sig <- vapply(nb, function(x) rlang::hash(sort(unique(x))), "")
twin <- match(sig, unique(sig))[match(seq_along(gnodes), as.integer(names(nb)))]
# map to genes table (twin NA = no called edge)
genes$twin <- twin[match(paste(genes$hog, genes$species, genes$gene, sep = "|"), gnodes)]
genes$n_called <- tabulate(
  c(
    match(paste(Ec$species1, Ec$gene1, sep = "|"), genes$node),
    match(paste(Ec$species2, Ec$gene2, sep = "|"), genes$node)
  ),
  nbins = nrow(genes)
)
bq <- c(Ec$q_value, Ec$q_value)
bn <- c(
  match(paste(Ec$species1, Ec$gene1, sep = "|"), genes$node),
  match(paste(Ec$species2, Ec$gene2, sep = "|"), genes$node)
)
genes$best_q <- Inf
mq <- fmin(bq, bn)
genes$best_q[as.integer(names(mq))] <- mq

# cell = (hog, species); copies
genes$cell <- paste(genes$hog, genes$species, sep = "|")
genes$copies <- as.vector(table(genes$cell)[genes$cell])

# hog tier (hog_class; one per hog?)
cls <- readRDS(file.path(Sys.getenv("ENC_DIR"), "cls_B_12000.rds"))
ht <- unique(cls[, c("hog", "hog_class")])
stopifnot(!anyDuplicated(ht$hog))
hog_tier <- setNames(ht$hog_class, ht$hog)

# same-species pair sets among gene-graph nodes of multi-copy cells
pair_sets <- function(seed = 1) {
  g <- genes[!is.na(genes$twin) & genes$copies > 1, ]
  sp <- split(seq_len(nrow(g)), g$cell)
  sp <- sp[lengths(sp) > 1]
  pr <- do.call(rbind, lapply(sp, function(ix) t(combn(ix, 2))))
  i <- match(g$node[pr[, 1]], genes$node)
  j <- match(g$node[pr[, 2]], genes$node)
  tw <- genes$twin[i] == genes$twin[j]
  set.seed(seed)
  rnd <- do.call(rbind, lapply(sp_names, function(s) {
    ix <- which(genes$species == s & !is.na(genes$twin))
    m <- cbind(sample(ix, 5000, TRUE), sample(ix, 5000, TRUE))
    m[m[, 1] != m[, 2], ]
  }))
  list(twin = cbind(i, j)[tw, ], nontwin = cbind(i, j)[!tw, ], random = rnd)
}

auroc <- function(x, pos) { # P(score_pos > score_neg), rank formula
  r <- rank(x)
  n1 <- sum(pos)
  n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
qtab <- function(v) {
  c(
    n = length(v), median = signif(median(v), 3),
    q90 = signif(quantile(v, .9, names = FALSE), 3)
  )
}
