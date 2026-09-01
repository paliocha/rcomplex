# Pure-R reference implementations for validating C++ output
# Direct port of RComPlEx.Rmd logic

#' Reference MR normalization (ascending ranks, raw mutual rank)
#' Matches RComPlEx.Rmd lines 167-172
reference_mr_raw <- function(net) {
  row_ranks <- t(apply(net, 1, rank))
  mr <- sqrt(row_ranks * t(row_ranks))
  diag(mr) <- 0
  mr
}

#' Reference CLR normalization
#' Matches RComPlEx.Rmd lines 157-165
reference_clr <- function(net) {
  z <- scale(net)
  z[z < 0] <- 0
  clr <- sqrt(t(z)^2 + z^2)
  diag(clr) <- 0
  clr
}

#' Reference density threshold
reference_density_threshold <- function(net, density) {
  vals <- sort(net[upper.tri(net, diag = FALSE)], decreasing = TRUE)
  vals[round(density * length(vals))]
}

#' Reference neighborhood comparison for a single pair
#' Direct port of RComPlEx.Rmd lines 217-276
reference_compare_pair <- function(net1, net2, thr1, thr2, ortho, g1, g2) {
  # Direction 1: sp1 -> sp2
  neigh <- net1[g1, ]
  neigh <- names(neigh[neigh >= thr1])

  ortho_neigh <- net2[g2, ]
  ortho_neigh <- names(ortho_neigh[ortho_neigh >= thr2])
  ortho_neigh <- unique(ortho$Species1[ortho$Species2 %in% ortho_neigh])

  n_genes <- nrow(net1)
  m <- length(neigh)
  k <- length(ortho_neigh)
  x <- length(intersect(neigh, ortho_neigh))
  p_val_con1 <- 1
  p_val_div1 <- 1
  effect1 <- 1
  if (x > 1) {
    p_val_con1 <- phyper(x - 1, m, n_genes - m, k, lower.tail = FALSE)
  }
  if (k > 0 && m > 0) {
    effect1 <- (x / k) / (m / n_genes)
    p_val_div1 <- phyper(x, m, n_genes - m, k, lower.tail = TRUE)
  }

  # Direction 2: sp2 -> sp1
  neigh2 <- net2[g2, ]
  neigh2 <- names(neigh2[neigh2 >= thr2])

  ortho_neigh2 <- net1[g1, ]
  ortho_neigh2 <- names(ortho_neigh2[ortho_neigh2 >= thr1])
  ortho_neigh2 <- unique(ortho$Species2[ortho$Species1 %in% ortho_neigh2])

  n_genes2 <- nrow(net2)
  m2 <- length(neigh2)
  k2 <- length(ortho_neigh2)
  x2 <- length(intersect(neigh2, ortho_neigh2))
  p_val_con2 <- 1
  p_val_div2 <- 1
  effect2 <- 1
  if (x2 > 1) {
    p_val_con2 <- phyper(x2 - 1, m2, n_genes2 - m2, k2, lower.tail = FALSE)
  }
  if (k2 > 0 && m2 > 0) {
    effect2 <- (x2 / k2) / (m2 / n_genes2)
    p_val_div2 <- phyper(x2, m2, n_genes2 - m2, k2, lower.tail = TRUE)
  }

  union1 <- m + k - x
  jaccard1 <- if (union1 > 0) x / union1 else 0

  union2 <- m2 + k2 - x2
  jaccard2 <- if (union2 > 0) x2 / union2 else 0

  data.frame(
    Species1.neigh = m,
    Species1.ortho.neigh = k,
    Species1.neigh.overlap = x,
    Species1.p.val.con = p_val_con1,
    Species1.p.val.div = p_val_div1,
    Species1.effect.size = effect1,
    Species1.jaccard = jaccard1,
    Species2.neigh = m2,
    Species2.ortho.neigh = k2,
    Species2.neigh.overlap = x2,
    Species2.p.val.con = p_val_con2,
    Species2.p.val.div = p_val_div2,
    Species2.effect.size = effect2,
    Species2.jaccard = jaccard2
  )
}


# ---- shared network fixtures for sparse (dgCMatrix) dispatch tests ----

#' Sparse (dgCMatrix) copy of a dense network object
#'
#' Converts at `thr` (default: the network's own threshold) and records it as
#' `store_threshold`, as compute_network(sparse = TRUE) does.
sparse_net <- function(net, thr = net$threshold) {
  modifyList(net, list(network = dense_to_dgc(net$network, thr),
                       store_threshold = thr))
}

#' Random compute_network() pair (50 x 10 and 40 x 10, density 0.1) with a
#' 1:1 ortholog table plus one paralog row
make_cmp_nets <- function() {
  set.seed(42)
  expr1 <- matrix(rnorm(500), nrow = 50, ncol = 10)
  rownames(expr1) <- paste0("A_", sprintf("%03d", 1:50))
  expr2 <- matrix(rnorm(400), nrow = 40, ncol = 10)
  rownames(expr2) <- paste0("B_", sprintf("%03d", 1:40))

  net1 <- compute_network(expr1, density = 0.1, mr_log_transform = FALSE)
  net2 <- compute_network(expr2, density = 0.1, mr_log_transform = FALSE)

  ortho <- data.frame(
    Species1 = c(paste0("A_", sprintf("%03d", 1:30)), "A_001"),
    Species2 = c(paste0("B_", sprintf("%03d", 1:30)), "B_031"),
    hog = c(1:30, 1)
  )
  list(net1 = net1, net2 = net2, ortho = ortho)
}

#' Graded-weight network pair for tighter-than-stored threshold tests
#'
#' Weight tiers 10 / 7 / 4: a store built at threshold 5 keeps tiers 10 and
#' 7, an analysis threshold of 8 keeps tier 10 only, so the cut provably
#' drops some stored, overlapping edges and keeps others independent of RNG.
#' 30 genes per species, 1:1 orthologs, 10 HOGs of 3. HOG1 (genes 1-3) hubs
#' have conserved neighbourhoods; HOG7 (genes 19-21) connect to different
#' genes in each species (zero overlap); genes 28-30 are isolated.
make_graded_nets <- function() {
  n <- 30
  build <- function(prefix, far) {
    m <- matrix(0, n, n)
    rownames(m) <- colnames(m) <- paste0(prefix, sprintf("%02d", 1:n))
    for (g in 1:3) {
      m[g, 4:9]   <- m[4:9, g]   <- 10   # kept at thr 8
      m[g, 10:15] <- m[10:15, g] <- 7    # stored at 5, dropped at 8
      m[g, 16:18] <- m[16:18, g] <- 4    # below the store threshold
    }
    for (g in 19:21) {
      m[g, far] <- m[far, g] <- 10
    }
    diag(m) <- 10
    m
  }
  net1 <- list(network = build("A", far = 22:24), threshold = 5)
  net2 <- list(network = build("B", far = 25:27), threshold = 5)
  ortho <- data.frame(
    Species1 = paste0("A", sprintf("%02d", 1:n)),
    Species2 = paste0("B", sprintf("%02d", 1:n)),
    hog = rep(paste0("HOG", 1:10), each = 3),
    stringsAsFactors = FALSE
  )
  list(net1 = net1, net2 = net2, ortho = ortho)
}
