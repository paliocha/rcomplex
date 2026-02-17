# Pure-R reference implementations for validating C++ output
# Direct port of RComPlEx.Rmd logic

#' Reference MR normalization (ascending ranks, raw mutual rank)
#' Matches RComPlEx.Rmd lines 167-172
reference_mr_raw <- function(net) {
  R <- t(apply(net, 1, rank))
  mr <- sqrt(R * t(R))
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

  N <- nrow(net1)
  m <- length(neigh)
  k <- length(ortho_neigh)
  x <- length(intersect(neigh, ortho_neigh))
  p_val_con1 <- 1
  p_val_div1 <- 1
  effect1 <- 1
  if (x > 1) {
    p_val_con1 <- phyper(x - 1, m, N - m, k, lower.tail = FALSE)
  }
  if (k > 0 && m > 0) {
    effect1 <- (x / k) / (m / N)
    p_val_div1 <- phyper(x, m, N - m, k, lower.tail = TRUE)
  }

  # Direction 2: sp2 -> sp1
  neigh2 <- net2[g2, ]
  neigh2 <- names(neigh2[neigh2 >= thr2])

  ortho_neigh2 <- net1[g1, ]
  ortho_neigh2 <- names(ortho_neigh2[ortho_neigh2 >= thr1])
  ortho_neigh2 <- unique(ortho$Species2[ortho$Species1 %in% ortho_neigh2])

  N2 <- nrow(net2)
  m2 <- length(neigh2)
  k2 <- length(ortho_neigh2)
  x2 <- length(intersect(neigh2, ortho_neigh2))
  p_val_con2 <- 1
  p_val_div2 <- 1
  effect2 <- 1
  if (x2 > 1) {
    p_val_con2 <- phyper(x2 - 1, m2, N2 - m2, k2, lower.tail = FALSE)
  }
  if (k2 > 0 && m2 > 0) {
    effect2 <- (x2 / k2) / (m2 / N2)
    p_val_div2 <- phyper(x2, m2, N2 - m2, k2, lower.tail = TRUE)
  }

  data.frame(
    Species1.neigh = m,
    Species1.ortho.neigh = k,
    Species1.neigh.overlap = x,
    Species1.p.val.con = p_val_con1,
    Species1.p.val.div = p_val_div1,
    Species1.effect.size = effect1,
    Species2.neigh = m2,
    Species2.ortho.neigh = k2,
    Species2.neigh.overlap = x2,
    Species2.p.val.con = p_val_con2,
    Species2.p.val.div = p_val_div2,
    Species2.effect.size = effect2
  )
}
