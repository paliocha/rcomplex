# Probe E3: species-blocked Bethe Hessian spectral modules.
# Not part of the package. Source into an environment whose parent is
# asNamespace("rcomplex"), after dev/probes/joint_louvain.R (jl_data()):
#   e <- new.env(parent = asNamespace("rcomplex"))
#   sys.source("dev/probes/joint_louvain.R", e)
#   sys.source("dev/probes/bethe.R", e)
#
# Graph W = co + lambda * hog from jl_data(). Bethe Hessian (Saade,
# Krzakala & Zdeborova 2014; weighted form, r = sqrt(rho), Dall'Amico,
# Couillet & Tremblay 2021):
#   H(r) = (r^2 - 1) I - r W + D,  D = diag(k),  rho = sum k^2 / sum k - 1
# global: one r for the graph. blocked: r_s per species from rho_s over
# that species' rows of W (co block plus its hog rows); H_ij uses
# sqrt(r_i r_j), the diagonal uses r_i. K = number of negative
# eigenvalues; labels = k-means (K centres, nstart 10) on the
# row-normalised eigenvectors of the negative eigenvalues.

bethe_hessian <- function(W, spi, variant = c("global", "blocked")) {
  variant <- match.arg(variant)
  k <- Matrix::rowSums(W)
  rho <- if (variant == "global") {
    sum(k^2) / sum(k) - 1
  } else {
    as.numeric(rowsum(k^2, spi, reorder = TRUE) /
      rowsum(k, spi, reorder = TRUE)) - 1
  }
  r <- sqrt(rho)
  ri <- if (variant == "global") rep(r, length(k)) else r[spi + 1L]
  sq <- Matrix::Diagonal(x = sqrt(ri))
  H <- Matrix::Diagonal(x = ri^2 - 1 + k) - sq %*% W %*% sq
  list(H = methods::as(H, "generalMatrix"), r = r, rho = rho)
}

# Negative eigenvalues of H; k doubles while every returned one is
# negative (up to k_max)
bethe_negative <- function(H, k0 = 40L, k_max = 1280L) {
  k <- min(k0, nrow(H) - 1L)
  repeat {
    e <- RSpectra::eigs_sym(H,
      k = k, which = "SA",
      opts = list(maxitr = 5000, retvec = TRUE)
    )
    neg <- e$values < 0
    if (!all(neg) || k >= k_max || k >= nrow(H) - 1L) break
    k <- min(2L * k, k_max, nrow(H) - 1L)
  }
  o <- order(e$values)
  list(
    values = e$values[o], vectors = e$vectors[, o, drop = FALSE],
    K = sum(neg), k = k, nconv = e$nconv, saturated = all(neg)
  )
}

bethe <- function(nets, orthologs, lambda,
                  variant = c("global", "blocked"), k0 = 40L, seed = NULL,
                  data = jl_data(nets, orthologs)) {
  variant <- match.arg(variant)
  d <- data
  t0 <- Sys.time()
  bh <- bethe_hessian(d$co + lambda * d$hog, d$spi, variant)
  ev <- bethe_negative(bh$H, k0)
  K <- ev$K
  if (K < 1L) stop("`lambda` = ", lambda, ": no negative eigenvalue.")
  X <- ev$vectors[, seq_len(K), drop = FALSE]
  nr <- sqrt(rowSums(X^2))
  X <- X / ifelse(nr > 0, nr, 1)
  .seed_scope(seed)
  lab <- if (K == 1L) {
    rep(1L, d$n)
  } else {
    stats::kmeans(X, K, nstart = 10, iter.max = 100)$cluster
  }
  lab <- match(lab, unique(lab))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  q <- do.call(jl_quality_cpp, c(
    list(lab = lab - 1L), .jl_args(d), list(lambda = lambda, gamma = 1)
  ))
  params <- list(
    lambda = lambda, variant = variant, seed = seed, K = K,
    r = stats::setNames(bh$r, if (variant == "blocked") d$species),
    eigenvalues = ev$values, k = ev$k, saturated = ev$saturated,
    secs = secs
  )
  modules <- lapply(
    stats::setNames(seq_along(d$species), d$species),
    function(k) {
      g <- d$genes[[k]]
      m <- stats::setNames(lab[d$off[k] + seq_along(g)], g)
      mg <- split(names(m), m)
      list(
        modules = m, module_genes = mg, n_modules = length(mg),
        method = paste0("bethe_", variant), params = params
      )
    }
  )
  list(
    modules = modules, hogs = .jl_hogs(lab, d$link), params = params,
    quality = q$quality, labels = lab
  )
}
