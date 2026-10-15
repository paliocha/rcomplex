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

# k-means (K centres, nstart 10) on the row-normalised embedding
bethe_kmeans <- function(ev, seed = NULL) {
  X <- ev$vectors[, seq_len(ev$K), drop = FALSE]
  nr <- sqrt(rowSums(X^2))
  X <- X / ifelse(nr > 0, nr, 1)
  .seed_scope(seed)
  if (ev$K == 1L) {
    return(rep(1L, nrow(X)))
  }
  lab <- stats::kmeans(X, ev$K, nstart = 10, iter.max = 100)$cluster
  match(lab, unique(lab))
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
  lab <- bethe_kmeans(ev, seed)
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

# ---- Hook for dev/probes/joint_halves.R ----
# bethe_global / bethe_blocked at lambda {0.25, 1, 2}, seeds 1-3 (k-means
# only). One eigen count per data set, half, variant and lambda: its
# eigenvalues go to <ds>_bethe_eig_<variant>_<lambda>_<half>.rds and a row
# to bethe_K.tsv; the vectors stay in memory for the other seeds. K at
# the cap (k_max 320) skips k-means: the run returns NULL.
BETHE_K_MAX <- 320L
bethe_hook <- function(variant) {
  list(lambdas = c(0.25, 1, 2), seeds = 1:3, run = function(ctx, lam, s) {
    g <- globalenv()
    key <- sprintf("%s_%g_%s", variant, lam, ctx$half)
    f <- file.path(g$OUT, sprintf("%s_bethe_eig_%s.rds", g$DS, key))
    if (is.null(ctx$bethe_ev)) ctx$bethe_ev <- list()
    ev <- ctx$bethe_ev[[key]]
    if (is.null(ev) && file.exists(f)) {
      ev <- readRDS(f)
      if (!ev$saturated) ev <- NULL # vectors are not on disk
    }
    if (is.null(ev)) {
      d <- ctx$dat()
      t0 <- Sys.time()
      bh <- bethe_hessian(d$co + lam * d$hog, d$spi, variant)
      ev <- bethe_negative(bh$H, k_max = BETHE_K_MAX)
      ev$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      ev$r <- bh$r
      saveRDS(ev[names(ev) != "vectors"], f)
      kf <- file.path(g$OUT, "bethe_K.tsv")
      utils::write.table(
        data.frame(
          ds = g$DS, half = ctx$half, variant = variant, lambda = lam,
          K = ev$K, saturated = ev$saturated, secs_eigs = ev$secs
        ), kf,
        sep = "\t", quote = FALSE, row.names = FALSE,
        append = file.exists(kf), col.names = !file.exists(kf)
      )
      g$log_t(
        "%s half %s: bethe %s lambda %g K %d%s, eigs %.0f s", g$DS,
        ctx$half, variant, lam, ev$K, if (ev$saturated) "+" else "",
        ev$secs
      )
      ctx$bethe_ev[[key]] <- ev
    }
    if (ev$saturated) {
      return(NULL)
    }
    t0 <- Sys.time()
    lab <- bethe_kmeans(ev, s)
    secs_km <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    list(
      lab = lab, secs = ev$secs + secs_km,
      params = list(
        variant = variant, K = ev$K, r = ev$r, eigenvalues = ev$values,
        secs_eigs = ev$secs, secs_kmeans = secs_km
      )
    )
  })
}
halves_engines <- list(
  bethe_global = bethe_hook("global"), bethe_blocked = bethe_hook("blocked")
)
