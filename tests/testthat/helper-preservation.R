# Shared module-preservation fixture (test-module-preservation.R,
# test-as-modules.R): two species sharing module structure, and a hand-built
# module map. Gene loadings on each module's latent factor are heavy-tailed
# and SHARED between species, so hub identity is conserved and cor.degree
# has signal. A fixture where every gene in a module is exchangeable (one
# factor, iid noise) correctly yields cor.degree ~ 0 even for preserved
# modules, and would look like a bug.

pres_expr <- function(seed, n, prefix, loadings, per, n_samp = 40) {
  set.seed(seed)
  n_mod <- length(loadings)
  e <- matrix(stats::rnorm(n * n_samp), n, n_samp)
  for (k in seq_len(n_mod)) {
    f <- stats::rnorm(n_samp)
    rows <- ((k - 1) * per + 1):(k * per)
    for (j in seq_along(rows)) {
      lam <- loadings[[k]][j]
      e[rows[j], ] <- lam * f +
        stats::rnorm(n_samp, sd = sqrt(max(1e-6, 1 - lam^2)))
    }
  }
  rownames(e) <- paste0(prefix, sprintf("%04d", seq_len(n)))
  e
}

pres_fixture <- function() {
  n_mod <- 4L
  per <- 40L
  set.seed(77)
  loadings <- lapply(seq_len(n_mod), function(k) {
    l <- stats::rlnorm(per, 0, 0.9)
    l / max(l)
  })
  eA <- pres_expr(31, 300, "A", loadings, per)  # nolint
  eB <- pres_expr(32, 500, "B", loadings, per)  # nolint

  # Map background genes too, not just module genes. The mappable universe is
  # what the permutation draws from, so an ortholog table covering only module
  # genes makes the null "genes from the other modules" -- and since every
  # module is dense, avg.weight then has almost no contrast. Real ortholog
  # tables include background, and so must the fixture.
  n_map <- 300L

  list(
    n_mod = n_mod, per = per,
    netA = compute_network(eA, density = 0.03, sparse = FALSE),
    netB = compute_network(eB, density = 0.03, sparse = FALSE),
    netB_sparse = compute_network(eB,
      density = 0.03, sparse = TRUE,
      store_density = 0.03
    ),
    ortho = data.frame(
      Species1 = paste0("A", sprintf("%04d", seq_len(n_map))),
      Species2 = paste0("B", sprintf("%04d", seq_len(n_map))),
      hog = paste0("H", seq_len(n_map)),
      stringsAsFactors = FALSE
    ),
    mods = lapply(
      seq_len(n_mod),
      function(k) ((k - 1) * per + 1):(k * per)
    )
  )
}

# Module labels matching the simulated block structure, in detect_modules()
# shape, so the tests do not depend on Leiden's partition.
true_modules <- function(net, mods) {
  genes <- rownames(net$network)
  membership <- stats::setNames(rep(NA_integer_, length(genes)), genes)
  for (k in seq_along(mods)) membership[mods[[k]]] <- k
  membership <- membership[!is.na(membership)]
  list(
    modules = membership,
    module_genes = split(names(membership), membership),
    n_modules = length(mods)
  )
}
