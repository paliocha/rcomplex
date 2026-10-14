# Probe E2: joint modules with a species-blocked objective.
# Not part of the package. Source into an environment whose parent is
# asNamespace("rcomplex"), as the probe scripts do:
#   e <- new.env(parent = asNamespace("rcomplex"))
#   sys.source("dev/probes/joint_louvain.R", e)
#
# Objective (per-species modularity plus a hog term with its own null):
#   Q = sum_s (1 / 2m_s) sum_{i,j in s} (A_ij - g k_i k_j / 2m_s) d(c_i, c_j)
#     + (l / 2m_bar) sum_{i in s, j in t, s != t}
#         (w_ij - h_i^t h_j^s / H_st) d(c_i, c_j)
# A: per-species edges of .joint_edges(), scaled by the species maximum.
# w_ij = 1 / (c_s(h) c_t(h)) for copies i, j of hog h in species s != t.
# h_i^t = sum of w_ij into species t; H_st = sum of h_i^t over i in s.

# JL_DIR: directory of this file (set it in the env before sys.source)
if (!exists("JL_DIR")) JL_DIR <- "dev/probes"
Rcpp::sourceCpp(file.path(JL_DIR, "joint_louvain.cpp"), env = environment())

# Node table, co-expression and hog matrices for joint_louvain()
jl_data <- function(nets, orthologs) {
  sp <- names(nets)
  long <- .driver_orthologs(orthologs, sp)
  genes <- lapply(nets, function(n) rownames(n$network))
  off <- cumsum(c(0L, lengths(genes)))[seq_along(sp)]
  n <- sum(lengths(genes))
  spi <- rep(seq_along(sp) - 1L, lengths(genes))
  e <- do.call(rbind, lapply(seq_along(sp), function(k) {
    e <- .joint_edges(nets[[k]], sp[k])
    e$w <- e$w / max(e$w)
    e$i <- e$i + off[k]
    e$j <- e$j + off[k]
    e
  }))
  co <- Matrix::sparseMatrix(
    i = c(e$i, e$j), j = c(e$j, e$i), x = c(e$w, e$w), dims = c(n, n)
  )
  rm(e)
  node <- rep(NA_integer_, nrow(long))
  for (k in seq_along(sp)) {
    r <- long$species == sp[k]
    node[r] <- off[k] + match(long$gene[r], genes[[k]])
  }
  link <- unique(data.frame(node = node, hog = long$hog)[!is.na(node), ])
  link$s <- spi[link$node]
  link$c_s <- stats::ave(link$node, link$hog, link$s, FUN = length)
  dt <- data.table::as.data.table(link)
  pr <- merge(dt, dt, by = "hog", allow.cartesian = TRUE)
  pr <- pr[pr$s.x != pr$s.y, ]
  w <- 1 / (pr$c_s.x * pr$c_s.y)
  hog <- Matrix::sparseMatrix(
    i = pr$node.x, j = pr$node.y, x = w, dims = c(n, n)
  )
  h <- as.matrix(Matrix::sparseMatrix(
    i = pr$node.x, j = pr$s.y + 1L, x = w, dims = c(n, length(sp))
  ))
  H <- rowsum(h, spi, reorder = TRUE)
  list(
    species = sp, genes = genes, off = off, n = n, spi = spi,
    co = co, hog = hog, h = h, H = unname(H),
    link = link[, c("node", "hog")]
  )
}

# The C++ arguments shared by the kernel and the quality
.jl_args <- function(d) {
  list(
    sp = d$spi, S = length(d$species),
    co_p = d$co@p, co_i = d$co@i, co_x = d$co@x,
    hog_p = d$hog@p, hog_i = d$hog@i, hog_x = d$hog@x,
    h = d$h, H = d$H
  )
}

# Joint Louvain. init: named list (species) of named label vectors over
# genes, with freeze = TRUE (two-step + match: per-species modules become
# supernodes; only supernodes move). Pass a prepared `data` from
# jl_data() to skip the build.
joint_louvain <- function(nets, orthologs, lambda, gamma = 1, seed = NULL,
                          max_passes = 50, tol = 1e-10, init = NULL,
                          freeze = FALSE, data = jl_data(nets, orthologs)) {
  d <- data
  if (freeze != !is.null(init)) {
    stop("`init` and `freeze = TRUE` go together.")
  }
  init0 <- integer(0)
  if (freeze) {
    key <- rep(NA_character_, d$n)
    for (k in seq_along(d$species)) {
      s <- d$species[k]
      m <- init[[s]][d$genes[[k]]]
      key[d$off[k] + seq_along(d$genes[[k]])] <- ifelse(
        is.na(m), NA, paste(s, m)
      )
    }
    solo <- is.na(key)
    key[solo] <- paste0("solo", which(solo))
    init0 <- match(key, unique(key)) - 1L
  }
  .seed_scope(seed)
  t0 <- Sys.time()
  r <- do.call(jl_louvain_cpp, c(.jl_args(d), list(
    lambda = lambda, gamma = gamma, init = init0, freeze = freeze,
    max_passes = as.integer(max_passes), tol = tol
  )))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  lab <- r$labels + 1L
  q <- do.call(jl_quality_cpp, c(
    list(lab = r$labels), .jl_args(d),
    list(lambda = lambda, gamma = gamma)
  ))
  params <- list(
    lambda = lambda, gamma = gamma, seed = seed, freeze = freeze,
    levels = r$levels, refine_moves = r$refine_moves, secs = secs
  )
  modules <- lapply(
    stats::setNames(seq_along(d$species), d$species),
    function(k) {
      g <- d$genes[[k]]
      m <- stats::setNames(lab[d$off[k] + seq_along(g)], g)
      mg <- split(names(m), m)
      list(
        modules = m, module_genes = mg, n_modules = length(mg),
        method = "joint_louvain", params = params
      )
    }
  )
  list(
    modules = modules, hogs = .jl_hogs(lab, d$link), params = params,
    quality = q$quality,
    quality_parts = c(stats::setNames(q$species, d$species),
      hog = q$hog, hog_raw = q$hog_raw
    ),
    labels = lab
  )
}

# Hog table in the joint_modules() shape: hogs with two or more copies
.jl_hogs <- function(lab, link) {
  hogs <- sort(unique(link$hog[duplicated(link$hog)]), method = "radix")
  link <- link[link$hog %in% hogs, , drop = FALSE]
  by_hog <- split(lab[link$node], factor(link$hog, levels = hogs))
  main <- vapply(by_hog, function(x) {
    u <- sort(unique(x))
    n <- tabulate(match(x, u))
    c(u[which.max(n)], max(n))
  }, integer(2))
  data.frame(
    hog = hogs, n_copies = lengths(by_hog, use.names = FALSE),
    n_modules = vapply(by_hog, function(x) length(unique(x)), 1L,
      USE.NAMES = FALSE
    ),
    module_main = main[1L, ], copies_main = main[2L, ],
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# The objective in plain R, from a fit or a label vector over d's nodes
joint_quality <- function(x, data, lambda = x$params$lambda,
                          gamma = if (is.list(x)) x$params$gamma else 1) {
  d <- data
  lab <- if (is.list(x)) x$labels else x
  S <- length(d$species)
  co <- Matrix::summary(d$co)
  k <- Matrix::colSums(d$co)
  same <- lab[co$i] == lab[co$j]
  parts <- vapply(seq_len(S), function(s) {
    v <- d$spi == s - 1L
    tm <- sum(k[v])
    vol <- tapply(k[v], lab[v], sum)
    r <- d$spi[co$i] == s - 1L
    sum(co$x[r & same]) / tm - gamma * sum(vol^2) / tm^2
  }, 0)
  tm_bar <- mean(vapply(seq_len(S), function(s) sum(k[d$spi == s - 1L]), 0))
  hg <- Matrix::summary(d$hog)
  hin <- sum(hg$x[lab[hg$i] == lab[hg$j]])
  nul <- 0
  for (s in seq_len(S)) {
    for (t in seq_len(S)) {
      if (s == t || d$H[s, t] == 0) next
      a <- rowsum(d$h[d$spi == s - 1L, t], lab[d$spi == s - 1L])
      b <- rowsum(d$h[d$spi == t - 1L, s], lab[d$spi == t - 1L])
      cm <- intersect(rownames(a), rownames(b))
      nul <- nul + sum(a[cm, 1] * b[cm, 1]) / d$H[s, t]
    }
  }
  hog_raw <- (hin - nul) / tm_bar
  list(
    quality = sum(parts) + lambda * hog_raw, species = parts,
    hog_raw = hog_raw
  )
}
