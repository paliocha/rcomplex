# Probe WP15b: species-blocked h-Leiden. Not part of the package. Source
# into an environment whose parent is asNamespace("rcomplex"):
#   e <- new.env(parent = asNamespace("rcomplex"))
#   e$JL_DIR <- "dev/probes"
#   sys.source("dev/probes/hleiden.R", e)
#
# Objective q(alpha) = (1 - alpha) q(0) + alpha q(1), see hleiden.cpp:
# q(0) is E2's objective (joint_quality()); q(1) swaps E2's hog-pair term
# for a strict hyperedge term (a hog scores when all its copies share a
# module, paralogs collapsed to 1 / c_s) with a species-blocked null.
# Search: Leiden (Traag 2019) with the alpha schedule of h-Louvain
# (Kaminski 2024). Reuses jl_data() from joint_louvain.R.

if (!exists("JL_DIR")) JL_DIR <- "dev/probes"
sys.source(file.path(JL_DIR, "joint_louvain.R"), environment())
Rcpp::sourceCpp(file.path(JL_DIR, "hleiden.cpp"), env = environment())

# Kernel arguments: E2's plus the (node, hog) link, both 0-based
.hl_args <- function(d) {
  c(.jl_args(d), list(
    ln = d$link$node - 1L,
    lh = match(d$link$hog, unique(d$link$hog)) - 1L
  ))
}

# q(alpha) of a label vector over d's nodes, from the kernel's structures
hl_quality <- function(lab, d, lambda, gamma = 1, alpha = 0) {
  r <- do.call(hl_quality_cpp, c(
    list(lab = match(lab, unique(lab)) - 1L), .hl_args(d),
    list(gamma = gamma)
  ))
  r$quality <- r$co + lambda * ((1 - alpha) * r$hog_raw +
    alpha * r$strict_raw)
  r
}

# Modules not connected in the 2-section (co-expression plus hog edges)
hl_disconnected <- function(lab, d) {
  ij <- function(m) {
    j <- rep.int(seq_len(ncol(m)), diff(m@p))
    i <- m@i + 1L
    k <- i < j & lab[i] == lab[j]
    cbind(i[k], j[k])
  }
  e <- rbind(ij(d$co), ij(d$hog))
  g <- igraph::make_graph(as.vector(t(e)), n = d$n, directed = FALSE)
  x <- unique(data.frame(lab = lab, comp = igraph::components(g)$membership))
  length(unique(x$lab[duplicated(x$lab)]))
}

# alpha: one value (fixed) or a schedule by aggregation level, last value
# held and used for every later iteration. Same return shape as
# joint_louvain(); quality is q(last alpha).
hleiden <- function(nets, orthologs, lambda, gamma = 1, alpha = c(0, 0.5, 1),
                    theta = 0.01, seed = NULL, max_iter = 20,
                    data = jl_data(nets, orthologs)) {
  d <- data
  .seed_scope(seed)
  t0 <- Sys.time()
  r <- do.call(hl_leiden_cpp, c(.hl_args(d), list(
    lambda = lambda, gamma = gamma, alpha = as.numeric(alpha),
    theta = theta, max_iter = as.integer(max_iter)
  )))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  lab <- r$labels + 1L
  a <- alpha[length(alpha)]
  q <- hl_quality(lab, d, lambda, gamma, a)
  sp <- do.call(jl_quality_cpp, c(
    list(lab = r$labels), .jl_args(d),
    list(lambda = lambda, gamma = gamma)
  ))
  params <- list(
    lambda = lambda, gamma = gamma, alpha = alpha, theta = theta,
    seed = seed, iterations = r$iterations, converged = r$converged,
    levels = r$levels, stalls = r$stalls, secs = secs
  )
  modules <- lapply(
    stats::setNames(seq_along(d$species), d$species),
    function(k) {
      g <- d$genes[[k]]
      m <- stats::setNames(lab[d$off[k] + seq_along(g)], g)
      mg <- split(names(m), m)
      list(
        modules = m, module_genes = mg, n_modules = length(mg),
        method = "hleiden", params = params
      )
    }
  )
  list(
    modules = modules, hogs = .jl_hogs(lab, d$link), params = params,
    quality = q$quality,
    quality_parts = c(stats::setNames(sp$species, d$species),
      hog = lambda * (1 - a) * q$hog_raw, hog_raw = q$hog_raw,
      strict = lambda * a * q$strict_raw, strict_raw = q$strict_raw
    ),
    labels = lab
  )
}
