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

# Connected component of every node in its module's 2-section subgraph
# (co-expression plus hog edges inside the module)
.hl_comp <- function(lab, d) {
  ij <- function(m) {
    j <- rep.int(seq_len(ncol(m)), diff(m@p))
    i <- m@i + 1L
    k <- i < j & lab[i] == lab[j]
    cbind(i[k], j[k])
  }
  e <- rbind(ij(d$co), ij(d$hog))
  g <- igraph::make_graph(as.vector(t(e)), n = d$n, directed = FALSE)
  igraph::components(g)$membership
}

# Modules not connected in the 2-section
hl_disconnected <- function(lab, d) {
  x <- unique(data.frame(lab = lab, comp = .hl_comp(lab, d)))
  length(unique(x$lab[duplicated(x$lab)]))
}

# Split every disconnected module into its components. The largest
# component keeps the module; n_moved counts genes in the others.
hl_split <- function(lab, d) {
  comp <- .hl_comp(lab, d)
  nc <- tapply(comp, lab, function(x) length(unique(x)))
  size <- stats::ave(comp, comp, FUN = length)
  keep <- tapply(seq_along(lab), lab, function(i) {
    comp[i][which.max(size[i])]
  })
  moved <- comp != keep[as.character(lab)]
  key <- ifelse(moved, paste0("c", comp), paste0("m", lab))
  list(
    labels = match(key, unique(key)), n_split = sum(nc > 1),
    n_moved = sum(moved)
  )
}

# alpha: one value (fixed) or a schedule by aggregation level, last value
# held and used for every later iteration. Modules disconnected in the
# 2-section are split afterwards (params n_split, n_moved); with the
# strict term (alpha > 0) Leiden cannot guarantee connectivity. Same
# return shape as joint_louvain(); quality is q(last alpha), after split.
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
  sp0 <- hl_split(r$labels + 1L, d)
  lab <- sp0$labels
  a <- alpha[length(alpha)]
  q <- hl_quality(lab, d, lambda, gamma, a)
  sp <- do.call(jl_quality_cpp, c(
    list(lab = lab - 1L), .jl_args(d),
    list(lambda = lambda, gamma = gamma)
  ))
  params <- list(
    lambda = lambda, gamma = gamma, alpha = alpha, theta = theta,
    seed = seed, iterations = r$iterations, converged = r$converged,
    levels = r$levels, stalls = r$stalls, n_split = sp0$n_split,
    n_moved = sp0$n_moved, secs = secs
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

# Hook for joint_halves.R: hleiden0 (alpha 0) and hleidenS (schedule).
# HL_ENGINES (comma list) restricts the set, so the two can run as
# separate processes.
halves_engines <- lapply(
  list(hleiden0 = 0, hleidenS = c(0, 0.5, 1)),
  function(a) {
    force(a)
    list(
      lambdas = c(0.25, 1, 4, 16, 64), seeds = 1:3,
      run = function(ctx, lambda, seed) {
        f <- hleiden(lambda = lambda, alpha = a, seed = seed,
                     data = ctx$dat())
        list(lab = f$labels, secs = f$params$secs, quality = f$quality,
             quality_parts = f$quality_parts, params = f$params)
      }
    )
  }
)
if (nzchar(Sys.getenv("HL_ENGINES"))) {
  halves_engines <- halves_engines[
    strsplit(Sys.getenv("HL_ENGINES"), ",", fixed = TRUE)[[1]]
  ]
}
