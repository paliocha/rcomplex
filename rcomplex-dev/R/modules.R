#' Detect co-expression modules in a network
#'
#' Applies the Leiden algorithm (Traag *et al.*, 2019) to a thresholded
#' co-expression network: modularity or CPM optimization with guaranteed
#' well-connected communities. The resolution parameter controls module
#' granularity.
#'
#' @section Reproducibility at small sample sizes:
#' Modules can be stable across seeds and still not replicate across
#' independent samples, and a correlation network from few samples has
#' cluster structure even on noise. Measured on Pooideae leaf networks
#' (2,000 genes, 20 samples) with \code{resolution = c(0.5, 1, 2)}; the
#' scripts and tables are in the package's source repository under
#' \file{dev/bench/}:
#' \itemize{
#'   \item \emph{Replication} (samples split into two halves of 10, two
#'     replicates per time point each): with
#'     \code{objective_function = "CPM"}, whether modules exist at
#'     all flipped between the halves -- in 4 of 8 species one half gave
#'     6--10 modules and the other 1, and two species gave 1 in both
#'     (median split-half adjusted Rand index 0). With
#'     \code{"modularity"} the modules replicated partially (median ARI
#'     0.36, range 0.05--0.44; about 0 on shuffled expression).
#'   \item \emph{Noise}: on random and per-gene shuffled expression (18 data
#'     sets, all 20 samples), CPM returned one module, but on 10-sample
#'     halves it gave 10 modules on shuffled expression in 3 of 8 species.
#'     Modularity returned 8--14 modules on noise.
#' }
#' These CPM measurements predate 0.3.2 and ran on raw mutual-rank edge
#' weights, of the order of the gene count, so a resolution of 0.5--2 was
#' effectively 0 and CPM collapsed to one module whenever the graph was
#' connected. Since 0.3.2 CPM weights are divided by the maximum edge
#' weight, so \code{resolution} is a density on the 0--1 weight scale,
#' and the default resolution is the network's edge density; the CPM
#' figures above describe the old scale, not the current one.
#'
#' The default objective is modularity since 0.3.2 (CPM before). On
#' 20,000-gene Pooideae leaf and wood networks, CPM at the density default
#' and modularity replicated equally between sample halves and conserved
#' their large modules equally across species; CPM only added a tail of
#' small communities that were not conserved. The 2,000-gene
#' \file{dev/bench/} numbers above predate the CPM scale fix.
#'
#' Treat modules as units only after they replicate on independent samples,
#' well above the same split on shuffled expression; \code{"modularity"}
#' replicated better here but needs that check most. To test gene sets from
#' other sources (pathways, regulons, another tool) use [as_modules()] with
#' [module_preservation()].
#'
#' @param net Network object from [compute_network()].
#' @param resolution Resolution parameter for Leiden. \code{NULL}
#'   (default) runs a single resolution: for CPM the edge density of the
#'   thresholded graph, \code{ecount / choose(vcount, 2)}, so modules must
#'   be denser than the network average; for modularity 1. The resolved
#'   value is returned in \code{params$resolution}. Pass a numeric vector
#'   (e.g., \code{seq(0.5, 2.0, by = 0.5)} for modularity) to run at
#'   multiple resolutions and produce a consensus partition via
#'   co-classification (Lancichinetti & Fortunato, 2012); consensus mode
#'   has no default. Under CPM it is read on
#'   the 0--1 weight scale (edge weights divided by their maximum, since
#'   0.3.2): a density each module must exceed, so values below 1 are the
#'   useful range (at 1 no edge is attractive and genes stay singletons).
#'   Earlier versions read it on raw mutual-rank weights with a default of
#'   1, where a resolution of 1 or more returned one module containing
#'   every gene. On 20,000-gene Pooideae and wood networks the density
#'   default gave 13--63 modules per species: a few large ones (over 1,000
#'   genes) that replicated between sample halves and were conserved
#'   across species as well as modularity's, plus a tail of small
#'   communities (under 300 genes) that were not conserved. Drop the tail
#'   by size before reading modules as units.
#' @param objective_function Leiden objective: `"modularity"` (default
#'   since 0.3.2; CPM before) or `"CPM"`. The
#'   two replicated and conserved large modules equally on 20,000-gene
#'   data, and modularity has no tail of small non-conserved communities;
#'   see the section on reproducibility. CPM sees edge weights
#'   rescaled to 0--1; modularity is scale-invariant and sees them
#'   unchanged. At small sample sizes
#'   `"modularity"` replicated better across independent samples but also
#'   finds modules in noise; see the section on reproducibility.
#' @param seed Random seed for reproducibility (default `NULL`).
#'
#'   With a seed, the result is reproducible and identical at any `n_cores`:
#'   every parallel task derives its own RNG stream from the seed and its task
#'   index, so the answer does not depend on how the work was distributed.
#'   Two limits are worth knowing. Reproducibility holds for one machine and
#'   one igraph build -- the partition is a function of what
#'   `igraph::cluster_leiden()` draws, so an igraph upgrade may move it. And
#'   the answer depends on `RNGkind` as well as on `seed`: a session that has
#'   set `RNGkind("L'Ecuyer-CMRG")` gets a different, still reproducible,
#'   partition. `RNGkind` is deliberately not pinned inside the function.
#'
#'   With an explicit seed the caller's stream is restored on exit, so the
#'   call is invisible to anything drawn afterwards: the clustering backend
#'   draws from a private stream started at `seed`. With `seed = NULL` the
#'   stream advances by exactly what was taken from it -- the one draw used
#'   to pick a root in consensus mode, whatever the backend consumed in
#'   single-resolution mode -- so consecutive unseeded calls still differ.
#' @param n_cores Number of parallel cores (default 1). Used for
#'   \code{mclapply} Leiden sweeps on Unix and OpenMP edge scans in C++.
#'   Uses fork-based parallelism; avoid combining with active CUDA
#'   contexts in the same session.
#' @param max_consensus_iter Maximum number of consensus iterations for
#'   consensus mode. Default 10.
#'   Iteration stops when the sweep reproduces its own input (a fixed
#'   point) or when all resolutions agree. On the eight Pooideae networks
#'   that stopping rule was measured on it fired at 5--12 iterations, so
#'   the default cap of 10 cuts two of those eight off before it fires:
#'   at the default the loop may return a partition it had not finished
#'   settling. Raise it (20 is above every measured stop) when the
#'   returned \code{params$n_consensus_iterations} equals
#'   \code{max_consensus_iter} -- a truncated run always does, though a
#'   run that settles on the last allowed iteration does too.
#'
#' @return A list with components:
#'   \describe{
#'     \item{modules}{Named integer vector of module assignments
#'       (gene -> module ID)}
#'     \item{module_genes}{Named list: module ID -> character vector of
#'       gene names}
#'     \item{n_modules}{Number of modules detected}
#'     \item{modularity}{Modularity score of the partition}
#'     \item{graph}{The igraph graph object used for community detection}
#'     \item{method}{Method used}
#'     \item{params}{List of parameters used}
#'   }
#'   When \code{resolution} is a vector, the output also includes
#'   \code{resolution_scan}, a data frame with columns \code{resolution},
#'   \code{n_modules}, \code{modularity}, \code{ari_next} (Adjusted
#'   Rand Index with the next resolution; NA for the last), and
#'   \code{expected_coclassification} (per-resolution expected scalar).
#'   The \code{params} list includes \code{n_consensus_iterations}
#'   (number of iterations until convergence).
#'
#' @details The adaptive consensus path uses sparse co-classification
#'   restricted to the original network's edge set, reducing memory from
#'   O(N^2) to O(|E|).
#'
#' @references
#' Traag, V. A., Waltman, L. & van Eck, N. J. (2019). From Louvain to Leiden:
#' guaranteeing well-connected communities. *Scientific Reports*, 9, 5233.
#' \doi{10.1038/s41598-019-41695-z}
#'
#' Lancichinetti, A. & Fortunato, S. (2012). Consensus clustering in complex
#' networks. *Scientific Reports*, 2, 336. \doi{10.1038/srep00336}
#'
#' Jeub, L. G. S., Sporns, O. & Fortunato, S. (2018). Multiresolution
#' consensus clustering in networks. *Scientific Reports*, 8, 3259.
#' \doi{10.1038/s41598-018-21352-7}
#'
#' Senbabaoglu, Y. et al. (2014). Critical limitations of consensus
#' clustering in class discovery. *Scientific Reports*, 4, 6207.
#' \doi{10.1038/srep06207}
#'
#' @examples
#' \dontrun{
#' # Single resolution (modularity at resolution 1)
#' mods <- detect_modules(net)
#' table(mods$modules) # module sizes
#'
#' # Multi-resolution consensus (Jeub et al. 2018)
#' mods_consensus <- detect_modules(net,
#'   resolution = c(0.5, 1.0, 2.0), objective_function = "modularity",
#'   n_cores = 4L
#' )
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
detect_modules <- function(net, ...) UseMethod("detect_modules")

#' @rdname detect_modules
#' @export
detect_modules.default <- function(net,
                                   resolution = NULL,
                                   objective_function = c("modularity", "CPM"),
                                   seed = NULL,
                                   n_cores = 1L,
                                   max_consensus_iter = 10L, ...) {
  objective_function <- match.arg(objective_function)

  # Consensus mode: vector resolution triggers multi-resolution + consensus
  if (length(resolution) > 1L) {
    return(detect_modules_consensus(
      net, resolution,
      objective_function, seed, as.integer(n_cores),
      as.integer(max_consensus_iter)
    ))
  }

  if (!is.list(net) || is.null(net$network)) {
    stop("net must be a network object from compute_network()")
  }

  mat <- .net_check(net, net$threshold)
  thr <- net$threshold

  # Build thresholded adjacency (sparse stays sparse; igraph accepts a
  # dgCMatrix directly)
  if (.net_is_sparse(net)) {
    adj <- mat
    adj@x[adj@x < thr] <- 0
    adj <- Matrix::drop0(adj)
    has_edges <- length(adj@x) > 0L
  } else {
    adj <- mat
    adj[adj < thr] <- 0
    has_edges <- any(adj[upper.tri(adj)] > 0)
  }

  if (!has_edges) {
    stop("No edges above threshold; cannot detect modules")
  }

  # A seeded call runs cluster_leiden() on a private stream and restores the
  # caller's, so a downstream set.seed()-free draw -- summarize_comparison()'s
  # randomized-p pi0, for one -- continues from the caller's own seed
  # rather than from wherever the clustering backend happened to stop.
  # See .seed_scope() in R/rng.R for the package-wide contract.
  .seed_scope(seed)

  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "upper", weighted = TRUE, diag = FALSE
  )

  if (is.null(resolution)) {
    # CPM: modules denser than the network average (unit weight scale)
    resolution <- if (objective_function == "CPM") {
      igraph::ecount(g) / choose(igraph::vcount(g), 2)
    } else {
      1
    }
  }
  comm <- igraph::cluster_leiden(
    g,
    resolution = resolution,
    objective_function = objective_function,
    weights = .leiden_weights(g, objective_function),
    n_iterations = 2L
  )
  params <- list(
    resolution = resolution,
    objective_function = objective_function,
    seed = seed
  )

  membership <- igraph::membership(comm)
  names(membership) <- igraph::V(g)$name
  module_genes <- split(names(membership), membership)

  list(
    modules = membership,
    module_genes = module_genes,
    n_modules = length(module_genes),
    modularity = igraph::modularity(g, membership),
    graph = g,
    method = "leiden",
    params = params
  )
}


#' Leiden edge weights on the CPM scale (internal)
#'
#' CPM quality is sum_ij (A_ij - gamma) delta(c_i, c_j), so gamma is read on
#' the weight scale. Raw mutual-rank weights are of the order of the gene
#' count, which made any gamma near 1 effectively 0 and every edge
#' attractive (one module). CPM weights are divided by the graph's maximum
#' so gamma is a density in [0, 1]. Modularity is scale-invariant: NULL
#' keeps igraph's default (the weight attribute), bit for bit.
#'
#' @noRd
.leiden_weights <- function(g, objective_function) {
  if (objective_function != "CPM") return(NULL)
  w <- igraph::E(g)$weight
  w / max(w)
}


#' Multi-resolution consensus module detection (internal)
#'
#' Runs Leiden at each resolution, builds a co-classification matrix,
#' subtracts per-pair expected co-classification (Jeub et al. 2018),
#' and iterates until the partition converges.
#'
#' @noRd
detect_modules_consensus <- function(net, resolutions,
                                     objective_function, seed,
                                     n_cores = 1L,
                                     max_consensus_iter = 10L) {
  if (!is.list(net) || is.null(net$network)) {
    stop("net must be a network object from compute_network()")
  }

  mat <- .net_check(net, net$threshold)
  genes <- rownames(mat)
  n_genes <- length(genes)

  # Build thresholded adjacency (sparse stays sparse)
  if (.net_is_sparse(net)) {
    adj <- mat
    adj@x[adj@x < net$threshold] <- 0
    adj <- Matrix::drop0(adj)
    has_edges <- length(adj@x) > 0L
  } else {
    adj <- mat
    adj[adj < net$threshold] <- 0
    has_edges <- any(adj[upper.tri(adj)] > 0)
  }

  if (!has_edges) {
    stop("No edges above threshold; cannot detect modules")
  }

  # The per-task set.seed() calls below run in the caller's session on the
  # serial path (no fork) but not under mclapply, so the scope below restores
  # the ambient stream on exit and leaves detect_modules() looking the same at
  # any n_cores. With seed = NULL the root is drawn from the ambient stream
  # first, so that one draw -- and only that one -- is what the caller sees,
  # and consecutive unseeded calls still differ.
  seed_root <- if (is.null(seed)) {
    sample.int(.Machine$integer.max, 1L)
  } else {
    as.integer(seed)
  }
  .seed_scope(seed_root)

  # Build original graph — then free the dense adjacency (~4.6 GB for N=24k)
  g <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "upper", weighted = TRUE, diag = FALSE
  )
  rm(adj)
  gc()

  resolutions <- sort(resolutions)
  n_res <- length(resolutions)

  # If only one resolution after dedup, fall back to single-resolution
  if (n_res == 1L) {
    return(detect_modules(net,
      resolution = resolutions,
      objective_function = objective_function, seed = seed
    ))
  }

  use_mc <- .can_fork(n_cores)

  # ---- Initial Leiden sweep on original graph ----
  run_initial <- function(ri) {
    set.seed(.task_seed(seed_root, 1L, ri))
    res <- resolutions[[ri]]
    comm <- igraph::cluster_leiden(
      g,
      resolution = res,
      objective_function = objective_function,
      weights = .leiden_weights(g, objective_function),
      n_iterations = 2L
    )
    mem <- igraph::membership(comm)
    names(mem) <- igraph::V(g)$name
    list(
      mem = mem,
      n_mod = length(unique(mem)),
      quality = igraph::modularity(g, mem)
    )
  }

  if (use_mc) {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
    Sys.setenv(OMP_NUM_THREADS = 1L)
    on.exit(
      {
        if (is.na(old_omp)) {
          Sys.unsetenv("OMP_NUM_THREADS")
        } else {
          Sys.setenv(OMP_NUM_THREADS = old_omp)
        }
      },
      add = TRUE
    )
    # one fork per task (not per core chunk), so a failure names its task
    results <- parallel::mclapply(seq_len(n_res), run_initial,
      mc.cores = n_cores, mc.preschedule = FALSE
    )
  } else {
    results <- lapply(seq_len(n_res), run_initial)
  }

  .check_fork_results(results, resolutions, "resolution")

  memberships <- lapply(results, `[[`, "mem")
  scan_n_modules <- vapply(results, `[[`, integer(1), "n_mod")
  scan_modularity <- vapply(results, `[[`, numeric(1), "quality")

  # ARI between consecutive resolutions
  scan_ari_next <- rep(NA_real_, n_res)
  for (r in seq_len(n_res - 1L)) {
    scan_ari_next[r] <- igraph::compare(
      memberships[[r]], memberships[[r + 1L]],
      method = "adjusted.rand"
    )
  }

  # ---- Extract edge list for sparse co-classification ----
  edge_list_0 <- igraph::as_edgelist(g, names = FALSE) - 1L
  storage.mode(edge_list_0) <- "integer"

  # ---- Build resolution scan (expected filled during consensus) ----
  resolution_scan <- data.frame(
    resolution = resolutions,
    n_modules = scan_n_modules,
    modularity = scan_modularity,
    ari_next = scan_ari_next,
    expected_coclassification = NA_real_
  )

  # ---- Consensus clustering ----
  n_consensus_iter <- 0L
  # Save initial sweep for fallback if consensus graph becomes empty
  initial_memberships <- memberships

  # Adaptive: iterate until convergence (Jeub et al. 2018, Algorithm 1)
  for (iter in seq_len(max_consensus_iter)) {
    n_consensus_iter <- iter

    cc_sparse <- build_sparse_coclassification_cpp(
      memberships, n_genes, edge_list_0, n_cores
    )

    if (iter == 1L) {
      resolution_scan$expected_coclassification <- cc_sparse$expected
    }

    keep <- cc_sparse$excess > 0
    if (!any(keep)) {
      memberships <- initial_memberships
      break
    }

    consensus_el <- edge_list_0[keep, , drop = FALSE] + 1L
    g_consensus <- igraph::make_empty_graph(n = n_genes, directed = FALSE)
    igraph::V(g_consensus)$name <- genes
    g_consensus <- igraph::add_edges(
      g_consensus,
      as.vector(t(consensus_el)),
      weight = cc_sparse$excess[keep]
    )

    new_memberships <- consensus_leiden_sweep(
      g_consensus, resolutions, n_cores,
      seed_root = seed_root, iter = iter
    )

    # Convergence: all K partitions are ~identical (vacuously true for
    # n_res == 1, but that case is caught by the early return above)
    converged <- all(vapply(seq_len(n_res - 1L), function(r) {
      igraph::compare(new_memberships[[r]], new_memberships[[r + 1L]],
        method = "adjusted.rand"
      ) > 0.999
    }, logical(1)))

    # Or the sweep has reproduced its own input. The criterion above asks
    # the K partitions to agree with EACH OTHER, which a stable
    # disagreement between the coarsest and the finest resolution can deny
    # forever: on 6 of 8 Pooideae networks the sweep stopped changing after
    # ~10 iterations and every further iteration rebuilt the same
    # co-classification, the same consensus graph and the same partitions
    # until max_consensus_iter.
    #
    # This is a stopping heuristic, not a proof of a fixed point. The
    # sweep is seeded per iteration (.task_seed(seed_root, 100L + iter,
    # ri) in consensus_leiden_sweep()), so the map applied at iteration
    # i + 1 is not the map applied at iteration i and one reproduction
    # does not entail the next. What is checked is that stopping here
    # agrees with running the full budget: on all 8 Pooideae networks the
    # break fires at 5-12 iterations and returns the same partition, with
    # the same n_modules, as max_consensus_iter = 1000.
    settled <- identical(
      lapply(new_memberships, .partition_id),
      lapply(memberships, .partition_id)
    )

    memberships <- new_memberships
    if (converged || settled) break
  }

  membership <- pick_best_partition(memberships, g)

  # Ensure integer membership with gene names
  membership <- stats::setNames(as.integer(membership), names(membership))
  module_genes <- split(names(membership), membership)

  result <- list(
    modules = membership,
    module_genes = module_genes,
    n_modules = length(module_genes),
    modularity = igraph::modularity(g, membership),
    graph = g,
    method = "leiden_consensus",
    params = list(
      resolutions = resolutions,
      objective_function = objective_function,
      n_resolutions = n_res,
      n_consensus_iterations = n_consensus_iter,
      seed = seed
    ),
    resolution_scan = resolution_scan
  )

  result
}


#' Run Leiden at all resolutions on a consensus graph (internal)
#'
#' Always uses modularity objective regardless of the original objective.
#' Consensus graph weights are between 0 and 1 (excess co-classification); CPM
#' with resolution > 1 would make all edges repulsive, collapsing to
#' singletons. Modularity is the correct objective per Jeub et al. (2018).
#' @noRd
consensus_leiden_sweep <- function(graph, resolutions,
                                   n_cores = 1L, seed_root = 0L, iter = 0L) {
  vertex_names <- igraph::V(graph)$name

  run_one <- function(ri) {
    set.seed(.task_seed(seed_root, 100L + iter, ri))
    res <- resolutions[[ri]]
    comm <- igraph::cluster_leiden(
      graph,
      resolution = res,
      objective_function = "modularity",
      n_iterations = 2L
    )
    mem <- igraph::membership(comm)
    names(mem) <- vertex_names
    mem
  }

  if (.can_fork(n_cores)) {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
    Sys.setenv(OMP_NUM_THREADS = 1L)
    on.exit(
      {
        if (is.na(old_omp)) {
          Sys.unsetenv("OMP_NUM_THREADS")
        } else {
          Sys.setenv(OMP_NUM_THREADS = old_omp)
        }
      },
      add = TRUE
    )
    .check_fork_results(
      parallel::mclapply(seq_along(resolutions), run_one,
        mc.cores = n_cores, mc.preschedule = FALSE
      ),
      resolutions, "resolution"
    )
  } else {
    lapply(seq_along(resolutions), run_one)
  }
}


#' Canonical form of a partition (internal)
#'
#' Relabels module IDs by first appearance, so two membership vectors that
#' describe the same grouping under different labels compare identical.
#' Leiden hands back arbitrary labels; comparing them raw would call a
#' fixed point of the consensus iteration a change.
#' @noRd
.partition_id <- function(mem) {
  mem <- as.integer(mem)
  match(mem, unique(mem))
}


#' Pick partition with best modularity on a reference graph (internal)
#' @noRd
pick_best_partition <- function(memberships, graph) {
  mods <- vapply(memberships, function(mem) {
    igraph::modularity(graph, mem)
  }, numeric(1))
  memberships[[which.max(mods)]]
}
