#' Detect co-expression modules in a network
#'
#' Applies community detection to a thresholded co-expression network using
#' the Leiden algorithm, Infomap, or a Stochastic Block Model (SBM).
#'
#' @section Methods:
#' \describe{
#'   \item{leiden}{Modularity or CPM optimization with guaranteed well-connected
#'     communities (Traag *et al.*, 2019). Resolution parameter controls module
#'     granularity.}
#'   \item{infomap}{Flow-based method that compresses the description of random
#'     walks on the network (Rosvall & Bergstrom, 2008). No resolution parameter;
#'     naturally handles weighted networks.}
#'   \item{sbm}{Gaussian Stochastic Block Model fit by variational EM
#'     (requires the \pkg{sbm} package). Number of blocks is selected
#'     automatically via the Integrated Classification Likelihood (ICL).
#'     Can detect both assortative and non-assortative structure.}
#' }
#'
#' @param net Network object from [compute_network()].
#' @param method Community detection method: `"leiden"` (default), `"infomap"`,
#'   or `"sbm"`.
#' @param resolution Resolution parameter for Leiden (default 1.0). Pass a
#'   numeric vector (e.g., \code{seq(0.5, 2.0, by = 0.5)}) to run at multiple
#'   resolutions and produce a consensus partition via co-classification
#'   (Lancichinetti & Fortunato, 2012). Ignored for other methods.
#' @param objective_function Leiden objective: `"CPM"` (default) or
#'   `"modularity"`. Ignored for other methods.
#' @param n_iterations Number of Leiden iterations (default 2). Ignored for
#'   other methods.
#' @param nb_trials Number of Infomap attempts; best result is kept
#'   (default 10). Ignored for other methods.
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
#'   With an explicit seed the global stream is left exactly where
#'   `set.seed(seed)` put it, so a seeded call does not displace the
#'   caller's stream by however much the clustering backend happened to
#'   consume. That is not the same as restoring the caller's pre-call
#'   state: `set.seed(seed)` has still happened, and anything drawn
#'   afterwards continues from there. With `seed = NULL` the stream
#'   advances -- by the one draw used to pick a root in consensus mode,
#'   and by whatever the backend consumed in single-resolution mode -- so
#'   consecutive unseeded calls still differ.
#' @param consensus_threshold Threshold for consensus mode. \code{NULL}
#'   (default) uses iterative adaptive thresholding per Jeub et al. (2018):
#'   subtracts the per-pair expected co-classification under random assignment
#'   and iterates until the partition converges. A numeric value in (0, 1)
#'   applies a fixed threshold without iteration (single pass).
#' @param n_cores Number of parallel cores (default 1). Used for
#'   \code{mclapply} Leiden sweeps on Unix and OpenMP edge scans in C++.
#'   Uses fork-based parallelism; avoid combining with active CUDA
#'   contexts in the same session.
#' @param max_consensus_iter Maximum number of consensus iterations for
#'   adaptive mode (\code{consensus_threshold = NULL}). Default 10.
#'   Typically converges in 2--5 iterations. Ignored when
#'   \code{consensus_threshold} is numeric.
#' @param test_k1 Logical. Test the null hypothesis K = 1 (no community
#'   structure) via permutation of the spectral norm of the excess
#'   co-classification matrix. Default \code{TRUE}. Only used in adaptive
#'   consensus mode.
#' @param n_perm_k1 Number of permutations for the K = 1 test.
#'   Default 100.
#' @param alpha_k1 Significance level for the K = 1 test.
#'   Default 0.05.
#'
#' @return A list with components:
#'   \describe{
#'     \item{modules}{Named integer vector of module assignments (gene -> module ID)}
#'     \item{module_genes}{Named list: module ID -> character vector of gene names}
#'     \item{n_modules}{Number of modules detected}
#'     \item{modularity}{Modularity score of the partition}
#'     \item{graph}{The igraph graph object used for community detection}
#'     \item{method}{Method used}
#'     \item{params}{List of parameters used}
#'   }
#'   When \code{resolution} is a vector, the output also includes:
#'   \describe{
#'     \item{resolution_scan}{Data frame with columns \code{resolution},
#'       \code{n_modules}, \code{modularity}, \code{ari_next} (Adjusted
#'       Rand Index with the next resolution; NA for the last), and
#'       \code{expected_coclassification} (per-resolution expected scalar).}
#'     \item{k1_test}{When \code{test_k1 = TRUE}, a list with components:
#'       \code{lambda_obs} (observed spectral norm), \code{lambda_null}
#'       (null distribution), \code{p_value}, \code{has_structure}, and
#'       \code{n_perm_completed} (actual permutations run; may be less
#'       than \code{n_perm_k1} due to early stopping).}
#'   }
#'   The \code{params} list includes \code{n_consensus_iterations}
#'   (number of iterations until convergence; 0 for fixed threshold).
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
#' Rosvall, M. & Bergstrom, C. T. (2008). Maps of random walks on complex
#' networks reveal community structure. *PNAS*, 105(4), 1118--1123.
#' \doi{10.1073/pnas.0706851105}
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
#' # Single resolution
#' mods <- detect_modules(net, method = "leiden", resolution = 1.0)
#' table(mods$modules)  # module sizes
#'
#' # Multi-resolution consensus (Jeub et al. 2018)
#' mods_consensus <- detect_modules(net, resolution = c(0.5, 1.0, 2.0),
#'                                  n_cores = 4L)
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
detect_modules <- function(net, ...) UseMethod("detect_modules")

#' @rdname detect_modules
#' @export
detect_modules.default <- function(net,
                           method = c("leiden", "infomap", "sbm"),
                           resolution = 1.0,
                           objective_function = c("CPM", "modularity"),
                           n_iterations = 2L,
                           nb_trials = 10L,
                           seed = NULL,
                           consensus_threshold = NULL,
                           n_cores = 1L,
                           max_consensus_iter = 10L,
                           test_k1 = TRUE,
                           n_perm_k1 = 100L,
                           alpha_k1 = 0.05, ...) {
  method <- match.arg(method)
  objective_function <- match.arg(objective_function)

  # Consensus mode: vector resolution triggers multi-resolution + consensus
  if (length(resolution) > 1L) {
    if (method != "leiden")
      stop("Consensus mode (vector resolution) only supported for method = \"leiden\"")
    return(detect_modules_consensus(
      net, resolution, consensus_threshold,
      objective_function, n_iterations, seed, as.integer(n_cores),
      as.integer(max_consensus_iter), test_k1, as.integer(n_perm_k1),
      alpha_k1))
  }

  n_iterations <- as.integer(n_iterations)
  nb_trials <- as.integer(nb_trials)

  if (!is.list(net) || is.null(net$network)) {
    stop("net must be a network object from compute_network()")
  }

  mat <- .net_check(net, net$threshold)
  thr <- net$threshold
  genes <- rownames(mat)

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

  if (!is.null(seed)) {
    set.seed(seed)
    # Same contract as consensus mode, which pins the stream at the
    # post-seed position (see detect_modules_consensus). Without this the
    # single-resolution path left the stream wherever cluster_leiden() /
    # cluster_infomap() / estimateSimpleSBM() stopped, which is
    # backend- and build-dependent, so a downstream set.seed()-free draw
    # -- summarize_comparison()'s randomized-p pi0, for one -- started
    # from an unpredictable position.
    old_rng <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_rng, envir = globalenv()),
            add = TRUE)
  }

  if (method == "sbm") {
    if (!requireNamespace("sbm", quietly = TRUE)) {
      stop("Package 'sbm' is required for method = \"sbm\". ",
           "Install it with install.packages(\"sbm\")")
    }

    if (.net_is_sparse(net)) {
      warning("SBM requires a dense matrix; densifying")
      adj <- as.matrix(adj)
    }

    fit <- sbm::estimateSimpleSBM(
      adj, model = "gaussian", directed = FALSE,
      estimOptions = list(verbosity = 0L, plot = FALSE)
    )

    membership <- stats::setNames(as.integer(fit$memberships), genes)
    module_genes <- split(names(membership), membership)

    g <- igraph::graph_from_adjacency_matrix(
      adj, mode = "upper", weighted = TRUE, diag = FALSE
    )

    return(list(
      modules = membership,
      module_genes = module_genes,
      n_modules = length(module_genes),
      modularity = igraph::modularity(g, membership),
      graph = g,
      method = method,
      params = list(n_blocks = fit$nbBlocks, ICL = fit$ICL, seed = seed)
    ))
  }

  # Graph-based methods (leiden, infomap)
  g <- igraph::graph_from_adjacency_matrix(
    adj, mode = "upper", weighted = TRUE, diag = FALSE
  )

  if (method == "leiden") {
    comm <- igraph::cluster_leiden(
      g,
      resolution = resolution,
      objective_function = objective_function,
      n_iterations = n_iterations
    )
    params <- list(
      resolution = resolution,
      objective_function = objective_function,
      n_iterations = n_iterations,
      seed = seed
    )
  } else {
    comm <- igraph::cluster_infomap(
      g, e.weights = igraph::E(g)$weight, nb.trials = nb_trials
    )
    params <- list(nb_trials = nb_trials, seed = seed)
  }

  membership <- igraph::membership(comm)
  names(membership) <- igraph::V(g)$name
  module_genes <- split(names(membership), membership)

  list(
    modules = membership,
    module_genes = module_genes,
    n_modules = length(module_genes),
    modularity = igraph::modularity(g, membership),
    graph = g,
    method = method,
    params = params
  )
}


#' Multi-resolution consensus module detection (internal)
#'
#' Runs Leiden at each resolution, builds a co-classification matrix,
#' subtracts per-pair expected co-classification (Jeub et al. 2018),
#' and iterates until the partition converges. For fixed thresholds,
#' performs a single pass without iteration.
#'
#' @noRd
detect_modules_consensus <- function(net, resolutions, consensus_threshold,
                                     objective_function, n_iterations, seed,
                                     n_cores = 1L,
                                     max_consensus_iter = 10L,
                                     test_k1 = TRUE,
                                     n_perm_k1 = 100L,
                                     alpha_k1 = 0.05) {
  # Validate threshold
  if (!is.null(consensus_threshold)) {
    if (!is.numeric(consensus_threshold) || consensus_threshold <= 0 ||
        consensus_threshold >= 1)
      stop("consensus_threshold must be NULL (adaptive) or numeric in (0, 1)")
  }

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

  if (is.null(seed)) {
    seed_root <- sample.int(.Machine$integer.max, 1L)
  } else {
    set.seed(seed)
    seed_root <- as.integer(seed)
  }
  # Per-task set.seed() below runs in the caller's session on the serial path
  # (no fork) but not under mclapply, so restore the ambient stream on exit and
  # leave detect_modules() looking the same at any n_cores. Snapshot taken
  # after the seed handling: with an explicit seed the stream is left exactly
  # where set.seed(seed) put it, and with seed = NULL it is left advanced by
  # the one draw above, so consecutive unseeded calls still differ.
  # Both branches above have seeded, so .Random.seed exists unconditionally
  # here -- unlike coexpressolog_null(), which snapshots before seeding.
  old_rng <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", old_rng, envir = globalenv()), add = TRUE)

  # Build original graph — then free the dense adjacency (~4.6 GB for N=24k)
  g <- igraph::graph_from_adjacency_matrix(
    adj, mode = "upper", weighted = TRUE, diag = FALSE
  )
  rm(adj); gc()

  resolutions <- sort(resolutions)
  n_res <- length(resolutions)

  # If only one resolution after dedup, fall back to single-resolution
  if (n_res == 1L) {
    return(detect_modules(net, method = "leiden", resolution = resolutions,
                          objective_function = objective_function,
                          n_iterations = n_iterations, seed = seed))
  }

  use_mc <- .Platform$OS.type == "unix" && n_cores > 1L

  # ---- Initial Leiden sweep on original graph ----
  run_initial <- function(ri) {
    set.seed(.task_seed(seed_root, 1L, ri))
    res <- resolutions[[ri]]
    comm <- igraph::cluster_leiden(
      g,
      resolution = res,
      objective_function = objective_function,
      n_iterations = as.integer(n_iterations)
    )
    mem <- igraph::membership(comm)
    names(mem) <- igraph::V(g)$name
    list(mem = mem,
         n_mod = length(unique(mem)),
         quality = igraph::modularity(g, mem))
  }

  if (use_mc) {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
    Sys.setenv(OMP_NUM_THREADS = 1L)
    on.exit({
      if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS")
      else Sys.setenv(OMP_NUM_THREADS = old_omp)
    }, add = TRUE)
    results <- parallel::mclapply(seq_len(n_res), run_initial,
                                  mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(n_res), run_initial)
  }

  errs <- which(vapply(results, inherits, logical(1), "try-error"))
  if (length(errs)) {
    e <- results[[errs[1L]]]
    stop(attr(e, "condition") %||% as.character(e))
  }
  failed <- vapply(results, is.null, logical(1))
  if (any(failed)) {
    stop("Parallel workers returned NULL at resolutions: ",
         paste(resolutions[failed], collapse = ", "))
  }

  memberships <- lapply(results, `[[`, "mem")
  scan_n_modules <- vapply(results, `[[`, integer(1), "n_mod")
  scan_modularity <- vapply(results, `[[`, numeric(1), "quality")

  # ARI between consecutive resolutions
  scan_ari_next <- rep(NA_real_, n_res)
  for (r in seq_len(n_res - 1L)) {
    scan_ari_next[r] <- igraph::compare(
      memberships[[r]], memberships[[r + 1L]], method = "adjusted.rand"
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
  k1_result <- NULL

  if (!is.null(consensus_threshold)) {
    # Fixed threshold: single pass, sparse co-classification
    cc_sparse <- build_sparse_coclassification_cpp(
      memberships, n_genes, edge_list_0, n_cores
    )
    resolution_scan$expected_coclassification <- cc_sparse$expected

    keep <- cc_sparse$coclassification >= consensus_threshold
    if (!any(keep)) {
      best_r <- which.max(scan_modularity)
      membership <- initial_memberships[[best_r]]
    } else {
      consensus_el <- edge_list_0[keep, , drop = FALSE] + 1L
      g_consensus <- igraph::make_empty_graph(n = n_genes, directed = FALSE)
      igraph::V(g_consensus)$name <- genes
      g_consensus <- igraph::add_edges(
        g_consensus,
        as.vector(t(consensus_el)),
        weight = cc_sparse$coclassification[keep]
      )

      consensus_mems <- consensus_leiden_sweep(
        g_consensus, resolutions, n_iterations, n_cores,
        seed_root = seed_root, iter = 0L
      )
      membership <- pick_best_partition(consensus_mems, g)
    }
  } else {
    # K = 1 test
    if (test_k1) {
      k1_result <- test_community_structure(
        g, genes, resolutions, objective_function, n_iterations,
        memberships, edge_list_0, n_perm_k1, n_cores, alpha_k1,
        seed_root = seed_root
      )
      if (!k1_result$has_structure) {
        membership <- stats::setNames(rep(1L, n_genes), genes)
        module_genes <- list(`1` = genes)
        return(list(
          modules = membership,
          module_genes = module_genes,
          n_modules = 1L,
          modularity = igraph::modularity(g, membership),
          graph = g,
          method = "leiden_consensus",
          params = list(
            resolutions = resolutions,
            consensus_threshold = consensus_threshold,
            objective_function = objective_function,
            n_resolutions = n_res,
            n_consensus_iterations = 0L,
            seed = seed
          ),
          resolution_scan = resolution_scan,
          k1_test = k1_result
        ))
      }
    }

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
        g_consensus, resolutions, n_iterations, n_cores,
        seed_root = seed_root, iter = iter
      )

      # Convergence: all K partitions are ~identical (vacuously true for
      # n_res == 1, but that case is caught by the early return above)
      converged <- all(vapply(seq_len(n_res - 1L), function(r) {
        igraph::compare(new_memberships[[r]], new_memberships[[r + 1L]],
                        method = "adjusted.rand") > 0.999
      }, logical(1)))

      memberships <- new_memberships
      if (converged) break
    }

    membership <- pick_best_partition(memberships, g)
  }

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
      consensus_threshold = consensus_threshold,
      objective_function = objective_function,
      n_resolutions = n_res,
      n_consensus_iterations = n_consensus_iter,
      seed = seed
    ),
    resolution_scan = resolution_scan
  )

  if (test_k1 && !is.null(k1_result)) {
    result$k1_test <- k1_result
  }

  result
}


#' Run Leiden at all resolutions on a consensus graph (internal)
#'
#' Always uses modularity objective regardless of the original objective.
#' Consensus graph weights are between 0 and 1 (excess co-classification); CPM
#' with resolution > 1 would make all edges repulsive, collapsing to
#' singletons. Modularity is the correct objective per Jeub et al. (2018).
#' @noRd
consensus_leiden_sweep <- function(graph, resolutions, n_iterations,
                                   n_cores = 1L, seed_root = 0L, iter = 0L) {
  vertex_names <- igraph::V(graph)$name

  run_one <- function(ri) {
    set.seed(.task_seed(seed_root, 100L + iter, ri))
    res <- resolutions[[ri]]
    comm <- igraph::cluster_leiden(
      graph,
      resolution = res,
      objective_function = "modularity",
      n_iterations = as.integer(n_iterations)
    )
    mem <- igraph::membership(comm)
    names(mem) <- vertex_names
    mem
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
    Sys.setenv(OMP_NUM_THREADS = 1L)
    on.exit({
      if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS")
      else Sys.setenv(OMP_NUM_THREADS = old_omp)
    }, add = TRUE)
    results <- parallel::mclapply(seq_along(resolutions), run_one,
                                  mc.cores = n_cores)
    errs <- which(vapply(results, inherits, logical(1), "try-error"))
    if (length(errs)) {
      e <- results[[errs[1L]]]
      stop(attr(e, "condition") %||% as.character(e))
    }
    results
  } else {
    lapply(seq_along(resolutions), run_one)
  }
}


#' Pick partition with best modularity on a reference graph (internal)
#' @noRd
pick_best_partition <- function(memberships, graph) {
  mods <- vapply(memberships, function(mem) {
    igraph::modularity(graph, mem)
  }, numeric(1))
  memberships[[which.max(mods)]]
}


#' Test for community structure (K = 1 null) via spectral norm permutation
#'
#' Compares the leading eigenvalue of the sparse excess co-classification
#' matrix against a null distribution from degree-preserving rewiring.
#' Uses batch-based early stopping: once enough permutations have been
#' completed without any exceedance (\code{ceil(1/alpha)} permutations
#' with lambda_null < lambda_obs), the test concludes that structure is
#' present without running all \code{n_perm} permutations.
#'
#' @section Performance:
#' Three optimizations reduce runtime vs naive implementation:
#' \enumerate{
#'   \item Rewiring uses 5 * |E| swap attempts (sufficient for mixing;
#'     Greenhill, 2015).
#'   \item Null Leiden sweeps use \code{n_iterations = 1} (partitions need
#'     not be optimal for the null distribution).
#'   \item Batch early stopping: permutations run in batches of
#'     \code{n_cores}. After each batch, if \code{ceil(1/alpha)}
#'     permutations have completed with zero exceedances, the test stops
#'     early (p < alpha is guaranteed). Similarly, if exceedances
#'     accumulate such that p > alpha is certain, the test stops.
#' }
#'
#' @noRd
test_community_structure <- function(g, genes, resolutions, objective_function,
                                      n_iterations, memberships_obs,
                                      edge_list_0, n_perm = 100L,
                                      n_cores = 1L, alpha = 0.05,
                                      seed_root = 0L) {
  n_genes <- length(genes)
  n_edges <- igraph::ecount(g)

  lambda_obs <- sparse_excess_spectral_norm_cpp(memberships_obs, n_genes,
                                                 edge_list_0, n_cores)

  run_one_perm <- function(b) {
    set.seed(.task_seed(seed_root, 2L, b))
    g_perm <- igraph::rewire(g, igraph::keeping_degseq(
      niter = 5L * n_edges))
    igraph::E(g_perm)$weight <- igraph::E(g)$weight[
      sample.int(n_edges)]

    mems_perm <- lapply(resolutions, function(res) {
      comm <- igraph::cluster_leiden(
        g_perm, resolution = res,
        objective_function = objective_function,
        n_iterations = 1L
      )
      mem <- igraph::membership(comm)
      names(mem) <- igraph::V(g_perm)$name
      mem
    })

    el_perm <- igraph::as_edgelist(g_perm, names = FALSE) - 1L
    storage.mode(el_perm) <- "integer"
    sparse_excess_spectral_norm_cpp(mems_perm, n_genes, el_perm)
  }

  # Minimum permutations before early stopping can trigger (ceil(1/alpha))
  min_for_sig <- as.integer(ceiling(1 / alpha))

  use_mc <- .Platform$OS.type == "unix" && n_cores > 1L
  # Batch on the significance grid, not on the core count: the early-stop rule
  # must be evaluated at the same points regardless of the machine. The batch
  # is still spread over mc.cores below, so on typical hardware concurrency is
  # unaffected -- but a batch of ceiling(1 / alpha) tasks cannot occupy more
  # than that many workers, so a large alpha_k1 on a many-core node leaves
  # cores idle during this test. Determinism is worth that.
  batch_size <- max(1L, min_for_sig)

  if (use_mc) {
    old_omp <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
    Sys.setenv(OMP_NUM_THREADS = 1L)
    on.exit({
      if (is.na(old_omp)) Sys.unsetenv("OMP_NUM_THREADS")
      else Sys.setenv(OMP_NUM_THREADS = old_omp)
    }, add = TRUE)
  }

  lambda_null <- numeric(n_perm)
  n_exceed <- 0L
  n_done <- 0L

  for (batch_start in seq(1L, n_perm, by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1L, n_perm)
    batch_idx <- seq.int(batch_start, batch_end)

    if (use_mc) {
      batch_vals <- unlist(parallel::mclapply(
        batch_idx, run_one_perm, mc.cores = n_cores
      ))
    } else {
      batch_vals <- vapply(batch_idx, run_one_perm, numeric(1))
    }

    lambda_null[batch_idx] <- batch_vals
    n_exceed <- n_exceed + sum(batch_vals >= lambda_obs)
    n_done <- batch_end

    # Early stop: clear structure — p = 1/(n_done+1) < alpha
    if (n_exceed == 0L && n_done >= min_for_sig) break
    # Early stop: no structure — p = (n_exceed+1)/(n_done+1) > alpha
    if (n_exceed > 0L && n_done >= min_for_sig &&
        (n_exceed + 1L) / (n_done + 1L) > alpha) break
  }

  lambda_null <- lambda_null[seq_len(n_done)]
  p_value <- (1 + n_exceed) / (1 + n_done)

  list(
    lambda_obs = lambda_obs,
    lambda_null = lambda_null,
    p_value = p_value,
    has_structure = p_value < alpha,
    n_perm_completed = n_done
  )
}


#' Identify hub genes within co-expression modules
#'
#' Computes within-module centrality for each gene and flags the top-ranked
#' genes as hubs.  Optionally maps genes to ortholog groups (HOGs) for
#' downstream conservation analysis with [classify_hub_conservation()].
#'
#' @section Centrality measures:
#' Centrality is computed on the **within-module subgraph** (edges between
#' genes in the same module only):
#' \describe{
#'   \item{degree}{Weighted degree (`igraph::strength`): sum of edge weights
#'     to other genes in the same module.}
#'   \item{betweenness}{Shortest-path betweenness using inverse edge weights
#'     as distances.  Identifies genes that bridge sub-clusters within a module.}
#'   \item{eigenvector}{Eigenvector centrality (`igraph::eigen_centrality`):
#'     high for genes connected to other high-centrality genes.}
#' }
#'
#' @section Tie-breaking cascade:
#' When genes share the same primary centrality score, hub selection uses a
#' biologically informed cascade (all available tiers evaluated):
#' \enumerate{
#'   \item Primary centrality (user-selected measure)
#'   \item Global weighted degree across the full network
#'   \item Alternative within-module centrality (betweenness if primary is
#'     degree; degree otherwise)
#'   \item Mean within-module edge weight (strength / degree)
#'   \item Per-gene conservation effect size (requires `comparison`)
#'   \item Per-HOG minimum q-value (requires `comparison`; lower = better)
#' }
#'
#' @param modules Output of [detect_modules()].
#' @param net Output of [compute_network()].
#' @param orthologs Optional data frame from [parse_orthologs()] with columns
#'   `Species1`, `Species2`, `hog`.  The function auto-detects which column
#'   matches the gene names in `modules`.  If `NULL`, the `hog` column in the
#'   result is all `NA`.
#' @param comparison Optional data frame: the `$results` element from
#'   [summarize_comparison()].  When provided, enables conservation-informed
#'   tie-breaking (tiers 5--6).  Must contain columns `Species1`, `Species2`,
#'   `hog`, `Species1.effect.size`, `Species2.effect.size`, plus at least one
#'   pair of q-value columns (`Species1.q.val.con`/`Species2.q.val.con` or
#'   the `.div` variants).
#' @param centrality Centrality measure: `"degree"` (default), `"betweenness"`,
#'   or `"eigenvector"`.
#' @param top_n Integer: flag the top N genes per module as hubs.  If `NULL`
#'   (default), uses `top_fraction` instead.
#' @param top_fraction Numeric in (0, 1): fraction of genes per module to flag
#'   as hubs (default 0.1).  Ignored when `top_n` is non-NULL.
#' @param min_module_size Integer: modules with fewer genes get
#'   `is_hub = FALSE` for all genes (default 3).
#'
#' @return A data frame with one row per gene, ordered by module then rank:
#'   \describe{
#'     \item{gene}{Gene identifier}
#'     \item{module}{Module ID (integer)}
#'     \item{degree}{Within-module weighted degree (`igraph::strength`)}
#'     \item{betweenness}{Within-module betweenness centrality}
#'     \item{eigenvector}{Within-module eigenvector centrality}
#'     \item{mean_edge_weight}{Mean weight of edges to other module members}
#'     \item{global_degree}{Weighted degree in the full (thresholded) network}
#'     \item{rank}{Rank within module by primary centrality
#'       (1 = highest; ties use `"min"`)}
#'     \item{is_hub}{`TRUE` if the gene is in the top slice after the
#'       6-tier tie-breaking cascade}
#'     \item{hog}{HOG identifier (`NA` if `orthologs` not provided or gene
#'       not in the ortholog table)}
#'   }
#'
#' @examples
#' \dontrun{
#' hubs <- identify_module_hubs(modules, net, orthologs,
#'                              comparison = summary$results)
#' hubs[hubs$is_hub, ]
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
identify_module_hubs <- function(modules, ...) {
  UseMethod("identify_module_hubs")
}

#' @rdname identify_module_hubs
#' @export
identify_module_hubs.default <- function(modules, net, orthologs = NULL,
                                 comparison = NULL,
                                 centrality = c("degree", "betweenness",
                                                "eigenvector"),
                                 top_n = NULL,
                                 top_fraction = 0.1,
                                 min_module_size = 3L, ...) {
  centrality <- match.arg(centrality)

  if (!is.list(modules) || is.null(modules$module_genes) ||
        is.null(modules$graph) || is.null(modules$modules)) {
    stop("modules must be output from detect_modules()")
  }
  if (!is.list(net) || is.null(net$network)) {
    stop("net must be output from compute_network()")
  }
  .net_check(net, net$threshold)
  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    if (top_n < 1L) stop("top_n must be >= 1")
  } else {
    if (top_fraction <= 0 || top_fraction >= 1) {
      stop("top_fraction must be in (0, 1)")
    }
  }
  min_module_size <- as.integer(min_module_size)

  g <- modules$graph

  # Pre-compute global weighted degree (tie-breaker tier 2)
  global_str <- igraph::strength(g)

  # Pre-compute conservation lookups if comparison provided (tiers 5-6)
  gene_conserv <- NULL
  hog_min_q <- NULL
  if (!is.null(comparison)) {
    if (!all(c("Species1", "Species2", "hog",
               "Species1.effect.size", "Species2.effect.size") %in%
             names(comparison))) {
      stop("comparison must be $results from summarize_comparison()")
    }
    # Auto-detect which column has our genes
    all_genes <- names(modules$modules)
    in_sp1 <- sum(all_genes %in% comparison$Species1)
    in_sp2 <- sum(all_genes %in% comparison$Species2)
    comp_col <- if (in_sp1 >= in_sp2) "Species1" else "Species2"

    # Per-row geometric mean of effect sizes
    geo_eff <- sqrt(comparison$Species1.effect.size *
                    comparison$Species2.effect.size)

    # Per-gene mean conservation effect (higher = more conserved)
    comp_genes <- comparison[[comp_col]]
    gene_conserv <- vapply(
      split(geo_eff, comp_genes), mean, numeric(1), na.rm = TRUE
    )

    # Per-HOG minimum q-value (lower = more conserved)
    q1_col <- if ("Species1.q.val.con" %in% names(comparison)) {
      "Species1.q.val.con"
    } else if ("Species1.q.val.div" %in% names(comparison)) {
      "Species1.q.val.div"
    } else {
      NULL
    }
    q2_col <- if ("Species2.q.val.con" %in% names(comparison)) {
      "Species2.q.val.con"
    } else if ("Species2.q.val.div" %in% names(comparison)) {
      "Species2.q.val.div"
    } else {
      NULL
    }
    if (!is.null(q1_col) && !is.null(q2_col)) {
      pair_q <- pmin(comparison[[q1_col]], comparison[[q2_col]], na.rm = TRUE)
      hog_min_q <- vapply(
        split(pair_q, comparison$hog), min, numeric(1), na.rm = TRUE
      )
    }
  }

  # HOG mapping (needed for tier-6 tie-breaking and output)
  hog_lookup <- NULL
  if (!is.null(orthologs)) {
    if (!all(c("Species1", "Species2", "hog") %in% names(orthologs))) {
      stop("orthologs must have columns: Species1, Species2, hog")
    }
    all_genes <- names(modules$modules)
    in_sp1 <- sum(all_genes %in% orthologs$Species1)
    in_sp2 <- sum(all_genes %in% orthologs$Species2)
    gene_col <- if (in_sp1 >= in_sp2) "Species1" else "Species2"

    gene_hog <- unique(orthologs[, c(gene_col, "hog"), drop = FALSE])
    gene_hog <- gene_hog[!duplicated(gene_hog[[gene_col]]), , drop = FALSE]
    hog_lookup <- stats::setNames(
      as.character(gene_hog$hog), gene_hog[[gene_col]]
    )
  }

  rows <- vector("list", length(modules$module_genes))

  for (i in seq_along(modules$module_genes)) {
    mod_id <- names(modules$module_genes)[i]
    genes <- modules$module_genes[[i]]
    n_genes <- length(genes)

    if (n_genes < min_module_size) {
      rows[[i]] <- data.frame(
        gene = genes, module = as.integer(mod_id),
        degree = NA_real_, betweenness = NA_real_, eigenvector = NA_real_,
        mean_edge_weight = NA_real_,
        global_degree = global_str[genes],
        rank = NA_integer_, is_hub = FALSE,
        stringsAsFactors = FALSE
      )
      next
    }

    sub <- igraph::induced_subgraph(g, genes)
    w <- igraph::E(sub)$weight
    inv_w <- if (!is.null(w)) 1 / w else NULL

    # Compute all three centrality measures once
    sub_str <- igraph::strength(sub)
    sub_btw <- igraph::betweenness(sub, weights = inv_w)
    sub_eig <- tryCatch(
      igraph::eigen_centrality(sub, weights = w)$vector,
      error = function(e) {
        warning("eigen_centrality failed for module ", mod_id, ": ",
                conditionMessage(e), "; using zero fallback")
        stats::setNames(rep(0, length(genes)), genes)
      }
    )

    # Primary centrality for ranking/tie-breaking (tier 1)
    cent_vals <- switch(centrality,
      degree = sub_str, betweenness = sub_btw, eigenvector = sub_eig
    )
    # Alternative centrality (tier 3): betweenness if primary is degree,
    # degree otherwise — the most complementary pair
    alt_cent <- if (centrality == "degree") sub_btw else sub_str

    # Mean within-module edge weight (tier 4): strength / degree
    sub_deg <- igraph::degree(sub)
    mean_ew <- ifelse(sub_deg > 0, sub_str / sub_deg, 0)

    rnk <- rank(-cent_vals, ties.method = "min")

    rows[[i]] <- data.frame(
      gene = names(cent_vals), module = as.integer(mod_id),
      degree = as.numeric(sub_str),
      betweenness = as.numeric(sub_btw),
      eigenvector = as.numeric(sub_eig),
      mean_edge_weight = as.numeric(mean_ew),
      global_degree = as.numeric(global_str[genes]),
      rank = as.integer(rnk),
      is_hub = FALSE,  # filled below
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  # HOG mapping (vectorized, once for all genes)
  result$hog <- NA_character_
  if (!is.null(hog_lookup)) {
    matched <- match(result$gene, names(hog_lookup))
    result$hog[!is.na(matched)] <- hog_lookup[matched[!is.na(matched)]]
  }

  # Conservation lookups (vectorized, once for all genes — tiers 5-6)
  result$conserv_eff <- 0
  result$hog_q <- 1
  if (!is.null(gene_conserv)) {
    matched <- match(result$gene, names(gene_conserv))
    result$conserv_eff[!is.na(matched)] <- gene_conserv[matched[!is.na(matched)]]
  }
  if (!is.null(hog_min_q) && !is.null(hog_lookup)) {
    matched <- match(result$hog, names(hog_min_q))
    result$hog_q[!is.na(matched)] <- hog_min_q[matched[!is.na(matched)]]
  }

  # Primary and alternative centrality column names for tie-breaking
  primary_col <- centrality  # "degree", "betweenness", or "eigenvector"
  alt_col <- if (centrality == "degree") "betweenness" else "degree"

  # Hub selection per module: tie-breaking cascade across all 6 tiers
  # Tiers 1-5 descending (higher = better), tier 6 ascending (lower = better)
  mod_ids <- unique(result$module[!is.na(result$degree)])
  for (m in mod_ids) {
    idx <- which(result$module == m & !is.na(result$degree))
    n_mod <- length(idx)
    hub_cutoff <- if (!is.null(top_n)) {
      min(top_n, n_mod)
    } else {
      max(1L, ceiling(top_fraction * n_mod))
    }
    ord <- order(-result[[primary_col]][idx], -result$global_degree[idx],
                 -result[[alt_col]][idx], -result$mean_edge_weight[idx],
                 -result$conserv_eff[idx], result$hog_q[idx])
    result$is_hub[idx[ord[seq_len(hub_cutoff)]]] <- TRUE
  }

  # Drop internal tie-breaking columns (conservation lookups)
  result$conserv_eff <- NULL
  result$hog_q <- NULL

  attr(result, "primary_centrality") <- centrality
  result
}


#' Classify hub gene conservation across species and traits
#'
#' Given per-species hub identification results (from
#' [identify_module_hubs()]), maps hub genes to HOGs and classifies each HOG
#' by its hub conservation pattern relative to a discrete trait (e.g.
#' annual / perennial).
#'
#' @section Classification waterfall:
#' For each HOG that appears in at least one species:
#' \describe{
#'   \item{conserved_hub}{Hub in multiple trait groups **and** the hub modules
#'     correspond across traits (checked via `module_comparisons`).}
#'   \item{rewired_hub}{Hub in multiple trait groups but in
#'     **non-corresponding** modules -- the gene kept its centrality but
#'     changed regulatory context.}
#'   \item{multi_trait_hub}{Hub in multiple trait groups; module correspondence
#'     unknown (`module_comparisons` not provided).}
#'   \item{\emph{trait}_specific_hub}{Hub in exactly one trait group (e.g.
#'     `"annual_specific_hub"`).}
#'   \item{sporadic_hub}{Hub in some species but does not reach
#'     `min_trait_fraction` in any trait group.}
#'   \item{non_hub}{Present in modules but not a hub in any species.}
#' }
#'
#' @param hub_results Named list keyed by species name.  Each element is the
#'   data frame output of [identify_module_hubs()] (with `orthologs`
#'   provided so the `hog` column is populated).
#' @param species_trait Named character or factor vector mapping species to
#'   trait groups, e.g. `c(SP_A = "annual", SP_B = "annual",
#'   SP_C = "perennial", SP_D = "perennial")`.
#' @param module_comparisons Optional named list of
#'   [module_correspondence()] outputs keyed by alphabetically sorted species
#'   pair (e.g. `"SP_A.SP_C"`). Required for the conserved_hub vs rewired_hub
#'   distinction. Two things the caller must now satisfy themselves: each
#'   element must be built with the alphabetically first species as
#'   `modules_ref`, and [module_correspondence()] needs a map from
#'   [resolve_ortholog_map()], so the networks are required for the gene
#'   universes. Pass `sp_ref` / `sp_test` to [module_correspondence()] and
#'   that orientation is checked here instead of taken on trust. A module
#'   pair absent from the table counts as not corresponding.
#' @param alpha Significance threshold for module correspondence (default 0.05).
#' @param jaccard_threshold Jaccard threshold for module correspondence
#'   (default 0.1). [module_correspondence()] computes this over the
#'   one-to-one paralog-resolved projection, whereas the retired gene-overlap
#'   engine used the paralog-expanded mappable set, so values now run
#'   systematically higher and the unchanged default is slightly more
#'   permissive.
#' @param min_trait_fraction Minimum fraction of species (within a trait group)
#'   where the HOG must be a hub for the group to count (default 0.5).
#' @param correspondence_threshold Fraction of cross-trait hub pairs that must
#'   have corresponding modules for the HOG to be classified as
#'   `conserved_hub` rather than `rewired_hub` (default 0.5).
#'
#' @return A data frame with one row per HOG:
#'   \describe{
#'     \item{hog}{HOG identifier}
#'     \item{classification}{Conservation category (see Classification
#'       waterfall)}
#'     \item{n_species_hub}{Number of species where the HOG is a hub}
#'     \item{n_species_present}{Number of species where the HOG has genes}
#'     \item{hub_trait_groups}{Comma-separated trait groups where it qualifies
#'       as hub (`NA` for non_hub)}
#'     \item{n_corresponding}{Cross-trait hub pairs with corresponding modules
#'       (`NA` without `module_comparisons`)}
#'     \item{n_cross_pairs}{Total cross-trait hub pairs checked (`NA` without
#'       `module_comparisons`)}
#'     \item{max_centrality}{Highest centrality score across species}
#'     \item{best_hub_species}{Species with highest centrality}
#'   }
#'
#' @examples
#' \dontrun{
#' hub_list <- list(
#'   SP_A = identify_module_hubs(mods_A, net_A, ortho_A),
#'   SP_B = identify_module_hubs(mods_B, net_B, ortho_B)
#' )
#' trait <- c(SP_A = "annual", SP_B = "perennial")
#' classify_hub_conservation(hub_list, trait)
#'
#' # With module correspondence, for the conserved_hub / rewired_hub split.
#' # The list key must be the alphabetically sorted species pair.
#' map <- resolve_ortholog_map(
#'   ortho_AB, rownames(net_A$network), rownames(net_B$network)
#' )
#' corr <- list(SP_A.SP_B = module_correspondence(
#'   mods_A, mods_B, map, sp_ref = "SP_A", sp_test = "SP_B"
#' ))
#' classify_hub_conservation(hub_list, trait, module_comparisons = corr)
#' }
#'
#' @param ... Additional arguments passed to the default method.
#' @export
classify_hub_conservation <- function(hub_results, ...) {
  UseMethod("classify_hub_conservation")
}

#' @rdname classify_hub_conservation
#' @export
classify_hub_conservation.default <- function(hub_results, species_trait,
                                      module_comparisons = NULL,
                                      alpha = 0.05,
                                      jaccard_threshold = 0.1,
                                      min_trait_fraction = 0.5,
                                      correspondence_threshold = 0.5, ...) {
  # --- Validation ---
  if (!is.list(hub_results) || is.null(names(hub_results))) {
    stop("hub_results must be a named list keyed by species")
  }
  if (!is.character(species_trait) && !is.factor(species_trait)) {
    stop("species_trait must be a named character or factor vector")
  }
  if (is.null(names(species_trait))) {
    stop("species_trait must be a named vector")
  }
  missing_sp <- setdiff(names(hub_results), names(species_trait))
  if (length(missing_sp) > 0) {
    stop("species_trait missing entries for: ",
         paste(missing_sp, collapse = ", "))
  }
  req_cols <- c("gene", "module", "is_hub", "hog", "degree")
  for (sp in names(hub_results)) {
    if (!is.data.frame(hub_results[[sp]]) ||
          !all(req_cols %in% names(hub_results[[sp]]))) {
      stop("hub_results[['", sp,
           "']] must be output from identify_module_hubs() with orthologs")
    }
  }

  trait_char <- as.character(species_trait[names(hub_results)])
  names(trait_char) <- names(hub_results)
  trait_levels <- unique(trait_char)
  species_by_trait <- split(names(trait_char), trait_char)

  # Determine which centrality column to use for max_centrality / hub_module.
  # Reads the attribute set by identify_module_hubs(); falls back to "degree".
  primary_col <- unique(vapply(hub_results, function(hr) {
    pc <- attr(hr, "primary_centrality")
    if (is.null(pc)) "degree" else pc
  }, character(1)))
  if (length(primary_col) != 1L) primary_col <- "degree"

  # --- Build HOG-level summary: stack all results, aggregate per (hog, sp) ---
  tagged <- lapply(names(hub_results), function(sp) {
    hr <- hub_results[[sp]]
    hr <- hr[!is.na(hr$hog), , drop = FALSE]
    if (nrow(hr) == 0L) return(NULL)
    hr$species <- sp
    hr
  })
  stacked <- do.call(rbind, tagged)

  # Empty result template
  empty <- data.frame(
    hog = character(0), classification = character(0),
    n_species_hub = integer(0), n_species_present = integer(0),
    hub_trait_groups = character(0),
    n_corresponding = integer(0), n_cross_pairs = integer(0),
    max_centrality = numeric(0), best_hub_species = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(stacked) || nrow(stacked) == 0L) return(empty)

  # One row per (hog, species): is_hub (OR), hub_module, max_centrality
  hog_df <- do.call(rbind, lapply(
    split(stacked, paste(stacked$hog, stacked$species, sep = "\x01")),
    function(df) {
      any_hub <- any(df$is_hub)
      hub_mod <- if (any_hub) {
        hub_rows <- df[df$is_hub, , drop = FALSE]
        hub_rows$module[which.max(hub_rows[[primary_col]])]
      } else {
        NA_integer_
      }
      cent <- df[[primary_col]][!is.na(df[[primary_col]])]
      data.frame(
        hog = df$hog[1], species = df$species[1],
        is_hub = any_hub, hub_module = hub_mod,
        max_centrality = if (length(cent) == 0L) NA_real_ else max(cent),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(hog_df) <- NULL

  # --- Pre-compute (hog x trait) hub fraction matrix ---
  hog_df$trait <- trait_char[hog_df$species]
  hub_frac <- tapply(hog_df$is_hub, list(hog_df$hog, hog_df$trait), mean)
  hub_frac[is.na(hub_frac)] <- 0
  is_hub_group <- hub_frac >= min_trait_fraction  # logical matrix

  # --- Pre-compute per-HOG aggregates ---
  hog_n_present <- tapply(hog_df$species, hog_df$hog, length)
  hog_n_hub <- tapply(hog_df$is_hub, hog_df$hog, sum)
  hog_max_cent <- tapply(hog_df$max_centrality, hog_df$hog, function(x) {
    cx <- x[!is.na(x)]
    if (length(cx) == 0L) NA_real_ else max(cx)
  })
  hog_best_sp <- tapply(
    seq_len(nrow(hog_df)), hog_df$hog,
    function(idx) {
      sub <- hog_df[idx, , drop = FALSE]
      cx <- sub$max_centrality
      if (all(is.na(cx))) sub$species[1] else sub$species[which.max(cx)]
    }
  )

  # --- Pre-build module correspondence lookup per species pair ---
  # Validate first: an element without a $pairs data frame (what
  # preservation_paired()$raw gives you) would otherwise leave is_match as
  # logical(0) and report every HOG as NA, indistinguishable from having
  # supplied no comparison at all.
  if (!is.null(module_comparisons)) {
    # Keys first. An unnamed list makes the loops below iterate over NULL and
    # a wrongly-ordered key never matches the sorted lookup, and both leave
    # every HOG at NA -- indistinguishable from supplying no comparison, which
    # is the failure this guard exists to prevent.
    nm <- names(module_comparisons)
    if (is.null(nm) || !all(nzchar(nm))) {
      stop("module_comparisons must be a named list keyed by ",
           "alphabetically sorted species pair (e.g. \"SP_A.SP_C\")")
    }
    known_sp <- names(species_trait)
    valid_keys <- if (length(known_sp) >= 2L) {
      apply(utils::combn(sort(known_sp), 2L), 2L, paste, collapse = ".")
    } else {
      character(0)
    }
    bad_keys <- setdiff(nm, valid_keys)
    if (length(bad_keys) > 0L) {
      stop("module_comparisons keys must be alphabetically sorted species ",
           "pairs drawn from species_trait; unusable: ",
           paste(bad_keys, collapse = ", "))
    }
    # Orientation, when the producer recorded it. A transposed call --
    # module_correspondence(mods_B, mods_A, ...) filed under "A.B" -- passes
    # the name and shape checks and then matches lookups with module_sp1 and
    # module_sp2 swapped, giving wrong verdicts rather than a detectable NA.
    for (k in nm) {
      ref <- module_comparisons[[k]]$sp_ref
      if (is.null(ref)) next
      tst <- module_comparisons[[k]]$sp_test
      # Rebuild the key from the recorded labels rather than splitting it.
      # Splitting on "." mangles species names that contain one, and
      # comparing the pair as a whole also catches a sp_test naming a third
      # species, which a first-element check would pass.
      rebuilt <- if (is.null(tst)) NULL else paste(c(ref, tst), collapse = ".")
      ok <- if (is.null(rebuilt)) {
        startsWith(k, paste0(ref, "."))
      } else {
        identical(k, rebuilt)
      }
      if (!ok) {
        stop("module_comparisons[[\"", k, "\"]] was built with sp_ref = \"",
             ref, "\"", if (!is.null(tst)) paste0(", sp_test = \"", tst, "\""),
             "; module_sp1 must belong to the first species of the key, so ",
             "the arguments or the key are wrong")
      }
    }

    req_corr <- c("module_sp1", "module_sp2", "jaccard", "q.value")
    for (k in nm) {
      pk <- module_comparisons[[k]]$pairs
      if (!is.data.frame(pk) || !all(req_corr %in% names(pk))) {
        stop("module_comparisons[[\"", k, "\"]] must be a ",
             "module_correspondence() result: a list with a `pairs` data ",
             "frame carrying ", paste(req_corr, collapse = ", "))
      }
    }
  }

  corresp_lookup <- list()  # keyed by "SP_A.SP_C", values = named logical
  if (!is.null(module_comparisons)) {
    for (pair_key in names(module_comparisons)) {
      pairs <- module_comparisons[[pair_key]]$pairs
      is_match <- pairs$q.value < alpha & pairs$jaccard >= jaccard_threshold
      keys <- paste(pairs$module_sp1, pairs$module_sp2, sep = "\x01")
      corresp_lookup[[pair_key]] <- stats::setNames(is_match, keys)
    }
  }

  # O(1) module correspondence check
  check_correspondence <- function(sp_a, mod_a, sp_b, mod_b) {
    if (length(corresp_lookup) == 0L) return(NA)
    pair_key <- paste(sort(c(sp_a, sp_b)), collapse = ".")
    lkp <- corresp_lookup[[pair_key]]
    if (is.null(lkp)) return(NA)
    sorted <- sort(c(sp_a, sp_b))
    mod_key <- if (sp_a == sorted[1]) {
      paste(mod_a, mod_b, sep = "\x01")
    } else {
      paste(mod_b, mod_a, sep = "\x01")
    }
    val <- lkp[mod_key]
    if (is.na(val)) FALSE else val
  }

  # --- Classify each HOG ---
  hog_groups <- split(hog_df, hog_df$hog)
  all_hogs <- names(hog_groups)

  out_rows <- lapply(all_hogs, function(h) {
    h_df <- hog_groups[[h]]
    n_present <- hog_n_present[[h]]
    n_hub <- hog_n_hub[[h]]
    max_cent <- hog_max_cent[[h]]
    best_sp <- hog_best_sp[[h]]

    # Trait-group hub status from pre-computed matrix (fix #6)
    hub_group_names <- colnames(is_hub_group)[is_hub_group[h, ]]

    n_corresponding <- NA_integer_
    n_cross_pairs <- NA_integer_

    if (n_hub == 0L) {
      classification <- "non_hub"
    } else if (length(hub_group_names) >= 2L) {
      # Hub in multiple trait groups -- check module correspondence
      hub_sp_by_group <- lapply(hub_group_names, function(g) {
        h_df$species[h_df$species %in% species_by_trait[[g]] & h_df$is_hub]
      })
      group_indices <- seq_along(hub_group_names)
      pair_mat <- if (length(group_indices) == 2L) {
        matrix(group_indices, nrow = 2)
      } else {
        utils::combn(group_indices, 2)
      }
      cross_pairs <- do.call(rbind, lapply(
        seq_len(ncol(pair_mat)), function(k) {
          expand.grid(sp_a = hub_sp_by_group[[pair_mat[1, k]]],
                      sp_b = hub_sp_by_group[[pair_mat[2, k]]],
                      stringsAsFactors = FALSE)
        }
      ))

      n_cross_pairs <- nrow(cross_pairs)

      corresp <- vapply(seq_len(n_cross_pairs), function(j) {
        mod_a <- h_df$hub_module[h_df$species == cross_pairs$sp_a[j]]
        mod_b <- h_df$hub_module[h_df$species == cross_pairs$sp_b[j]]
        check_correspondence(cross_pairs$sp_a[j], mod_a,
                             cross_pairs$sp_b[j], mod_b)
      }, logical(1))

      if (all(is.na(corresp))) {
        classification <- "multi_trait_hub"
        n_corresponding <- NA_integer_
      } else {
        n_corresponding <- sum(corresp, na.rm = TRUE)
        n_available <- sum(!is.na(corresp))
        classification <- if (n_corresponding / n_available >=
                              correspondence_threshold) {
          "conserved_hub"
        } else {
          "rewired_hub"
        }
      }
    } else if (length(hub_group_names) == 1L) {
      classification <- paste0(hub_group_names, "_specific_hub")
    } else {
      classification <- "sporadic_hub"
    }

    data.frame(
      hog = h,
      classification = classification,
      n_species_hub = as.integer(n_hub),
      n_species_present = as.integer(n_present),
      hub_trait_groups = if (length(hub_group_names) > 0L) {
        paste(sort(hub_group_names), collapse = ",")
      } else {
        NA_character_
      },
      n_corresponding = n_corresponding,
      n_cross_pairs = n_cross_pairs,
      max_centrality = max_cent,
      best_hub_species = best_sp,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, out_rows)
  rownames(result) <- NULL
  result
}


#' Characterize hub genes by regulatory potential
#'
#' Enriches the output of \code{\link{identify_module_hubs}} with
#' metrics that help distinguish regulatory hubs (transcription
#' factors, signalling genes) from downstream effectors.
#'
#' @section Metrics:
#' \describe{
#'   \item{bridge_fraction}{Fraction of a gene's weighted co-expression
#'     connections that reach genes in \emph{other} modules. Values
#'     near 0 indicate a gene embedded within its own module;
#'     values above ~0.3 suggest a regulatory coordinator.
#'     Computed as \code{1 - degree / global_degree}.}
#'   \item{bt_degree_ratio}{Betweenness centrality divided by
#'     within-module weighted degree. A high ratio indicates
#'     importance for information flow (path bridging) relative to
#'     raw connectivity --- a signature of regulatory hubs. A low
#'     ratio indicates a well-connected but non-bridging gene,
#'     typical of housekeeping or effector roles.}
#'   \item{cv}{Coefficient of variation of expression across samples.
#'     Regulators often show higher expression variability than their
#'     targets because they \emph{drive} transcriptional changes.
#'     Only computed when \code{expr} is provided.}
#' }
#'
#' @param hub_result Data frame from \code{\link{identify_module_hubs}}.
#'   Must contain columns \code{gene}, \code{module}, \code{degree},
#'   \code{global_degree}, \code{betweenness}.
#' @param modules Optional output of \code{\link{detect_modules}}.
#'   Currently unused; reserved for future module-aware metrics.
#' @param expr Optional expression matrix (genes as rows, samples as
#'   columns) with rownames matching gene identifiers. When provided,
#'   a \code{cv} column is appended.
#' @param annotations Optional data frame with a \code{gene} column
#'   and any additional annotation columns (e.g., \code{is_tf},
#'   \code{domain}, \code{family}). Left-joined onto the result.
#'
#' @return The input data frame with additional columns:
#'   \code{bridge_fraction}, \code{bt_degree_ratio}, and optionally
#'   \code{cv} and annotation columns.
#'
#' @examples
#' \dontrun{
#' hubs <- identify_module_hubs(mods, net, orthologs)
#' hubs <- characterize_hubs(hubs, mods)
#'
#' # With expression variability
#' hubs <- characterize_hubs(hubs, mods, expr = expr_matrix)
#'
#' # With TF annotations
#' tf_db <- data.frame(gene = c("AT1G01010", "AT2G02020"),
#'                     is_tf = c(TRUE, TRUE),
#'                     family = c("MYB", "WRKY"))
#' hubs <- characterize_hubs(hubs, mods, annotations = tf_db)
#' }
#'
#' @export
characterize_hubs <- function(hub_result, modules = NULL,
                              expr = NULL, annotations = NULL) {
  # --- Validation ---
  if (!is.data.frame(hub_result)) {
    stop("hub_result must be a data frame from identify_module_hubs()")
  }
  req_cols <- c("gene", "module", "degree", "global_degree", "betweenness")
  missing <- setdiff(req_cols, names(hub_result))
  if (length(missing) > 0L) {
    stop("hub_result missing required columns: ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(modules) && (!is.list(modules) || is.null(modules$modules))) {
    stop("modules must be output of detect_modules() or NULL")
  }
  if (!is.null(expr)) {
    if (!is.matrix(expr) || is.null(rownames(expr))) {
      stop("expr must be a matrix with gene names as rownames")
    }
  }
  if (!is.null(annotations)) {
    if (!is.data.frame(annotations) || !"gene" %in% names(annotations)) {
      stop("annotations must be a data frame with a 'gene' column")
    }
  }

  n <- nrow(hub_result)

  # --- Bridge fraction ---
  gd <- hub_result$global_degree
  wd <- hub_result$degree
  bridge_fraction <- ifelse(gd > 0, 1 - wd / gd, NA_real_)
  bridge_fraction <- pmin(pmax(bridge_fraction, 0), 1)
  hub_result$bridge_fraction <- bridge_fraction

  # --- Betweenness / degree ratio ---
  hub_result$bt_degree_ratio <- ifelse(
    wd > 0, hub_result$betweenness / wd, NA_real_
  )

  # --- Expression variability (CV) ---
  if (!is.null(expr)) {
    cv_vec <- rep(NA_real_, n)
    matched <- hub_result$gene %in% rownames(expr)
    if (any(matched)) {
      matched_genes <- hub_result$gene[matched]
      expr_sub <- expr[matched_genes, , drop = FALSE]
      cv_vals <- apply(expr_sub, 1L, function(x) {
        m <- mean(x)
        if (abs(m) < .Machine$double.eps) NA_real_ else stats::sd(x) / abs(m)
      })
      cv_vec[matched] <- cv_vals
    }
    hub_result$cv <- cv_vec
  }

  # --- Annotation join ---
  if (!is.null(annotations)) {
    dup_genes <- duplicated(annotations$gene)
    if (any(dup_genes)) {
      warning("annotations contains duplicate gene entries; ",
              "keeping first occurrence for each gene")
      annotations <- annotations[!dup_genes, , drop = FALSE]
    }
    orig_order <- hub_result$gene
    hub_result <- merge(hub_result, annotations, by = "gene",
                        all.x = TRUE, sort = FALSE)
    # Restore original row order
    hub_result <- hub_result[match(orig_order, hub_result$gene), ,
                             drop = FALSE]
    rownames(hub_result) <- NULL
  }

  hub_result
}


#' Deterministic per-task RNG seed (internal)
#'
#' Maps (root, stream, index) to a legal R seed, so every parallel task seeds
#' itself from its own identity rather than inheriting one. That is what makes
#' the result independent of how mclapply chunks the work: with the default
#' RNGkind a forked child runs parallel:::mc.set.stream(), whose non-L'Ecuyer
#' branch deletes .Random.seed, and the child then re-seeds from clock and PID
#' at its first draw. igraph::cluster_leiden() consumes the R stream, so
#' set.seed() in the parent reaches no worker at n_cores > 1.
#'
#' The arithmetic is done in doubles and folded modulo 2^31 - 1 so it cannot
#' overflow integer range: a naive root + offset form gives NA for a seed near
#' the limit, and set.seed(NA) errors. Products stay exact in double
#' (2654435761 * 1e6 is well under 2^53).
#'
#' Streams: 1 = initial sweep, 2 = K=1 permutations, 100 + iter = consensus
#' sweep at that iteration. The consensus stream must vary with iter because
#' the same resolutions are re-run on a different graph each round.
#'
#' @noRd
.task_seed <- function(root, stream, index) {
  as.integer((root + 2654435761 * stream + 40503 * index) %% 2147483647)
}
