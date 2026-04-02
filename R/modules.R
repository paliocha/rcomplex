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
#' @export
detect_modules <- function(net,
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
                           alpha_k1 = 0.05) {
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

  mat <- net$network
  thr <- net$threshold
  genes <- rownames(mat)

  # Build thresholded adjacency
  adj <- mat
  adj[adj < thr] <- 0

  if (!any(adj[upper.tri(adj)] > 0)) {
    stop("No edges above threshold; cannot detect modules")
  }

  if (!is.null(seed)) set.seed(seed)

  if (method == "sbm") {
    if (!requireNamespace("sbm", quietly = TRUE)) {
      stop("Package 'sbm' is required for method = \"sbm\". ",
           "Install it with install.packages(\"sbm\")")
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

  genes <- rownames(net$network)
  n_genes <- length(genes)

  # Build thresholded adjacency
  adj <- net$network
  adj[adj < net$threshold] <- 0

  if (!any(adj[upper.tri(adj)] > 0)) {
    stop("No edges above threshold; cannot detect modules")
  }

  if (!is.null(seed)) set.seed(seed)

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
  run_initial <- function(res) {
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
    results <- parallel::mclapply(resolutions, run_initial, mc.cores = n_cores)
  } else {
    results <- lapply(resolutions, run_initial)
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
        g_consensus, resolutions, n_iterations, n_cores
      )
      membership <- pick_best_partition(consensus_mems, g)
    }
  } else {
    # K = 1 test
    if (test_k1) {
      k1_result <- test_community_structure(
        g, genes, resolutions, objective_function, n_iterations,
        memberships, edge_list_0, n_perm_k1, n_cores, alpha_k1
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
        g_consensus, resolutions, n_iterations, n_cores
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
#' Consensus graph weights are in [0, 1] (excess co-classification); CPM
#' with resolution > 1 would make all edges repulsive, collapsing to
#' singletons. Modularity is the correct objective per Jeub et al. (2018).
#' @noRd
consensus_leiden_sweep <- function(graph, resolutions, n_iterations,
                                   n_cores = 1L) {
  vertex_names <- igraph::V(graph)$name

  run_one <- function(res) {
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
    results <- parallel::mclapply(resolutions, run_one, mc.cores = n_cores)
    errs <- which(vapply(results, inherits, logical(1), "try-error"))
    if (length(errs)) {
      e <- results[[errs[1L]]]
      stop(attr(e, "condition") %||% as.character(e))
    }
    results
  } else {
    lapply(resolutions, run_one)
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
                                      n_cores = 1L, alpha = 0.05) {
  n_genes <- length(genes)
  n_edges <- igraph::ecount(g)

  lambda_obs <- sparse_excess_spectral_norm_cpp(memberships_obs, n_genes,
                                                 edge_list_0, n_cores)

  run_one_perm <- function(b) {
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
  batch_size <- if (use_mc) n_cores else max(1L, min_for_sig)

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


#' Compare co-expression modules across species
#'
#' Tests all pairs of modules between two species for significant overlap of
#' ortholog-mapped genes, using either Jaccard + permutation or hypergeometric
#' tests.
#'
#' @section Jaccard + permutation method (recommended):
#' Computes observed Jaccard index between ortholog-mapped module_i genes and
#' module_j genes (restricted to ortholog-mappable genes). The null permutes
#' the ortholog mapping (which sp1 genes map to which sp2 genes) using
#' Besag-Clifford adaptive stopping. Q-values use the Liang (2016) discrete
#' method via [DiscreteQvalue::DQ()].
#'
#' @section Hypergeometric method:
#' For each (module_i, module_j) pair, maps module_i genes to species 2 via
#' orthologs, then tests overlap with module_j using [stats::phyper()].
#' Q-values are computed via [qvalue::qvalue()]. The hypergeometric assumes
#' independent sampling, which is violated when orthologs have multi-copy
#' HOGs (one sp1 gene mapping to multiple sp2 genes). This inflates
#' overlap relative to the hypergeometric expectation, making the test
#' **anti-conservative** for conserved modules. The Jaccard permutation
#' null avoids this by preserving module structure and ortholog count
#' while shuffling gene identity.
#'
#' @param modules1 Module detection result for species 1
#'   (output of [detect_modules()]).
#' @param modules2 Module detection result for species 2
#'   (output of [detect_modules()]).
#' @param orthologs Data frame with columns `Species1`, `Species2`, and
#'   `hog` (output of [parse_orthologs()]).
#' @param method Test method: `"jaccard"` (default) or `"hypergeometric"`.
#' @param min_exceedances Besag-Clifford stopping parameter for Jaccard method
#'   (default 50).
#' @param max_permutations Maximum permutations for Jaccard method
#'   (default 10000).
#' @param n_cores Number of threads for Jaccard permutation (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{pairs}{Data frame with one row per module pair, including
#'       `module_sp1`, `module_sp2`, `size_sp1`, `size_sp2`, `overlap`,
#'       `jaccard`, `p.value`, `q.value`}
#'     \item{best_matches}{Data frame with best match per module
#'       (from both species directions)}
#'   }
#'
#' @examples
#' \dontrun{
#' mod_comp <- compare_modules(modules_A, modules_B, orthologs)
#' mod_comp$pairs[mod_comp$pairs$q.value < 0.05, ]
#' }
#'
#' @export
compare_modules <- function(modules1, modules2, orthologs,
                            method = c("jaccard", "hypergeometric"),
                            min_exceedances = 50L,
                            max_permutations = 10000L,
                            n_cores = 1L) {
  method <- match.arg(method)
  min_exceedances <- as.integer(min_exceedances)
  max_permutations <- as.integer(max_permutations)
  if (!is.list(modules1) || is.null(modules1$module_genes)) {
    stop("modules1 must be output from detect_modules()")
  }
  if (!is.list(modules2) || is.null(modules2$module_genes)) {
    stop("modules2 must be output from detect_modules()")
  }
  if (!all(c("Species1", "Species2") %in% names(orthologs))) {
    stop("orthologs must have columns: Species1, Species2")
  }

  mg1 <- modules1$module_genes
  mg2 <- modules2$module_genes

  # Filter orthologs to genes present in both module sets
  ortho <- orthologs[orthologs$Species1 %in% names(modules1$modules) &
                       orthologs$Species2 %in% names(modules2$modules), ,
                     drop = FALSE]
  ortho <- unique(ortho[, c("Species1", "Species2"), drop = FALSE])

  if (nrow(ortho) == 0L) {
    stop("No orthologs found in both module sets")
  }

  sp2_mappable <- unique(ortho$Species2)
  sp1_to_sp2 <- split(ortho$Species2, ortho$Species1)

  if (method == "hypergeometric") {
    pairs <- compare_modules_hypergeometric(
      mg1, mg2, sp1_to_sp2, sp2_mappable
    )
  } else {
    pairs <- compare_modules_jaccard(
      mg1, mg2, ortho, sp2_mappable, sp1_to_sp2,
      min_exceedances, max_permutations, n_cores
    )
  }

  best <- compute_best_matches(pairs)
  list(pairs = pairs, best_matches = best)
}


#' Hypergeometric module comparison (internal)
#' @noRd
compare_modules_hypergeometric <- function(mg1, mg2, sp1_to_sp2, sp2_mappable) {
  mod_names1 <- names(mg1)
  mod_names2 <- names(mg2)
  n1 <- length(mod_names1)
  n2 <- length(mod_names2)
  N <- length(sp2_mappable)

  # Binary membership matrices: rows = sp2_mappable genes, cols = modules
  # M1: sp1 modules mapped to sp2 space; M2: sp2 modules in mappable space
  sp2_idx <- stats::setNames(seq_along(sp2_mappable), sp2_mappable)
  M1 <- matrix(0L, nrow = N, ncol = n1)
  M2 <- matrix(0L, nrow = N, ncol = n2)

  k_vec <- integer(n1)  # mapped set size per sp1 module
  for (i in seq_len(n1)) {
    mapped <- unique(unlist(sp1_to_sp2[mg1[[i]]], use.names = FALSE))
    if (!is.null(mapped)) {
      hits <- sp2_idx[mapped]
      hits <- hits[!is.na(hits)]
      M1[hits, i] <- 1L
      k_vec[i] <- length(hits)
    }
  }

  m_vec <- integer(n2)  # mappable set size per sp2 module
  for (j in seq_len(n2)) {
    hits <- sp2_idx[mg2[[j]]]
    hits <- hits[!is.na(hits)]
    M2[hits, j] <- 1L
    m_vec[j] <- length(hits)
  }

  # All overlaps at once via BLAS (n1 x n2 matrix)
  overlap_mat <- crossprod(M1, M2)

  # Expand to pair vectors (column-major: i varies fastest)
  idx_i <- rep(seq_len(n1), n2)
  idx_j <- rep(seq_len(n2), each = n1)
  ol <- as.integer(overlap_mat)
  k <- k_vec[idx_i]
  m <- m_vec[idx_j]
  un <- k + m - ol

  pairs <- data.frame(
    module_sp1 = mod_names1[idx_i],
    module_sp2 = mod_names2[idx_j],
    size_sp1 = vapply(mg1, length, integer(1))[idx_i],
    size_sp2 = vapply(mg2, length, integer(1))[idx_j],
    overlap = ol,
    jaccard = ifelse(un > 0L, ol / un, 0),
    p.value = ifelse(ol > 0L & k > 0L & m > 0L,
                     stats::phyper(ol - 1L, m, N - m, k, lower.tail = FALSE),
                     1)
  )

  n_pairs <- nrow(pairs)
  pairs$q.value <- if (n_pairs < 2L) pairs$p.value else {
    compute_qvalues(pairs$p.value)
  }
  pairs
}


#' Jaccard + permutation module comparison (internal)
#' @noRd
compare_modules_jaccard <- function(mg1, mg2, ortho, sp2_mappable,
                                    sp1_to_sp2,
                                    min_exceedances, max_permutations,
                                    n_cores) {
  mod_names1 <- names(mg1)
  mod_names2 <- names(mg2)

  # Integer-index mapping for sp2 genes (mappable only)
  sp2_idx_map <- stats::setNames(seq_along(sp2_mappable) - 1L, sp2_mappable)
  n_sp2 <- length(sp2_mappable)

  # Unique sp1 genes and their indices
  sp1_unique <- unique(ortho$Species1)
  sp1_idx_map <- stats::setNames(seq_along(sp1_unique) - 1L, sp1_unique)

  # Per ortholog row: sp1 unique gene index and sp2 gene index
  ortho_sp1_int <- as.integer(sp1_idx_map[ortho$Species1])
  ortho_sp2_int <- as.integer(sp2_idx_map[ortho$Species2])

  # Per sp2 module: mappable gene indices
  mod_sp2_sets <- lapply(mg2, function(genes) {
    as.integer(sp2_idx_map[intersect(genes, sp2_mappable)])
  })

  # Per sp1 module: unique sp1 gene indices (for C++ permutation)
  mod1_sp1_genes <- lapply(mg1, function(genes) {
    as.integer(sp1_idx_map[intersect(genes, sp1_unique)])
  })

  # Compute observed Jaccard for all pairs via crossprod
  n1 <- length(mod_names1)
  n2 <- length(mod_names2)
  n_pairs <- n1 * n2

  # Binary membership matrices (0-based sp2 indices -> 1-based rows)
  M1 <- matrix(0L, nrow = n_sp2, ncol = n1)
  k_vec <- integer(n1)
  for (i in seq_len(n1)) {
    mapped <- unique(unlist(sp1_to_sp2[mg1[[i]]], use.names = FALSE))
    if (!is.null(mapped)) {
      hits <- sp2_idx_map[mapped]
      hits <- hits[!is.na(hits)] + 1L  # 0-based -> 1-based row index
      M1[hits, i] <- 1L
      k_vec[i] <- length(hits)
    }
  }

  M2 <- matrix(0L, nrow = n_sp2, ncol = n2)
  m_vec <- integer(n2)
  for (j in seq_len(n2)) {
    genes_j <- mod_sp2_sets[[j]]
    if (length(genes_j) > 0L) {
      M2[genes_j + 1L, j] <- 1L  # 0-based -> 1-based row index
      m_vec[j] <- length(genes_j)
    }
  }

  overlap_mat <- crossprod(M1, M2)  # n1 x n2

  # Expand to pair vectors (column-major: i varies fastest)
  idx_i <- rep(seq_len(n1), n2)
  idx_j <- rep(seq_len(n2), each = n1)
  obs_overlap <- as.integer(overlap_mat)
  k <- k_vec[idx_i]
  m <- m_vec[idx_j]
  un <- k + m - obs_overlap
  obs_jaccard <- ifelse(un > 0L, obs_overlap / un, 0)

  perm_result <- module_jaccard_permutation_cpp(
    ortho_sp1_gene = ortho_sp1_int,
    ortho_sp2_gene = ortho_sp2_int,
    n_sp1_unique = length(sp1_unique),
    n_sp2_universe = n_sp2,
    mod1_sp1_genes = mod1_sp1_genes,
    mod_sp2_sets = mod_sp2_sets,
    mod_i_idx = idx_i - 1L,
    mod_j_idx = idx_j - 1L,
    obs_jaccard = obs_jaccard,
    min_exceedances = min_exceedances,
    max_permutations = max_permutations,
    n_cores = n_cores
  )

  pairs <- data.frame(
    module_sp1 = mod_names1[idx_i],
    module_sp2 = mod_names2[idx_j],
    size_sp1 = vapply(mg1, length, integer(1))[idx_i],
    size_sp2 = vapply(mg2, length, integer(1))[idx_j],
    overlap = obs_overlap, jaccard = obs_jaccard,
    p.value = perm_result$p_value
  )

  # Discrete q-values (Liang 2016) for Besag-Clifford p-values
  bc_support <- bc_pvalue_support(min_exceedances, max_permutations)

  pairs$q.value <- if (n_pairs < 2L) pairs$p.value else {
    DiscreteQvalue::DQ(
      pairs$p.value, ss = bc_support, method = "Liang"
    )$q.values
  }
  pairs
}


#' Extract best match per module in one direction (internal)
#' @noRd
best_match_direction <- function(pairs, from_col, to_col, species_label) {
  do.call(rbind, lapply(
    split(pairs, pairs[[from_col]]),
    function(df) {
      best <- df[which.min(df$p.value), , drop = FALSE]
      data.frame(
        module = best[[from_col]],
        species = species_label,
        best_match = best[[to_col]],
        overlap = best$overlap,
        jaccard = best$jaccard,
        p.value = best$p.value,
        q.value = best$q.value,
        n_significant = sum(df$q.value < 0.05)
      )
    }
  ))
}


#' Compute best matches for each module (internal)
#' @noRd
compute_best_matches <- function(pairs) {
  if (nrow(pairs) == 0L) {
    return(data.frame(
      module = character(0), species = character(0),
      best_match = character(0), overlap = integer(0),
      jaccard = numeric(0), p.value = numeric(0),
      q.value = numeric(0), n_significant = integer(0)
    ))
  }

  result <- rbind(
    best_match_direction(pairs, "module_sp1", "module_sp2", "sp1"),
    best_match_direction(pairs, "module_sp2", "module_sp1", "sp2")
  )
  rownames(result) <- NULL
  result
}


#' Classify modules by conservation status
#'
#' Classifies each module as conserved, partially conserved, or
#' species-specific based on cross-species module comparison results.
#'
#' @section Classification criteria:
#' \describe{
#'   \item{Conserved}{Best-match q-value < `alpha` AND Jaccard >= `jaccard_threshold`}
#'   \item{Partially conserved}{Best-match q-value < `alpha` AND
#'     Jaccard < `jaccard_threshold` (module split or merged across species)}
#'   \item{Species-specific}{No significant match (best-match q-value >= `alpha`)}
#' }
#'
#' @param module_comparison Output of [compare_modules()].
#' @param alpha Significance threshold for q-values (default 0.05).
#' @param jaccard_threshold Minimum Jaccard index to call a module conserved
#'   (default 0.1).
#'
#' @return A data frame with one row per module (from both species), including:
#'   \describe{
#'     \item{module}{Module identifier}
#'     \item{species}{Which species this module belongs to (`"sp1"` or `"sp2"`)}
#'     \item{classification}{One of `"conserved"`, `"partially_conserved"`,
#'       or `"species_specific"`}
#'     \item{best_match}{Best matching module in the other species}
#'     \item{best_jaccard}{Jaccard index with best match}
#'     \item{best_q}{Q-value for best match}
#'     \item{n_significant}{Number of significant matches in the other species}
#'   }
#'
#' @examples
#' \dontrun{
#' classes <- classify_modules(mod_comp)
#' table(classes$classification)  # conserved / partially / species_specific
#' }
#'
#' @export
classify_modules <- function(module_comparison,
                             alpha = 0.05,
                             jaccard_threshold = 0.1) {
  if (!is.list(module_comparison) ||
        !all(c("pairs", "best_matches") %in% names(module_comparison))) {
    stop("module_comparison must be output from compare_modules()")
  }

  best <- module_comparison$best_matches

  if (nrow(best) == 0L) {
    return(data.frame(
      module = character(0), species = character(0),
      classification = character(0), best_match = character(0),
      best_jaccard = numeric(0), best_q = numeric(0),
      n_significant = integer(0)
    ))
  }

  significant <- best$q.value < alpha
  classification <- ifelse(
    !significant, "species_specific",
    ifelse(best$jaccard >= jaccard_threshold,
           "conserved", "partially_conserved")
  )

  data.frame(
    module = best$module,
    species = best$species,
    classification = classification,
    best_match = best$best_match,
    best_jaccard = best$jaccard,
    best_q = best$q.value,
    n_significant = best$n_significant
  )
}


#' Coarsen a module partition to fewer modules
#'
#' Merges modules in a fine-grained partition to a target count by
#' hierarchical clustering on inter-module edge density. This produces
#' matched-granularity partitions for fair cross-species comparison when
#' species have intrinsically different modular structure.
#'
#' @section Algorithm:
#' \enumerate{
#'   \item Compute inter-module edge density for all module pairs:
#'     \code{n_edges(i,j) / (size_i * size_j)}.
#'   \item Hierarchical clustering (complete linkage) on distance
#'     \code{1 - D / max(D)}.
#'   \item Cut the dendrogram at \code{target_n_modules}.
#'   \item Relabel module assignments contiguously and recompute
#'     module gene lists, modularity, and sizes.
#' }
#'
#' @param modules Output of [detect_modules()].
#' @param target_n_modules Integer, desired number of modules after
#'   coarsening. Must be at least 1 and less than the current number
#'   of modules.
#'
#' @return A list compatible with [detect_modules()] output (usable
#'   directly with [compare_modules()]), with additional components:
#'   \describe{
#'     \item{merge_map}{Data frame with columns \code{original_module}
#'       and \code{coarsened_module} mapping old to new module IDs.}
#'     \item{merge_dendrogram}{The [hclust] object used for merging.}
#'   }
#'
#' @examples
#' \dontrun{
#' mods <- detect_modules(net, resolution = seq(0.25, 2.5, by = 0.25))
#' mods$n_modules  # 52
#' coarse <- coarsen_modules(mods, target_n_modules = 15)
#' coarse$n_modules  # 15
#' }
#'
#' @export
coarsen_modules <- function(modules, target_n_modules) {
  if (!is.list(modules) || is.null(modules$modules) ||
      is.null(modules$module_genes) || is.null(modules$graph)) {
    stop("modules must be output from detect_modules()")
  }

  target_n_modules <- as.integer(target_n_modules)
  if (is.na(target_n_modules) || target_n_modules < 1L)
    stop("target_n_modules must be a finite integer >= 1")

  n_current <- modules$n_modules
  if (target_n_modules >= n_current) return(modules)

  g <- modules$graph
  mem <- modules$modules
  mg <- modules$module_genes
  mod_ids <- names(mg)
  n_mod <- length(mod_ids)

  # Module sizes
  mod_sizes <- vapply(mg, length, integer(1))

  # Inter-module edge density via vectorized tabulation
  el <- igraph::as_edgelist(g, names = TRUE)
  mod_idx <- stats::setNames(seq_len(n_mod), mod_ids)
  i_mod <- mod_idx[as.character(mem[el[, 1]])]
  j_mod <- mod_idx[as.character(mem[el[, 2]])]

  # Keep only inter-module edges, count by (i, j) pair.
  # Row-major linear index; matrix() fills column-major so D is
  # transposed, but D + t(D) symmetrizes immediately after.
  inter <- i_mod != j_mod
  counts <- tabulate(
    (i_mod[inter] - 1L) * n_mod + j_mod[inter],
    nbins = n_mod * n_mod
  )
  D <- matrix(counts, n_mod, n_mod)
  D <- D + t(D)  # symmetrize (each undirected edge counted once)

  # Normalize: edge density = n_edges / (size_i * size_j)
  size_outer <- outer(mod_sizes, mod_sizes)
  diag(size_outer) <- 1  # avoid division by zero on diagonal
  D <- D / size_outer
  diag(D) <- 0

  # Hierarchical clustering (complete linkage on 1 - normalized density)
  max_d <- max(D)
  if (max_d == 0) {
    warning("No inter-module edges; merging by module order")
    dist_mat <- stats::dist(seq_len(n_mod))
  } else {
    dist_mat <- stats::as.dist(1 - D / max_d)
  }
  hc <- stats::hclust(dist_mat, method = "complete")
  new_groups <- stats::cutree(hc, k = target_n_modules)

  # Build merge map
  merge_map <- data.frame(
    original_module = mod_ids,
    coarsened_module = as.integer(new_groups),
    stringsAsFactors = FALSE
  )

  # Relabel gene-level membership
  group_lookup <- stats::setNames(new_groups, mod_ids)
  new_mem <- group_lookup[as.character(mem)]
  names(new_mem) <- names(mem)
  new_mem <- as.integer(new_mem)
  names(new_mem) <- names(mem)

  new_module_genes <- split(names(new_mem), new_mem)

  list(
    modules = new_mem,
    module_genes = new_module_genes,
    n_modules = length(new_module_genes),
    modularity = igraph::modularity(g, new_mem),
    graph = g,
    method = paste0(modules$method, "_coarsened"),
    params = c(modules$params, list(
      coarsened_from = n_current,
      target_n_modules = target_n_modules
    )),
    merge_map = merge_map,
    merge_dendrogram = hc
  )
}


#' Compare modules across contrasting conditions
#'
#' Wraps [compare_modules()], [classify_modules()], and optionally
#' [coarsen_modules()] into a single call that processes all contrasts
#' (species pairs, tissue pairs, condition pairs, etc.), tags group-
#' specificity, and performs matched-scale comparison when module
#' counts differ substantially.
#'
#' @param modules Named list of [detect_modules()] outputs. Names are
#'   sample/species/condition identifiers.
#' @param orthologs Data frame with columns `Species1`, `Species2`,
#'   and `hog`.
#' @param pairs Data frame defining contrasts. Must have columns `sp1`
#'   and `sp2` (identifiers matching names of `modules`). An optional
#'   `pair_name` column is used for labelling; if absent,
#'   `"sp1.sp2"` is generated.
#' @param group Optional named character vector mapping identifiers to
#'   group labels (e.g., `c(BDIS = "annual", BSYL = "perennial")` for
#'   species, or `c(leaf = "source", root = "sink")` for tissues).
#'   When provided, species-specific modules are tagged with the group
#'   label of their side. All identifiers in `pairs` must have entries
#'   in this vector.
#' @param method Comparison method passed to [compare_modules()]:
#'   `"jaccard"` (default) or `"hypergeometric"`.
#' @param matched_scale Logical. If `TRUE` (default), contrasts whose
#'   module-count ratio exceeds `coarsen_ratio` are re-compared at
#'   matched granularity via [coarsen_modules()].
#' @param coarsen_ratio Numeric threshold for triggering coarsening
#'   (default 2). Only used when `matched_scale = TRUE`.
#' @param alpha Significance level for [classify_modules()] (default
#'   0.05).
#' @param jaccard_threshold Minimum Jaccard for a module to be called
#'   conserved (default 0.1).
#' @param min_exceedances Besag-Clifford parameter for Jaccard method
#'   (default 50).
#' @param max_permutations Maximum permutations for Jaccard method
#'   (default 10000).
#' @param n_cores Number of threads (default 1).
#'
#' @return A list with components:
#'   \describe{
#'     \item{classification}{Data frame with columns `pair_name`,
#'       `module`, `species`, `classification`, `group` (if `group`
#'       provided), `best_match`, `best_jaccard`, `best_q`,
#'       `n_significant`.}
#'     \item{matched_classification}{Same structure at matched scale,
#'       or `NULL` if no contrast exceeds `coarsen_ratio`.}
#'     \item{summary}{Group-aggregated counts per contrast and scale.}
#'     \item{module_counts}{Data frame with `pair_name`, `n_sp1`,
#'       `n_sp2`, `ratio`.}
#'     \item{raw}{Named list of per-contrast [compare_modules()]
#'       outputs, keyed as `"sp1.sp2"`. Compatible with
#'       [classify_hub_conservation(module_comparisons=)].}
#'   }
#'
#' @examples
#' \dontrun{
#' # Species-pair contrasts with trait groups
#' mod_results <- compare_modules_paired(
#'   modules, orthologs,
#'   pairs = data.frame(sp1 = c("BDIS", "HVUL"),
#'                      sp2 = c("BSYL", "HJUB")),
#'   group = c(BDIS = "annual", BSYL = "perennial",
#'             HVUL = "annual", HJUB = "perennial"),
#'   method = "jaccard", n_cores = 4L
#' )
#'
#' # Tissue contrasts (no trait grouping)
#' mod_results <- compare_modules_paired(
#'   modules, orthologs,
#'   pairs = data.frame(sp1 = "leaf", sp2 = "root")
#' )
#' }
#'
#' @export
compare_modules_paired <- function(modules, orthologs, pairs,
                                   group = NULL,
                                   method = c("jaccard", "hypergeometric"),
                                   matched_scale = TRUE,
                                   coarsen_ratio = 2,
                                   alpha = 0.05,
                                   jaccard_threshold = 0.1,
                                   min_exceedances = 50L,
                                   max_permutations = 10000L,
                                   n_cores = 1L) {
  method <- match.arg(method)
  if (!is.list(modules) || is.null(names(modules)))
    stop("modules must be a named list keyed by identifier")
  if (!all(c("sp1", "sp2") %in% names(pairs)))
    stop("pairs must have columns 'sp1' and 'sp2'")
  missing_sp <- setdiff(c(pairs$sp1, pairs$sp2), names(modules))
  if (length(missing_sp))
    stop("identifiers not found in modules: ",
         paste(missing_sp, collapse = ", "))

  if (!"pair_name" %in% names(pairs))
    pairs$pair_name <- paste(pairs$sp1, pairs$sp2, sep = ".")

  has_group <- !is.null(group)
  if (has_group) {
    missing_grp <- setdiff(c(pairs$sp1, pairs$sp2), names(group))
    if (length(missing_grp))
      stop("group missing entries for: ",
           paste(missing_grp, collapse = ", "))
  }

  # Tag one contrast's classification with group info
  tag_classification <- function(cls, s1, s2, pair_name) {
    cls$pair_name <- pair_name
    if (has_group) {
      cls$group <- ifelse(
        cls$classification != "species_specific", "conserved",
        ifelse(cls$species == "sp1", group[s1], group[s2])
      )
    }
    cls
  }

  # Process each contrast
  n_pairs <- nrow(pairs)
  raw_list <- vector("list", n_pairs)
  class_list <- vector("list", n_pairs)
  matched_list <- vector("list", n_pairs)
  n_sp1 <- integer(n_pairs)
  n_sp2 <- integer(n_pairs)

  for (p in seq_len(n_pairs)) {
    s1 <- pairs$sp1[p]
    s2 <- pairs$sp2[p]
    pn <- pairs$pair_name[p]
    m1 <- modules[[s1]]
    m2 <- modules[[s2]]
    n_sp1[p] <- m1$n_modules
    n_sp2[p] <- m2$n_modules

    comp <- compare_modules(m1, m2, orthologs, method = method,
                            min_exceedances = min_exceedances,
                            max_permutations = max_permutations,
                            n_cores = n_cores)
    raw_list[[p]] <- comp
    class_list[[p]] <- tag_classification(
      classify_modules(comp, alpha = alpha,
                       jaccard_threshold = jaccard_threshold),
      s1, s2, pn
    )

    ratio <- max(n_sp1[p], n_sp2[p]) / max(1L, min(n_sp1[p], n_sp2[p]))
    if (matched_scale && ratio > coarsen_ratio) {
      target <- min(n_sp1[p], n_sp2[p])
      cm1 <- if (m1$n_modules > target) coarsen_modules(m1, target) else m1
      cm2 <- if (m2$n_modules > target) coarsen_modules(m2, target) else m2
      comp_m <- compare_modules(cm1, cm2, orthologs, method = method,
                                min_exceedances = min_exceedances,
                                max_permutations = max_permutations,
                                n_cores = n_cores)
      matched_list[[p]] <- tag_classification(
        classify_modules(comp_m, alpha = alpha,
                         jaccard_threshold = jaccard_threshold),
        s1, s2, pn
      )
    }
  }

  names(raw_list) <- paste(pairs$sp1, pairs$sp2, sep = ".")
  classification <- do.call(rbind, class_list)

  matched_notnull <- !vapply(matched_list, is.null, logical(1))
  matched_classification <- if (any(matched_notnull)) {
    do.call(rbind, matched_list[matched_notnull])
  }

  # Summary table (bind_rows-safe: union column sets, fill missing with 0)
  ratios <- pmax(n_sp1, n_sp2) / pmax(1L, pmin(n_sp1, n_sp2))
  module_counts <- data.frame(
    pair_name = pairs$pair_name,
    n_sp1 = n_sp1, n_sp2 = n_sp2,
    ratio = round(ratios, 1),
    stringsAsFactors = FALSE
  )

  build_summary <- function(cls, scale_label) {
    if (is.null(cls) || nrow(cls) == 0L) return(NULL)
    count_col <- if (has_group) "group" else "classification"
    agg <- stats::aggregate(
      stats::as.formula(paste("module ~", "pair_name +", count_col)),
      data = cls, FUN = length
    )
    names(agg)[3] <- "n"
    wide <- stats::reshape(agg, idvar = "pair_name",
                           timevar = count_col, direction = "wide")
    names(wide) <- sub("^n\\.", "", names(wide))
    wide[is.na(wide)] <- 0L
    wide$scale <- scale_label
    wide
  }

  summary_parts <- Filter(Negate(is.null), list(
    build_summary(classification, "natural"),
    build_summary(matched_classification, "matched")
  ))
  if (length(summary_parts) > 0L) {
    # Union columns across scales (matched may lack some group levels)
    all_cols <- unique(unlist(lapply(summary_parts, names)))
    summary_parts <- lapply(summary_parts, function(df) {
      missing <- setdiff(all_cols, names(df))
      for (col in missing) df[[col]] <- 0L
      df[all_cols]
    })
    summary_df <- do.call(rbind, summary_parts)
    rownames(summary_df) <- NULL
  } else {
    summary_df <- NULL
  }

  list(
    classification = classification,
    matched_classification = matched_classification,
    summary = summary_df,
    module_counts = module_counts,
    raw = raw_list
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
#' @export
identify_module_hubs <- function(modules, net, orthologs = NULL,
                                 comparison = NULL,
                                 centrality = c("degree", "betweenness",
                                                "eigenvector"),
                                 top_n = NULL,
                                 top_fraction = 0.1,
                                 min_module_size = 3L) {
  centrality <- match.arg(centrality)

  if (!is.list(modules) || is.null(modules$module_genes) ||
        is.null(modules$graph) || is.null(modules$modules)) {
    stop("modules must be output from detect_modules()")
  }
  if (!is.list(net) || is.null(net$network)) {
    stop("net must be output from compute_network()")
  }
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
#' @param module_comparisons Optional named list of [compare_modules()] outputs
#'   keyed by alphabetically sorted species pair (e.g. `"SP_A.SP_C"`).
#'   Required for the conserved_hub vs rewired_hub distinction.
#' @param alpha Significance threshold for module correspondence (default 0.05).
#' @param jaccard_threshold Jaccard threshold for module correspondence
#'   (default 0.1).
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
#' }
#'
#' @export
classify_hub_conservation <- function(hub_results, species_trait,
                                      module_comparisons = NULL,
                                      alpha = 0.05,
                                      jaccard_threshold = 0.1,
                                      min_trait_fraction = 0.5,
                                      correspondence_threshold = 0.5) {
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
