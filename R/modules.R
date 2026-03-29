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
#' @param n_cores Number of parallel cores for consensus mode (default 1).
#'   Uses a PSOCK cluster (\code{parallel::parLapply()}) for fork-safe
#'   parallelism (compatible with prior torch/CUDA initialization).
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
#'       (null distribution), \code{p_value}, and \code{has_structure}.}
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

  # Create PSOCK cluster if needed (used for initial sweep and iterations)
  cl <- NULL
  if (n_cores > 1L) {
    cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }

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

  if (!is.null(cl)) {
    results <- parallel::parLapply(cl, resolutions, run_initial)
  } else {
    results <- lapply(resolutions, run_initial)
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
    expected_coclassification = NA_real_,
    stringsAsFactors = FALSE
  )

  # ---- Consensus clustering ----
  n_consensus_iter <- 0L
  # Save initial sweep for fallback if consensus graph becomes empty
  initial_memberships <- memberships
  k1_result <- NULL

  if (!is.null(consensus_threshold)) {
    # Fixed threshold: single pass, no iteration (dense co-classification)
    cc_result <- build_coclassification_cpp(memberships, n_genes, FALSE)
    resolution_scan$expected_coclassification <- cc_result$expected

    consensus_adj <- cc_result$coclassification
    rm(cc_result)  # drop ref before diag<- triggers COW
    consensus_adj[consensus_adj < consensus_threshold] <- 0
    diag(consensus_adj) <- 0
    dimnames(consensus_adj) <- list(genes, genes)

    g_consensus <- igraph::graph_from_adjacency_matrix(
      consensus_adj, mode = "upper", weighted = TRUE, diag = FALSE
    )
    rm(consensus_adj)

    if (igraph::ecount(g_consensus) == 0L) {
      best_r <- which.max(scan_modularity)
      membership <- initial_memberships[[best_r]]
    } else {
      consensus_mems <- consensus_leiden_sweep(
        g_consensus, resolutions, n_iterations, cl
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
        memberships, n_genes, edge_list_0
      )

      if (iter == 1L) {
        resolution_scan$expected_coclassification <- cc_sparse$expected
      }

      keep <- cc_sparse$excess > 0
      if (!any(keep)) {
        memberships <- initial_memberships
        break
      }

      consensus_el <- cbind(cc_sparse$from[keep] + 1L,
                            cc_sparse$to[keep] + 1L)
      g_consensus <- igraph::make_empty_graph(n = n_genes, directed = FALSE)
      igraph::V(g_consensus)$name <- genes
      g_consensus <- igraph::add_edges(
        g_consensus,
        as.vector(t(consensus_el)),
        weight = cc_sparse$excess[keep]
      )

      new_memberships <- consensus_leiden_sweep(
        g_consensus, resolutions, n_iterations, cl
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
                                   cl = NULL) {
  run_one <- function(res) {
    comm <- igraph::cluster_leiden(
      graph,
      resolution = res,
      objective_function = "modularity",
      n_iterations = as.integer(n_iterations)
    )
    mem <- igraph::membership(comm)
    names(mem) <- igraph::V(graph)$name
    mem
  }

  if (!is.null(cl)) {
    parallel::parLapply(cl, resolutions, run_one)
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
#' @noRd
test_community_structure <- function(g, genes, resolutions, objective_function,
                                      n_iterations, memberships_obs,
                                      edge_list_0, n_perm = 100L,
                                      n_cores = 1L, alpha = 0.05) {
  n_genes <- length(genes)

  lambda_obs <- sparse_excess_spectral_norm_cpp(memberships_obs, n_genes,
                                                 edge_list_0)

  run_one_perm <- function(b) {
    g_perm <- igraph::rewire(g, igraph::keeping_degseq(
      niter = 10L * igraph::ecount(g)))
    igraph::E(g_perm)$weight <- igraph::E(g)$weight[
      sample.int(igraph::ecount(g))]

    mems_perm <- lapply(resolutions, function(res) {
      comm <- igraph::cluster_leiden(
        g_perm, resolution = res,
        objective_function = objective_function,
        n_iterations = as.integer(n_iterations)
      )
      mem <- igraph::membership(comm)
      names(mem) <- igraph::V(g_perm)$name
      mem
    })

    el_perm <- igraph::as_edgelist(g_perm, names = FALSE) - 1L
    storage.mode(el_perm) <- "integer"
    sparse_excess_spectral_norm_cpp(mems_perm, n_genes, el_perm)
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    lambda_null <- unlist(parallel::mclapply(
      seq_len(n_perm), run_one_perm, mc.cores = n_cores
    ))
  } else {
    lambda_null <- vapply(seq_len(n_perm), run_one_perm, numeric(1))
  }

  p_value <- (1 + sum(lambda_null >= lambda_obs)) / (1 + n_perm)

  list(
    lambda_obs = lambda_obs,
    lambda_null = lambda_null,
    p_value = p_value,
    has_structure = p_value < alpha
  )
}


#' Compare co-expression modules across species
#'
#' Tests all pairs of modules between two species for significant overlap of
#' ortholog-mapped genes, using either hypergeometric or Jaccard + permutation
#' tests.
#'
#' @section Hypergeometric method:
#' For each (module_i, module_j) pair, maps module_i genes to species 2 via
#' orthologs, then tests overlap with module_j using [stats::phyper()].
#' Q-values are computed via [qvalue::qvalue()].
#'
#' @section Jaccard + permutation method:
#' Computes observed Jaccard index between ortholog-mapped module_i genes and
#' module_j genes (restricted to ortholog-mappable genes). The null permutes
#' the ortholog mapping (which sp1 genes map to which sp2 genes) using
#' Besag-Clifford adaptive stopping. Q-values use the Liang (2016) discrete
#' method via [DiscreteQvalue::DQ()].
#'
#' @param modules1 Module detection result for species 1
#'   (output of [detect_modules()]).
#' @param modules2 Module detection result for species 2
#'   (output of [detect_modules()]).
#' @param orthologs Data frame with columns `Species1`, `Species2`, and
#'   `OrthoGroup` (output of [parse_orthologs()]).
#' @param method Test method: `"hypergeometric"` (default) or `"jaccard"`.
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
#' @export
compare_modules <- function(modules1, modules2, orthologs,
                            method = c("hypergeometric", "jaccard"),
                            min_exceedances = 50L,
                            max_permutations = 10000L,
                            n_cores = 1L) {
  method <- match.arg(method)
  min_exceedances <- as.integer(min_exceedances)
  max_permutations <- as.integer(max_permutations)
  n_cores <- as.integer(n_cores)

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
                     1),
    stringsAsFactors = FALSE
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
    p.value = perm_result$p_value,
    stringsAsFactors = FALSE
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
        n_significant = sum(df$q.value < 0.05),
        stringsAsFactors = FALSE
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
      q.value = numeric(0), n_significant = integer(0),
      stringsAsFactors = FALSE
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
      n_significant = integer(0),
      stringsAsFactors = FALSE
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
    n_significant = best$n_significant,
    stringsAsFactors = FALSE
  )
}
