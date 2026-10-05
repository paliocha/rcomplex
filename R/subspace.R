#' Partition-free species-pair agreement of Laplacian subspaces
#'
#' Scores how much of one species' large-scale co-expression structure
#' reappears in another's, without detecting modules in either. For every
#' ordered species pair (A, B) the statistic is the mean squared cosine of
#' the principal angles between A's top-K Laplacian subspace and B's,
#' mapped through orthology.
#'
#' @details
#' Per species the network is read at `net$threshold` as a binary
#' adjacency W, restricted to its largest connected component. The
#' normalised Laplacian is \eqn{L = I - D^{-1/2} W D^{-1/2}}; \eqn{U}
#' (genes x K) holds the eigenvectors of its K smallest non-trivial
#' eigenvalues (the trivial one, eigenvalue 0, is dropped). For an
#' ordered pair, \eqn{P_{AB}} is the gene(A) x gene(B) indicator of a
#' shared ortholog group, each row scaled to sum 1 (genes of A without a
#' B ortholog are zero rows), and
#' \deqn{S_{AB} = \|U_A^\top P_{AB} U_B\|_F^2 / K.}
#' With \eqn{P} the identity, this is the mean squared cosine of the
#' principal angles between the two subspaces: 1 when they coincide, about
#' \eqn{K / n} for random ones. It is the agreement term that co-regularised
#' multi-view spectral clustering (Kumar, Rai & Daume 2011) maximises, read
#' here as a statistic rather than optimised.
#'
#' The statistic is partition-free: it compares smooth bases of the whole
#' network, so no module is detected and per-species module instability
#' never enters (on the n = 20 Pooideae leaf data, Leiden modules did not
#' replicate between sample halves). It does not say *which* genes agree.
#'
#' The null permutes the ortholog-group labels among A's genes (rows of
#' \eqn{P_{AB}}) and keeps both subspaces fixed: `n_null` draws give
#' `S_null_mean`, `S_null_sd` and `z = (S - S_null_mean) / S_null_sd`. On
#' the Pooideae and wood data this null sat at the level of the two
#' shuffled-expression nulls, which need a fresh eigendecomposition per
#' draw, and the real S was above all three on every pair
#' (`dev/design-notes/module-engine-borrow-scope.md`, Section 10.5).
#'
#' K is a resolution knob, to be calibrated against the null on the data
#' at hand: K = 100 separated real from null best on the Pooideae networks
#' (about 20,000 genes), K = 10 on the EVOTREE wood networks. Binary and
#' MR-weighted adjacency gave the same S there, so only the binary
#' adjacency is used. S itself is small (about 0.02 on real pairs); read it
#' against a ceiling, which is the same statistic on the networks of the
#' two sample halves of one species with an identity map (`hog = gene`):
#' 0.06-0.25 on Pooideae, 0.18-0.53 on wood. [as_preservation_matrix()]
#' passes the symmetrised z to [preservation_matrix_test()] as a trait
#' readout.
#'
#' @param nets Named list of at least two network objects from
#'   [compute_network()], one per species.
#' @param hog_map Data frame with columns `species`, `gene` and `hog`, one
#'   row per gene copy (many-to-many), as in [recurrence_graph()]. Genes
#'   outside a species' largest component are ignored.
#' @param K Number of non-trivial eigenvectors per species. Must be below
#'   the size of every species' largest connected component.
#' @param n_null Number of label permutations per ordered pair.
#' @param n_cores Species eigendecompositions run in parallel (forked) on
#'   Unix when this exceeds 1 and the BLAS is fork-safe. Draws are made in
#'   the parent, so the result does not depend on it.
#' @param seed Integer seed, or `NULL` (default) to draw from the ambient
#'   stream and leave it advanced; a seed draws from a private stream and
#'   restores the caller's on exit (see [detect_modules()]).
#'
#' @return An object of class `subspace_preservation`: `pairs` (data frame,
#'   one row per ordered pair: `species1`, `species2`, `S`, `S_null_mean`,
#'   `S_null_sd`, `z`, `K`, `n_genes1`, `n_genes2`, the last two the sizes
#'   of the largest components), `S` and `z` (species x species matrices,
#'   row = species1, `NA` on the diagonal) and `params`.
#'
#' @references
#' Kumar, A., Rai, P. & Daume III, H. (2011). Co-regularized multi-view
#' spectral clustering. \emph{Advances in Neural Information Processing
#' Systems} 24, 1413--1421.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200 * 12), 200)
#' x[1:40, ] <- x[1:40, ] + rnorm(12)
#' rownames(x) <- paste0("g", 1:200)
#' net <- compute_network(x, density = 0.1)
#' map <- data.frame(
#'   species = rep(c("A", "B"), each = 200),
#'   gene = rep(rownames(x), 2), hog = rep(rownames(x), 2)
#' )
#' sp <- subspace_preservation(list(A = net, B = net), map,
#'   K = 5L, n_null = 10L, seed = 1L
#' )
#' sp$S
#'
#' @seealso [as_preservation_matrix()], [recurrence_graph()]
#' @export
subspace_preservation <- function(nets, hog_map, K = 50L, # nolint
                                  n_null = 50L, n_cores = 1L,
                                  seed = NULL) {
  .seed_scope(seed)
  sp <- names(nets)
  if (!is.list(nets) || length(nets) < 2L || is.null(sp) ||
        any(!nzchar(sp)) || anyDuplicated(sp)) {
    stop("nets must be a list of at least two networks with unique names")
  }
  if (!is.data.frame(hog_map) ||
        !all(c("species", "gene", "hog") %in% names(hog_map))) {
    stop("hog_map must be a data frame with columns species, gene, hog")
  }
  ok_int <- function(v, lo) {
    is.numeric(v) && length(v) == 1L && !is.na(v) && v == round(v) &&
      v >= lo
  }
  if (!ok_int(K, 1)) stop("K must be a single positive whole number")
  if (!ok_int(n_null, 2)) stop("n_null must be a whole number >= 2")
  k <- as.integer(K)
  n_null <- as.integer(n_null)
  hog_map <- hog_map[!is.na(hog_map$hog), , drop = FALSE]
  hog_sp <- as.character(hog_map$species)
  hog_gene <- as.character(hog_map$gene)
  hogs <- unique(as.character(hog_map$hog))
  hog_id <- match(as.character(hog_map$hog), hogs)

  one <- function(s) .subspace_basis(nets[[s]], k)
  basis <- if (.can_fork(n_cores) && .blas_fork_safe()) {
    .check_fork_results(
      parallel::mclapply(sp, one, mc.cores = n_cores),
      sp, "species"
    )
  } else {
    lapply(sp, one)
  }
  names(basis) <- sp
  # gene x group incidence over the component genes of each species
  inc <- lapply(sp, function(s) {
    g <- rownames(basis[[s]])
    r <- hog_sp == s & hog_gene %in% g
    Matrix::sparseMatrix(
      i = match(hog_gene[r], g), j = hog_id[r], x = 1,
      dims = c(length(g), length(hogs))
    )
  })
  names(inc) <- sp

  ss <- function(ua, y) sum(crossprod(ua, y)^2) / k
  grid <- expand.grid(b = seq_along(sp), a = seq_along(sp))
  grid <- grid[grid$a != grid$b, c("a", "b")]
  rows <- lapply(seq_len(nrow(grid)), function(r) {
    a <- sp[grid$a[r]]
    b <- sp[grid$b[r]]
    pm <- inc[[a]] %*% Matrix::t(inc[[b]])
    rs <- Matrix::rowSums(pm)
    rs[rs > 0] <- 1 / rs[rs > 0]
    # P U_B with rows of P scaled to sum 1; permuting the rows of this is
    # permuting the group labels among A's genes
    y <- as.matrix(Matrix::Diagonal(x = rs) %*% (pm %*% basis[[b]]))
    na <- nrow(y)
    null <- vapply(seq_len(n_null), function(i) {
      ss(basis[[a]], y[sample.int(na), , drop = FALSE])
    }, numeric(1))
    s_obs <- ss(basis[[a]], y)
    data.frame(
      species1 = a, species2 = b, S = s_obs,
      S_null_mean = mean(null), S_null_sd = stats::sd(null),
      z = (s_obs - mean(null)) / stats::sd(null), K = k,
      n_genes1 = na, n_genes2 = nrow(basis[[b]]),
      stringsAsFactors = FALSE
    )
  })
  pairs <- do.call(rbind, rows)
  mat <- function(v) {
    m <- matrix(NA_real_, length(sp), length(sp), dimnames = list(sp, sp))
    m[cbind(match(pairs$species1, sp), match(pairs$species2, sp))] <- v
    m
  }
  structure(list(
    pairs = pairs, S = mat(pairs$S), z = mat(pairs$z),
    params = list(K = k, n_null = n_null, species = sp)
  ), class = "subspace_preservation")
}


#' Top-K non-trivial normalised-Laplacian eigenvectors of one network
#'
#' Binary adjacency at `net$threshold` through the sparse choke point,
#' largest connected component, then the K + 1 largest eigenvectors of
#' \eqn{D^{-1/2} W D^{-1/2}} (the K + 1 smallest of the Laplacian) with
#' the trivial one dropped. Rows are named by gene.
#' @noRd
.subspace_basis <- function(net, k) {
  net <- .net_as_sparse(net)
  a <- .net_cpp_args(net, net$threshold)
  g <- rownames(net$network)
  n <- length(g)
  col <- rep.int(seq_len(n), diff(a$p))
  keep <- a$x >= a$thr
  w <- Matrix::sparseMatrix(
    i = a$i[keep] + 1L, j = col[keep], x = 1, dims = c(n, n)
  )
  cc <- igraph::components(
    igraph::graph_from_adjacency_matrix(w, mode = "undirected")
  )
  lcc <- cc$membership == which.max(cc$csize)
  if (k + 1L > sum(lcc)) {
    stop(
      "K = ", k, " needs at least ", k + 1L, " genes in the largest ",
      "connected component; it has ", sum(lcc)
    )
  }
  w <- w[lcc, lcc]
  dm <- Matrix::Diagonal(x = 1 / sqrt(Matrix::colSums(w)))
  nw <- methods::as(dm %*% w %*% dm, "generalMatrix")
  e <- top_eigs_sym_cpp(nw@p, nw@i, nw@x, nrow(nw), k + 1L)
  u <- e$vectors[, -1L, drop = FALSE]
  rownames(u) <- g[lcc]
  u
}


#' @export
print.subspace_preservation <- function(x, ...) {
  p <- x$pairs
  cat(
    "Subspace preservation over", length(x$params$species), "species,",
    "K =", x$params$K, ", n_null =", x$params$n_null, "\n"
  )
  cat(
    "  S: median", signif(stats::median(p$S), 3), "range",
    paste(signif(range(p$S), 3), collapse = "-"), "\n"
  )
  cat(
    "  z: median", signif(stats::median(p$z), 3), "range",
    paste(signif(range(p$z), 3), collapse = " to "), "\n"
  )
  invisible(x)
}


#' Symmetrised subspace z in the shape preservation_matrix_test() reads
#'
#' [preservation_matrix_test()] reads `reference`, `test` and
#' `Zsummary_std` from a classification table. This returns one row per
#' unordered species pair with `Zsummary_std` set to the mean of the two
#' directed z of [subspace_preservation()], so the trait test runs on the
#' partition-free statistic unchanged; the unit is then the species pair,
#' not the module-direction.
#'
#' @param x A `subspace_preservation` object.
#' @return Data frame with columns `reference`, `test` and `Zsummary_std`.
#'
#' @examples
#' \dontrun{
#' sp <- subspace_preservation(nets, hog_map, K = 100L)
#' preservation_matrix_test(as_preservation_matrix(sp), group, block = genus)
#' }
#' @seealso [subspace_preservation()], [preservation_matrix_test()]
#' @export
as_preservation_matrix <- function(x) {
  if (!inherits(x, "subspace_preservation")) {
    stop("x must be a subspace_preservation object")
  }
  zs <- (x$z + t(x$z)) / 2
  idx <- which(upper.tri(zs), arr.ind = TRUE)
  sp <- rownames(zs)
  data.frame(
    reference = sp[idx[, 1L]], test = sp[idx[, 2L]], Zsummary_std = zs[idx],
    stringsAsFactors = FALSE
  )
}
