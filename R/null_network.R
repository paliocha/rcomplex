#' Shuffled-partner null network
#'
#' Rebuilds a network with the parameters recorded in `net` after
#' permuting every gene's expression values across samples independently,
#' which keeps each gene's marginal distribution and the gene universe
#' but destroys all co-expression. It is the partner network for the
#' empirical calibration of `method = "rank"`: p-values from a
#' comparison against it are the null draws that
#' [summarize_specificity()] turns into empirical p-values.
#'
#' With `block`, each gene's samples are permuted independently *within*
#' each block (one permutation per gene per block). The permutations
#' must differ between genes: one permutation applied to every gene only
#' reorders the samples and leaves the correlation matrix unchanged.
#'
#' **Which null.** The full shuffle (`block = NULL`) is the null of no
#' co-expression at all. The within-block shuffle keeps each gene's block
#' means, so the shared time course (or tree effect) survives, and is the
#' null of no co-expression beyond the shared axis. That is the null a
#' cross-species statistic on raw networks should be read against: on the
#' Pooideae time series, cross-species neighbour AUROC is 0.54 for
#' networks shuffled within time point against 0.50 for the full shuffle,
#' and that 0.04 is the shared time axis, not conserved wiring. The block
#' null leaves the random-geometric-graph geometry of the shared axis in
#' place, so within-species module statistics stay uninterpretable under
#' it. After regressing time out, one gene's residuals correlate -1/3
#' within a time point of four samples, which is a second reason to
#' permute within time point rather than across it.
#'
#' @param x Expression matrix (genes x samples, row names are gene
#'   identifiers) or a `SummarizedExperiment`, whose first assay is used.
#'   The same object `net` was built from.
#' @param net Network object from [compute_network()] built from `x`; its
#'   `params` and `store_density` are reused so the null is built the same
#'   way.
#' @param block `NULL` (default) for the full per-gene shuffle, or a
#'   per-sample grouping (vector or factor of length `ncol(x)`, no `NA`),
#'   e.g. time point or tree, within which each gene is permuted.
#' @param seed Integer seed for the permutation, or \code{NULL} (default)
#'   to draw from the global RNG. A seeded call draws from a private
#'   stream and restores the caller's on exit, so it does not displace the
#'   caller by however many draws the permutation consumed; with
#'   \code{seed = NULL} the draws come from the ambient stream and leave
#'   it advanced. Same contract as \code{\link{summarize_comparison}}.
#' @param n_cores,use_torch,block_size Passed to [compute_network()];
#'   pass the `block_size` `net` was built with to keep the null build
#'   blockwise (it is not recorded in `params`).
#' @return A sparse network object, as [compute_network()] returns, on the
#'   same genes as `net`. `params$log10_perms_per_gene` is the log10 size
#'   of one gene's permutation space (`lgamma(n + 1) / log(10)` for the
#'   full shuffle, the sum over blocks otherwise); with `block`,
#'   `params$block` holds the block sizes, named by level.
#' @export
null_network <- function(x, net, seed = NULL, n_cores = 1L,
                         use_torch = FALSE, block_size = NULL,
                         block = NULL) {
  .seed_scope(seed)
  if (methods::is(x, "SummarizedExperiment")) {
    x <- SummarizedExperiment::assay(x, 1L)
  }
  if (!is.null(block)) {
    if (length(block) != ncol(x)) {
      stop("block must have one entry per sample (ncol(x) = ", ncol(x), ")")
    }
    if (anyNA(block)) stop("block must not contain NA")
  }
  xp <- .shuffle_genes(x, block)
  tab <- if (is.null(block)) ncol(x) else table(block)
  p <- net$params
  out <- compute_network(xp,
    cor_method = p$cor_method, norm_method = p$norm_method,
    density = p$density, abs_cor = p$abs_cor,
    mr_log_transform = p$mr_log_transform, min_var = p$min_var,
    sparse = TRUE, store_density = net$store_density,
    n_cores = n_cores, use_torch = use_torch, block_size = block_size
  )
  if (!is.null(block)) {
    out$params$block <- stats::setNames(as.integer(tab), names(tab))
  }
  out$params$log10_perms_per_gene <- sum(lgamma(tab + 1)) / log(10)
  out
}


#' Permute each gene's samples independently (internal)
#'
#' Across all samples when `block` is `NULL`, otherwise within each block.
#' @noRd
.shuffle_genes <- function(x, block = NULL) {
  if (is.null(block)) {
    xp <- t(apply(x, 1L, sample))
  } else {
    xp <- x
    for (j in split(seq_len(ncol(x)), block)) {
      if (length(j) > 1L) {
        xp[, j] <- t(apply(x[, j, drop = FALSE], 1L, function(v) {
          v[sample.int(length(v))]
        }))
      }
    }
  }
  dimnames(xp) <- dimnames(x)
  xp
}


#' Validate and normalise a set of null networks (internal)
#'
#' `NULL` passes through. Otherwise `null_networks` must name every
#' species in `species`, each entry being a network object or a list of
#' them, built by [null_network()] on the same genes as the matching
#' entry of `networks`. The result is always `list(sp = list(net, ...))`.
#'
#' @param null_networks `NULL`, or a named list of network objects or
#'   lists of network objects.
#' @param networks Named list of the observed network objects.
#' @param species Character vector of species that need a null.
#' @return `NULL`, or the normalised list, one list of networks per
#'   species in `species`.
#' @noRd
.check_null_networks <- function(null_networks, networks, species) {
  if (is.null(null_networks)) {
    return(NULL)
  }
  missing <- setdiff(species, names(null_networks))
  if (length(missing) > 0L) {
    stop(
      "null_networks must name every species (missing: ",
      paste(missing, collapse = ", "), "); build them with null_network()"
    )
  }
  is_net <- function(n) {
    is.list(n) && all(c("network", "threshold") %in% names(n))
  }
  out <- lapply(species, function(sp) {
    nulls <- null_networks[[sp]]
    if (is_net(nulls)) {
      nulls <- list(nulls)
    }
    ok <- is.list(nulls) && length(nulls) > 0L &&
      all(vapply(nulls, is_net, logical(1)))
    if (!ok) {
      stop(
        "null_networks[['", sp, "']] must be a network object or a ",
        "list of them from null_network()"
      )
    }
    genes <- rownames(networks[[sp]]$network)
    for (n in nulls) {
      .net_check(n, n$threshold)
      if (!identical(rownames(n$network), genes)) {
        stop(
          "null network for ", sp, " is not on the genes of networks[['",
          sp, "']]; rebuild it with null_network() from the same ",
          "expression matrix"
        )
      }
    }
    nulls
  })
  names(out) <- species
  out
}
