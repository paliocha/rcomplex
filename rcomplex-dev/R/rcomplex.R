#' Compare co-expression networks across species in one call
#'
#' `rcomplex()` runs the whole analysis. It builds one network per
#' species, tests every ortholog pair in every species pair, and
#' classifies the gene cliques. Give it expression or networks, and
#' orthologs.
#'
#' @param expr Named list of expression data, one entry per species: a
#'   genes x samples matrix or a `SummarizedExperiment`. Give `expr` or
#'   `networks`, not both.
#' @param orthologs The long ortholog table (`species`, `gene`, `hog`),
#'   or a path to an ortholog file for [read_orthologs()].
#' @param networks Named list of network objects from [compute_network()]
#'   or [as_network()], one entry per species.
#' @param block Named list of per-sample factors, one per species in
#'   `expr`. Each species is split with [split_layers()]. The network is
#'   built from the wiring layer.
#' @param clades Named list of species vectors, one per clade, as for
#'   [classify_gene_cliques()]. `NULL` skips the clade tiers.
#' @param density Fraction of gene pairs kept as edges in each network.
#' @param sign `"positive"` or `"negative"`: the correlation sign that
#'   makes an edge. See [compute_network()].
#' @param method `"hypergeometric"` or `"rank"`: the test of each
#'   ortholog pair. See [find_coexpressologs()].
#' @param alpha The q-value threshold for a call.
#' @param modules If `TRUE`, detect modules in each species and test
#'   their preservation in every species pair.
#' @param null If `TRUE`, repeat the comparison on shuffled networks.
#'   The null calls estimate the false-call rate.
#' @param n_cores Number of threads.
#' @param seed Integer seed, or `NULL` to draw from the global stream.
#'   A seed makes the run reproducible and leaves your stream unchanged.
#'
#' @return A list of class `rcomplex`:
#'   \describe{
#'     \item{networks}{The network objects, by species.}
#'     \item{edges}{The edge table from [find_coexpressologs()].}
#'     \item{cliques}{The cliques from [gene_clique_graph()].}
#'     \item{classification}{One row per clique, from
#'       [classify_gene_cliques()].}
#'     \item{modules, preservation}{With `modules = TRUE`: the
#'       [detect_modules()] results by species, and the
#'       [preservation_paired()] result.}
#'     \item{edges_null}{With `null = TRUE`: the edge table of the
#'       shuffled networks.}
#'     \item{call}{The call.}
#'   }
#'   `print()` shows the species, edges and tier counts. `summary()`
#'   prints the tier table and the null table. `as.data.frame()` returns
#'   the classification. [write_rcomplex()] writes the tables.
#'
#' @details
#' The steps are [compute_network()] per species, [find_coexpressologs()]
#' over all species pairs, [gene_clique_graph()] at `q_value < alpha`,
#' and [classify_gene_cliques()] with `alpha_call = alpha`. Rank edges
#' are classified at `min_power = 0.9`, hypergeometric edges at 0.8.
#' With two species a clique is one edge, so cliques need at least
#' `min(3, number of species)` genes.
#'
#' `density` and `sign` apply only to networks built from `expr`.
#' `method = "rank"`, `null = TRUE` and `block` need `expr`, because
#' the shuffled networks are built from expression with
#' [null_network()]. With `block`, the shuffle stays within each block.
#'
#' With `null = TRUE`, `summary()` reports `calls` and `calls_null` per
#' species pair at `q_value < alpha`. Their ratio is the empirical
#' false-call rate.
#'
#' @examples
#' f <- function(x) system.file("extdata", x, package = "rcomplex")
#' read <- function(x) as.matrix(read.delim(f(x), row.names = 1L))
#' expr <- list(
#'   SpA = read("expr_sp1_small.txt"),
#'   SpB = read("expr_sp2_small.txt")
#' )
#' res <- rcomplex(expr, f("orthologs_small.txt"), density = 0.1, seed = 1)
#' res
#' summary(res)
#' write_rcomplex(res, tempfile())
#' @export
rcomplex <- function(expr = NULL, orthologs, networks = NULL, block = NULL,
                     clades = NULL, density = 0.03,
                     sign = c("positive", "negative"),
                     method = c("hypergeometric", "rank"), alpha = 0.1,
                     modules = FALSE, null = FALSE, n_cores = 1L,
                     seed = NULL) {
  cl <- match.call()
  sign <- match.arg(sign)
  method <- match.arg(method)
  if (is.null(expr) == is.null(networks)) {
    stop("give exactly one of expr and networks")
  }
  input <- if (is.null(expr)) networks else expr
  sp <- names(input)
  bad <- !is.list(input) || length(sp) < 2L || anyDuplicated(sp) > 0L
  if (bad || any(sp == "")) {
    stop("expr or networks must be a named list of at least two species")
  }
  if (is.null(expr) && (null || method == "rank" || !is.null(block))) {
    stop("null = TRUE, method = \"rank\" and block each needs expr")
  }
  ortho <- .driver_orthologs(orthologs, sp)
  .seed_scope(seed)

  xs <- NULL
  if (!is.null(expr)) {
    xs <- if (is.null(block)) expr else .driver_wiring(expr, block)
    networks <- lapply(xs, compute_network,
      density = density, sign = sign, n_cores = n_cores
    )
    for (s in sp) networks[[s]]$params$n_samples <- ncol(expr[[s]])
  }
  edges <- .driver_edges(networks, xs, ortho, method, block, n_cores)
  cliques <- gene_clique_graph(edges,
    min_size = min(3L, length(sp)), alpha_graph = alpha
  )
  res <- list(
    networks = networks, edges = edges, cliques = cliques,
    classification = classify_gene_cliques(cliques, edges, sp,
      clades = clades, alpha_call = alpha,
      min_power = if (method == "rank") 0.9 else 0.8
    )
  )
  if (modules) {
    res$modules <- lapply(networks, detect_modules, n_cores = n_cores)
    pairs <- t(utils::combn(sp, 2L))
    res$preservation <- preservation_paired(res$modules, networks, ortho,
      data.frame(species1 = pairs[, 1L], species2 = pairs[, 2L]),
      edges = edges
    )
  }
  if (null) {
    nulls <- .driver_nulls(xs, networks, block, n_cores)
    res$edges_null <- .driver_edges(nulls, xs, ortho, method, block, n_cores)
  }
  res$call <- cl
  structure(res, class = "rcomplex")
}


#' Pairwise ortholog table for the driver's species, in their order
#' @noRd
.driver_orthologs <- function(orthologs, sp) {
  long <- if (is.character(orthologs)) read_orthologs(orthologs) else orthologs
  cols <- c("species", "gene", "hog")
  if (!is.data.frame(long) || !all(cols %in% names(long))) {
    stop("orthologs must be a file or a data frame with species, gene, hog")
  }
  miss <- setdiff(sp, long$species)
  if (length(miss) > 0L) {
    stop("species not in orthologs: ", paste(miss, collapse = ", "))
  }
  long <- long[long$species %in% sp, , drop = FALSE]
  # gene1 must belong to the first species of each pair
  prepare_orthologs(long[order(match(long$species, sp)), , drop = FALSE])
}


#' Wiring layer per species, with a message for missing block levels
#' @noRd
.driver_wiring <- function(expr, block) {
  if (!is.list(block) || !setequal(names(block), names(expr))) {
    stop("block must be a named list with one factor per species in expr")
  }
  lv <- lapply(block, function(b) unique(as.character(b)))
  all_lv <- unique(unlist(lv))
  for (s in names(expr)) {
    miss <- setdiff(all_lv, lv[[s]])
    if (length(miss) > 0L) {
      message(
        s, " lacks block ", ngettext(length(miss), "level", "levels"),
        ": ", paste(miss, collapse = ", ")
      )
    }
  }
  sapply(names(expr), function(s) split_layers(expr[[s]], block[[s]])$wiring,
    simplify = FALSE
  )
}


#' Shuffled-expression null network per species
#' @noRd
.driver_nulls <- function(xs, networks, block, n_cores) {
  sapply(names(networks), function(s) {
    null_network(xs[[s]], networks[[s]], n_cores = n_cores, block = block[[s]])
  }, simplify = FALSE)
}


#' find_coexpressologs() with the rank test's null networks when needed
#' @noRd
.driver_edges <- function(networks, xs, ortho, method, block, n_cores) {
  nulls <- if (method == "rank") .driver_nulls(xs, networks, block, n_cores)
  find_coexpressologs(networks, ortho,
    method = method, n_cores = n_cores, null_networks = nulls
  )
}


#' @export
print.rcomplex <- function(x, ...) {
  cls <- x$classification
  tiers <- table(cls$classification)
  sp <- names(x$networks)
  big <- function(n) format(n, big.mark = ",")
  cat(
    "rcomplex: ", length(sp), " species, ", big(length(unique(x$edges$hog))),
    " hogs, ", length(tiers), " tiers\n",
    sep = ""
  )
  genes <- vapply(x$networks, function(n) big(n$n_genes), character(1))
  samples <- vapply(x$networks, function(n) {
    if (is.null(n$params$n_samples)) "?" else big(n$params$n_samples)
  }, character(1))
  cat(paste0(
    "  ", format(sp), "  ", format(genes, justify = "right"), " genes  ",
    format(samples, justify = "right"), " samples  density ",
    vapply(x$networks, function(n) format(n$params$density), character(1)),
    "\n"
  ), sep = "")
  alpha <- attr(cls, "alpha_call")
  cat(
    "  edges ", big(nrow(x$edges)), " tested, ",
    big(sum(x$edges$q_value < alpha, na.rm = TRUE)), " called at q < ",
    alpha, "\n",
    sep = ""
  )
  cat("  cliques ", big(nrow(cls)), sep = "")
  if (length(tiers) > 0L) {
    cat(":", paste(names(tiers), big(as.vector(tiers)), collapse = ", "))
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.rcomplex <- function(object, ...) {
  cls <- object$classification
  alpha <- attr(cls, "alpha_call")
  tiers <- as.data.frame(table(tier = cls$classification),
    responseName = "cliques", stringsAsFactors = FALSE
  )
  null_tab <- NULL
  if (!is.null(object$edges_null)) {
    pr <- t(utils::combn(names(object$networks), 2L))
    calls <- function(e) {
      hit <- !is.na(e$q_value) & e$q_value < alpha
      key <- paste(e$species1, e$species2)
      unname(vapply(paste(pr[, 1L], pr[, 2L]), function(k) {
        sum(hit & key == k)
      }, 1L))
    }
    n <- calls(object$edges)
    n0 <- calls(object$edges_null)
    null_tab <- data.frame(
      species1 = pr[, 1L], species2 = pr[, 2L], calls = n, calls_null = n0,
      false_call_rate = ifelse(n > 0L, n0 / n, NA_real_)
    )
  }
  print(tiers, row.names = FALSE)
  if (is.null(null_tab)) {
    cat("null: not run\n")
  } else {
    cat("null: calls at q <", alpha, "\n")
    print(null_tab, row.names = FALSE)
  }
  invisible(list(tiers = tiers, null = null_tab))
}


#' @export
as.data.frame.rcomplex <- function(x, ...) x$classification


#' Write an rcomplex result as tab-separated files
#'
#' `write_rcomplex()` writes each table of the result to its own file,
#' named after the element: `edges.tsv`, `cliques.tsv`,
#' `classification.tsv` and, with `null = TRUE`, `edges_null.tsv`.
#'
#' @param x A result of [rcomplex()].
#' @param dir Output directory. It is created when it does not exist.
#' @return The file paths, invisibly.
#' @seealso The example of [rcomplex()].
#' @export
write_rcomplex <- function(x, dir) {
  if (!inherits(x, "rcomplex")) stop("x must be a result of rcomplex()")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tabs <- Filter(is.data.frame, unclass(x))
  paths <- file.path(dir, paste0(names(tabs), ".tsv"))
  for (i in seq_along(tabs)) {
    utils::write.table(tabs[[i]], paths[[i]],
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  invisible(paths)
}
