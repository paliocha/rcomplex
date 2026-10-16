#' Compare co-expression networks across species in one call
#'
#' `rcomplex()` runs the whole analysis. It builds one network per
#' species, tests every ortholog pair in every species pair, and
#' classifies the gene cliques. Give it expression or networks, plus
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
#' @param null If `TRUE`, repeat the run with hog labels shuffled among
#'   genes of the same expression level. The null calls and tiers estimate
#'   what abundance and network structure alone produce.
#' @param n_cores Number of threads.
#' @param block_size `NULL` (default) builds each species network from the
#'   dense gene x gene matrix. A positive whole number builds it in blocks
#'   of that many genes. The full matrix then never forms.
#' @param r_min `NULL` (default), a number, or `"null"`: the correlation
#'   floor of each species network. See [compute_network()].
#' @param seed Integer seed, or `NULL` to draw from the global stream.
#'   A seed makes the run reproducible and leaves your stream unchanged.
#'
#' @return A list of class `rcomplex`:
#'   \describe{
#'     \item{networks}{The network objects, by species.}
#'     \item{edges}{The edge table from [find_coexpressologs()], plus a
#'       `sign` column: the sign of every network.}
#'     \item{cliques, members, twins}{The three tables of
#'       [gene_clique_graph()]: one row per clique, one row per clique
#'       member, and one row per collapsed twin gene. The tables of the
#'       three graphs are stacked. `run` numbers the graph: 1 strict, 2
#'       loose, 3 all tested pairs. `alpha_graph` holds its cut: `alpha`,
#'       0.9 or `Inf`. `clique_id` counts within a run, so `run` and
#'       `clique_id` together name a clique.}
#'     \item{capped}{One row per hog and run over the clique budget of
#'       [gene_clique_graph()], with `run` and `alpha_graph`. Its cliques
#'       are species blocks. Usually empty.}
#'     \item{classification}{One row per clique, from
#'       [classify_gene_cliques()]. Its `run` and `clique_id` join
#'       `cliques` and `members`.}
#'     \item{modules, preservation}{With `modules = TRUE`: the
#'       [detect_modules()] results by species, and the
#'       [preservation_paired()] result.}
#'     \item{edges_null, classification_null}{With `null = TRUE`: the
#'       edge table and the classification of the shuffled ortholog
#'       table.}
#'     \item{call}{The call.}
#'   }
#'   `as.data.frame()` returns the classification.
#'
#' @details
#' The driver lists cliques on three gene graphs. The strict graph holds
#' the edges with `q_value < alpha`. It feeds `complete_conserved`,
#' `lineage_specific`, `partial_present` and `trait_specific`. The loose
#' graph holds the edges with `q_value < 0.9`. It feeds
#' `partial_significant`. The third graph holds every tested pair. It
#' covers only the hogs with no complete or partial tier on the first two
#' graphs. It feeds `differentiated`, as does the loose graph where every
#' cross pair has `q_value < 0.9`. [classify_gene_cliques()] classifies
#' the three graphs together, so `hog_class` reads all of them. A clique
#' maximal on two graphs appears once for each. The third graph is built
#' only with `clades`, since only the clade tiers read it.
#'
#' Every clique has at least `min(3, number of species)` genes. Rank
#' edges are classified at `min_power = 0.9`, hypergeometric edges at 0.8.
#'
#' `density` and `sign` apply only to networks built from `expr`.
#' `method = "rank"` and `block` need `expr`. The rank test calibrates
#' against [null_network()] partners, built within each block.
#'
#' `block_size` passes to [compute_network()] and needs `expr`. The
#' blockwise network equals the dense one except at floating-point
#' near-ties at the density cut. 512 is a good value. A species with more
#' than 40,000 genes gets a message that suggests it.
#'
#' `summary()` counts each tier in `cliques`, in `hogs` and in
#' `patterns`. A pattern is one HOG with one clique species set. The
#' paralog copies of a large HOG can form millions of cliques, but they
#' share a few patterns. Paralogs with the same partners are twins: one
#' clique stands for all their combinations. The `multiplicity` column
#' counts those combinations, and `twins` lists the collapsed genes.
#'
#' With `null = TRUE`, the driver shuffles the `gene` column of
#' `orthologs` within each species. A gene moves only among genes of
#' the same expression decile, by `gene_mean` from [compute_network()].
#' A network without `gene_mean`, such as one from [as_network()], forms
#' one decile. The driver then reruns edges, cliques and tiers on the
#' real networks. Both runs share the rank calibration networks.
#'
#' With `null = TRUE`, `summary()` reports `calls` and `calls_null` per
#' species pair at `q_value < alpha`. Their ratio is the empirical
#' false-call rate.
#' The rate is `NA` for a species pair without calls. The tier table
#' gains `cliques_null`: the null cliques per tier.
#'
#' `print()` counts the cliques of each graph. It shows, per species, the
#' genes, the samples, the density and `r_threshold` from
#' [compute_network()]: the weakest correlation that
#' passed the density threshold. The same density needs a stronger
#' correlation at few samples. A network from [as_network()] has no
#' correlation and shows \verb{r ?}, as does a network without edges.
#' With `r_min`, the line also shows the floor and the edges it removed.
#'
#' @examples
#' f <- function(x) system.file("extdata", x, package = "rcomplex")
#' read <- function(x) as.matrix(read.delim(f(x), row.names = 1L))
#' expr <- lapply(c(SpA = "expr_sp1_small.txt", SpB = "expr_sp2_small.txt"),
#'   read)
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
                     block_size = NULL, r_min = NULL, seed = NULL) {
  sign <- match.arg(sign)
  method <- match.arg(method)
  if (is.null(expr) == is.null(networks)) {
    stop("`expr` and `networks`: give exactly one, not both or neither.")
  }
  input <- if (is.null(expr)) networks else expr
  sp <- names(input)
  bad <- !is.list(input) || length(sp) < 2L || anyDuplicated(sp) > 0L
  if (bad || any(sp == "")) {
    stop("`expr` or `networks` must be a named list of at least two species.")
  }
  if (is.null(expr) && (method == "rank" || !is.null(block))) {
    stop("method = \"rank\" and block need expr; you gave networks")
  }
  if (is.null(expr) && !is.null(block_size)) {
    stop("block_size needs expr; you gave networks")
  }
  if (is.null(expr) && !is.null(r_min)) {
    stop("r_min needs expr; you gave networks")
  }
  long <- .driver_orthologs(orthologs, sp)
  ortho <- prepare_orthologs(long)
  .seed_scope(seed)

  xs <- NULL
  if (!is.null(expr)) {
    .driver_size_note(expr, block_size)
    xs <- if (is.null(block)) expr else .driver_wiring(expr, block)
    networks <- lapply(xs, compute_network,
      density = density, sign = sign, n_cores = n_cores,
      block_size = block_size, r_min = r_min
    )
    for (s in sp) {
      x <- expr[[s]]
      if (methods::is(x, "SummarizedExperiment")) {
        x <- SummarizedExperiment::assay(x, 1L)
      }
      networks[[s]]$params$n_samples <- ncol(x)
      # with block the network is built from the centred wiring layer
      networks[[s]]$gene_mean <- rowMeans(x)[rownames(networks[[s]]$network)]
    }
  } else {
    signs <- vapply(networks, function(n) n$params$sign %||% "positive", "")
    if (length(unique(signs)) > 1L) {
      stop("networks differ in sign: ", paste(sp, signs, collapse = ", "))
    }
    sign <- signs[[1L]]
  }
  nulls <- if (method == "rank") .driver_nulls(xs, networks, block, n_cores)
  edges <- .driver_edges(networks, ortho, method, nulls, n_cores, sign)
  cliques <- function(e) {
    .driver_cliques(e, sp, clades, alpha,
      min_power = if (method == "rank") 0.9 else 0.8
    )
  }
  g <- cliques(edges)
  res <- list(
    networks = networks, edges = edges, cliques = g$cliques,
    members = g$members, twins = g$twins, capped = g$capped,
    classification = g$classification
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
    ortho0 <- prepare_orthologs(.driver_shuffle(long, networks))
    res$edges_null <- .driver_edges(networks, ortho0, method, nulls, n_cores,
      sign
    )
    res$classification_null <- cliques(res$edges_null)$classification
  }
  res$call <- match.call()
  structure(res, class = "rcomplex")
}


#' Cliques of three gene graphs and their classification
#'
#' The strict graph holds edges with `q_value < alpha`, the loose graph
#' edges with `q_value < 0.9`. With `clades`, the third graph holds every
#' tested pair, but only of hogs that the first two leave without a
#' complete or partial tier. The tables of the runs are stacked with
#' `run` and `alpha_graph` columns.
#' @noRd
.driver_cliques <- function(e, sp, clades, alpha, min_power) {
  k <- min(3L, length(sp))
  classify <- function(runs) {
    classify_gene_cliques(runs, e, sp,
      clades = clades, alpha_call = alpha, min_power = min_power
    )
  }
  runs <- list(
    gene_clique_graph(e, min_size = k, alpha_graph = alpha),
    gene_clique_graph(e, min_size = k, alpha_graph = 0.9)
  )
  # Without clades no tier reads the all-tested graph.
  if (!is.null(clades)) {
    first <- classify(runs)
    claim <- c("complete_conserved", "partial_significant", "partial_present")
    claimed <- first$hog[first$hog_class %in% claim]
    rm(first)
    runs[[3L]] <- gene_clique_graph(e[!e$hog %in% claimed, , drop = FALSE],
      min_size = k, alpha_graph = Inf
    )
  }
  a_graph <- c(alpha, 0.9, Inf)
  stack <- function(f) {
    do.call(rbind, lapply(seq_along(runs), function(i) {
      t <- runs[[i]][[f]]
      cbind(
        run = rep(i, nrow(t)), alpha_graph = rep(a_graph[i], nrow(t)), t
      )
    }))
  }
  out <- list(classification = classify(runs))
  # Each table leaves the runs as it is stacked, to bound the peak memory.
  for (f in c("cliques", "members", "twins", "capped")) {
    out[[f]] <- stack(f)
    for (i in seq_along(runs)) runs[[i]][[f]] <- NULL
  }
  attr(out$cliques, "alpha_graph") <- a_graph[seq_along(runs)]
  out
}


#' Long ortholog table of the driver's species, sorted in their order
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
  long[order(match(long$species, sp)), , drop = FALSE]
}


#' Shuffle each species' genes among rows of the same expression decile
#'
#' Rows whose gene is in the network swap genes within deciles of
#' `gene_mean`. Hog sizes per species stay; orthology breaks.
#' @noRd
.driver_shuffle <- function(long, networks) {
  for (s in names(networks)) {
    net <- networks[[s]]
    i <- which(long$species == s & long$gene %in% rownames(net$network))
    strata <- 1L
    if (is.null(net$gene_mean)) {
      message(
        "null: ", s, " has no gene_mean, so the shuffle ignores expression"
      )
    } else {
      m <- net$gene_mean[long$gene[i]]
      br <- unique(stats::quantile(m, 0:10 / 10, names = FALSE))
      if (length(br) > 1L) strata <- cut(m, br, include.lowest = TRUE)
    }
    for (k in split(i, strata)) {
      long$gene[k] <- long$gene[k][sample.int(length(k))]
    }
  }
  long
}


#' One message per species too large for the dense build
#' @noRd
.driver_size_note <- function(expr_list, block_size) {
  if (!is.null(block_size)) return(invisible())
  for (s in names(expr_list)) {
    n <- nrow(expr_list[[s]])
    if (n > 40000L) {
      message(
        s, ": ", n, " genes need a ", round(8 * n^2 / 2^30, 1),
        " GB dense matrix; block_size = 512 builds it in pieces"
      )
    }
  }
  invisible()
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


#' Shuffled-expression null network per species, for the rank test
#' @noRd
.driver_nulls <- function(xs, networks, block, n_cores) {
  sapply(names(networks), function(s) {
    null_network(xs[[s]], networks[[s]], n_cores = n_cores, block = block[[s]])
  }, simplify = FALSE)
}


#' find_coexpressologs() with the rank nulls, plus `sign`
#' @noRd
.driver_edges <- function(networks, ortho, method, nulls, n_cores, sign) {
  e <- find_coexpressologs(networks, ortho,
    method = method, n_cores = n_cores, null_networks = nulls
  )
  e$sign <- rep(sign, nrow(e))
  e
}


#' The weakest edge correlation of a network, as print() shows it
#' @noRd
.print_r_threshold <- function(n) {
  r <- n$params$r_threshold
  if (is.null(r) || is.na(r)) {
    return("r ?")
  }
  neg <- identical(n$params$sign, "negative")
  op <- if (neg) "<=" else ">="
  out <- sprintf("r %s %.2f", op, r)
  n_fl <- n$params$n_floored
  if (!is.null(n_fl) && n_fl > 0L) {
    out <- sprintf("%s  floor r %s %.2f: %s edges removed",
      out, op, if (neg) -n$params$r_min else n$params$r_min,
      format(n_fl, big.mark = ",")
    )
  }
  out
}


#' @export
print.rcomplex <- function(x, ...) {
  cls <- x$classification
  tiers <- table(cls$classification)
  sp <- names(x$networks)
  big <- function(n) format(n, big.mark = ",")
  cat(
    "rcomplex: ", length(sp), " species, ", big(length(unique(x$edges$hog))),
    " hogs, ", length(tiers), ngettext(length(tiers), " tier", " tiers"),
    "\n",
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
    "  ", vapply(x$networks, .print_r_threshold, character(1)), "\n"
  ), sep = "")
  alpha <- attr(cls, "alpha_call")
  cat(
    "  edges ", big(nrow(x$edges)), " tested, ",
    big(sum(x$edges$q_value < alpha, na.rm = TRUE)), " called at q < ",
    alpha, "\n",
    sep = ""
  )
  ag <- attr(x$cliques, "alpha_graph") %||%
    c(alpha, 0.9, Inf)[seq_len(max(2L, x$cliques$run))]
  n_run <- format(tabulate(x$cliques$run, length(ag)),
    big.mark = ",", trim = TRUE
  )
  lab <- ifelse(is.finite(ag), paste("at q <", ag), "on all tested pairs")
  lab[1L] <- paste("cliques", lab[1L])
  cat("  graphs: ", paste(n_run, lab, collapse = ", "), "\n", sep = "")
  cat("  cliques ", big(nrow(cls)), sep = "")
  n_all <- sum(as.numeric(cls$multiplicity))
  if (n_all > nrow(cls)) {
    cat(" distinct (", big(n_all), " with paralog copies)", sep = "")
  }
  if (length(tiers) > 0L) {
    cat(":", paste(names(tiers), big(as.vector(tiers)), collapse = ", "))
  }
  n_cap <- length(unique(x$capped$hog))
  if (n_cap > 0L) {
    cat("; ", big(n_cap), " hogs over the clique budget, reported as blocks",
      sep = ""
    )
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.rcomplex <- function(object, ...) {
  cls <- object$classification
  alpha <- attr(cls, "alpha_call")
  cls0 <- object$classification_null
  tier <- factor(cls$classification,
    levels = sort(unique(as.character(
      c(cls$classification, cls0$classification)
    )))
  )
  # A species pattern is a HOG with one clique species set.
  pattern <- collapse::group(list(cls$hog, cls$missing_species))
  n_distinct <- function(x) {
    unname(vapply(split(x, tier), function(v) length(unique(v)), 1L))
  }
  tiers <- data.frame(
    tier = levels(tier), cliques = as.vector(table(tier)),
    cliques_null = if (is.null(cls0)) {
      rep(NA_integer_, nlevels(tier))
    } else {
      as.vector(table(factor(cls0$classification, levels(tier))))
    },
    multiplicity = unname(vapply(
      split(as.numeric(cls$multiplicity), tier), sum, 1
    )),
    hogs = n_distinct(cls$hog), patterns = n_distinct(pattern)
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
#' `write_rcomplex()` writes each table of the result to its own file.
#' Each file takes the name of its element: `edges.tsv`, `cliques.tsv`,
#' `members.tsv`, `twins.tsv`, `capped.tsv`, `classification.tsv` and,
#' with `null = TRUE`, `edges_null.tsv` and `classification_null.tsv`.
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
