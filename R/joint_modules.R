#' Find modules jointly across species
#'
#' Clusters the networks of all species and their hogs as one graph.
#' Labels are shared: the same label in two species is the same module.
#' The hog table shows which hogs keep their copies in one module.
#'
#' @details The graph is a star expansion. Each species adds its genes and
#'   the edges at or above its `threshold`. Edge weights are divided by
#'   that species' largest weight, so every species lies on the same
#'   scale. No edge crosses species. Each hog with two or more copies
#'   among the genes adds one hog node, linked to each copy with weight
#'   `weight`. Hog nodes carry no other edge. Leiden clusters this graph;
#'   the hog nodes are then dropped from the modules. `weight` is the
#'   coupling between the co-expression edges and the hogs. A small
#'   `weight` gives per-species modules; a large one pulls the copies of a
#'   hog together. The default, the median co-expression weight, is a
#'   guess until a sweep over `weight` tests it.
#'
#' @param nets Named list of at least two [compute_network()] results,
#'   dense or sparse. The names are the species.
#' @param orthologs Long table with columns `species`, `gene`, `hog`, or a
#'   file for [read_orthologs()].
#' @param weight Weight of the edge from a copy to its hog. `NULL` uses the
#'   median of all scaled co-expression edge weights.
#' @param objective Leiden objective: `"modularity"` or `"CPM"`.
#' @param resolution Leiden resolution. `NULL` uses 1 for modularity. For
#'   CPM it uses the edge density of the joint graph.
#' @param seed Integer seed, or `NULL` to use the current stream.
#' @param engine Clustering engine. Only `"star"` exists.
#'
#' @return A list with:
#'   \describe{
#'     \item{modules}{Named list with one module assignment per species, in
#'       the order of `nets`. Each has the [as_modules()] shape: `modules`
#'       (named integer vector, gene -> label), `module_genes`, `n_modules`,
#'       `method = "joint"` and `params`.}
#'     \item{hogs}{Data frame with one row per hog node, sorted by `hog`:
#'       `hog`, `n_copies` (genes linked), `n_modules` (distinct labels
#'       among the copies), `module_main` (the label with most copies; a tie
#'       goes to the smallest label) and `copies_main` (copies in it).}
#'     \item{params}{List of `weight` (the value used), `objective`,
#'       `resolution`, `engine` and `seed`.}
#'     \item{modularity}{Modularity of the joint graph, hog nodes included.}
#'   }
#'
#' @seealso [detect_modules()] for one species,
#'   [module_preservation()] to test a module in another species.
#' @examples
#' \dontrun{
#' res <- joint_modules(list(A = net_a, B = net_b), "orthologs.tsv",
#'   seed = 1
#' )
#' res$hogs[res$hogs$n_modules > 1, ]
#' pres <- module_preservation(res$modules$A, net_a, net_b, ortho)
#' }
#' @export
joint_modules <- function(nets, orthologs, weight = NULL,
                          objective = c("modularity", "CPM"),
                          resolution = NULL, seed = NULL, engine = "star") {
  objective <- match.arg(objective)
  if (!isTRUE(engine == "star")) {
    stop("`engine` must be \"star\"; got ", deparse(engine))
  }
  sp <- names(nets)
  bad <- !is.list(nets) || length(sp) < 2L || anyDuplicated(sp) > 0L
  if (bad || any(sp == "")) {
    stop("`nets` must be a named list of at least two species.")
  }
  ok <- is.numeric(weight) && length(weight) == 1L && isTRUE(weight > 0)
  if (!is.null(weight) && !ok) {
    stop("`weight` must be one positive number; got ", deparse(weight))
  }
  long <- .driver_orthologs(orthologs, sp)

  # Gene nodes 1..n_gene in species order, hog nodes after them
  genes <- lapply(nets, function(n) rownames(n$network))
  off <- stats::setNames(cumsum(c(0L, lengths(genes)))[seq_along(sp)], sp)
  n_gene <- sum(lengths(genes))
  edges <- do.call(rbind, lapply(sp, function(s) {
    e <- .joint_edges(nets[[s]], s)
    e$w <- e$w / max(e$w)
    e$i <- e$i + off[[s]]
    e$j <- e$j + off[[s]]
    e
  }))

  node <- rep(NA_integer_, nrow(long))
  for (s in sp) {
    r <- long$species == s
    node[r] <- off[[s]] + match(long$gene[r], genes[[s]])
  }
  link <- unique(data.frame(node = node, hog = long$hog)[!is.na(node), ])
  hogs <- sort(unique(link$hog[duplicated(link$hog)]), method = "radix")
  link <- link[link$hog %in% hogs, , drop = FALSE]

  if (is.null(weight)) weight <- stats::median(edges$w)
  edges <- rbind(edges, data.frame(
    i = link$node, j = n_gene + match(link$hog, hogs), w = weight
  ))
  g <- igraph::make_graph(rbind(edges$i, edges$j),
    n = n_gene + length(hogs), directed = FALSE
  )
  igraph::E(g)$weight <- edges$w
  if (is.null(resolution)) {
    resolution <- if (objective == "CPM") {
      igraph::ecount(g) / choose(igraph::vcount(g), 2)
    } else {
      1
    }
  }

  .seed_scope(seed)
  comm <- igraph::cluster_leiden(g,
    resolution = resolution, objective_function = objective,
    weights = .leiden_weights(g, objective), n_iterations = 2L
  )
  memb <- as.integer(igraph::membership(comm))
  lab <- memb[seq_len(n_gene)]
  lab <- match(lab, unique(lab))

  params <- list(
    weight = weight, objective = objective, resolution = resolution,
    engine = engine, seed = seed
  )
  modules <- lapply(stats::setNames(sp, sp), function(s) {
    m <- stats::setNames(lab[off[[s]] + seq_along(genes[[s]])], genes[[s]])
    module_genes <- split(names(m), m)
    list(
      modules = m, module_genes = module_genes,
      n_modules = length(module_genes), method = "joint", params = params
    )
  })

  by_hog <- split(lab[link$node], factor(link$hog, levels = hogs))
  main <- vapply(by_hog, function(x) {
    u <- sort(unique(x))
    n <- tabulate(match(x, u))
    c(u[which.max(n)], max(n))
  }, integer(2))
  list(
    modules = modules,
    hogs = data.frame(
      hog = hogs, n_copies = lengths(by_hog, use.names = FALSE),
      n_modules = vapply(by_hog, function(x) length(unique(x)), 1L,
        USE.NAMES = FALSE
      ),
      module_main = main[1L, ], copies_main = main[2L, ],
      row.names = NULL, stringsAsFactors = FALSE
    ),
    params = params,
    modularity = igraph::modularity(g, memb)
  )
}


#' Upper-triangle edges of one species at or above its threshold
#' @noRd
.joint_edges <- function(net, s) {
  m <- .net_check(net, net$threshold)
  thr <- net$threshold
  if (.net_is_sparse(net)) {
    e <- Matrix::summary(m)
    e <- e[e$i < e$j & e$x >= thr & e$x > 0, , drop = FALSE]
    e <- data.frame(i = e$i, j = e$j, w = e$x)
  } else {
    k <- which(upper.tri(m) & m >= thr & m > 0, arr.ind = TRUE)
    e <- data.frame(i = k[, 1L], j = k[, 2L], w = m[k])
  }
  if (nrow(e) == 0L) {
    stop(
      "`nets$", s, "` has no edge above its threshold; ",
      "there is nothing to cluster."
    )
  }
  e
}
