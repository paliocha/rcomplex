#' Line up per-species modules across species
#'
#' Lines up the modules of each species into shared labels. The same label
#' in two species is the same module. Modules come from [detect_modules()],
#' one species at a time.
#'
#' @details Each module of each species is one node. Louvain joins these
#'   nodes into shared labels. It never moves a single gene, so labels are
#'   species modules lined up, never genes clustered again. A shared label
#'   may hold several modules of one species. It may hold modules of one
#'   species only, when nothing lines up.
#'
#'   The objective has two terms. Within a species, it is modularity with
#'   that species' own null. Across species, it is the hog overlap of two
#'   modules against the overlap expected for that species pair. A hog with
#'   `c_s` copies in species s and `c_t` in species t gives each copy pair
#'   the weight `1 / (c_s c_t)`. So hogs with many copies count less. The
#'   edges are those at or above each species' `threshold`. Each species'
#'   edge weights are divided by its largest weight.
#'
#'   On wood and leaf, every `lambda` from 0.25 to 16 gave the same labels.
#'   For why family-steered clustering of genes was withdrawn, see
#'   `dev/design-notes/hypergraph-community-literature.md`, sections 10-12.
#'
#' @param nets Named list of at least two [compute_network()] results,
#'   dense or sparse. The names are the species.
#' @param orthologs Long table with columns `species`, `gene`, `hog`, or a
#'   file for [read_orthologs()].
#' @param modules Named list with one [detect_modules()] or [as_modules()]
#'   result per species of `nets`, in any order. Each must label every gene
#'   of its network. `NULL` runs [detect_modules()] on each species.
#' @param lambda Weight of the hog term against the co-expression term.
#'   One number of 0 or more. At 0 no label spans two species.
#' @param seed Integer seed, or `NULL` to use the current stream.
#'
#' @return A list with:
#'   \describe{
#'     \item{modules}{Named list with one module assignment per species, in
#'       the order of `nets`. Each has the [as_modules()] shape: `modules`
#'       (named integer vector, gene -> shared label), `module_genes`,
#'       `n_modules`, `method = "joint"` and `params`.}
#'     \item{members}{Data frame with one row per species module: `species`,
#'       `module` (the label in `modules`), `label` (the shared label) and
#'       `n_genes`. Rows are sorted by `label`, species and `module`.}
#'     \item{hogs}{Data frame with one row per hog with two or more copies
#'       among the network genes, sorted by `hog`: `hog`, `n_copies`,
#'       `n_modules` (distinct shared labels among the copies),
#'       `module_main` (the label with most copies; a tie goes to the
#'       smallest label) and `copies_main` (copies in it).}
#'     \item{params}{List of `lambda`, `gamma` (1), `seed`, `modules_given`
#'       and `quality` (the objective of the final labels).}
#'   }
#'
#' @seealso [detect_modules()] for one species,
#'   [module_correspondence()] for pairs of species with p-values,
#'   [module_preservation()] to test a module in another species.
#' @examples
#' \dontrun{
#' res <- joint_modules(list(A = net_a, B = net_b), "orthologs.tsv",
#'   seed = 1
#' )
#' res$members
#' res$hogs[res$hogs$n_modules > 1, ]
#' }
#' @export
joint_modules <- function(nets, orthologs, modules = NULL, lambda = 1,
                          seed = NULL) {
  sp <- names(nets)
  bad <- !is.list(nets) || length(sp) < 2L || anyDuplicated(sp) > 0L
  if (bad || any(sp == "")) {
    stop(
      "`nets` must be a named list of at least two species; got ",
      if (is.null(sp)) "no names" else paste("names", toString(sp))
    )
  }
  ok <- is.numeric(lambda) && length(lambda) == 1L
  if (!ok || !isTRUE(is.finite(lambda) && lambda >= 0)) {
    stop("`lambda` must be one number of 0 or more; got ", deparse(lambda))
  }
  given <- !is.null(modules)
  if (given && !(is.list(modules) && setequal(names(modules), sp))) {
    stop(
      "`modules` must be named by the species of `nets` (", toString(sp),
      "); got ", if (is.null(names(modules))) {
        "no names"
      } else {
        toString(names(modules))
      }
    )
  }
  long <- .driver_orthologs(orthologs, sp)
  genes <- lapply(nets, function(n) rownames(n$network))
  off <- cumsum(c(0L, lengths(genes)))
  # link: gene (species concatenated, network genes only) to hog
  gene <- rep(NA_integer_, nrow(long))
  for (k in seq_along(sp)) {
    r <- long$species == sp[k]
    gene[r] <- off[k] + match(long$gene[r], genes[[k]])
  }
  link <- unique(data.frame(gene = gene, hog = long$hog)[!is.na(gene), ])

  .seed_scope(seed)
  if (!given) {
    root <- if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
    modules <- lapply(seq_along(sp), function(k) {
      detect_modules(nets[[k]], seed = .task_seed(root, 1L, k))
    })
    names(modules) <- sp
  }

  # Module nodes: species in nets order, labels in order of first gene
  mods <- lapply(sp, function(s) {
    m <- as_modules(modules[[s]])$modules
    miss <- setdiff(genes[[s]], names(m))
    if (length(miss) > 0L) {
      stop(
        "`modules$", s, "` has no label for ", length(miss),
        " gene(s) of its network, e.g. ", miss[1L],
        "; it needs one label per gene."
      )
    }
    m[genes[[s]]]
  })
  lv <- lapply(mods, unique)
  nsp <- rep(seq_along(sp), lengths(lv))
  node <- unlist(lapply(seq_along(sp), function(k) {
    sum(lengths(lv)[seq_len(k - 1L)]) + match(mods[[k]], lv[[k]])
  }), use.names = FALSE)

  g <- .joint_gain(nets, off, node, nsp, link, lambda)
  lab <- .joint_louvain(g)
  members <- data.frame(
    species = sp[nsp], module = unlist(lv, use.names = FALSE), label = lab,
    n_genes = tabulate(node, length(lab)), stringsAsFactors = FALSE
  )
  o <- order(members$label, nsp, members$module, method = "radix")
  members <- members[o, , drop = FALSE]
  rownames(members) <- NULL

  params <- list(
    lambda = lambda, gamma = 1, seed = seed, modules_given = given,
    quality = sum(g[outer(lab, lab, "==")])
  )
  glab <- lab[node]
  out <- lapply(stats::setNames(seq_along(sp), sp), function(k) {
    m <- stats::setNames(glab[off[k] + seq_along(genes[[k]])], genes[[k]])
    module_genes <- split(names(m), m)
    list(
      modules = m, module_genes = module_genes,
      n_modules = length(module_genes), method = "joint", params = params
    )
  })
  list(
    modules = out, members = members,
    hogs = .joint_hogs(glab, link), params = params
  )
}


#' Gain matrix of the module graph
#'
#' Entry (C, D) is what the pair adds to the objective when C and D share
#' a label. `node` maps each gene (species concatenated) to its module.
#' @noRd
.joint_gain <- function(nets, off, node, nsp, link, lambda) {
  n <- length(nsp)
  a <- matrix(0, n, n)
  for (k in seq_along(nets)) {
    e <- .joint_edges(nets[[k]], names(nets)[k])
    i <- node[off[k] + e$i]
    j <- node[off[k] + e$j]
    v <- rowsum(e$w / max(e$w), (j - 1L) * n + i)
    cell <- as.integer(rownames(v))
    a[cell] <- a[cell] + v[, 1L]
  }
  a <- a + t(a) # a[C, C] is twice the internal weight
  deg <- rowSums(a)
  two_m <- as.vector(rowsum(deg, nsp))
  tm <- two_m[nsp]

  # Hog weight 1 / (c_s c_t) per copy pair, summed per module pair
  gsp <- nsp[node[link$gene]]
  copies <- stats::ave(link$gene, link$hog, gsp, FUN = length)
  uh <- unique(link$hog)
  x <- Matrix::sparseMatrix(
    i = node[link$gene], j = match(link$hog, uh), x = 1 / copies,
    dims = c(n, length(uh))
  )
  same <- outer(nsp, nsp, "==")
  w <- as.matrix(Matrix::tcrossprod(x))
  w[same] <- 0
  h <- t(rowsum(t(w), nsp)) # h[C, t]: hog mass from C into species t
  hh <- rowsum(h, nsp) # hh[s, t]: total hog mass between s and t
  p <- h[, nsp, drop = FALSE]
  hn <- hh[nsp, nsp]
  null <- ifelse(hn > 0, p * t(p) / hn, 0)
  cross <- lambda * (w - null) / mean(two_m)
  within <- (a - outer(deg, deg) / tm) / tm
  ifelse(same, within, cross)
}


#' Louvain on a gain matrix: levels of local moves, then one refinement
#' pass on the original nodes. Labels 1..K in order of first appearance.
#' @noRd
.joint_louvain <- function(g) {
  part <- seq_len(nrow(g))
  gl <- g
  repeat {
    r <- .joint_moves(gl, seq_len(nrow(gl)))
    if (r$moves == 0L) break
    comm <- match(r$part, unique(r$part))
    part <- comm[part]
    gl <- rowsum(t(rowsum(gl, comm)), comm)
  }
  part <- .joint_moves(g, part)$part
  match(part, unique(part))
}


#' Local moves: each node joins the part with the largest gain, an empty
#' part counting 0. Passes repeat until none moves, at most 50.
#' @noRd
.joint_moves <- function(g, part) {
  n <- nrow(g)
  np <- matrix(0, n, n) # np[v, p]: summed gain of v to the nodes of part p
  np[, sort(unique(part))] <- t(rowsum(g, part))
  size <- tabulate(part, n)
  moves <- 0L
  for (pass in seq_len(50L)) {
    moved <- 0L
    for (v in sample.int(n)) {
      d <- part[v]
      size[d] <- size[d] - 1L
      gain <- np[v, ]
      gain[d] <- gain[d] - g[v, v]
      gain[size == 0L] <- 0
      b <- which.max(gain)
      if (gain[b] > gain[d]) {
        part[v] <- b
        np[, d] <- np[, d] - g[, v]
        np[, b] <- np[, b] + g[, v]
        moved <- moved + 1L
      } else {
        b <- d
      }
      size[b] <- size[b] + 1L
    }
    moves <- moves + moved
    if (moved == 0L) break
  }
  list(part = part, moves = moves)
}


#' Hog table on the shared gene labels: hogs with two or more copies
#' @noRd
.joint_hogs <- function(glab, link) {
  hogs <- sort(unique(link$hog[duplicated(link$hog)]), method = "radix")
  link <- link[link$hog %in% hogs, , drop = FALSE]
  by_hog <- split(glab[link$gene], factor(link$hog, levels = hogs))
  main <- vapply(by_hog, function(x) {
    u <- sort(unique(x))
    n <- tabulate(match(x, u))
    c(u[which.max(n)], max(n))
  }, integer(2))
  data.frame(
    hog = hogs, n_copies = lengths(by_hog, use.names = FALSE),
    n_modules = vapply(by_hog, function(x) length(unique(x)), 1L,
      USE.NAMES = FALSE
    ),
    module_main = main[1L, ], copies_main = main[2L, ],
    row.names = NULL, stringsAsFactors = FALSE
  )
}


#' Upper-triangle edges of one species at or above its threshold
#' @noRd
.joint_edges <- function(net, s) {
  m <- .net_check(net, net$threshold)
  thr <- net$threshold
  if (.net_is_sparse(net)) {
    i <- m@i + 1L
    j <- rep.int(seq_len(ncol(m)), diff(m@p))
    k <- i < j & m@x >= thr & m@x > 0
    e <- list(i = i[k], j = j[k], w = m@x[k])
  } else {
    k <- which(upper.tri(m) & m >= thr & m > 0, arr.ind = TRUE)
    e <- list(i = k[, 1L], j = k[, 2L], w = m[k])
  }
  if (length(e$w) == 0L) {
    stop(
      "`nets$", s, "` has no edge above its threshold; ",
      "there is nothing to cluster."
    )
  }
  e
}
