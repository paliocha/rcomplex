# Structify-Net-style toy networks on a species tree: one shared rank
# rule, a noise knob, and planted gains, losses, switches and
# duplications of co-expression context, inherited down the tree.
# Model: Cazabet, Citraro & Rossetti 2023, Peer Community Journal 3:e103,
# doi:10.24072/pcjournal.335. Mixture instead of the Bezier curve, same
# ends. Paralogs only through `dup` / `dup_switch` events (two copies).
#
# Rank rules, ported from structify_net.zoo:
#   blocks  same block -> 0, else 1 (sort_blocks_assortative)
#   spatial Euclidean distance on (x, y) (sort_distances; the zoo
#           default is 1 dimension, this probe uses 2)
#   core    |x_u - .5| * |x_v - .5| * |x_u - x_v|: the file's product of
#           three distances, core at 0.5 (sort_core_distance), not x_u+x_v
#   nested  x_u + x_v, ascending (sort_nestedness ranks u + v on node
#           index; uniform x stands in for the index), not max(x_u, x_v)
# A gene with none == TRUE has d = Inf with every partner.
suppressPackageStartupMessages(library(rcomplex))

STRUCTURES <- c("blocks", "spatial", "core", "nested")
N_BLOCKS <- 20L

# d for every unordered pair, as a vector in combn(n, 2) order
# (= lower triangle of the n x n matrix, column-major).
structify_ranks <- function(attr, structure) {
  structure <- match.arg(structure, STRUCTURES)
  d <- switch(structure,
    blocks = 1 * outer(attr$block, attr$block, "!="),
    spatial = as.matrix(stats::dist(cbind(attr$x, attr$y))),
    core = {
      c0 <- abs(attr$x - 0.5)
      outer(c0, c0) * abs(outer(attr$x, attr$x, "-"))
    },
    nested = outer(attr$x, attr$x, "+")
  )
  d[attr$none, ] <- Inf
  d[, attr$none] <- Inf
  d[lower.tri(d)]
}

# Exactly m edges with the marginals of P(r) = (1 - eps) [r <= m] +
# eps m / L. P has two levels, so draw K1 ~ Binom(m, p1) pairs from the
# m top ranks and m - K1 from the rest, uniformly within each level.
structify_edges <- function(d, m, epsilon, genes) {
  L <- length(d)
  stopifnot(m >= 1, m < L, epsilon >= 0, epsilon <= 1)
  r <- integer(L)
  r[order(d, stats::runif(L))] <- seq_len(L) # random tie-break
  p1 <- (1 - epsilon) + epsilon * m / L
  k1 <- min(stats::rbinom(1L, m, p1), m)
  k1 <- max(k1, m - (L - m))
  top <- which(r <= m)
  rest <- which(r > m)
  keep <- c(top[sample.int(m, k1)], rest[sample.int(L - m, m - k1)])
  # combn index -> (i, j): lower triangle column-major, row j > col i
  n <- round((1 + sqrt(1 + 8 * L)) / 2)
  ij <- which(lower.tri(matrix(0, n, n)), arr.ind = TRUE)[keep, ]
  data.frame(
    gene1 = genes[ij[, 2L]], gene2 = genes[ij[, 1L]],
    weight = (L - r[keep] + 1) / L, stringsAsFactors = FALSE
  )
}

# A fresh attribute for k genes; for blocks, never the current block.
.fresh <- function(structure, k, block = NULL) {
  if (structure == "blocks") {
    b <- sample.int(N_BLOCKS, k, replace = TRUE)
    if (!is.null(block)) {
      same <- !is.na(block) & b == block
      while (any(same)) {
        b[same] <- sample.int(N_BLOCKS, sum(same), replace = TRUE)
        same <- !is.na(block) & b == block
      }
    }
    list(block = b)
  } else if (structure == "spatial") {
    list(x = stats::runif(k), y = stats::runif(k))
  } else {
    list(x = stats::runif(k))
  }
}

# Root attributes: one row per gene; a fraction p_none has no structure.
structify_root <- function(n, structure, p_none = 0.1) {
  gene <- sprintf("g%04d", seq_len(n))
  a <- if (structure == "blocks") {
    list(block = sample(rep_len(seq_len(N_BLOCKS), n)))
  } else {
    .fresh(structure, n)
  }
  none <- seq_len(n) %in% sample.int(n, round(p_none * n))
  data.frame(
    copy = gene, hog = gene, a, none = none,
    stringsAsFactors = FALSE
  )
}

.branches <- function(tree) {
  nm <- rcomplex:::.node_names(tree)
  nm[-(length(tree$tip.label) + 1L)]
}

# Events drawn uniformly over genes and non-root branches; gain only from
# none genes, every other kind only from structured ones; one per gene.
# n_root_dup genes get a dup on "root" (two copies in every species).
structify_events <- function(tree, genes, n_gain, n_loss, n_switch,
                             n_dup = 0L, n_dup_switch = 0L,
                             n_root_dup = 0L, seed = 1L) {
  withr::with_seed(seed, {
    br <- .branches(tree)
    none <- genes$hog[genes$none]
    some <- genes$hog[!genes$none]
    kinds <- c("loss", "switch", "dup", "dup_switch", "root_dup")
    n_k <- c(n_loss, n_switch, n_dup, n_dup_switch, n_root_dup)
    stopifnot(n_gain <= length(none), sum(n_k) <= length(some))
    g_s <- sample(some, sum(n_k))
    kind <- rep(kinds, n_k)
    ev <- data.frame(
      gene = c(sample(none, n_gain), g_s),
      kind = c(rep("gain", n_gain), sub("root_dup", "dup", kind)),
      branch = c(
        sample(br, n_gain, replace = TRUE),
        ifelse(kind == "root_dup", "root",
          sample(br, length(kind), replace = TRUE)
        )
      ),
      stringsAsFactors = FALSE
    )
    ev[order(ev$gene), ]
  })
}

# Apply the events on one node to the table inherited from its parent.
.apply_events <- function(tab, ev, structure) {
  for (k in seq_len(nrow(ev))) {
    h <- ev$gene[k]
    i <- which(tab$hog == h)
    fr <- if (structure == "blocks") tab$block[i] else NULL
    kind <- ev$kind[k]
    if (kind %in% c("gain", "switch")) {
      tab[i, names(.fresh(structure, 1))] <- .fresh(structure, length(i), fr)
      tab$none[i] <- FALSE
    } else if (kind == "loss") {
      tab$none[i] <- TRUE
    } else { # dup, dup_switch: twin copies <hog>.1 and <hog>.2
      two <- tab[c(i, i), ]
      two$copy <- paste0(h, c(".1", ".2"))
      if (kind == "dup_switch") {
        two[2L, names(.fresh(structure, 1))] <- .fresh(structure, 1, fr)
      }
      tab <- rbind(tab[-i, ], two)
    }
  }
  tab
}

# One network per tip; see the header for the model.
structify_tree <- function(tree, n, structure, m = 0.03, epsilon = 0,
                           events = NULL, seed = 1L, p_none = 0.1) {
  structure <- match.arg(structure, STRUCTURES)
  withr::with_seed(seed, {
    if (is.null(tree$node.label)) tree$node.label <- character(tree$Nnode)
    nm <- rcomplex:::.node_names(tree)
    tree$node.label <- nm[-seq_along(tree$tip.label)]
    tips <- tree$tip.label
    root_id <- length(tips) + 1L
    root <- if (is.data.frame(n)) n else structify_root(n, structure, p_none)
    if (is.null(events)) {
      events <- data.frame(
        gene = character(), kind = character(),
        branch = character()
      )
    }
    stopifnot(
      !anyDuplicated(events$gene), events$branch %in% nm,
      events$branch != "root" | events$kind == "dup"
    )
    on <- function(b) events[events$branch == b, , drop = FALSE]
    attr <- vector("list", length(nm))
    names(attr) <- nm
    attr[["root"]] <- .apply_events(root, on(nm[root_id]), structure)
    ed <- ape::reorder.phylo(tree, "cladewise")$edge # parents first
    for (k in seq_len(nrow(ed))) {
      ch <- nm[ed[k, 2L]]
      attr[[ch]] <- .apply_events(attr[[nm[ed[k, 1L]]]], on(ch), structure)
    }
    if (length(epsilon) == 1L) epsilon <- rep(epsilon, length(tips))
    names(epsilon) <- tips
    networks <- lapply(stats::setNames(tips, tips), function(sp) {
      a <- attr[[sp]][order(attr[[sp]]$copy), ]
      d <- structify_ranks(a, structure)
      mm <- floor(m * length(d) + 0.5) # as_network()'s k
      el <- structify_edges(d, mm, epsilon[[sp]], a$copy)
      as_network(el, density = m, genes = a$copy)
    })
    orthologs <- do.call(rbind, lapply(tips, function(sp) {
      data.frame(
        species = sp, gene = attr[[sp]]$copy,
        hog = attr[[sp]]$hog, stringsAsFactors = FALSE
      )
    }))
    orthologs <- orthologs[order(orthologs$species, orthologs$gene), ]
    rownames(orthologs) <- NULL
    hogs <- root$hog
    ev <- events[match(hogs, events$gene), ]
    per_tip <- function(f) {
      vapply(hogs, function(h) {
        paste(vapply(tips, function(sp) f(attr[[sp]], h), ""),
          collapse = ","
        )
      }, "")
    }
    truth <- data.frame(
      hog = hogs, kind = ifelse(is.na(ev$kind), "none", ev$kind),
      branch = ev$branch,
      root_none = root$none,
      structured_in = vapply(hogs, function(h) {
        paste(tips[vapply(tips, function(sp) {
          any(!attr[[sp]]$none[attr[[sp]]$hog == h])
        }, NA)], collapse = ",")
      }, ""),
      copies = per_tip(function(a, h) as.character(sum(a$hog == h))),
      kept = ifelse(ev$kind %in% "dup_switch", ".1", NA_character_),
      stringsAsFactors = FALSE, row.names = NULL
    )
    list(
      networks = networks, orthologs = orthologs, events = events,
      truth = truth, attr = attr, tree = tree,
      params = list(
        structure = structure, n = nrow(root), m = m,
        epsilon = epsilon, seed = seed, p_none = p_none
      )
    )
  })
}
