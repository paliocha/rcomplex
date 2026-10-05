#' Cross-species conservation of gene modules by neighbour voting
#'
#' @description
#' Tests whether each module of species 1 is still a co-expressed set in
#' species 2. The module's ortholog groups are translated to every
#' species-2 copy, and the translated set is scored by how well
#' species-2 network neighbours recognise held-out parts of it (a
#' cross-validated AUROC). The score is compared against random sets of
#' ortholog groups matched on size, copy number and degree.
#'
#' @details
#' **Statistic.** For module \eqn{S}, \eqn{T} is the set of species-2
#' genes orthologous to any gene of \eqn{S}, through every copy. The
#' adjacency is the binary species-2 network at its analysis threshold,
#' with edges inside one ortholog group dropped when `drop_within_hog`.
#' The ortholog groups of \eqn{S} are split into `n_fold` folds (by group,
#' so paralogs of a held-out gene are never training genes). With fold
#' \eqn{f} held out, every gene \eqn{j} outside the training part scores
#' \eqn{v(j) = } (edges from \eqn{j} to the training part of \eqn{T}) /
#' degree(\eqn{j}), and the fold's AUROC is that of the held-out genes of
#' \eqn{T} against all genes outside \eqn{T} (analytic rank sum, ties
#' counted one half). `auroc` is the mean over folds. This is the
#' neighbour-voting score of EGAD (Ballouz et al. 2017), applied to
#' a translated module as in CoCoCoNet (Lee et al. 2020): the
#' [compare_specificity()] score lifted from one anchor's neighbourhood
#' to a set. `degree_auroc` ranks species-2 genes by degree alone against
#' \eqn{T}; near 0.5 means the set is not a hub set and `auroc` is not
#' a degree artefact.
#'
#' **Null.** A random set of ortholog groups drawn from the species-1
#' groups that map to species 2, matched to the module group by group on
#' a stratum: species-1 copy class (1, 2, 3, 4+), species-2 copy class and
#' species-2 degree decile (mean degree of the group's species-2 copies).
#' The species-1 copy class is necessary: a module of genes picks groups
#' in proportion to their species-1 copy number, and without the stratum
#' modules built on shuffled expression scored z of about +2. The null set
#' goes through the identical statistic, folds included, so the p-value is
#' a conditional Monte Carlo p-value given the strata. A stratum with no
#' groups to spare (the module holds all of it) is taken whole;
#' `n_exhausted` counts the module's groups in such strata, where the null
#' cannot differ from the module.
#'
#' **Sequential stopping.** Every module gets `n_null` draws, then
#' batches of `batch` until `h` null draws reach the observed AUROC or
#' `max_draws` are made (Besag & Clifford 1991). Stopping at the `h`-th
#' exceedance makes `h / n_draws` the unbiased estimate; the reported
#' `p.val = (1 + n_exceed) / (1 + n_draws)` is the conservative form.
#' `p.val.gt` and `p.val.eq` split it into strictly greater and tied
#' draws (the observed value counted as a tie), so `p.val.gt + U *
#' p.val.eq` is the randomized p-value on which the Storey pi0 of `q.val`
#' is estimated, as in [summarize_comparison()].
#'
#' **No within-species p-value.** The test asks whether species 2 keeps a
#' set that species 1 defines. Whether the module is a module in species
#' 1 is a separate question (replication between sample halves, scored
#' with this same function on the other half's network), and a p-value
#' for it would not change this one.
#'
#' @section Calibration:
#' Check the null on the data before reading `z` or comparing module
#' engines. Shuffle species-1 expression per gene ([null_network()]),
#' detect modules on that network ([detect_modules()], or any partition
#' through [as_modules()]), and score them with `module_auroc()` against
#' the real species-2 network. These modules carry no co-expression, so
#' `z` should have mean about 0 and SD about 1, and about 5 % of `p.val`
#' should fall below 0.05:
#'
#' ```r
#' nn <- null_network(x1, net1, seed = 1)
#' m0 <- detect_modules(nn, objective_function = "modularity", seed = 1)
#' r0 <- module_auroc(m0, nn, net2, orthologs, seed = 1)
#' c(mean(r0$z), sd(r0$z), mean(r0$p.val < 0.05))
#' ```
#'
#' A departure means the null misses a stratum the module engine selects
#' on (the module benchmark saw mean 0.9 and SD 1.4 on one tissue). Then
#' recalibrate the per-species `z` against these shuffled-module scores,
#' \eqn{(z - \bar z_0) / s_0}, before comparing engines or species.
#'
#' @param modules Module map over species-1 genes, from [as_modules()]
#'   (or [detect_modules()]).
#' @param net1,net2 Network objects from [compute_network()] for species 1
#'   and 2. `net1` only restricts the ortholog table to its genes, as in
#'   [compare_specificity()]; `net2` is scored at `net2$threshold`.
#' @param orthologs Ortholog pair table with columns `Species1`,
#'   `Species2` and `hog` (e.g. from [prepare_orthologs()]); many-to-many.
#' @param n_null Draws made for every module before stopping is
#'   considered; `null_mean`, `null_sd` and `z` rest on at least this many.
#' @param max_draws Most draws per module.
#' @param h Exceedances at which a module stops drawing.
#' @param batch Draws per module per round after the first.
#' @param n_fold Folds over the module's ortholog groups (at least 2).
#' @param drop_within_hog Drop species-2 edges inside one ortholog group.
#' @param n_cores Threads for the scoring kernel (OpenMP); the result does
#'   not depend on it.
#' @param seed Seed for the folds, the null draws and the pi0 draws, or
#'   `NULL` for the ambient stream (package RNG contract).
#'
#' @return A data frame, one row per module with at least
#'   `max(3, n_fold)` mappable ortholog groups (others are dropped with a
#'   message), in module order: `module`, `n_hogs`, `n_genes_2` (size of
#'   \eqn{T}), `auroc`, `degree_auroc`, `null_mean`, `null_sd`, `z =
#'   (auroc - null_mean) / null_sd`, `n_draws`, `n_exceed` (draws >=
#'   `auroc`), `p.val`, `p.val.gt`, `p.val.eq`, `q.val` (Storey over the
#'   modules of the call), `n_exhausted`. Attributes `resolution`
#'   ([pvalue_resolution()] of `p.val` at `max_draws`) and `params`.
#'
#' @references
#' Ballouz, S., Weber, M., Pavlidis, P. & Gillis, J. (2017). EGAD:
#' ultra-fast functional analysis of gene networks. \emph{Bioinformatics},
#' 33(4), 612--614. \doi{10.1093/bioinformatics/btw695}
#'
#' Besag, J. & Clifford, P. (1991). Sequential Monte Carlo p-values.
#' \emph{Biometrika}, 78(2), 301--304. \doi{10.1093/biomet/78.2.301}
#'
#' Lee, J., Shah, M., Ballouz, S., Crow, M. & Gillis, J. (2020).
#' CoCoCoNet: conserved and comparative co-expression across a diverse
#' set of species. \emph{Nucleic Acids Research}, 48(W1), W566--W571.
#' \doi{10.1093/nar/gkaa348}
#'
#' @export
module_auroc <- function(modules, net1, net2, orthologs, n_null = 300L,
                         max_draws = 5000L, h = 10L, batch = 300L,
                         n_fold = 3L, drop_within_hog = TRUE, n_cores = 1L,
                         seed = NULL) {
  modules <- as_modules(modules)
  for (a in c("n_null", "max_draws", "h", "batch", "n_cores")) {
    v <- get(a)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 1) {
      stop(a, " must be a single number >= 1")
    }
  }
  if (!is.numeric(n_fold) || length(n_fold) != 1L || is.na(n_fold) ||
        n_fold < 2) {
    stop("n_fold must be a single number >= 2")
  }
  if (max_draws < n_null) stop("max_draws must be at least n_null")
  .seed_scope(seed)

  op <- .ortholog_pair_index(net1, net2, orthologs)
  ort <- op$orthologs
  tg <- .module_auroc_target(net2, ort, op$sp2_idx, drop_within_hog)

  # HOG universe: groups with a copy in both networks, and their strata
  hog <- as.character(ort$hog)
  copies2 <- split(op$sp2_idx, hog)
  copies2 <- lapply(copies2, unique)
  univ <- names(copies2)
  n1 <- lengths(lapply(split(ort$Species1, hog), unique))[univ]
  n2 <- lengths(copies2)
  hdeg <- vapply(copies2, function(g) mean(tg$deg[g + 1L]), numeric(1))
  dec <- ceiling(10 * rank(hdeg, ties.method = "first") / length(hdeg))
  stratum <- stats::setNames(
    paste(pmin(n1, 4L), pmin(n2, 4L), dec), univ
  )
  pools <- split(univ, stratum)

  min_hogs <- max(3L, n_fold)
  mod_hogs <- lapply(modules$module_genes, function(g) {
    unique(hog[ort$Species1 %in% g])
  })
  keep <- lengths(mod_hogs) >= min_hogs
  if (!all(keep)) {
    message(
      sum(!keep), " module(s) with fewer than ", min_hogs,
      " mappable ortholog groups dropped"
    )
  }
  mod_hogs <- mod_hogs[keep]
  nm <- length(mod_hogs)
  if (nm == 0L) stop("no module has ", min_hogs, " mappable ortholog groups")

  score <- function(sets) {
    gs <- lapply(sets, function(hs) {
      f <- sample(rep_len(seq_len(n_fold), length(hs)))
      g <- copies2[hs]
      list(g = unlist(g, use.names = FALSE), f = rep.int(f, lengths(g)))
    })
    module_auroc_cpp(
      tg$p, tg$i,
      c(0L, cumsum(vapply(gs, function(x) length(x$g), integer(1)))),
      unlist(lapply(gs, `[[`, "g"), use.names = FALSE),
      unlist(lapply(gs, `[[`, "f"), use.names = FALSE),
      as.integer(n_fold), as.integer(n_cores)
    )
  }
  obs <- score(mod_hogs)

  counts <- lapply(mod_hogs, function(hs) table(stratum[hs]))
  n_exhausted <- vapply(counts, function(cn) {
    sum(cn[lengths(pools[names(cn)]) == cn])
  }, numeric(1))
  draw <- function(cn) {
    unlist(lapply(names(cn), function(s) {
      pool <- pools[[s]]
      pool[sample.int(length(pool), cn[[s]])]
    }), use.names = FALSE)
  }

  null <- vector("list", nm)
  n_draws <- n_ge <- integer(nm)
  active <- seq_len(nm)
  size <- n_null
  while (length(active)) {
    k <- pmin(size, max_draws - n_draws[active])
    sets <- unlist(lapply(seq_along(active), function(a) {
      replicate(k[a], draw(counts[[active[a]]]), simplify = FALSE)
    }), recursive = FALSE)
    au <- split(score(sets), rep.int(seq_along(active), k))
    for (a in seq_along(active)) {
      m <- active[a]
      null[[m]] <- c(null[[m]], au[[a]])
      n_draws[m] <- length(null[[m]])
      n_ge[m] <- sum(null[[m]] >= obs[m])
    }
    active <- active[n_ge[active] < h & n_draws[active] < max_draws]
    size <- batch
  }

  n_gt <- vapply(seq_len(nm), function(m) sum(null[[m]] > obs[m]), 0)
  null_mean <- vapply(null, mean, 0)
  null_sd <- vapply(null, stats::sd, 0)
  p_gt <- n_gt / (1 + n_draws)
  p_eq <- (1 + n_ge - n_gt) / (1 + n_draws)
  p_val <- (1 + n_ge) / (1 + n_draws)
  q_val <- compute_qvalues(p_val, function() {
    p_gt + stats::runif(nm) * p_eq
  })$qvalues

  n_genes_2 <- vapply(mod_hogs, function(hs) {
    length(unlist(copies2[hs], use.names = FALSE))
  }, integer(1))
  r <- rank(tg$deg)
  degree_auroc <- vapply(mod_hogs, function(hs) {
    pos <- logical(length(tg$deg))
    pos[unlist(copies2[hs], use.names = FALSE) + 1L] <- TRUE
    n_pos <- sum(pos)
    (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * sum(!pos))
  }, numeric(1))

  res <- data.frame(
    module = names(mod_hogs), n_hogs = lengths(mod_hogs),
    n_genes_2 = n_genes_2, auroc = obs, degree_auroc = degree_auroc,
    null_mean = null_mean, null_sd = null_sd,
    z = (obs - null_mean) / null_sd, n_draws = n_draws,
    n_exceed = n_ge, p.val = p_val, p.val.gt = p_gt, p.val.eq = p_eq,
    q.val = q_val, n_exhausted = n_exhausted,
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(res, "resolution") <- pvalue_resolution(p_val, n_perm = max_draws)
  attr(res, "params") <- list(
    n_null = n_null, max_draws = max_draws, h = h, batch = batch,
    n_fold = n_fold, drop_within_hog = drop_within_hog
  )
  res
}


#' Reciprocal cross-species module conservation
#'
#' @description
#' Runs [module_auroc()] in both directions, species 1 modules in the
#' species-2 network and species 2 modules in the species-1 network,
#' pairs the modules of the two species by reciprocal best hit on shared
#' ortholog groups, and combines the two directions' p-values per pair.
#'
#' @details
#' **Matching.** Each module is reduced to the ortholog groups (`hog`)
#' of its genes. Module \eqn{a} of species 1 and \eqn{b} of species 2
#' pair when \eqn{b} has the highest Jaccard index over groups with
#' \eqn{a} among species-2 modules and \eqn{a} the highest with \eqn{b}
#' among species-1 modules (ties go to the first module), the Jaccard is
#' above 0, and both were tested by [module_auroc()].
#' [module_correspondence()] is not used: its Jaccard is over genes
#' projected through a [resolve_ortholog_map()] map, one copy per
#' group, while a module here is scored through every copy of its groups,
#' so the overlap that matches the test is the overlap of groups.
#'
#' **Combination.** `pval_combine = "max"` (default) takes `pmax` of the
#' two directions: a pair counts as conserved only if each species'
#' module holds in the other's network, the reciprocal criterion
#' [comparison_to_edges()] applies to gene pairs (Netotea et al. 2014).
#' `pmax` is a valid p-value for that intersection-union null. `"min"`
#' accepts either direction and is not corrected for taking the smaller
#' of two. The Storey pi0 behind `q.val` is estimated on the same
#' combination of each direction's randomized p-value `p.val.gt + U *
#' p.val.eq`, as in [module_auroc()].
#'
#' @param modules1,modules2 Module maps over species-1 and species-2
#'   genes ([as_modules()] or [detect_modules()]).
#' @param net1,net2 Network objects from [compute_network()].
#' @param orthologs Ortholog pair table with `Species1` (species-1 genes),
#'   `Species2` and `hog`; swapped internally for the reverse direction.
#' @param pval_combine `"max"` (both directions) or `"min"` (either).
#' @param ... Further arguments to [module_auroc()], used in both
#'   directions (`n_null`, `max_draws`, `n_cores`, ...).
#' @param seed Seed for both directions and the pi0 draws, or `NULL` for
#'   the ambient stream (package RNG contract).
#'
#' @return A data frame, one row per matched pair in species-1 module
#'   order: `module1`, `module2`, `jaccard`, `p.val.1to2`, `p.val.2to1`,
#'   `p.val`, `q.val`, `z.1to2`, `z.2to1`. Attributes `unmatched` (a list
#'   with `species1` and `species2`: modules with no reciprocal partner,
#'   untested modules included) and `resolution` ([pvalue_resolution()]
#'   of `p.val`).
#'
#' @references
#' Netotea, S., Sundell, D., Street, N. R. & Hvidsten, T. R. (2014).
#' Evolution of the plant co-expression network. \emph{BMC Genomics},
#' 15, 106. \doi{10.1186/1471-2164-15-106}
#'
#' @export
module_auroc_reciprocal <- function(modules1, modules2, net1, net2,
                                    orthologs,
                                    pval_combine = c("max", "min"), ...,
                                    seed = NULL) {
  pval_combine <- match.arg(pval_combine)
  comb <- if (pval_combine == "max") pmax else pmin
  modules1 <- as_modules(modules1)
  modules2 <- as_modules(modules2)
  .seed_scope(seed)

  rev_ort <- orthologs
  rev_ort$Species1 <- orthologs$Species2
  rev_ort$Species2 <- orthologs$Species1
  r12 <- module_auroc(modules1, net1, net2, orthologs, ..., seed = NULL)
  r21 <- module_auroc(modules2, net2, net1, rev_ort, ..., seed = NULL)

  hog <- as.character(orthologs$hog)
  hogs <- function(mods, genes) {
    lapply(mods$module_genes, function(g) unique(hog[genes %in% g]))
  }
  h1 <- hogs(modules1, orthologs$Species1)[r12$module]
  h2 <- hogs(modules2, orthologs$Species2)[r21$module]
  inc <- function(h) {
    Matrix::sparseMatrix(
      i = match(unlist(h, use.names = FALSE), unique(hog)),
      j = rep.int(seq_along(h), lengths(h)),
      x = 1, dims = c(length(unique(hog)), length(h))
    )
  }
  inter <- as.matrix(Matrix::crossprod(inc(h1), inc(h2)))
  jac <- inter / (outer(lengths(h1), lengths(h2), "+") - inter)
  best2 <- max.col(jac, ties.method = "first")
  best1 <- max.col(t(jac), ties.method = "first")
  i <- which(best1[best2] == seq_along(best2))
  j <- best2[i]
  i <- i[jac[cbind(i, j)] > 0]
  j <- best2[i]
  if (!length(i)) stop("no module pair is a reciprocal best hit")

  a <- r12[i, ]
  b <- r21[j, ]
  p_val <- comb(a$p.val, b$p.val)
  q_val <- compute_qvalues(p_val, function() {
    comb(
      a$p.val.gt + stats::runif(length(i)) * a$p.val.eq,
      b$p.val.gt + stats::runif(length(i)) * b$p.val.eq
    )
  })$qvalues

  res <- data.frame(
    module1 = a$module, module2 = b$module, jaccard = jac[cbind(i, j)],
    p.val.1to2 = a$p.val, p.val.2to1 = b$p.val, p.val = p_val,
    q.val = q_val, z.1to2 = a$z, z.2to1 = b$z,
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(res, "unmatched") <- list(
    species1 = setdiff(names(modules1$module_genes), res$module1),
    species2 = setdiff(names(modules2$module_genes), res$module2)
  )
  attr(res, "resolution") <- pvalue_resolution(
    p_val,
    n_perm = attr(r12, "params")$max_draws
  )
  res
}


#' Binary species-2 adjacency for module_auroc()
#'
#' Entries at or above `net$threshold` from the sparse store (dense
#' networks converted, every nonzero stored), minus edges whose two genes
#' share an ortholog group when `drop_within_hog`. A gene in several groups
#' keeps its first. Returns the CSC slots and the degree.
#' @noRd
.module_auroc_target <- function(net, ort, sp2_idx, drop_within_hog) {
  a <- .net_cpp_args(.net_as_sparse(net), net$threshold)
  n <- length(a$p) - 1L
  col <- rep.int(seq_len(n), diff(a$p))
  keep <- a$x >= a$thr
  if (drop_within_hog) {
    hid <- rep(NA_integer_, n)
    first <- !duplicated(sp2_idx)
    hid[sp2_idx[first] + 1L] <- match(ort$hog, unique(ort$hog))[first]
    same <- hid[a$i + 1L] == hid[col]
    keep <- keep & (is.na(same) | !same)
  }
  deg <- tabulate(col[keep], n)
  list(p = c(0L, cumsum(deg)), i = a$i[keep], deg = deg)
}
