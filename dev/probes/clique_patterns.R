# was rcomplex::conservation_pattern_table conservation_lattice bicm_species_z until 0.4.0; see git history (main @ b245ac0)

# Clique x species conservation patterns: the 3-valued pattern table,
# its closed-set (iceberg) lattice, and the BiCM species-pair z.

# Columns of conservation_pattern_table() that are not species states.
.cpt_fixed <- c(
  "clique_id", "hog", "classification", "pattern", "n_plus",
  "n_minus", "pairs_sig", "pairs_ns", "pairs_untested"
)


#' Clique x species conservation pattern table
#'
#' One row per gene clique from [gene_clique_graph()], with a 3-valued
#' state per species: `"+"` the species is a member of the clique; `"-"`
#' it is not, but was compared against the clique and rejected at
#' adequate power (`missing_reason == "tested_ns"` in
#' [classify_gene_cliques()]) -- negative evidence; `"?"` otherwise
#' (`absent`, `untested`, `underpowered` or `extendable`) -- absent
#' evidence. The states are read off [classify_gene_cliques()]'s own
#' bookkeeping, so `min_power` and the other gates act exactly as they
#' do there.
#'
#' Rows are never merged over paralogs: a HOG whose copy 1 is conserved
#' in one trait group and copy 2 in the other stays two rows, where an
#' OR over copies would read "conserved in all species".
#'
#' @param cliques,edges,species As for [classify_gene_cliques()]: the
#'   gene-graph cliques, the **unfiltered** edge table, and every species
#'   in the analysis.
#' @param min_power Passed to [classify_gene_cliques()]; `NULL` (default)
#'   keeps its default.
#' @param ... Further arguments to [classify_gene_cliques()] (`lineage`,
#'   `alpha_call`, ...).
#' @return A data frame, one row per clique: `clique_id`, `hog`,
#'   `classification` (the tier of [classify_gene_cliques()]), `pattern`
#'   (the states concatenated in `species` order), `n_plus`, `n_minus`,
#'   and the member species pairs split into `pairs_sig` (significant at
#'   `alpha_call`), `pairs_ns` (tested, not significant) and
#'   `pairs_untested` (no row in `edges`); then one state column per
#'   species, named by the species.
#' @examples
#' edges <- data.frame(
#'   gene1 = c("a1", "a1", "b1", "a1"), gene2 = c("b1", "c1", "c1", "d1"),
#'   species1 = c("SP_A", "SP_A", "SP_B", "SP_A"),
#'   species2 = c("SP_B", "SP_C", "SP_C", "SP_D"),
#'   hog = "HOG1", q.value = c(0.01, 0.02, 0.03, 0.9)
#' )
#' cl <- gene_clique_graph(edges, alpha_graph = 0.1)
#' conservation_pattern_table(cl, edges, c("SP_A", "SP_B", "SP_C", "SP_D"))
#' @seealso [conservation_lattice()], [bicm_species_z()]
#' @export
conservation_pattern_table <- function(cliques, edges, species,
                                       min_power = NULL, ...) {
  species <- as.character(species)
  clash <- intersect(species, .cpt_fixed)
  if (length(clash) > 0L) {
    stop("species names clash with output columns: ", toString(clash))
  }
  args <- list(cliques, edges, species, ...)
  if (!is.null(min_power)) args$min_power <- min_power
  cls <- do.call(classify_gene_cliques, args)

  st <- matrix("?", nrow(cls), length(species),
    dimnames = list(NULL, species)
  )
  cl_id <- as.character(cliques$clique_id)
  mem <- split(as.character(cliques$species), cl_id)[cls$clique_id]
  for (i in seq_len(nrow(cls))) {
    st[i, intersect(mem[[i]], species)] <- "+"
    if (cls$n_missing[i] > 0L) {
      gone <- strsplit(cls$missing_species[i], ",", fixed = TRUE)[[1L]]
      why <- strsplit(cls$missing_reason[i], ",", fixed = TRUE)[[1L]]
      st[i, gone[why == "tested_ns"]] <- "-"
    }
  }
  out <- data.frame(
    clique_id = cls$clique_id, hog = cls$hog,
    classification = cls$classification,
    pattern = apply(st, 1L, paste, collapse = ""),
    n_plus = rowSums(st == "+"), n_minus = rowSums(st == "-"),
    pairs_sig = cls$n_sig, pairs_ns = cls$n_present - cls$n_sig,
    pairs_untested = cls$n_pairs - cls$n_present,
    stringsAsFactors = FALSE
  )
  if (nrow(out) == 0L) out$pattern <- character(0)
  out <- cbind(out, as.data.frame(st, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}


# Species columns and the logical "+" matrix of a pattern table.
.cpt_plus <- function(pattern_table) {
  ok <- is.data.frame(pattern_table) &&
    all(.cpt_fixed %in% names(pattern_table))
  if (!ok) {
    stop("pattern_table must come from conservation_pattern_table()")
  }
  sp <- setdiff(names(pattern_table), .cpt_fixed)
  if (length(sp) < 2L) stop("pattern_table needs at least 2 species")
  m <- as.matrix(pattern_table[sp]) == "+"
  dimnames(m) <- list(NULL, sp)
  m
}


#' Closed species sets (iceberg concept lattice) of a pattern table
#'
#' The intents of the formal context cliques x species under `"+"`: every
#' species set that is the intersection of the `+`-sets of some cliques
#' (closure by intersection). An intent is the largest species set shared
#' by exactly the cliques that contain it, so the lattice lists every
#' distinct conservation pattern and every pattern implied by a group of
#' them, with no species set repeated under a different support. Pure
#' base R over bitmasks, so the species count is limited to 30; with
#' `S` species there are at most `2^S` intents.
#'
#' With `trait`, intents equal to all species are labelled `"complete"`
#' and intents equal to one trait level's whole species set carry that
#' level's name, so the species-level tiers of [classify_gene_cliques()]
#' read off as named intents: `complete_conserved` cliques sit at
#' `"complete"`; `lineage_specific` and `trait_specific` cliques at the
#' trait level's intent, told apart by `n_minus_outside` (0 for a
#' lineage-specific clique, whose outside species are `"?"`, positive
#' for a trait-specific one, whose rejections are `"-"`). Pair-level
#' tiers (`partial_significant`, `differentiated`) are not species sets
#' and are not intents.
#'
#' @param pattern_table Output of [conservation_pattern_table()].
#' @param min_support Keep intents contained in at least this many
#'   cliques' `+`-sets (default 1).
#' @param trait Optional named vector mapping species to trait level, used
#'   only for the `label` column.
#' @return A list with `intents`, a data frame with `intent` (species
#'   comma-separated, in table order), `size`, `support` (cliques whose
#'   `+`-set contains the intent), `support_exact` (cliques whose `+`-set
#'   equals it), `n_minus_outside` (of those exact cliques, how many have
#'   at least one `"-"` species) and `label` (`NA` without `trait`),
#'   ordered by decreasing size, then decreasing support; and
#'   `containment`, a data frame of every strict inclusion `subset` <
#'   `superset` between kept intents.
#' @examples
#' pt <- data.frame(
#'   clique_id = c("c1", "c2", "c3"), hog = "H", classification = "x",
#'   pattern = c("+++-", "+++?", "++++"), n_plus = c(3, 3, 4),
#'   n_minus = c(1, 0, 0), pairs_sig = 3, pairs_ns = 0, pairs_untested = 0,
#'   A = "+", B = "+", C = "+", D = c("-", "?", "+")
#' )
#' conservation_lattice(pt, trait = c(A = "x", B = "x", C = "x", D = "y"))
#' @seealso [conservation_pattern_table()]
#' @export
conservation_lattice <- function(pattern_table, min_support = 1L,
                                 trait = NULL) {
  m <- .cpt_plus(pattern_table)
  sp <- colnames(m)
  if (length(sp) > 30L) stop("at most 30 species are supported")
  bits <- 2^(seq_along(sp) - 1L)
  mask <- as.integer(m %*% bits)
  cnt <- table(mask)
  rows <- as.integer(names(cnt))
  # Closure by intersection: add pairwise ANDs until nothing new appears.
  # ponytail: quadratic in the number of intents (<= 2^S), fine for S <= 10.
  int <- rows
  repeat {
    nxt <- unique(c(int, as.vector(outer(int, int, bitwAnd))))
    if (length(nxt) == length(int)) break
    int <- nxt
  }
  int <- int[int != 0L]
  supp <- vapply(int, function(i) {
    sum(cnt[bitwAnd(rows, i) == i])
  }, numeric(1))
  int <- int[supp >= min_support]
  supp <- supp[supp >= min_support]
  has_minus <- pattern_table$n_minus > 0
  n_ex <- vapply(int, function(i) sum(mask == i), numeric(1))
  n_mo <- vapply(int, function(i) sum(mask == i & has_minus), numeric(1))
  name_of <- function(i) paste(sp[bitwAnd(i, bits) > 0], collapse = ",")
  size <- vapply(int, function(i) sum(bitwAnd(i, bits) > 0), integer(1))

  label <- rep(NA_character_, length(int))
  if (!is.null(trait)) {
    if (is.null(names(trait)) || !all(sp %in% names(trait))) {
      stop("trait must be a named vector covering every species")
    }
    tr <- as.character(trait[sp])
    for (lv in unique(tr)) {
      label[int == sum(bits[tr == lv])] <- lv
    }
    label[int == sum(bits)] <- "complete"
  }
  ord <- order(-size, -supp)
  int <- int[ord]
  intents <- data.frame(
    intent = vapply(int, name_of, character(1)), size = size[ord],
    support = as.integer(supp[ord]), support_exact = as.integer(n_ex[ord]),
    n_minus_outside = as.integer(n_mo[ord]), label = label[ord],
    stringsAsFactors = FALSE
  )
  sub <- outer(int, int, function(a, b) bitwAnd(a, b) == a & a != b)
  ij <- which(sub, arr.ind = TRUE)
  containment <- data.frame(
    subset = intents$intent[ij[, 1L]], superset = intents$intent[ij[, 2L]],
    stringsAsFactors = FALSE
  )
  list(intents = intents, containment = containment)
}


# Fit the BiCM to a binary matrix: p_is = x_i y_s / (1 + x_i y_s) with
# expected row and column sums equal to the observed ones. Cells forced
# by the margins (a row or column whose residual degree is 0 or equal to
# its free cells) are fixed to 0 / 1 first, repeatedly; what remains is a
# free rectangle, fitted by the fixed point
#   x_k = r_k / sum_s y_s / (1 + x_k y_s)
#   y_s = c_s / sum_k n_k x_k / (1 + x_k y_s)
# over the distinct residual row degrees k (n_k rows each) and species s:
# the score equations of the BiCM likelihood, solved for x and y in turn
# (Vallarano et al. 2021).
.bicm_fit <- function(m, tol = 1e-10, max_iter = 100000L) {
  p <- matrix(NA_real_, nrow(m), ncol(m))
  k <- rowSums(m)
  d <- colSums(m)
  repeat {
    fr <- rowSums(is.na(p))
    rr <- k - rowSums(p, na.rm = TRUE)
    fc <- colSums(is.na(p))
    rc <- d - colSums(p, na.rm = TRUE)
    row_fix <- fr > 0 & (rr == 0 | rr == fr)
    col_fix <- fc > 0 & (rc == 0 | rc == fc)
    if (!any(row_fix) && !any(col_fix)) break
    for (i in which(row_fix)) p[i, is.na(p[i, ])] <- as.numeric(rr[i] > 0)
    for (s in which(col_fix)) p[is.na(p[, s]), s] <- as.numeric(rc[s] > 0)
  }
  fi <- which(fr > 0)
  fs <- which(fc > 0)
  if (length(fi) == 0L) {
    return(p)
  }
  cls <- table(rr[fi])
  r <- as.numeric(names(cls))
  n <- as.numeric(cls)
  cs <- rc[fs]
  tot <- sum(cs)
  x <- r / sqrt(tot)
  y <- cs / sqrt(tot)
  for (it in seq_len(max_iter)) {
    x <- r / colSums(y / (1 + outer(y, x)))
    q <- outer(x, y)
    y <- cs / colSums(n * x / (1 + q))
    q <- outer(x, y)
    pk <- q / (1 + q)
    err <- max(abs(rowSums(pk) - r), abs(colSums(n * pk) - cs))
    if (err < tol) break
  }
  if (err >= tol) {
    warning("BiCM fixed point did not converge (max degree error ",
      signif(err, 3), ")",
      call. = FALSE
    )
  }
  p[fi, fs] <- pk[match(rr[fi], r), , drop = FALSE]
  p
}


#' Margin-corrected species-pair co-membership z (BiCM)
#'
#' Fits the bipartite configuration model (BiCM; Saracco et al. 2017) to
#' the clique x species `"+"` matrix of a pattern table -- independent
#' cells with `p_is = x_i y_s / (1 + x_i y_s)`, the maximum-entropy model
#' whose expected clique sizes and species degrees equal the observed ones
#' -- and standardises each species pair's V-motif count `V_ab` (cliques
#' containing both) by its Poisson-binomial mean `sum_i p_ia p_ib` and
#' variance `sum_i p_ia p_ib (1 - p_ia p_ib)`. A positive z means two
#' species share more cliques than their own conservation rates and the
#' clique sizes imply.
#'
#' BiCM is a *standardiser*, not a test: the cells are dependent (a
#' clique needs at least three species, edges are called pairwise, and
#' phylogeny ties species together), so the z values are not standard
#' normal under any null of interest and within-genus pairs all come out
#' high. Ask the trait question with the relabelling machinery of
#' [preservation_matrix_test()], whose `classification` input this
#' `pairs` table fills once z is named as its statistic:
#' `preservation_matrix_test(transform(res$pairs, Zsummary_std = z),
#' group, block)`.
#'
#' The fitted parameters are solved over the distinct clique sizes (at
#' most `S + 1` classes) and the `S` species by the fixed point
#' `x_k = k / sum_s y_s / (1 + x_k y_s)`,
#' `y_s = d_s / sum_k n_k x_k / (1 + x_k y_s)`, after fixing cells the
#' margins force (complete cliques, species in every clique or none) to
#' 1 or 0.
#'
#' @param pattern_table Output of [conservation_pattern_table()].
#' @return A list: `z`, the species x species z matrix (diagonal `NA`;
#'   `NA` where the variance is 0); `pairs`, a data frame with
#'   `reference`, `test`, `V`, `mean`, `sd` and `z`, one row per species
#'   pair; and `p`, the fitted clique x species probability matrix.
#' @references
#' Saracco F, Straka MJ, Di Clemente R, Gabrielli A, Caldarelli G,
#' Squartini T (2017). Inferring monopartite projections of bipartite
#' networks: an entropy-based approach. \emph{New Journal of Physics}
#' 19:053022. \doi{10.1088/1367-2630/aa6b38}
#' @examples
#' pt <- data.frame(
#'   clique_id = paste0("c", 1:4), hog = "H", classification = "x",
#'   pattern = "", n_plus = 3, n_minus = 0, pairs_sig = 3, pairs_ns = 0,
#'   pairs_untested = 0,
#'   A = c("+", "+", "+", "?"), B = c("+", "+", "?", "+"),
#'   C = c("+", "?", "+", "+"), D = c("?", "+", "+", "+")
#' )
#' bicm_species_z(pt)$pairs
#' @seealso [conservation_pattern_table()], [preservation_matrix_test()]
#' @export
bicm_species_z <- function(pattern_table) {
  m <- .cpt_plus(pattern_table)
  sp <- colnames(m)
  p <- .bicm_fit(m)
  dimnames(p) <- dimnames(m)
  cmb <- utils::combn(length(sp), 2L)
  a <- cmb[1L, ]
  b <- cmb[2L, ]
  v <- colSums(m[, a, drop = FALSE] & m[, b, drop = FALSE])
  q <- p[, a, drop = FALSE] * p[, b, drop = FALSE]
  mu <- colSums(q)
  sdv <- sqrt(colSums(q * (1 - q)))
  z <- ifelse(sdv > 0, (v - mu) / sdv, NA_real_)
  zm <- matrix(NA_real_, length(sp), length(sp), dimnames = list(sp, sp))
  zm[cbind(a, b)] <- z
  zm[cbind(b, a)] <- z
  pairs <- data.frame(
    reference = sp[a], test = sp[b], V = as.integer(v), mean = mu,
    sd = sdv, z = z, stringsAsFactors = FALSE
  )
  list(z = zm, pairs = pairs, p = p)
}
