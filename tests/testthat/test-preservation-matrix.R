# --- Fixture: a Pooideae-shaped all-pairs preservation matrix ---
# n_gen genera, each contributing one annual and one perennial, every
# choose(n, 2) contrast run in both directions with n_mod modules each.
# Zsummary_std carries the signal; q.value is deliberately saturated the way
# the real data is (many rows tied at a floor, many at exactly 1) so tests
# can prove the statistic never reads it.
make_pmt_fixture <- function(n_gen = 4L, n_mod = 4L, effect = 2,
                             sd = 0.5, seed = 42L) {
  set.seed(seed)
  gen <- LETTERS[seq_len(n_gen)]
  sp <- c(paste0(gen, "a"), paste0(gen, "p"))
  group <- stats::setNames(rep(c("annual", "perennial"), each = n_gen), sp)
  block <- stats::setNames(rep(gen, 2L), sp)
  pairs <- rcomplex::all_species_pairs(sp)

  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    both <- list(
      c(pairs$sp1[i], pairs$sp2[i]),
      c(pairs$sp2[i], pairs$sp1[i])
    )
    do.call(rbind, lapply(both, function(d) {
      disc <- group[[d[1L]]] != group[[d[2L]]]
      data.frame(
        pair_name = pairs$pair_name[i],
        module = as.character(seq_len(n_mod)),
        species = d[1L],
        reference = d[1L],
        test = d[2L],
        classification = "conserved",
        Zsummary = stats::rnorm(n_mod, 12, sd),
        Zsummary_std = stats::rnorm(n_mod, 10 - effect * disc, sd),
        # Saturated on purpose, and anti-aligned with the effect: the
        # discordant rows carry the *smaller* q-values, so a q-driven
        # statistic would report the opposite sign.
        q.value = if (disc) rep(0.00071, n_mod) else rep(1, n_mod),
        size = 20L,
        size_mapped = 15L,
        stringsAsFactors = FALSE
      )
    }))
  })
  cls <- do.call(rbind, rows)
  rownames(cls) <- NULL
  list(
    classification = cls, group = group, block = block, pairs = pairs,
    species = sp
  )
}

# The statistic, written out longhand for the two-level case.
pmt_manual_diff <- function(cls, group, block = NULL, exclude = TRUE) {
  if (!is.null(block) && exclude) {
    keep <- block[cls$reference] != block[cls$test]
    cls <- cls[keep, , drop = FALSE]
  }
  conc <- group[cls$reference] == group[cls$test]
  mean(cls$Zsummary_std[conc]) - mean(cls$Zsummary_std[!conc])
}


# A three-level trait on nine species: enough labellings (9!/(3!^3) = 1680)
# that the g! renamings of each labelling, which the dispersion cannot tell
# apart, still leave a floor below alpha.
make_pmt_multilevel <- function(n_per = 3L, n_mod = 2L, z_of = NULL) {
  # a-vs-c pairs are much less preserved; nothing else differs. A constant
  # effect would make the dispersion 0 under every labelling, which is a
  # degenerate null rather than a null with no signal.
  if (is.null(z_of)) {
    z_of <- function(a, b) if (identical(sort(c(a, b)), c("a", "c"))) 2 else 10
  }
  lev <- c("a", "b", "c")
  sp <- paste0(rep(lev, each = n_per), seq_len(n_per))
  grp <- stats::setNames(rep(lev, each = n_per), sp)
  pairs <- rcomplex::all_species_pairs(sp)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    both <- list(
      c(pairs$sp1[i], pairs$sp2[i]),
      c(pairs$sp2[i], pairs$sp1[i])
    )
    do.call(rbind, lapply(both, function(d) {
      data.frame(
        pair_name = pairs$pair_name[i],
        module = as.character(seq_len(n_mod)),
        reference = d[1L],
        test = d[2L],
        Zsummary = 10,
        Zsummary_std = z_of(grp[[d[1L]]], grp[[d[2L]]]),
        q.value = 1,
        stringsAsFactors = FALSE
      )
    }))
  })
  cls <- do.call(rbind, rows)
  rownames(cls) <- NULL
  list(classification = cls, group = grp, species = sp)
}


# A three-level trait with ragged module counts per contrast and a distinct
# mean for every class, the three concordant ones included. Both the row
# weighting and the pooling of the concordant classes change the statistic
# here; on a balanced fixture with one concordant mean neither does, so a
# statistic that dropped either would go unnoticed.
make_pmt_ragged <- function(seed = 5L) {
  lev <- c("a", "b", "c")
  sp <- paste0(rep(lev, each = 3L), 1:3)
  grp <- stats::setNames(rep(lev, each = 3L), sp)
  z_of <- c(aa = 10, bb = 11, cc = 9, ab = 7, ac = 3, bc = 8)
  pairs <- rcomplex::all_species_pairs(sp)
  set.seed(seed)
  n_mod <- sample.int(8L, nrow(pairs), replace = TRUE)
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    both <- list(
      c(pairs$sp1[i], pairs$sp2[i]),
      c(pairs$sp2[i], pairs$sp1[i])
    )
    do.call(rbind, lapply(both, function(d) {
      k <- paste(sort(c(grp[[d[1L]]], grp[[d[2L]]])), collapse = "")
      data.frame(
        reference = d[1L], test = d[2L],
        module = as.character(seq_len(n_mod[i])),
        Zsummary_std = z_of[[k]], q.value = 1,
        stringsAsFactors = FALSE
      )
    }))
  })
  list(classification = do.call(rbind, rows), group = grp, species = sp)
}


# Four species over a shared gene namespace, so one 1:1 ortholog table serves
# every contrast. Small and real: this is the actual preservation pipeline,
# used to prove preservation_matrix_test() consumes what preservation_paired()
# emits rather than what the fixtures above imitate.
make_pmt_real <- function(n_mod = 3L, per = 30L, n_gene = 120L,
                          n_samp = 30L) {
  gene <- sprintf("G%04d", seq_len(n_gene))
  set.seed(77)
  loadings <- lapply(seq_len(n_mod), function(k) {
    l <- stats::rlnorm(per, 0, 0.9)
    l / max(l)
  })
  expr <- function(seed) {
    set.seed(seed)
    e <- matrix(stats::rnorm(n_gene * n_samp), n_gene, n_samp)
    for (k in seq_len(n_mod)) {
      f <- stats::rnorm(n_samp)
      rows <- ((k - 1L) * per + 1L):(k * per)
      for (j in seq_along(rows)) {
        lam <- loadings[[k]][j]
        e[rows[j], ] <- lam * f +
          stats::rnorm(n_samp, sd = sqrt(max(1e-6, 1 - lam^2)))
      }
    }
    rownames(e) <- gene
    e
  }
  sp <- c("A", "B", "C", "D")
  nets <- lapply(seq_along(sp), function(i) {
    compute_network(expr(100L + i), density = 0.05, sparse = FALSE)
  })
  names(nets) <- sp
  mod_genes <- lapply(seq_len(n_mod), function(k) {
    gene[((k - 1L) * per + 1L):(k * per)]
  })
  mods <- lapply(nets, function(net) {
    mem <- stats::setNames(rep(NA_integer_, n_gene), gene)
    for (k in seq_along(mod_genes)) mem[mod_genes[[k]]] <- k
    mem <- mem[!is.na(mem)]
    list(
      modules = mem, module_genes = split(names(mem), mem),
      n_modules = length(mod_genes)
    )
  })
  list(
    modules = mods, networks = nets, species = sp,
    orthologs = data.frame(
      Species1 = gene, Species2 = gene,
      hog = paste0("H", seq_len(n_gene)), stringsAsFactors = FALSE
    ),
    group = stats::setNames(
      c("annual", "annual", "perennial", "perennial"), sp
    ),
    block = stats::setNames(c("g1", "g2", "g1", "g2"), sp)
  )
}


test_that("all_species_pairs builds every unordered contrast", {
  p <- rcomplex::all_species_pairs(c("BDIS", "BSYL", "HVUL", "HJUB"))

  expect_s3_class(p, "data.frame")
  expect_equal(nrow(p), choose(4L, 2L))
  expect_equal(names(p), c("sp1", "sp2", "pair_name"))
  expect_equal(p$pair_name, paste(p$sp1, p$sp2, sep = "."))
  # Unordered: no contrast appears twice in either orientation.
  key <- paste(pmin(p$sp1, p$sp2), pmax(p$sp1, p$sp2))
  expect_equal(anyDuplicated(key), 0L)
  expect_true(all(p$sp1 != p$sp2))
})


test_that("all_species_pairs reads a named vector as its names", {
  grp <- c(BDIS = "annual", BSYL = "perennial", HVUL = "annual")

  expect_equal(
    rcomplex::all_species_pairs(grp),
    rcomplex::all_species_pairs(names(grp))
  )
  expect_error(rcomplex::all_species_pairs("BDIS"), "at least two")
  expect_error(rcomplex::all_species_pairs(c("A", "A", "B")), "unique")
  expect_error(rcomplex::all_species_pairs(c("A", NA)), "NA or empty")
})


test_that("preservation_matrix_test returns the documented structure", {
  fix <- make_pmt_fixture()

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  expect_type(res, "list")
  expect_true(all(c(
    "observed", "statistic", "form", "class_means",
    "rows_per_pair", "free", "blocked", "p_free",
    "p_blocked", "saturation", "n_rows", "n_pairs",
    "n_excluded"
  ) %in% names(res)))
  expect_equal(res$form, "difference")
  expect_equal(res$statistic, "zsummary")
  # 4 genera give 8 species and choose(8, 2) = 28 contrasts, of which the 4
  # within-genus ones are excluded.
  expect_equal(res$n_pairs, 24L)
  expect_equal(res$n_excluded, 4L * 2L * 4L)
  expect_equal(res$n_rows, 24L * 2L * 4L)
  expect_s3_class(res$class_means, "data.frame")
  expect_equal(
    res$class_means$class,
    c("concordant", "annual vs perennial")
  )
})


test_that("the binary statistic is exactly the difference of means", {
  fix <- make_pmt_fixture()

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  expect_equal(
    res$observed,
    pmt_manual_diff(fix$classification, fix$group, fix$block)
  )
  # Positive means the discordant pairs are the less preserved ones, which
  # is the sign the biological hypothesis predicts.
  expect_gt(res$observed, 0)
  conc <- res$class_means$class == "concordant"
  expect_equal(
    res$observed,
    res$class_means$mean_z[conc] - res$class_means$mean_z[!conc]
  )
})


test_that("the free and within-block label spaces are the right size", {
  fix <- make_pmt_fixture()

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  # choose(8, 4) = 70 free labellings against 2^4 = 16 within-genus ones:
  # the whole reason the all-pairs matrix is worth running, since a floor of
  # 1/70 = 0.0143 can reach alpha = 0.05 and 1/16 = 0.0625 cannot.
  expect_true(res$free$exact)
  expect_equal(res$free$n_labellings, 70)
  expect_equal(res$free$n_scored, 70L)
  expect_equal(res$free$p_min, 1 / 70)
  expect_true(res$blocked$exact)
  expect_equal(res$blocked$n_labellings, 16)
  expect_equal(res$blocked$p_min, 1 / 16)
  # A binary trait's global flip leaves every pair's concordance alone, so
  # it reproduces the statistic exactly and the floor is never 1/n.
  expect_gte(res$free$n_tied_max, 2L)
  expect_equal(res$free$p_attainable, res$free$n_tied_max / 70)
})


test_that("the observed labelling is a point of both enumerated nulls", {
  fix <- make_pmt_fixture()

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  # Exactness rests on this: the enumeration contains the truth, so the
  # p-value needs no +1 and can never be 0.
  expect_true(any(abs(res$free$null_distribution - res$observed) < 1e-9))
  expect_true(any(abs(res$blocked$null_distribution - res$observed) < 1e-9))
  expect_gte(res$free$p_value, res$free$p_min)
  expect_lte(res$free$p_value, 1)
  expect_equal(res$p_free, res$free$p_value)
  expect_equal(res$p_blocked, res$blocked$p_value)
  at_or_above <- sum(res$free$null_distribution >= res$observed - 1e-9)
  expect_equal(res$free$rank, at_or_above)
})


test_that("a planted signal ranks near the top of both nulls", {
  fix <- make_pmt_fixture(effect = 4, sd = 0.3)

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  # Only the global flip should tie the truth, so the observed labelling is
  # at the maximum of both spaces.
  expect_equal(res$free$p_value, res$free$p_attainable)
  expect_equal(res$blocked$p_value, res$blocked$p_attainable)
  expect_lt(res$p_free, 0.05)
  # The two nulls agreeing is the diagnostic that the larger label space is
  # not buying resolution by discarding phylogenetic control.
  expect_lt(res$p_blocked, 0.2)
})


test_that("the statistic reads Zsummary_std and never the q-values", {
  fix <- make_pmt_fixture()
  cls <- fix$classification

  base <- suppressWarnings(preservation_matrix_test(
    cls, fix$group,
    block = fix$block
  ))
  # q.value is anti-aligned with the effect in the fixture, so a q-driven
  # statistic would flip the sign. Scrambling it, then dropping it, must
  # leave everything but $saturation untouched.
  set.seed(7)
  scrambled <- cls
  scrambled$q.value <- sample(scrambled$q.value)
  res_scr <- suppressWarnings(preservation_matrix_test(
    scrambled, fix$group,
    block = fix$block
  ))
  no_q <- cls[, setdiff(names(cls), "q.value"), drop = FALSE]
  res_noq <- suppressWarnings(preservation_matrix_test(
    no_q, fix$group,
    block = fix$block
  ))

  expect_equal(res_scr$observed, base$observed)
  expect_equal(res_noq$observed, base$observed)
  expect_equal(res_noq$p_free, base$p_free)
  expect_equal(res_noq$free$null_distribution, base$free$null_distribution)
  expect_equal(res_noq$saturation$n_tests, 0L)
})


test_that("saturation reports the resolution of the supplied q-values", {
  fix <- make_pmt_fixture()

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  q <- fix$classification$q.value
  expect_equal(res$saturation$n_tests, length(q))
  expect_equal(res$saturation$n_distinct, 2L)
  expect_equal(res$saturation$q_floor, min(q))
  expect_equal(res$saturation$n_at_floor, sum(q == min(q)))
  expect_equal(res$saturation$n_at_one, sum(q >= 1))
  # The multiplicity population is every test the matrix ran, not the subset
  # this statistic averages: n_tests counts the rows as supplied, while
  # n_rows is what survived the NA and within-block drops.
  expect_gt(res$saturation$n_tests, res$n_rows)
})


test_that("saturation reports the global BH floor the matrix would need", {
  fix <- make_pmt_fixture()
  n_q <- nrow(fix$classification)

  res <- preservation_matrix_test(fix$classification, fix$group,
    n_perm_pres = 20000
  )

  # q.value from preservation_paired() is corrected per contrast, over a
  # dozen or two modules. Read as one all-pairs analysis the population is
  # every module-direction, and BH over n_tests of them cannot reach below
  # n_tests / (n_perm + 1) whatever the data says.
  expect_equal(res$saturation$p_min_pres, 1 / 20001)
  expect_equal(res$saturation$q_floor_global, n_q / 20001)
  expect_lt(res$saturation$q_floor_global, 0.05)

  # At a low n_perm no global correction is available at all, and that is
  # worth saying out loud rather than leaving in a list element.
  w <- capture_warnings(
    low <- preservation_matrix_test(fix$classification, fix$group,
      n_perm_pres = 2000
    )
  )
  expect_true(any(grepl("global Benjamini-Hochberg", w)))
  expect_equal(low$saturation$q_floor_global, n_q / 2001)
  expect_gt(low$saturation$q_floor_global, 0.05)

  # Unsupplied, the field is NA rather than a guess, and nothing warns.
  plain <- preservation_matrix_test(fix$classification, fix$group)
  expect_true(is.na(plain$saturation$p_min_pres))
  expect_true(is.na(plain$saturation$q_floor_global))
  expect_error(
    preservation_matrix_test(fix$classification, fix$group,
      n_perm_pres = 0
    ),
    "n_perm_pres must be"
  )
})


test_that("statistic = zsummary_raw switches the effect column", {
  fix <- make_pmt_fixture()
  cls <- fix$classification
  # Give the raw column the opposite signal so the two cannot coincide.
  disc <- fix$group[cls$reference] != fix$group[cls$test]
  cls$Zsummary <- 10 + 3 * disc

  res <- suppressWarnings(preservation_matrix_test(
    cls, fix$group,
    block = fix$block, statistic = "zsummary_raw"
  ))

  expect_equal(res$statistic, "zsummary_raw")
  expect_equal(res$observed, -3)
  expect_error(
    preservation_matrix_test(cls[, setdiff(names(cls), "Zsummary")],
      fix$group,
      statistic = "zsummary_raw"
    ),
    "missing columns"
  )
})


test_that("within-block rows are excluded and the exclusion is fixed", {
  fix <- make_pmt_fixture()

  kept <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block,
    exclude_within_block = FALSE
  ))
  dropped <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))

  expect_equal(kept$n_excluded, 0L)
  expect_equal(kept$n_pairs, 28L)
  expect_equal(
    kept$observed,
    pmt_manual_diff(fix$classification, fix$group, fix$block,
      exclude = FALSE
    )
  )
  # Every within-genus pair is trait-discordant in this design, so keeping
  # them changes the discordant mean and hence the statistic.
  expect_false(isTRUE(all.equal(kept$observed, dropped$observed)))
  # The label space is a property of the species, not of which rows survive.
  expect_equal(kept$free$n_labellings, dropped$free$n_labellings)
  expect_true(all(!dropped$rows_per_pair$same_block))
})


test_that("block = NULL runs the free null only", {
  fix <- make_pmt_fixture()

  res <- preservation_matrix_test(fix$classification, fix$group)

  expect_null(res$blocked)
  expect_true(is.na(res$p_blocked))
  expect_null(res$block)
  expect_true(all(is.na(res$rows_per_pair$same_block)))
  # Without a block nothing is excluded, whatever exclude_within_block says.
  expect_equal(res$n_excluded, 0L)
  expect_equal(res$n_pairs, 28L)
})


test_that("more than two trait levels give a dispersion statistic", {
  fix <- make_pmt_multilevel()

  res <- preservation_matrix_test(fix$classification, fix$group)

  expect_equal(res$form, "dispersion")
  # One pooled concordant class plus one per unordered pair of levels.
  expect_equal(
    sort(res$class_means$class),
    sort(c("concordant", "a vs b", "a vs c", "b vs c"))
  )
  # The row-weighted standard deviation of the class means, longhand.
  w <- res$class_means$n
  m <- res$class_means$mean_z
  mbar <- sum(w * m) / sum(w)
  expect_equal(res$observed, sqrt(sum(w * (m - mbar)^2) / sum(w)))
  expect_gte(res$observed, 0)
  # 9!/(3!^3) = 1680 labellings.
  expect_equal(res$free$n_labellings, 1680)
})


test_that("dispersion is one-sided upward and detects planted structure", {
  fix <- make_pmt_multilevel()

  res <- preservation_matrix_test(fix$classification, fix$group)

  expect_gt(res$observed, 0)
  expect_lt(res$p_free, 0.05)
  # The dispersion cannot tell the level names apart, so every one of the
  # 3! = 6 renamings of a labelling scores identically and the attainable
  # floor is 6 / 1680, not 1 / 1680. A three-level design on six species
  # would have 90 labellings and a floor of 0.067 -- above alpha, with no
  # signal strong enough to fix it.
  expect_gte(res$free$n_tied_max, 6L)
  expect_equal(res$free$p_attainable, res$free$n_tied_max / 1680)
})


test_that("rows with no measured effect are dropped", {
  fix <- make_pmt_fixture()
  cls <- fix$classification
  extra <- cls[seq_len(6L), , drop = FALSE]
  extra$classification <- "untested"
  extra$Zsummary_std <- NA_real_
  extra$module <- paste0("x", seq_len(6L))

  base <- suppressWarnings(preservation_matrix_test(
    cls, fix$group,
    block = fix$block
  ))
  with_na <- suppressWarnings(preservation_matrix_test(
    rbind(cls, extra), fix$group,
    block = fix$block
  ))

  expect_equal(with_na$n_rows, base$n_rows)
  expect_equal(with_na$observed, base$observed)
})


test_that("a species outside group is dropped with a warning", {
  fix <- make_pmt_fixture()
  cls <- fix$classification
  cls$reference[cls$reference == "Ap"] <- "ZZZ"
  cls$test[cls$test == "Ap"] <- "ZZZ"

  expect_warning(
    res <- preservation_matrix_test(cls, fix$group),
    "no group entry"
  )
  expect_false("ZZZ" %in% res$species)
  expect_equal(length(res$species), 7L)
})


test_that("a space too large to enumerate is sampled reproducibly", {
  fix <- make_pmt_fixture()

  set.seed(11)
  a <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block,
    n_perm = 200L, enum_max = 5L
  ))
  set.seed(11)
  b <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block,
    n_perm = 200L, enum_max = 5L
  ))

  expect_false(a$free$exact)
  expect_equal(a$free$n_scored, 200L)
  expect_equal(a$free$p_min, 1 / 201)
  # Sampled nulls do not contain the observed labelling, hence (r + 1)/(m + 1)
  # and a p-value that can never be 0.
  expect_gte(a$free$p_value, 1 / 201)
  expect_equal(a$free$null_distribution, b$free$null_distribution)
  expect_equal(a$p_free, b$p_free)
  # The label space is still reported at its true size, not at the draw count.
  expect_equal(a$free$n_labellings, 70)
})


test_that("n_perm is ignored, with a message, when a null is enumerated", {
  fix <- make_pmt_fixture()

  expect_message(
    suppressWarnings(preservation_matrix_test(
      fix$classification, fix$group,
      n_perm = 500L
    )),
    "enumerated exactly"
  )
  # The default leaves n_perm unspoken, so an enumerated null says nothing.
  expect_silent(preservation_matrix_test(fix$classification, fix$group))
})


test_that("an unreachable p-value floor is warned about", {
  fix <- make_pmt_fixture()

  # Four genera give a within-block space of 16 with the global flip tied,
  # so the blocked floor is 2/16 = 0.125 and alpha = 0.05 is unreachable.
  expect_warning(
    preservation_matrix_test(fix$classification, fix$group,
      block = fix$block
    ),
    "smallest attainable p-value"
  )
  # The free space of 70 reaches it, so that null must not warn.
  expect_silent(preservation_matrix_test(fix$classification, fix$group))
})


test_that("preservation_matrix_test validates its inputs", {
  fix <- make_pmt_fixture()
  cls <- fix$classification

  expect_error(
    preservation_matrix_test(as.list(cls), fix$group),
    "must be a data frame"
  )
  expect_error(
    preservation_matrix_test(
      cls[, setdiff(names(cls), "Zsummary_std")],
      fix$group
    ),
    "missing columns"
  )
  expect_error(
    preservation_matrix_test(cls, unname(fix$group)),
    "named vector"
  )
  expect_error(
    preservation_matrix_test(cls, fix$group, block = unname(fix$block)),
    "named vector"
  )
  expect_error(
    preservation_matrix_test(cls, fix$group,
      block = fix$block[1:3]
    ),
    "block missing entries"
  )
  expect_error(
    preservation_matrix_test(cls, fix$group, n_perm = 0),
    "positive number"
  )
  expect_error(
    preservation_matrix_test(cls, fix$group, enum_max = -1),
    "positive number"
  )
  # A trait taking one value over the tested species has no contrast at all.
  flat <- stats::setNames(
    rep("annual", length(fix$group)),
    names(fix$group)
  )
  expect_error(preservation_matrix_test(cls, flat), "single value")
  # Nothing measured anywhere.
  none <- cls
  none$Zsummary_std <- NA_real_
  expect_error(preservation_matrix_test(none, fix$group), "no usable rows")
})


test_that("uneven blocks give the product of per-block label spaces", {
  sp <- paste0("s", 1:8)
  grp <- stats::setNames(c("a", "a", "p", "p", "a", "a", "p", "p"), sp)
  # Blocks of 1, 3 and 4, with the trait split unevenly across them: the
  # restricted space is the product of per-block multinomials, not 2^blocks,
  # and a singleton block contributes exactly one labelling.
  blk <- stats::setNames(c("B1", rep("B2", 3L), rep("B3", 4L)), sp)
  pairs <- rcomplex::all_species_pairs(sp)
  cls <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    both <- list(
      c(pairs$sp1[i], pairs$sp2[i]),
      c(pairs$sp2[i], pairs$sp1[i])
    )
    do.call(rbind, lapply(both, function(d) {
      data.frame(
        reference = d[1L], test = d[2L], module = "1",
        Zsummary_std = 10 + (grp[[d[1L]]] == grp[[d[2L]]]),
        stringsAsFactors = FALSE
      )
    }))
  }))

  res <- suppressWarnings(preservation_matrix_test(cls, grp, block = blk))

  expect_equal(res$free$n_labellings, choose(8, 4))
  expect_equal(res$blocked$n_labellings, 1 * 3 * 6)
  expect_true(any(abs(res$blocked$null_distribution - res$observed) < 1e-9))
  # 3 + 6 within-block pairs, both directions, one module each.
  expect_equal(res$n_excluded, 18L)
  expect_equal(res$n_pairs, choose(8L, 2L) - 9L)
})


test_that("the sampled within-block null still permutes within blocks", {
  fix <- make_pmt_fixture()

  exact <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))
  set.seed(4)
  samp <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block,
    n_perm = 300L, enum_max = 5L
  ))

  expect_false(samp$blocked$exact)
  # Every draw must be one of the 16 labellings the blocks allow, so its
  # statistic must be one the exact null already holds. A free permutation
  # in the sampled branch -- the documented route in, via enum_max -- would
  # leave that set at once and silently swap the phylogenetically
  # controlled null for the free one.
  in_space <- vapply(samp$blocked$null_distribution, function(v) {
    any(abs(v - exact$blocked$null_distribution) < 1e-8)
  }, logical(1))
  expect_true(all(in_space))
  # Non-vacuous: the free space does hold values the blocks cannot reach.
  free_only <- vapply(exact$free$null_distribution, function(v) {
    !any(abs(v - exact$blocked$null_distribution) < 1e-8)
  }, logical(1))
  expect_gt(sum(free_only), 0L)
})


test_that("the dispersion is row-weighted, not a plain sd of the means", {
  fix <- make_pmt_ragged()

  res <- preservation_matrix_test(fix$classification, fix$group)

  w <- res$class_means$n
  m <- res$class_means$mean_z
  mbar <- sum(w * m) / sum(w)
  weighted <- sqrt(sum(w * (m - mbar)^2) / sum(w))
  unweighted <- sqrt(mean((m - mean(m))^2))
  expect_equal(res$observed, weighted)
  # Ragged module counts are what makes the two differ at all: on a
  # fixture with a constant n_mod every class carries equal weight and the
  # weighting the documentation promises is untestable.
  expect_false(isTRUE(all.equal(weighted, unweighted)))
  expect_false(isTRUE(all.equal(res$observed, unweighted)))
})


test_that("concordant classes are pooled inside the statistic itself", {
  fix <- make_pmt_ragged()
  cls <- fix$classification
  grp <- fix$group

  res <- preservation_matrix_test(cls, grp)

  conc <- grp[cls$reference] == grp[cls$test]
  key <- ifelse(conc, paste0("C", grp[cls$reference]),
    paste(
      pmin(grp[cls$reference], grp[cls$test]),
      pmax(grp[cls$reference], grp[cls$test])
    )
  )
  parts <- split(cls$Zsummary_std, key)
  m <- vapply(parts, mean, numeric(1))
  w <- lengths(parts)
  split_stat <- sqrt(sum(w * (m - sum(w * m) / sum(w))^2) / sum(w))
  # The three concordant means differ here (10, 11, 9), so splitting them
  # into one class per level is a different number -- which is what makes
  # the pooling testable rather than a comment.
  expect_false(isTRUE(all.equal(res$observed, split_stat)))
  expect_equal(nrow(res$class_means), 4L)
  expect_equal(
    res$class_means$n[res$class_means$class == "concordant"], sum(conc)
  )
  expect_equal(
    res$class_means$mean_z[res$class_means$class == "concordant"],
    mean(cls$Zsummary_std[conc])
  )
})


test_that("rows_per_pair accounts for every tested row exactly once", {
  fix <- make_pmt_fixture()
  rag <- make_pmt_ragged()

  res <- suppressWarnings(preservation_matrix_test(
    fix$classification, fix$group,
    block = fix$block
  ))
  ragged <- preservation_matrix_test(rag$classification, rag$group)

  # The diagnostic for "one hugely partitioned species dominates" is only
  # a diagnostic if n really counts that pair's rows.
  expect_equal(sum(res$rows_per_pair$n), res$n_rows)
  expect_equal(nrow(res$rows_per_pair), res$n_pairs)
  expect_true(all(res$rows_per_pair$n == 8L))
  expect_equal(sum(ragged$rows_per_pair$n), ragged$n_rows)
  expect_gt(length(unique(ragged$rows_per_pair$n)), 1L)
})


test_that("a self-contrast row is dropped, not averaged in", {
  fix <- make_pmt_fixture()
  cls <- fix$classification
  self <- cls[seq_len(4L), , drop = FALSE]
  self$test <- self$reference
  # Trait-concordant by construction and wildly off scale, so keeping it
  # would move the concordant mean and the statistic with it.
  self$Zsummary_std <- 1000

  # Run without a block: a species shares its own block, so the within-block
  # exclusion would drop these rows for the wrong reason and hide whether
  # the self-contrast guard does anything.
  base <- preservation_matrix_test(cls, fix$group)
  with_self <- preservation_matrix_test(rbind(cls, self), fix$group)

  expect_equal(with_self$n_rows, base$n_rows)
  expect_equal(with_self$observed, base$observed)
  pairs_seen <- with_self$rows_per_pair
  expect_false(any(pairs_seen$sp1 == pairs_seen$sp2))
})


test_that("a species named twice in group or block is refused", {
  fix <- make_pmt_fixture()
  cls <- fix$classification

  # group[species] takes the first match, so a duplicate silently assigns a
  # species the wrong trait and the null then permutes the wrong design.
  dup_group <- c(fix$group, Aa = "perennial")
  expect_error(
    preservation_matrix_test(cls, dup_group),
    "group names a species more than once"
  )
  expect_error(
    preservation_matrix_test(cls, fix$group,
      block = c(fix$block, Aa = "Z")
    ),
    "block names a species more than once"
  )
})


test_that("an unbalanced multilevel design ties fewer than g! labellings", {
  # Counts 2/2/1 over five species: 5! / (2! 2! 1!) = 30 labellings, of
  # which only the two renamings that swap the two size-2 levels stay in
  # the space. A renaming has to map each level onto one of the same size,
  # so the tie count is prod(factorial(table(table(labels)))) -- g! = 6
  # only when the design is balanced.
  sp <- c("a1", "a2", "b1", "b2", "c1")
  grp <- stats::setNames(c("a", "a", "b", "b", "c"), sp)
  z_of <- c(aa = 10, bb = 11, cc = 9, ab = 7, ac = 3, bc = 8)
  pairs <- rcomplex::all_species_pairs(sp)
  cls <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    both <- list(
      c(pairs$sp1[i], pairs$sp2[i]),
      c(pairs$sp2[i], pairs$sp1[i])
    )
    do.call(rbind, lapply(both, function(d) {
      k <- paste(sort(c(grp[[d[1L]]], grp[[d[2L]]])), collapse = "")
      data.frame(
        reference = d[1L], test = d[2L], module = "1",
        Zsummary_std = z_of[[k]], stringsAsFactors = FALSE
      )
    }))
  }))

  res <- suppressWarnings(preservation_matrix_test(cls, grp))

  expect_equal(res$free$n_labellings, 30)
  expect_equal(res$free$n_tied_max, 2L)
  expect_equal(res$free$p_attainable, 2 / 30)
  expect_equal(
    res$free$n_tied_max,
    as.integer(prod(factorial(table(table(grp)))))
  )
})


test_that("a sampled null reports its floor over draws, not labellings", {
  fix <- make_pmt_fixture()

  set.seed(11)
  sampled <- capture_warnings(
    samp <- preservation_matrix_test(fix$classification, fix$group,
      block = fix$block, n_perm = 200L, enum_max = 5L
    )
  )
  enumerated <- capture_warnings(
    preservation_matrix_test(fix$classification, fix$group,
      block = fix$block
    )
  )

  # A drawn labelling can repeat, so the tie count and p_min of a sampled
  # null describe the draws and not the space: reporting them as
  # labellings would claim a design property the run never measured.
  expect_true(any(grepl("draws share the maximum", sampled)))
  expect_false(any(grepl("labellings share the maximum", sampled)))
  expect_true(any(grepl("labellings share the maximum", enumerated)))
  expect_gt(samp$blocked$n_tied_max, 2L)
  expect_equal(samp$blocked$p_min, 1 / 201)
  # The space itself is still reported at its true size.
  expect_equal(samp$blocked$n_labellings, 16)
})


test_that("preservation_matrix_test consumes a real preservation_paired run", {
  fx <- make_pmt_real()
  pairs <- rcomplex::all_species_pairs(fx$species)

  paired <- suppressWarnings(preservation_paired(
    fx$modules, fx$networks, fx$orthologs, pairs,
    n_perm = 50L, seed = 1
  ))
  res <- suppressWarnings(preservation_matrix_test(
    paired$classification, fx$group,
    block = fx$block, n_perm_pres = 50L
  ))

  # The column names are the contract: reference, test, Zsummary_std and
  # q.value as preservation_paired() writes them, not as a fixture
  # imitates them.
  expect_equal(res$form, "difference")
  expect_equal(res$statistic, "zsummary")
  expect_true(is.finite(res$observed))
  # 4 species, 6 contrasts x 2 directions x 3 modules, of which the two
  # within-block contrasts are excluded.
  expect_equal(nrow(paired$classification), 36L)
  expect_equal(res$n_pairs, 4L)
  expect_equal(res$n_rows, 24L)
  expect_equal(res$free$n_labellings, choose(4, 2))
  expect_equal(res$blocked$n_labellings, 4)
  expect_equal(
    sort(res$class_means$class),
    sort(c("concordant", "annual vs perennial"))
  )
  expect_gte(res$p_free, res$free$p_attainable)
  expect_lte(res$p_free, 1)
  expect_equal(res$saturation$n_tests, sum(is.finite(
    paired$classification$q.value
  )))
})
