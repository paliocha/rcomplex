# --- Synthetic fixtures for tag_permutation ---

make_tag_perm_fixtures <- function() {
  # 3 pairs, 6 species, 2 traits (3 annual, 3 perennial)
  group <- c(A1 = "annual", P1 = "perennial",
             A2 = "annual", P2 = "perennial",
             A3 = "annual", P3 = "perennial")

  pairs <- data.frame(
    sp1 = c("A1", "A2", "A3"),
    sp2 = c("P1", "P2", "P3"),
    pair_name = c("pair1", "pair2", "pair3"),
    stringsAsFactors = FALSE
  )

  # Modules: each species gets 2 modules (1, 2)
  make_modules <- function(sp, genes) {
    membership <- stats::setNames(rep(c(1L, 2L), each = length(genes) / 2),
                                  genes)
    list(
      modules = membership,
      module_genes = split(names(membership), membership),
      n_modules = 2L,
      modularity = 0.3,
      graph = NULL,
      method = "leiden",
      params = list()
    )
  }

  modules <- list(
    A1 = make_modules("A1", paste0("A1_g", 1:6)),
    P1 = make_modules("P1", paste0("P1_g", 1:6)),
    A2 = make_modules("A2", paste0("A2_g", 1:6)),
    P2 = make_modules("P2", paste0("P2_g", 1:6)),
    A3 = make_modules("A3", paste0("A3_g", 1:6)),
    P3 = make_modules("P3", paste0("P3_g", 1:6))
  )

  # Orthologs: map genes to HOGs (consistent across species)
  orthologs <- data.frame(
    Species1 = c(paste0("A1_g", 1:6), paste0("A2_g", 1:6), paste0("A3_g", 1:6)),
    Species2 = c(paste0("P1_g", 1:6), paste0("P2_g", 1:6), paste0("P3_g", 1:6)),
    hog = rep(paste0("HOG", 1:6), 3),
    stringsAsFactors = FALSE
  )

  # Classification in the preservation vocabulary: one row per module per
  # DIRECTION. Module 1 of each annual is diverged in its perennial partner;
  # module 2 of each perennial is diverged in its annual partner. One
  # "moderate" and one "untested" row prove neither contributes.
  annuals <- c("A1", "A2", "A3")
  perennials <- c("P1", "P2", "P3")
  classification <- data.frame(
    pair_name = rep(c("pair1", "pair2", "pair3"), each = 4),
    module = rep(c("1", "2", "1", "2"), 3),
    reference = rep(annuals, each = 4),
    test = rep(perennials, each = 4),
    classification = rep(c("diverged", "conserved",
                           "conserved", "diverged"), 3),
    stringsAsFactors = FALSE
  )
  # The reverse direction: rows 3 and 4 of each pair describe the perennial's
  # own modules tested in the annual.
  rev_rows <- rep(c(FALSE, FALSE, TRUE, TRUE), 3)
  classification$reference[rev_rows] <- rep(perennials, each = 2)
  classification$test[rev_rows] <- rep(annuals, each = 2)
  # Neither of these may contribute.
  classification$classification[2] <- "moderate"
  classification$classification[3] <- "untested"

  list(classification = classification, modules = modules,
       orthologs = orthologs, pairs = pairs, group = group)
}


test_that("tag_permutation returns correct structure", {
  fix <- make_tag_perm_fixtures()

  result <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 100L, min_recurrence = 2L
  ))

  expect_type(result, "list")
  expect_true(all(c("observed", "null_distribution", "p_value",
                     "p_min", "exact", "n_swappable",
                     "recurrence_table", "target_group",
                     "min_recurrence", "n_perm") %in% names(result)))
  expect_type(result$observed, "integer")
  # 3 swappable pairs -> the null is 2^3 labellings, enumerated, and
  # n_perm = 100 is ignored rather than sampled.
  expect_true(result$exact)
  expect_equal(result$n_swappable, 3L)
  expect_length(result$null_distribution, 8L)
  expect_equal(result$n_perm, 8L)
  expect_equal(result$p_min, 1 / 8)
  expect_true(result$p_value >= result$p_min && result$p_value <= 1)
  expect_s3_class(result$recurrence_table, "data.frame")
  expect_equal(result$target_group, "annual")
  expect_equal(result$min_recurrence, 2L)
})


test_that("tag_permutation detects known parallel signal", {
  fix <- make_tag_perm_fixtures()

  # All 3 annual species have module 1 as species-specific.
  # Module 1 genes are A*_g1:g3, mapping to HOG1:HOG3.
  # So HOG1, HOG2, HOG3 each recur in 3 pairs.
  result <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 100L, min_recurrence = 2L
  ))

  expect_equal(result$observed, 3L)  # HOG1, HOG2, HOG3

  # Recurrence table should have 3 HOGs with n_pairs = 3
  recurring <- result$recurrence_table[result$recurrence_table$n_pairs >= 2, ]
  expect_equal(nrow(recurring), 3L)
})


test_that("tag_permutation handles min_recurrence thresholds", {
  fix <- make_tag_perm_fixtures()

  # min_recurrence = 3: all 3 pairs must have the HOG
  r3 <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 50L, min_recurrence = 3L
  ))
  expect_equal(r3$observed, 3L)

  # min_recurrence = 4: impossible with 3 pairs
  r4 <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 50L, min_recurrence = 4L
  ))
  expect_equal(r4$observed, 0L)
})


test_that("tag_permutation pair exclusion: both sides same trait", {
  fix <- make_tag_perm_fixtures()

  # If we make all species "annual", no pair has exactly one annual
  all_annual <- stats::setNames(rep("annual", 6),
                                names(fix$group))

  result <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    all_annual, target_group = "annual",
    n_perm = 50L, min_recurrence = 2L
  ))

  expect_equal(result$observed, 0L)
  # No pair has two different labels, so nothing is swappable and the
  # null is the single observed labelling.
  expect_equal(result$n_swappable, 0L)
  expect_length(result$null_distribution, 1L)
  expect_equal(result$p_min, 1)
})


test_that("tag_permutation works with > 2 trait values", {
  fix <- make_tag_perm_fixtures()

  # 3 trait values, unequal frequencies
  group3 <- c(A1 = "annual", P1 = "perennial",
              A2 = "annual", P2 = "biennial",
              A3 = "biennial", P3 = "perennial")

  result <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    group3, target_group = "annual",
    n_perm = 100L, min_recurrence = 2L
  ))

  expect_type(result, "list")
  # Only pair1 and pair2 have exactly one "annual" species
  # pair3 has biennial + perennial, so it doesn't contribute
  # HOG1:3 from pair1, HOG1:3 from pair2 -> 3 HOGs in 2 pairs
  expect_equal(result$observed, 3L)
})


test_that("tag_permutation validates inputs", {
  fix <- make_tag_perm_fixtures()

  expect_error(
    tag_permutation(fix$classification[, -1], fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "annual"),
    "classification missing columns"
  )

  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "nonexistent"),
    "target_group.*not found"
  )

  # An old gene-overlap table must error, not return observed = 0. Only the
  # diverged rows are relabelled, so the table keeps its "conserved" rows --
  # a guard that merely looked for known levels would pass this.
  old_cls <- fix$classification
  old_cls$classification[old_cls$classification == "diverged"] <-
    "species_specific"
  expect_error(
    tag_permutation(old_cls, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "annual"),
    "retired gene-overlap vocabulary"
  )

  # A foreign vocabulary must error too, even mixed with a known level --
  # relabelled rows or another tool's output would otherwise filter to zero
  # diverged rows and return observed = 0.
  foreign <- fix$classification
  foreign$classification <- ifelse(
    foreign$classification == "diverged", "Diverged", "Conserved")
  foreign$classification[1] <- "conserved"
  expect_error(
    tag_permutation(foreign, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "annual"),
    "levels outside"
  )

  # A pair_name the classification does not carry must error rather than
  # silently yielding observed = 0.
  bad_pn <- fix$pairs
  bad_pn$pair_name <- c("nope1", "nope2", "nope3")
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    bad_pn, fix$group, "annual"),
    "no rows for pair_name"
  )

  no_side <- fix$classification
  no_side$reference <- NULL
  expect_error(
    tag_permutation(no_side, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, "annual"),
    "classification missing columns"
  )

  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, c("annual", "perennial"), "annual"),
    "named character"
  )

  bad_group <- fix$group[-1]
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, bad_group, "annual"),
    "group missing entries"
  )
})


test_that("the null stays inside the observed design", {
  # This is the property the unconditional shuffle broke. Every pair here
  # contributes a non-empty HOG set under the observed labelling, and a
  # within-pair swap cannot change that: each pair keeps one annual and
  # one perennial in every draw, so no draw can collapse to a design with
  # fewer contributing pairs. Under the retired null a draw could make a
  # pair trait-concordant, which zeroes its contribution: enumerating its
  # C(6, 3) = 20 labellings on this fixture gives a statistic of exactly
  # 0 for 12 of them, because only the 2^3 = 8 all-discordant labellings
  # leave enough pairs to meet min_recurrence = 2.
  fix <- make_tag_perm_fixtures()

  result <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 1000L, min_recurrence = 2L
  ))

  # Every labelling selects one side of all 3 pairs, and every side
  # carries 3 HOGs shared across pairs, so every draw recurs 3 HOGs.
  expect_true(all(result$null_distribution == 3L))
  expect_false(any(result$null_distribution == 0L))
})


test_that("the enumerated null is the exact 2^k label space", {
  # The default fixture gives every labelling the same statistic, which
  # cannot distinguish a correct enumeration from a constant. Break the
  # symmetry so the 8 labellings produce a spread: strip P3's diverged
  # module down to a single HOG, so selecting the P3 side of pair3 leaves
  # too few shared HOGs to meet min_recurrence.
  fix <- make_tag_perm_fixtures()
  fix$modules$P3$module_genes[["2"]] <- "P3_g4"

  result <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 8L, min_recurrence = 2L
  ))
  expect_true(result$exact)
  expect_length(result$null_distribution, 8L)
  # A real spread, so the assertions below have something to bite on.
  expect_gt(length(unique(result$null_distribution)), 1L)

  # The observed labelling is one of the enumerated points, so the exact
  # p-value needs no +1 and can never fall below 1 / 2^k.
  expect_gte(result$p_value, result$p_min)
  expect_equal(result$p_min, 1 / 8)
  expect_true(result$observed %in% result$null_distribution)
  expect_equal(result$p_value, mean(result$null_distribution >=
                                      result$observed))

  # n_perm below the label space is the one case that does NOT enumerate:
  # it falls back to sampling that many independent swap vectors, and the
  # +1 correction returns because the observed labelling is not
  # guaranteed to be among the draws.
  small <- suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 4L, min_recurrence = 2L
  ))
  expect_false(small$exact)
  expect_length(small$null_distribution, 4L)
  expect_equal(small$p_min, 1 / 5)
})


test_that("tag_permutation warns when significance is unreachable", {
  fix <- make_tag_perm_fixtures()

  expect_warning(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, target_group = "annual",
                    n_perm = 100L, min_recurrence = 2L),
    "smallest attainable p-value"
  )
})


test_that("tag_permutation requires a disjoint pairing", {
  fix <- make_tag_perm_fixtures()

  # A1 in two pairs: swapping pair1 would change pair3's labels too, so
  # the swaps are not independent and 2^k is not the support.
  bad <- fix$pairs
  bad$sp1[3] <- "A1"
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    bad, fix$group, target_group = "annual"),
    "disjoint"
  )

  same <- fix$pairs
  same$sp2[1] <- same$sp1[1]
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    same, fix$group, target_group = "annual"),
    "two distinct species"
  )
})


test_that("complementary trait values share one null distribution", {
  # Swapping every pair maps the annual statistic onto the perennial one,
  # so the two runs read the same 2^k numbers. They are not independent
  # evidence and must not be corrected as two tests.
  # Needs a fixture whose null is not constant: on the default one both
  # runs return eight copies of the same value and the multiset equality
  # holds for any implementation, symmetry or not.
  fix <- make_tag_perm_fixtures()
  fix$modules$P3$module_genes[["2"]] <- "P3_g4"

  run_for <- function(tg) {
    suppressMessages(suppressWarnings(tag_permutation(
      fix$classification,
      fix$modules,
      fix$orthologs,
      fix$pairs,
      fix$group,
      target_group = tg,
      min_recurrence = 2L
    )))
  }
  ann <- run_for("annual")
  per <- run_for("perennial")

  expect_gt(length(unique(ann$null_distribution)), 1L)
  expect_equal(sort(ann$null_distribution), sort(per$null_distribution))
  expect_true(per$observed %in% ann$null_distribution)
  # The complement of the observed labelling is the all-swapped one, and
  # its annual statistic is the observed perennial statistic.
  expect_equal(ann$null_distribution[length(ann$null_distribution)],
               per$observed)
})


test_that("the conditional null is calibrated under H0", {
  # Type I error must not exceed nominal when no trait effect exists.
  # The fixture plants a nuisance the retired null could not survive: a
  # "hot" HOG block that sits in diverged modules of every species
  # regardless of trait.
  set.seed(11)
  k <- 6L                       # 2^6 = 64 labellings, p_min = 0.0156
  sp1 <- paste0("A", seq_len(k))
  sp2 <- paste0("P", seq_len(k))
  group <- stats::setNames(rep(c("annual", "perennial"), each = k),
                           c(sp1, sp2))
  pairs <- data.frame(sp1 = sp1, sp2 = sp2,
                      pair_name = paste0("pair", seq_len(k)),
                      stringsAsFactors = FALSE)
  hogs <- paste0("HOG", seq_len(40))

  one_rep <- function() {
    # Each side of each pair draws its diverged HOG set from the same
    # distribution -- no trait effect -- with the hot block favoured.
    wt <- c(rep(6, 8), rep(1, 32))
    draw <- function() sample(hogs, 12L, prob = wt)
    sides <- lapply(seq_len(k), function(i) {
      list(annual = draw(), perennial = draw())
    })

    modules <- list()
    orth <- list()
    cls <- list()
    for (i in seq_len(k)) {
      for (side in c("annual", "perennial")) {
        sp <- if (side == "annual") sp1[i] else sp2[i]
        g <- paste0(sp, "_g", seq_along(sides[[i]][[side]]))
        memb <- stats::setNames(rep(1L, length(g)), g)
        modules[[sp]] <- list(
          modules = memb,
          module_genes = list(`1` = g),
          n_modules = 1L,
          modularity = 0.3,
          graph = NULL,
          method = "leiden",
          params = list()
        )
        orth[[sp]] <- data.frame(Species1 = g, Species2 = g,
                                 hog = sides[[i]][[side]],
                                 stringsAsFactors = FALSE)
      }
      cls[[i]] <- data.frame(
        pair_name = rep(pairs$pair_name[i], 2),
        module = c("1", "1"),
        reference = c(sp1[i], sp2[i]),
        test = c(sp2[i], sp1[i]),
        classification = c("diverged", "diverged"),
        stringsAsFactors = FALSE
      )
    }
    suppressMessages(suppressWarnings(tag_permutation(
      do.call(rbind, cls),
      modules,
      do.call(rbind, orth),
      pairs,
      group,
      target_group = "annual",
      min_recurrence = 2L
    )))
  }

  ps <- vapply(seq_len(60), function(i) one_rep()$p_value, numeric(1))
  # Exact tests are conservative on a discrete lattice, never
  # anti-conservative. The retired null measured 0.17 here.
  expect_lte(mean(ps <= 0.05), 0.05 + 3 * sqrt(0.05 * 0.95 / 60))
  expect_lte(mean(ps <= 0.10), 0.10 + 3 * sqrt(0.10 * 0.90 / 60))
})


test_that("pair_sizes exposes the exchangeability condition", {
  # The within-pair swap is exchangeable only if the target side is not
  # systematically the larger one. Simulation puts the false-positive
  # rate at 0.74 for a 13% systematic size excess with no recurrence
  # signal, so the sizes have to be visible and a clean sweep has to be
  # called out.
  fix <- make_tag_perm_fixtures()

  balanced <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual", min_recurrence = 2L
  )))
  ps <- balanced$pair_sizes
  expect_s3_class(ps, "data.frame")
  expect_equal(nrow(ps), nrow(fix$pairs))
  expect_true(all(ps$swappable))
  # Every annual side carries HOG1:3 and every perennial side HOG4:6, so
  # the sides are the same size and no asymmetry warning is due.
  expect_equal(ps$n_hogs_target, ps$n_hogs_partner)

  # Now make the annual side the larger one in all three pairs.
  skewed <- fix
  for (sp in c("P1", "P2", "P3")) {
    skewed$modules[[sp]]$module_genes[["2"]] <-
      skewed$modules[[sp]]$module_genes[["2"]][1]
  }
  # Two warnings fire here (the p_min floor and the size sweep), so
  # collect both rather than letting expect_warning swallow the first.
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(tag_permutation(
      skewed$classification, skewed$modules, skewed$orthologs,
      skewed$pairs, skewed$group, target_group = "annual",
      min_recurrence = 2L
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("larger HOG set in all 3 swappable pairs", warns)))
  expect_true(all(res$pair_sizes$n_hogs_target >
                    res$pair_sizes$n_hogs_partner))

  # A non-swappable pair has no target side, so its sizes are NA rather
  # than silently counted toward the sweep.
  group3 <- fix$group
  group3["A3"] <- "biennial"
  mixed <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    group3, target_group = "annual", min_recurrence = 2L
  )))
  expect_equal(sum(mixed$pair_sizes$swappable), 2L)
  unswappable <- mixed$pair_sizes[!mixed$pair_sizes$swappable, ]
  expect_true(all(is.na(unswappable$n_hogs_target)))
  expect_true(all(is.na(unswappable$n_hogs_partner)))
})
