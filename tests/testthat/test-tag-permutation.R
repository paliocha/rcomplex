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

  # min_recurrence = 4 is impossible with 3 contributing pairs. That used
  # to return observed = 0 with p = 1 in silence, which reads as a
  # measured absence rather than an unsatisfiable request.
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, target_group = "annual",
                    n_perm = 50L, min_recurrence = 4L),
    "exceeds the number of pairs contributing"
  )
})


test_that("tag_permutation pair exclusion: both sides same trait", {
  fix <- make_tag_perm_fixtures()

  # If we make all species "annual", no pair has exactly one annual
  all_annual <- stats::setNames(rep("annual", 6),
                                names(fix$group))

  # No pair has exactly one annual, so nothing contributes and no HOG can
  # reach two pairs -- an unsatisfiable request, now an error rather than
  # a silent observed = 0 that reads as a measured absence.
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, all_annual, target_group = "annual",
                    n_perm = 50L, min_recurrence = 2L),
    "exceeds the number of pairs contributing"
  )

  # With nothing contributing, every min_recurrence is unsatisfiable, and
  # the message says so by naming zero contributing pairs.
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, all_annual, target_group = "annual",
                    n_perm = 50L, min_recurrence = 1L),
    "contributing to the statistic [(]0[)]"
  )
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

  # statistic = "count" so the null holds recurrence counts and the
  # assertion reads directly; the default "excess" subtracts a size
  # expectation and is exercised separately.
  result <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 1000L, min_recurrence = 2L, statistic = "count"
  )))

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
  # null_distribution holds the statistic, so compare against
  # statistic_observed rather than the raw recurrence count.
  expect_true(result$statistic_observed %in% result$null_distribution)
  expect_equal(result$p_value,
               mean(result$null_distribution >=
                      result$statistic_observed))

  # Enumeration is decided on cost (at most 20 swappable pairs), not on
  # n_perm, so lowering n_perm below the label space no longer demotes
  # the null to sampling. This assertion is the reverse of what it was:
  # keying the decision on n_perm meant raising n_perm for precision
  # could flip an exact null to a sampled one.
  small <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 4L, min_recurrence = 2L
  )))
  expect_true(small$exact)
  expect_length(small$null_distribution, 8L)
  expect_equal(small$p_min, 1 / 8)
  expect_equal(small$p_value, result$p_value)

  # p_attainable is the floor after ties, which this fixture has.
  expect_equal(result$p_attainable,
               mean(result$null_distribution >=
                      max(result$null_distribution)))
  expect_gte(result$p_value, result$p_attainable)
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


test_that("contrasts sharing a species are coupled, not refused", {
  # A species in more than one contrast makes the relabellings dependent:
  # flipping one contrast changes the other. That used to be refused
  # outright, which restricted the test to a disjoint matched-pairs
  # design. The unit of independence is the connected component of the
  # species-by-contrast graph, so the coupling is modelled rather than
  # excluded, and a disjoint design is the special case where every
  # component is one contrast.
  fix <- make_tag_perm_fixtures()
  chain <- data.frame(
    sp1 = c("A1", "P1"), sp2 = c("P1", "A2"),
    pair_name = c("pair1", "pair2"), stringsAsFactors = FALSE
  )
  bl <- .tp_blocks(chain, fix$group)
  # One component over A1, P1, A2; fixing A1 forces the other two, so
  # there are 2 labellings, not the 4 a disjoint reading would claim.
  expect_equal(length(bl$labellings), 1L)
  expect_equal(unname(bl$n), 2L)
  lab <- bl$labellings[[1]]
  expect_true(all(vapply(lab, function(l) {
    l[["A1"]] != l[["P1"]] && l[["P1"]] != l[["A2"]]
  }, logical(1))))

  # A disjoint design still gives one component per contrast and 2^k.
  bl2 <- .tp_blocks(fix$pairs, fix$group)
  expect_equal(length(bl2$labellings), 3L)
  expect_equal(prod(bl2$n), 8)

  # A contrast whose two species share a label is pinned, not swappable.
  same <- fix$group
  same["P1"] <- "annual"
  bl3 <- .tp_blocks(fix$pairs[1, , drop = FALSE], same)
  expect_equal(unname(bl3$n), 1L)

  # Self-pairs and duplicate pair names remain errors.
  bad <- fix$pairs
  bad$sp2[1] <- bad$sp1[1]
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    bad, fix$group, target_group = "annual"),
    "two distinct species"
  )
  dup <- fix$pairs
  dup$pair_name[2] <- dup$pair_name[1]
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    dup, fix$group, target_group = "annual"),
    "pair_name must be unique"
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
      min_recurrence = 2L,
      statistic = "count"
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
  res <- suppressMessages(suppressWarnings(tag_permutation(
    skewed$classification, skewed$modules, skewed$orthologs,
    skewed$pairs, skewed$group, target_group = "annual",
    min_recurrence = 2L
  )))
  # A clean sweep of three pairs is a sign-test p of 0.125, which is above
  # the advisory 0.10 threshold: the diagnostic is honestly underpowered
  # at small k, exactly as p_min is. The value is still reported.
  expect_equal(res$size_asymmetry_p,
               stats::binom.test(3L, 3L, 0.5,
                                 alternative = "greater")$p.value)
  expect_gt(res$size_asymmetry_p, 0.10)
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


# Fixture generalised to k pairs, for the resolution and floor tests.
make_tag_perm_fixtures_k <- function(k) {
  ann <- paste0("A", seq_len(k))
  per <- paste0("P", seq_len(k))
  group <- stats::setNames(rep(c("annual", "perennial"), each = k),
                           c(ann, per))
  pairs <- data.frame(
    sp1 = ann, sp2 = per,
    pair_name = paste0("pair", seq_len(k)),
    stringsAsFactors = FALSE
  )
  mk <- function(sp) {
    g <- paste0(sp, "_g", 1:6)
    memb <- stats::setNames(rep(c(1L, 2L), each = 3L), g)
    list(modules = memb, module_genes = split(names(memb), memb),
         n_modules = 2L, modularity = 0.3, graph = NULL,
         method = "leiden", params = list())
  }
  modules <- stats::setNames(lapply(c(ann, per), mk), c(ann, per))
  orthologs <- data.frame(
    Species1 = unlist(lapply(ann, function(s) paste0(s, "_g", 1:6))),
    Species2 = unlist(lapply(per, function(s) paste0(s, "_g", 1:6))),
    hog = rep(paste0("HOG", 1:6), k),
    stringsAsFactors = FALSE
  )
  # Module 1 of each annual and module 2 of each perennial are diverged.
  cls <- do.call(rbind, lapply(seq_len(k), function(i) {
    data.frame(
      pair_name = rep(pairs$pair_name[i], 4),
      module = c("1", "2", "1", "2"),
      reference = c(ann[i], ann[i], per[i], per[i]),
      test = c(per[i], per[i], ann[i], ann[i]),
      classification = c("diverged", "conserved",
                         "conserved", "diverged"),
      stringsAsFactors = FALSE
    )
  }))
  list(classification = cls, modules = modules, orthologs = orthologs,
       pairs = pairs, group = group)
}


test_that("p_min counts only pairs whose swap changes the statistic", {
  # A pair with no diverged module on either side carries the same
  # (empty) HOG set both ways, so its bit duplicates every labelling
  # without adding a point of resolution. Counting it reported
  # p_min = 1/32 at five pairs when nothing could score below 1/16 --
  # a false resolution claim, and silent, which is exactly what this
  # guard exists to refuse.
  fix <- make_tag_perm_fixtures_k(5)
  degen <- fix$classification$pair_name == "pair5"
  fix$classification$classification[degen] <- "conserved"

  res <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual", min_recurrence = 2L
  )))

  expect_equal(res$n_swappable, 4L)      # not 5
  expect_equal(res$n_contributing, 5L)   # it still feeds the statistic
  expect_equal(res$p_min, 1 / 16)        # not 1/32
  expect_length(res$null_distribution, 16L)
  expect_gte(res$p_value, res$p_attainable)

  # And the warning fires, naming the degenerate pair. Under the old
  # definition p_min was 0.031 and no warning was raised at all.
  warns <- character(0)
  withCallingHandlers(
    suppressMessages(tag_permutation(
      fix$classification, fix$modules, fix$orthologs, fix$pairs,
      fix$group, target_group = "annual", min_recurrence = 2L
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("carry the same HOG set on both sides", warns)))

  # With all five pairs diverged the fifth bit is real again.
  full <- make_tag_perm_fixtures_k(5)
  res5 <- suppressMessages(suppressWarnings(tag_permutation(
    full$classification, full$modules, full$orthologs, full$pairs,
    full$group, target_group = "annual", min_recurrence = 2L
  )))
  expect_equal(res5$n_swappable, 5L)
  expect_equal(res5$p_min, 1 / 32)
  expect_length(res5$null_distribution, 32L)
})


test_that("enumeration is decided on cost, not on n_perm", {
  # Ten pairs is 1024 labellings: previously the default n_perm = 1000
  # sampled 1000 points from that space -- inexact, and slower than
  # walking all of it.
  fix <- make_tag_perm_fixtures_k(10)
  res <- suppressMessages(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 1000L, min_recurrence = 2L
  ))
  expect_true(res$exact)
  expect_length(res$null_distribution, 1024L)
  expect_equal(res$p_min, 1 / 1024)
  expect_equal(res$n_swappable, 10L)
})


test_that("tag_permutation rejects unusable n_perm and NA traits", {
  fix <- make_tag_perm_fixtures()
  call_with <- function(...) {
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, target_group = "annual", ...)
  }
  # n_perm = 0 used to return p = 1 with p_min = 1 in silence; 3e9
  # overflowed as.integer() to NA and died on an unrelated `if`.
  expect_error(call_with(n_perm = 0L), "n_perm must be")
  expect_error(call_with(n_perm = -5L), "n_perm must be")
  expect_error(call_with(n_perm = c(10L, 20L)), "n_perm must be")
  expect_error(call_with(min_recurrence = 0L), "min_recurrence must be")
  # 3e9 is above .Machine$integer.max: clamped rather than overflowed to
  # NA, and the null is enumerated anyway. (The p_min warning still
  # fires -- three pairs.)
  big <- suppressMessages(suppressWarnings(call_with(n_perm = 3e9)))
  expect_true(big$exact)
  expect_equal(big$n_perm, 8L)

  bad <- fix$group
  bad["P3"] <- NA_character_
  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, bad, target_group = "annual"),
    "NA trait values"
  )
})


test_that("the sampled branch runs when enumeration is capped", {
  # Moving the enumeration decision off n_perm made this branch
  # unreachable in practice (21 disjoint pairs = 42 species), so nothing
  # exercised the runif() draws, the +1-corrected p-value or its p_min.
  # enum_max is the hook that lets a test reach it.
  fix <- make_tag_perm_fixtures_k(4)

  set.seed(3)
  res <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual",
    n_perm = 50L, min_recurrence = 2L, enum_max = 2L
  )))
  expect_false(res$exact)
  expect_length(res$null_distribution, 50L)
  # The +1 correction returns: the observed labelling is not guaranteed
  # to be among the draws.
  expect_equal(res$p_min, 1 / 51)
  expect_equal(res$p_attainable, res$p_min)
  expect_equal(res$p_value,
               (sum(res$null_distribution >= res$statistic_observed) + 1L) /
                 51L)
  expect_equal(res$n_swappable, 4L)

  # The enumerated answer on the same data, for comparison.
  ex <- suppressMessages(suppressWarnings(tag_permutation(
    fix$classification, fix$modules, fix$orthologs, fix$pairs,
    fix$group, target_group = "annual", min_recurrence = 2L
  )))
  expect_true(ex$exact)
  expect_length(ex$null_distribution, 16L)

  expect_error(
    tag_permutation(fix$classification, fix$modules, fix$orthologs,
                    fix$pairs, fix$group, target_group = "annual",
                    enum_max = 31L),
    "enum_max must be"
  )
})


test_that("the size-asymmetry warning fires when the sign test can see it", {
  # Four pairs is the smallest clean sweep the advisory 0.10 threshold
  # can reach (binom.test(4, 4) = 0.0625). Nothing previously exercised
  # the warning branch at all.
  fix <- make_tag_perm_fixtures_k(4)
  for (sp in paste0("P", 1:4)) {
    fix$modules[[sp]]$module_genes[["2"]] <-
      fix$modules[[sp]]$module_genes[["2"]][1]
  }

  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(tag_permutation(
      fix$classification, fix$modules, fix$orthologs, fix$pairs,
      fix$group, target_group = "annual", min_recurrence = 2L
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(res$size_asymmetry_p,
               stats::binom.test(4L, 4L, 0.5,
                                 alternative = "greater")$p.value)
  expect_lte(res$size_asymmetry_p, 0.10)
  expect_true(any(grepl("larger HOG set in 4 of 4", warns)))
  expect_true(any(grepl("reading set size rather than recurrence", warns)))
})


test_that("the unreachable-significance warning names the right cause", {
  # With enough labellings but a statistic that cannot separate them, the
  # floor is set by ties, not by k. Blaming k there would prescribe more
  # pairs for a problem more pairs do not solve.
  fix <- make_tag_perm_fixtures_k(6)
  warns <- character(0)
  res <- withCallingHandlers(
    suppressMessages(tag_permutation(
      fix$classification, fix$modules, fix$orthologs, fix$pairs,
      fix$group, target_group = "annual", min_recurrence = 2L
    )),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # Every side of this fixture carries the same three HOGs, so the null
  # is heavily tied: p_min is 1/64 but the maximum is shared by most
  # labellings, putting the realised floor far above it.
  expect_equal(res$p_min, 1 / 64)
  expect_gt(res$p_attainable, 0.05)
  expect_gt(res$p_attainable, res$p_min)
  expect_true(any(grepl("^ties in the null", warns)))
  expect_false(any(grepl("At least 5 swappable pairs are needed", warns)))
})


# Build a k-contrast design whose sides are random subsets of a HOG pool,
# with one contrast deliberately lopsided so set size has something to
# dominate. Returns the pieces tag_permutation() consumes plus the sizes,
# so a test can regress the null on total selected size.
lopsided_design <- function(n_hog, sizes, seed) {
  set.seed(seed)
  k <- length(sizes)
  hogs <- paste0("HOG", seq_len(n_hog))
  sp1 <- paste0("A", seq_len(k))
  sp2 <- paste0("P", seq_len(k))
  group <- stats::setNames(rep(c("annual", "perennial"), each = k),
                           c(sp1, sp2))
  pairs <- data.frame(sp1 = sp1, sp2 = sp2,
                      pair_name = paste0("pair", seq_len(k)),
                      stringsAsFactors = FALSE)
  modules <- list()
  orth <- list()
  cls <- list()
  for (i in seq_len(k)) {
    for (side in 1:2) {
      sp <- if (side == 1L) sp1[i] else sp2[i]
      hs <- sample(hogs, sizes[[i]][side])
      g <- paste0(sp, "_g", seq_along(hs))
      modules[[sp]] <- list(
        modules = stats::setNames(rep(1L, length(g)), g),
        module_genes = list(`1` = g), n_modules = 1L,
        modularity = 0.3, graph = NULL, method = "leiden",
        params = list()
      )
      orth[[sp]] <- data.frame(Species1 = g, Species2 = g, hog = hs,
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
  list(cl = do.call(rbind, cls), og = do.call(rbind, orth),
       modules = modules, pairs = pairs, group = group,
       k = k, sizes = sizes)
}

# R^2 of the null on the total size of the selected sides, walked in the
# same mixed-radix order the kernel uses.
size_dependence <- function(fx, ...) {
  res <- suppressMessages(suppressWarnings(tag_permutation(
    fx$cl, fx$modules, fx$og, fx$pairs, fx$group,
    target_group = "annual", min_recurrence = 2L, ...
  )))
  tot <- vapply(seq_len(2^fx$k) - 1L, function(idx) {
    bits <- as.logical(bitwAnd(idx, bitwShiftL(1L, seq_len(fx$k) - 1L)))
    sum(vapply(seq_len(fx$k), function(i) {
      fx$sizes[[i]][if (bits[i]) 2L else 1L]
    }, numeric(1)))
  }, numeric(1))
  # cor^2 rather than lm(), which is the same quantity for a single
  # predictor and keeps the variables out of a formula.
  stats::cor(res$null_distribution, tot)^2
}


test_that("the raw count null is dominated by selected set size", {
  # This is the motivation for offering "excess" at all, and it is a
  # property of the statistic rather than of any one dataset: a contrast
  # whose two sides differ greatly in size decides the ordering.
  fx <- lopsided_design(6000, list(c(60, 60), c(60, 60), c(60, 60),
                                   c(300, 8)), seed = 21)
  # Measured 0.71 here and 0.98 on the eight-species Pooideae set; assert
  # that size is the dominant term rather than either exact level, since
  # how dominant it is depends on the pool the sides are drawn from.
  expect_gt(size_dependence(fx, statistic = "count"), 0.6)
})


test_that("the excess correction helps or hurts with the universe", {
  # The correction subtracts a Poisson-binomial expectation over
  # `universe`, which the data cannot identify. Assert both directions so
  # the regime dependence is pinned rather than assumed away.
  fx <- lopsided_design(6000, list(c(60, 60), c(60, 60), c(60, 60),
                                   c(300, 8)), seed = 21)
  raw <- size_dependence(fx, statistic = "count")

  # Sides small relative to the universe: the correction is well posed
  # and reduces the size dependence.
  big <- size_dependence(fx, statistic = "excess", universe = 20000)
  expect_lt(big, raw)

  # A universe near the union of the selected sides: the independence
  # model predicts more overlap than disjoint sides can deliver, so the
  # statistic goes systematically negative and tracks size again.
  small <- suppressMessages(suppressWarnings(tag_permutation(
    fx$cl, fx$modules, fx$og, fx$pairs, fx$group,
    target_group = "annual", min_recurrence = 2L,
    statistic = "excess", universe = 500
  )))
  expect_lt(small$statistic_observed, 0)
  expect_gt(size_dependence(fx, statistic = "excess", universe = 500),
            raw)

  # Inference stays exact under both: same label space, same floor.
  a <- suppressMessages(suppressWarnings(tag_permutation(
    fx$cl, fx$modules, fx$og, fx$pairs, fx$group,
    target_group = "annual", min_recurrence = 2L, statistic = "count"
  )))
  expect_equal(small$n_labellings, a$n_labellings)
  expect_equal(small$p_min, a$p_min)
})


test_that("min_recurrence scales with the number of contrasts", {
  # A fixed threshold describes one design size. The chance a HOG reaches
  # two of k sides unaided grows steeply with k, so with min_recurrence
  # pinned at 2 the statistic saturates on chance recurrence as contrasts
  # are added and power stops rising with the design.
  for (k in c(4L, 6L, 10L)) {
    fix <- make_tag_perm_fixtures_k(k)
    res <- suppressMessages(suppressWarnings(tag_permutation(
      fix$classification, fix$modules, fix$orthologs, fix$pairs,
      fix$group, target_group = "annual"
    )))
    expect_equal(res$n_contributing, k)
    expect_equal(res$min_recurrence, max(2L, as.integer(round(k / 2))))
  }

  # Four contrasts resolve to 2, the previous fixed default, so existing
  # small designs are unaffected by the change.
  fix4 <- make_tag_perm_fixtures_k(4)
  auto <- suppressMessages(suppressWarnings(tag_permutation(
    fix4$classification, fix4$modules, fix4$orthologs, fix4$pairs,
    fix4$group, target_group = "annual"
  )))
  fixed <- suppressMessages(suppressWarnings(tag_permutation(
    fix4$classification, fix4$modules, fix4$orthologs, fix4$pairs,
    fix4$group, target_group = "annual", min_recurrence = 2L
  )))
  expect_equal(auto$min_recurrence, 2L)
  expect_equal(auto$observed, fixed$observed)
  expect_equal(auto$p_value, fixed$p_value)

  # An explicit value still wins, and still cannot exceed the design.
  expl <- suppressMessages(suppressWarnings(tag_permutation(
    fix4$classification, fix4$modules, fix4$orthologs, fix4$pairs,
    fix4$group, target_group = "annual", min_recurrence = 3L
  )))
  expect_equal(expl$min_recurrence, 3L)
  expect_error(
    tag_permutation(fix4$classification, fix4$modules, fix4$orthologs,
                    fix4$pairs, fix4$group, target_group = "annual",
                    min_recurrence = 9L),
    "exceeds the number of pairs contributing"
  )
  expect_error(
    tag_permutation(fix4$classification, fix4$modules, fix4$orthologs,
                    fix4$pairs, fix4$group, target_group = "annual",
                    min_recurrence = 0L),
    "min_recurrence must be"
  )
})
