# Fixtures for the package-wide RNG contract test (test-rng-contract.R).
#
# Every case is a self-contained thunk of the form function(seed), so the
# contract assertions can be written once and looped over the table.

#' Two 40-gene block networks plus a 1:1 ortholog table
#'
#' Four disjoint blocks of ten genes, edge weights drawn once from a fixed
#' seed so the fixture itself is deterministic. Ten genes per block is
#' exactly `module_preservation()`'s default `min_module_size`, so every
#' module is tested and no coverage message is emitted.
rng_module_fixture <- function() {
  n <- 40L
  blk <- rep(seq_len(4L), each = 10L)
  same <- outer(blk, blk, "==")
  build <- function(prefix, seed) {
    w <- withr::with_seed(seed, {
      m <- matrix(stats::runif(n * n, 0.4, 1), n, n)
      (m + t(m)) / 2
    })
    w[!same] <- 0
    diag(w) <- 0
    g <- paste0(prefix, seq_len(n))
    dimnames(w) <- list(g, g)
    list(network = w, threshold = 0.2)
  }
  net_a <- build("A", 11L)
  net_b <- build("B", 12L)
  ortho <- data.frame(
    Species1 = rownames(net_a$network),
    Species2 = rownames(net_b$network),
    hog = paste0("HOG", seq_len(n)),
    stringsAsFactors = FALSE
  )
  mods <- function(net) {
    detect_modules(net,
      resolution = 1.0, objective_function = "modularity",
      seed = 1L
    )
  }
  list(
    net_a = net_a, net_b = net_b, ortho = ortho,
    mods_a = mods(net_a), mods_b = mods(net_b),
    map = resolve_ortholog_map(
      ortho, rownames(net_a$network), rownames(net_b$network)
    )
  )
}


#' Synthetic all-pairs preservation classification over four species
#'
#' `preservation_matrix_test()` reads only `reference`, `test` and the
#' effect column, so the table is written directly rather than run through
#' a preservation pipeline that would add minutes for no extra coverage.
rng_matrix_classification <- function() {
  sp <- c("A1", "A2", "P1", "P2")
  grid <- expand.grid(
    reference = sp, test = sp, module = c("1", "2"),
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$reference != grid$test, , drop = FALSE]
  trait <- c(A1 = "annual", A2 = "annual", P1 = "peren", P2 = "peren")
  concordant <- trait[grid$reference] == trait[grid$test]
  grid$Zsummary_std <- ifelse(concordant, 8, 2)
  grid$Zsummary <- grid$Zsummary_std
  grid$classification <- "conserved"
  rownames(grid) <- NULL
  list(classification = grid, group = trait)
}


#' Six-species trait-recurrence fixture for `tag_permutation()`
#'
#' Three annual/perennial contrasts, two modules per species, six HOGs.
#' The same shape as the fixture in test-tag-permutation.R, kept separate
#' so a change there cannot silently alter what the contract test runs.
rng_tag_fixture <- function() {
  annuals <- c("A1", "A2", "A3")
  perennials <- c("P1", "P2", "P3")
  group <- stats::setNames(
    rep(c("annual", "peren"), each = 3L), c(annuals, perennials)
  )
  pairs <- data.frame(
    sp1 = annuals, sp2 = perennials,
    pair_name = paste0("pair", 1:3), stringsAsFactors = FALSE
  )
  one_module_set <- function(genes) {
    membership <- stats::setNames(
      rep(c(1L, 2L), each = length(genes) / 2L), genes
    )
    list(
      modules = membership,
      module_genes = split(names(membership), membership),
      n_modules = 2L, modularity = 0.3, graph = NULL,
      method = "leiden", params = list()
    )
  }
  modules <- stats::setNames(
    lapply(c(annuals, perennials), function(sp) {
      one_module_set(paste0(sp, "_g", 1:6))
    }),
    c(annuals, perennials)
  )
  orthologs <- data.frame(
    Species1 = as.vector(vapply(annuals, function(sp) {
      paste0(sp, "_g", 1:6)
    }, character(6))),
    Species2 = as.vector(vapply(perennials, function(sp) {
      paste0(sp, "_g", 1:6)
    }, character(6))),
    hog = rep(paste0("HOG", 1:6), 3),
    stringsAsFactors = FALSE
  )
  classification <- data.frame(
    pair_name = rep(pairs$pair_name, each = 4L),
    module = rep(c("1", "2", "1", "2"), 3L),
    reference = rep(c(rbind(annuals, annuals, perennials, perennials)), 1L),
    test = rep(c(rbind(perennials, perennials, annuals, annuals)), 1L),
    classification = rep(
      c("diverged", "conserved", "conserved", "diverged"), 3L
    ),
    stringsAsFactors = FALSE
  )
  list(
    classification = classification, modules = modules,
    orthologs = orthologs, pairs = pairs, group = group
  )
}


#' Every exported entry point that takes a `seed`
#'
#' Generics forward through `...`, so the `seed` formal lives on the
#' `.default` method; both are searched. A new seeded export appears here
#' automatically, and the contract test fails until the table covers it.
rng_seeded_exports <- function() {
  ns <- asNamespace("rcomplex")
  out <- character(0)
  for (nm in sort(getNamespaceExports("rcomplex"))) {
    f <- get(nm, envir = ns)
    if (!is.function(f)) next
    cand <- nm
    default <- paste0(nm, ".default")
    if (exists(default, envir = ns, inherits = FALSE)) {
      cand <- c(cand, default)
    }
    for (cn in cand) {
      g <- get(cn, envir = ns)
      if (is.function(g) && "seed" %in% names(formals(g))) {
        out <- c(out, cn)
      }
    }
  }
  unique(out)
}


#' The contract table: one entry per seeded entry point
#'
#' `call` runs the function at the given seed (possibly `NULL`) and returns
#' something comparable with `identical()`. Everything is sized for speed;
#' the assertions are about the RNG stream, not about the statistics.
#'
#' @param fx The fixture bundle built at the top of test-rng-contract.R.
rng_contract_cases <- function(fx) {
  td <- fx$td
  nets <- fx$nets
  cmp <- fx$cmp
  sparse_nets <- fx$sparse_nets
  mf <- fx$mf
  cf <- fx$cf
  mx <- fx$mx
  tf <- fx$tf

  list(
    list(
      name = "summarize_comparison",
      call = function(seed) summarize_comparison(cmp, seed = seed)$results
    ),
    list(
      name = "find_coexpressologs.default",
      call = function(seed) find_coexpressologs(nets, td$ortho, seed = seed)
    ),
    list(
      name = "density_sweep.default",
      call = function(seed) {
        suppressMessages(density_sweep(nets, td$ortho,
          multipliers = 1, method = "analytical", seed = seed
        ))$edges
      }
    ),
    list(
      name = "permutation_hog_test",
      call = function(seed) {
        permutation_hog_test(td$net1, td$net2, cmp,
          min_exceedances = 3L, max_permutations = 40L, seed = seed
        )
      }
    ),
    list(
      name = "coexpressolog_null",
      variant = "pi0 none",
      # n_perm = 2 is far below the 19 permutations p < 0.05 needs, and an
      # unseeded call announces the seed it drew; neither is what the
      # contract is about, and the assertions are unaffected by both.
      call = function(seed) {
        suppressMessages(suppressWarnings(
          coexpressolog_null(sparse_nets, td$ortho,
            n_perm = 2L, swap_factor = 1L, seed = seed,
            pi0_method = "none", pval_combine = "max"
          )
        ))
      }
    ),
    list(
      name = "coexpressolog_null",
      variant = "pi0 randomized",
      # The scope covers the observed run, not just the permutation loop.
      # Under pi0_method = "none" the observed run draws nothing, so that
      # is the only case that would notice the scope sliding back below it.
      call = function(seed) {
        suppressMessages(suppressWarnings(
          coexpressolog_null(sparse_nets, td$ortho,
            n_perm = 2L, swap_factor = 1L, seed = seed,
            pi0_method = "randomized", pval_combine = "max"
          )
        ))
      }
    ),
    list(
      name = "detect_modules.default",
      variant = "single resolution",
      call = function(seed) {
        detect_modules(mf$net_a,
          resolution = 1.0, objective_function = "modularity", seed = seed
        )$modules
      }
    ),
    list(
      name = "detect_modules.default",
      variant = "consensus",
      # The consensus branch returns before the .seed_scope() in
      # detect_modules.default() and carries its own inside
      # detect_modules_consensus(), so the source grep cannot see it and
      # the single-resolution case never reaches it.
      call = function(seed) {
        detect_modules(mf$net_a,
          resolution = c(0.8, 1.0), objective_function = "modularity",
          n_iterations = 1L, max_consensus_iter = 1L, test_k1 = FALSE,
          seed = seed
        )$modules
      }
    ),
    list(
      name = "module_preservation",
      call = function(seed) {
        module_preservation(mf$mods_a, mf$net_a, mf$net_b,
          map = mf$map, n_perm = 20L, seed = seed
        )$preservation
      }
    ),
    list(
      name = "module_correspondence",
      call = function(seed) {
        module_correspondence(mf$mods_a, mf$mods_b, mf$map,
          seed = seed
        )$pairs
      }
    ),
    list(
      name = "preservation_paired.default",
      call = function(seed) {
        preservation_paired(
          list(A = mf$mods_a, B = mf$mods_b),
          list(A = mf$net_a, B = mf$net_b),
          mf$ortho,
          data.frame(sp1 = "A", sp2 = "B", stringsAsFactors = FALSE),
          n_perm = 20L, seed = seed
        )$classification
      }
    ),
    list(
      name = "preservation_matrix_test",
      call = function(seed) {
        # 30 draws is far too few for a usable p-value and the function
        # says so; the contract is about the stream, not the inference.
        suppressWarnings(preservation_matrix_test(
          mx$classification, mx$group,
          n_perm = 30L, enum_max = 1L, seed = seed
        ))$free$null_distribution
      }
    ),
    list(
      name = "tag_permutation",
      call = function(seed) {
        suppressWarnings(suppressMessages(tag_permutation(
          tf$classification, tf$modules, tf$orthologs, tf$pairs,
          tf$group,
          target_group = "annual", n_perm = 30L,
          min_recurrence = 2L, enum_max = 1L, seed = seed
        )))$null_distribution
      }
    ),
    list(
      name = "clique_perturbation_test.default",
      call = function(seed) {
        clique_perturbation_test(cf$cliques, cf$target_species,
          cf$networks, cf$orthologs,
          n_boot = 2L, seed = seed, pi0_method = "none"
        )
      }
    ),
    list(
      name = "clique_intensity_test.default",
      call = function(seed) {
        clique_intensity_test(cf$cliques, cf$target_species,
          cf$networks, cf$orthologs,
          n_perm = 2L, edges = cf$edges, seed = seed,
          pi0_method = "none"
        )
      }
    )
  )
}
