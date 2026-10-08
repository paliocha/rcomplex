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
    gene1 = rownames(net_a$network),
    gene2 = rownames(net_b$network),
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


#' Synthetic all-pairs preservation classification over twenty species
#'
#' `preservation_matrix_test()` reads only `reference`, `test` and the
#' effect column, so the table is written directly rather than run through
#' a preservation pipeline that would add minutes for no extra coverage.
#' Ten species per trait give 184756 free labellings, above the 50000 the
#' test enumerates, so the null is sampled and draws from the stream.
rng_matrix_classification <- function() {
  sp <- c(paste0("A", 1:10), paste0("P", 1:10))
  grid <- expand.grid(
    reference = sp, test = sp, module = "1",
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$reference != grid$test, , drop = FALSE]
  trait <- stats::setNames(rep(c("annual", "peren"), each = 10L), sp)
  concordant <- trait[grid$reference] == trait[grid$test]
  grid$Zsummary_std <- ifelse(concordant, 8, 2)
  grid$Zsummary <- grid$Zsummary_std
  grid$classification <- "conserved"
  rownames(grid) <- NULL
  list(classification = grid, group = trait)
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
  mx <- fx$mx
  list(
    list(
      name = "summarize_comparison",
      internal = TRUE,
      call = function(seed) {
        rcomplex:::summarize_comparison(cmp, seed = seed)$results
      }
    ),
    list(
      name = "find_coexpressologs.default",
      call = function(seed) find_coexpressologs(nets, td$ortho, seed = seed)
    ),
    list(
      name = "density_sweep.default",
      call = function(seed) {
        suppressMessages(density_sweep(nets, td$ortho,
          multipliers = 1, method = "hypergeometric", seed = seed
        ))$edges
      }
    ),
    list(
      name = "permutation_hog_test",
      internal = TRUE,
      call = function(seed) {
        rcomplex:::permutation_hog_test(td$net1, td$net2, cmp,
          min_exceedances = 3L, max_permutations = 40L, seed = seed
        )
      }
    ),
    list(
      name = "coexpressolog_null",
      # An unseeded call announces the seed it drew; that is not what the
      # contract is about, and the assertions are unaffected by it. The
      # scope covers the observed run, not just the permutation loop: its
      # randomized pi0 draws.
      call = function(seed) {
        suppressMessages(suppressWarnings(
          coexpressolog_null(sparse_nets, td$ortho,
            swap_factor = 1L, seed = seed, pval_combine = "max"
          )
        ))
      }
    ),
    list(
      name = "null_network",
      call = function(seed) {
        null_network(fx$null_x, fx$null_net, seed = seed)$network
      }
    ),
    list(
      name = "null_network",
      variant = "block",
      call = function(seed) {
        null_network(fx$null_x, fx$null_net,
          seed = seed,
          block = rep_len(1:3, ncol(fx$null_x))
        )$network
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
          max_consensus_iter = 1L, seed = seed
        )$modules
      }
    ),
    list(
      name = "module_preservation",
      call = function(seed) {
        module_preservation(mf$mods_a, mf$net_a, mf$net_b, mf$ortho,
          n_perm = 20L, seed = seed
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
          data.frame(species1 = "A", species2 = "B", stringsAsFactors = FALSE),
          n_perm = 20L, seed = seed
        )$classification
      }
    ),
    list(
      name = "rcomplex",
      call = function(seed) {
        rcomplex(fx$drv_expr, fx$drv_ortho,
          density = 0.1, null = TRUE, seed = seed
        )$edges_null
      }
    ),
    list(
      name = "preservation_matrix_test",
      call = function(seed) {
        suppressWarnings(preservation_matrix_test(
          mx$classification, mx$group,
          seed = seed
        ))$free$null_distribution
      }
    )
  )
}
