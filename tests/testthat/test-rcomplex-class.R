# Tests for S3 rcomplex container class

# --- Shared fixtures ---

make_rcx_fixtures <- function(n_sp = 3) {
  sp <- paste0("SP_", LETTERS[seq_len(n_sp)])
  traits <- stats::setNames(
    rep(c("annual", "perennial"), length.out = n_sp), sp
  )

  set.seed(42)
  networks <- stats::setNames(lapply(sp, function(s) {
    expr <- matrix(rnorm(200), nrow = 20, ncol = 10)
    rownames(expr) <- paste0(s, "_G", 1:20)
    compute_network(expr, density = 0.1, mr_log_transform = FALSE)
  }), sp)

  # Orthologs: first 15 genes of each pair share HOGs
  pairs <- utils::combn(sp, 2, simplify = FALSE)
  ortho_list <- lapply(pairs, function(p) {
    data.frame(
      Species1 = paste0(p[1], "_G", 1:15),
      Species2 = paste0(p[2], "_G", 1:15),
      hog = paste0("HOG", 1:15)
    )
  })
  orthologs <- do.call(rbind, ortho_list)

  list(species = sp, traits = traits, networks = networks,
       orthologs = orthologs)
}


# --- Constructor tests ---

test_that("rcomplex constructor creates valid object", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)

  expect_s3_class(rcx, "rcomplex")
  expect_equal(rcx$species, fix$species)
  expect_equal(rcx$traits, fix$traits)
  expect_equal(length(rcx$species_pairs), 3)  # combn(3,2)
  expect_null(rcx$phylo_pairs)
  # All result slots are NULL
  expect_null(rcx$edges)
  expect_null(rcx$modules)
  expect_null(rcx$cliques)
  expect_null(rcx$classification)
})


test_that("rcomplex constructor validates species", {
  fix <- make_rcx_fixtures()
  expect_error(rcomplex("A", fix$traits, fix$networks, fix$orthologs),
               "at least 2")
  expect_error(rcomplex(c("A", "A"), fix$traits, fix$networks, fix$orthologs),
               "duplicates")
  expect_error(rcomplex(1:3, fix$traits, fix$networks, fix$orthologs),
               "character vector")
})


test_that("rcomplex constructor validates traits", {
  fix <- make_rcx_fixtures()
  expect_error(rcomplex(fix$species, c("a", "b", "c"),
                         fix$networks, fix$orthologs),
               "named")
  expect_error(rcomplex(fix$species, c(SP_A = "annual"),
                         fix$networks, fix$orthologs),
               "traits missing species")
})


test_that("rcomplex constructor validates networks", {
  fix <- make_rcx_fixtures()
  bad_nets <- fix$networks
  names(bad_nets)[1] <- "WRONG"
  expect_error(rcomplex(fix$species, fix$traits, bad_nets, fix$orthologs),
               "network names must match")
  bad_nets2 <- fix$networks
  bad_nets2[[1]]$threshold <- NULL
  expect_error(rcomplex(fix$species, fix$traits, bad_nets2, fix$orthologs),
               "threshold")
})


test_that("rcomplex constructor validates orthologs", {
  fix <- make_rcx_fixtures()
  expect_error(rcomplex(fix$species, fix$traits, fix$networks,
                         data.frame(x = 1)),
               "Species1, Species2, hog")
})


test_that("rcomplex constructor validates species_pairs", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
             species_pairs = list(c("SP_A", "BOGUS"))),
    "unknown species"
  )
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
             species_pairs = list("SP_A")),
    "length-2"
  )
})


test_that("rcomplex constructor validates phylo_pairs", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
             phylo_pairs = data.frame(x = 1)),
    "sp1, sp2"
  )
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
             phylo_pairs = data.frame(sp1 = "BOGUS", sp2 = "SP_A")),
    "unknown species"
  )
  # Auto-generates pair_name
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
                   phylo_pairs = data.frame(sp1 = "SP_A", sp2 = "SP_B"))
  expect_equal(rcx$phylo_pairs$pair_name, "SP_A.SP_B")
})


# --- Print / summary tests ---

test_that("print.rcomplex shows pipeline status", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  out <- capture.output(print(rcx))
  expect_true(any(grepl("rcomplex:", out)))
  expect_true(any(grepl("\\[pending\\]", out)))
})


test_that("summary.rcomplex returns structured list", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  s <- summary(rcx)
  expect_s3_class(s, "summary.rcomplex")
  expect_equal(s$n_species, 3)
  expect_equal(s$n_orthologs, nrow(fix$orthologs))
})


# --- S3 dispatch tests ---

test_that("find_coexpressologs dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- find_coexpressologs(rcx)
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$edges))
  expect_true(nrow(rcx$edges) > 0)
})


test_that("find_coexpressologs.default still works with raw args", {
  fix <- make_rcx_fixtures()
  result <- find_coexpressologs(fix$networks, fix$orthologs)
  expect_true(is.data.frame(result))
  expect_false(inherits(result, "rcomplex"))
})


test_that("find_cliques dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- find_coexpressologs(rcx)
  rcx <- find_cliques(rcx, min_species = 2L)
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$cliques))
})


test_that("find_cliques errors when edges not yet computed", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(find_cliques(rcx), "run find_coexpressologs")
})


test_that("classify_cliques dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- find_coexpressologs(rcx)
  rcx <- classify_cliques(rcx, min_species = 2L)
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$classification))
})


test_that("detect_modules dispatches on rcomplex (loops species)", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- detect_modules(rcx, method = "leiden",
                         objective_function = "modularity")
  expect_s3_class(rcx, "rcomplex")
  expect_equal(length(rcx$modules), length(fix$species))
  expect_equal(names(rcx$modules), fix$species)
  for (sp in fix$species) {
    expect_true(!is.null(rcx$modules[[sp]]$modules))
  }
})


test_that("identify_module_hubs dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- detect_modules(rcx, method = "leiden",
                         objective_function = "modularity")
  rcx <- identify_module_hubs(rcx)
  expect_s3_class(rcx, "rcomplex")
  expect_equal(length(rcx$hubs), length(fix$species))
  for (sp in fix$species) {
    expect_true(is.data.frame(rcx$hubs[[sp]]))
  }
})


test_that("identify_module_hubs errors without modules", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(identify_module_hubs(rcx), "run detect_modules")
})


test_that("density_sweep dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- suppressMessages(density_sweep(rcx, multipliers = c(0.95, 1.0, 1.05)))
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$sweep))
  expect_equal(nrow(rcx$sweep), 3L)
})


test_that("clique_stability dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs) |>
    find_coexpressologs() |>
    find_cliques(min_species = 2L)
  rcx <- clique_stability(rcx, min_species = 2L, max_k = 1L)
  expect_s3_class(rcx, "rcomplex")
  expect_true(!is.null(rcx$stability))
})


test_that("clique_stability errors without edges", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(clique_stability(rcx), "run find_coexpressologs")
})


test_that("classify_cliques errors without edges", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(classify_cliques(rcx), "run find_coexpressologs")
})


test_that("compare_modules_paired dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  phylo <- data.frame(sp1 = "SP_A", sp2 = "SP_B")
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
                   phylo_pairs = phylo)
  rcx <- detect_modules(rcx, method = "leiden",
                          objective_function = "modularity")
  rcx <- compare_modules_paired(rcx)
  expect_s3_class(rcx, "rcomplex")
  expect_true(!is.null(rcx$module_comparisons))
})


test_that("compare_modules_paired errors without modules", {
  fix <- make_rcx_fixtures()
  phylo <- data.frame(sp1 = "SP_A", sp2 = "SP_B")
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
                   phylo_pairs = phylo)
  expect_error(compare_modules_paired(rcx), "run detect_modules")
})


test_that("compare_modules_paired errors without phylo_pairs", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- detect_modules(rcx, method = "leiden",
                          objective_function = "modularity")
  expect_error(compare_modules_paired(rcx), "phylo_pairs not set")
})


test_that("classify_hub_conservation dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- detect_modules(rcx, method = "leiden",
                          objective_function = "modularity")
  rcx <- identify_module_hubs(rcx)
  rcx <- classify_hub_conservation(rcx)
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$hub_classification))
})


test_that("classify_hub_conservation errors without hubs", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(classify_hub_conservation(rcx), "run identify_module_hubs")
})


test_that("... passthrough works for find_coexpressologs.rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- find_coexpressologs(rcx, alpha = 0.01)
  expect_true(nrow(rcx$edges) > 0)
})


# --- Pipeline chaining ---

test_that("pipe chaining works through gene-level pipeline", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs) |>
    find_coexpressologs() |>
    find_cliques(min_species = 2L) |>
    classify_cliques(min_species = 2L)

  expect_s3_class(rcx, "rcomplex")
  expect_true(!is.null(rcx$edges))
  expect_true(!is.null(rcx$cliques))
  expect_true(!is.null(rcx$classification))
})


test_that("run_pairwise_comparisons dispatches correctly as wrapper", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- run_pairwise_comparisons(rcx)
  expect_s3_class(rcx, "rcomplex")
  expect_true(!is.null(rcx$edges))
})
