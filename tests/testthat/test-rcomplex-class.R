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
    compute_network(expr,
      density = 0.1, mr_log_transform = FALSE,
      store_density = 0.2
    )
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

  list(
    species = sp, traits = traits, networks = networks,
    orthologs = orthologs
  )
}


# --- Constructor tests ---

test_that("rcomplex constructor creates valid object", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)

  expect_s3_class(rcx, "rcomplex")
  expect_equal(rcx$species, fix$species)
  expect_equal(rcx$traits, fix$traits)
  expect_equal(length(rcx$species_pairs), 3) # combn(3,2)  # nolint
  expect_null(rcx$phylo_pairs)
  # All result slots are NULL
  expect_null(rcx$edges)
  expect_null(rcx$modules)
  expect_null(rcx$preservation)
  expect_null(rcx$correspondence)
  expect_null(rcx$cliques)
  expect_null(rcx$classification)
})


test_that("rcomplex constructor validates species", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex("A", fix$traits, fix$networks, fix$orthologs),
    "at least 2"
  )
  expect_error(
    rcomplex(c("A", "A"), fix$traits, fix$networks, fix$orthologs),
    "duplicates"
  )
  expect_error(
    rcomplex(1:3, fix$traits, fix$networks, fix$orthologs),
    "character vector"
  )
})


test_that("rcomplex constructor validates traits", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex(
      fix$species, c("a", "b", "c"),
      fix$networks, fix$orthologs
    ),
    "named"
  )
  expect_error(
    rcomplex(
      fix$species, c(SP_A = "annual"),
      fix$networks, fix$orthologs
    ),
    "traits missing species"
  )
})


test_that("rcomplex constructor validates networks", {
  fix <- make_rcx_fixtures()
  bad_nets <- fix$networks
  names(bad_nets)[1] <- "WRONG"
  expect_error(
    rcomplex(fix$species, fix$traits, bad_nets, fix$orthologs),
    "network names must match"
  )
  bad_nets2 <- fix$networks
  bad_nets2[[1]]$threshold <- NULL
  expect_error(
    rcomplex(fix$species, fix$traits, bad_nets2, fix$orthologs),
    "threshold"
  )
})


test_that("rcomplex constructor validates orthologs", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex(
      fix$species, fix$traits, fix$networks,
      data.frame(x = 1)
    ),
    "Species1, Species2, hog"
  )
})


test_that("rcomplex constructor validates species_pairs", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
      species_pairs = list(c("SP_A", "BOGUS"))
    ),
    "unknown species"
  )
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
      species_pairs = list("SP_A")
    ),
    "length-2"
  )
})


test_that("rcomplex constructor validates phylo_pairs", {
  fix <- make_rcx_fixtures()
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
      phylo_pairs = data.frame(x = 1)
    ),
    "sp1, sp2"
  )
  expect_error(
    rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
      phylo_pairs = data.frame(sp1 = "BOGUS", sp2 = "SP_A")
    ),
    "unknown species"
  )
  # Auto-generates pair_name
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
    phylo_pairs = data.frame(sp1 = "SP_A", sp2 = "SP_B")
  )
  expect_equal(rcx$phylo_pairs$pair_name, "SP_A.SP_B")
})


# --- Print / summary tests ---

test_that("print.rcomplex shows pipeline status", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  out <- capture.output(print(rcx))
  expect_true(any(grepl("rcomplex:", out)))
  expect_true(any(grepl("\\[pending\\]", out)))
  expect_true(any(grepl("Preservation:", out)))
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
  rcx <- detect_modules(rcx,
    method = "leiden",
    objective_function = "modularity"
  )
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
  rcx <- detect_modules(rcx,
    method = "leiden",
    objective_function = "modularity"
  )
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


test_that("preservation_paired dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  phylo <- data.frame(sp1 = "SP_A", sp2 = "SP_B")
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
    phylo_pairs = phylo
  )
  rcx <- detect_modules(rcx,
    method = "leiden",
    objective_function = "modularity"
  )
  # classify_preservation() warns that some modules could not be tested on a
  # fixture this small.
  rcx <- suppressWarnings(
    preservation_paired(rcx, min_module_size = 3L, n_perm = 99L, seed = 1)
  )
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$preservation$classification))
  expect_setequal(
    names(rcx$preservation),
    c("classification", "summary", "raw")
  )
  # Preservation is directional: both orientations always run.
  expect_setequal(
    names(rcx$preservation$raw),
    c("SP_A.SP_B", "SP_B.SP_A")
  )
})


test_that("preservation_paired errors without modules", {
  fix <- make_rcx_fixtures()
  phylo <- data.frame(sp1 = "SP_A", sp2 = "SP_B")
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
    phylo_pairs = phylo
  )
  expect_error(preservation_paired(rcx), "run detect_modules")
})


test_that("preservation_paired errors without phylo_pairs", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- detect_modules(rcx,
    method = "leiden",
    objective_function = "modularity"
  )
  expect_error(preservation_paired(rcx), "phylo_pairs not set")
})


test_that("classify_hub_conservation dispatches on rcomplex", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx <- detect_modules(rcx,
    method = "leiden",
    objective_function = "modularity"
  )
  rcx <- identify_module_hubs(rcx)
  rcx <- classify_hub_conservation(rcx)
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$hub_classification))
})


test_that("the container builds a usable module correspondence", {
  fix <- make_rcx_fixtures()
  phylo <- data.frame(sp1 = "SP_A", sp2 = "SP_B")
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs,
    phylo_pairs = phylo
  )
  rcx <- detect_modules(rcx,
    method = "leiden",
    objective_function = "modularity"
  )
  rcx <- suppressWarnings(
    preservation_paired(rcx, min_module_size = 3L, n_perm = 99L, seed = 1)
  )
  rcx <- identify_module_hubs(rcx)
  rcx <- classify_hub_conservation(rcx)

  expect_false(is.null(rcx$correspondence))
  # Keys must be the ALPHABETICALLY SORTED pair -- that is what
  # classify_hub_conservation() looks up.
  expect_equal(names(rcx$correspondence), "SP_A.SP_B")
  expect_true(all(c("module_sp1", "module_sp2", "jaccard", "q.value") %in%
                    names(rcx$correspondence[["SP_A.SP_B"]]$pairs)))
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


# --- Gene-clique taxonomy ---

# The random 20-gene fixture yields 4 co-expressolog edges in total, far
# too few for any 3-node gene clique, so the taxonomy methods get a
# hand-built edge table carried on a real container. HOG1 is a complete
# triangle at alpha_call, HOG2 only enters the graph at alpha_graph 0.9.
make_gcg_edges <- function() {
  data.frame(
    gene1 = c(
      "SP_A_G1", "SP_A_G1", "SP_B_G1",
      "SP_A_G2", "SP_A_G2", "SP_B_G2"
    ),
    gene2 = c(
      "SP_B_G1", "SP_C_G1", "SP_C_G1",
      "SP_B_G2", "SP_C_G2", "SP_C_G2"
    ),
    species1 = c("SP_A", "SP_A", "SP_B", "SP_A", "SP_A", "SP_B"),
    species2 = c("SP_B", "SP_C", "SP_C", "SP_B", "SP_C", "SP_C"),
    hog = c(rep("HOG1", 3), rep("HOG2", 3)),
    q.value = c(0.01, 0.02, 0.03, 0.5, 0.6, 0.7),
    effect_size = c(3, 3, 3, 1, 1, 1),
    stringsAsFactors = FALSE
  )
}


make_gcg_rcx <- function() {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  rcx$edges <- make_gcg_edges()
  rcx
}


test_that("gene_clique_graph dispatches on rcomplex", {
  rcx <- gene_clique_graph(make_gcg_rcx())
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$gene_cliques))
  expect_equal(unique(rcx$gene_cliques$clique_id), "HOG1_1")
  expect_equal(nrow(rcx$gene_cliques), 3L)
})


test_that("gene_clique_graph.rcomplex passes ... to the default", {
  # alpha_graph = 0.9 exceeds the fixture's largest q.value (0.7) on
  # purpose, to build the looser HOG2 clique; the ceiling warning this
  # triggers is expected here, not a regression.
  rcx <- suppressWarnings(gene_clique_graph(make_gcg_rcx(),
    alpha_graph = 0.9, id_prefix = "loose_"
  ))
  expect_setequal(
    unique(rcx$gene_cliques$clique_id),
    c("loose_HOG1_1", "loose_HOG2_1")
  )
})


test_that("gene_clique_graph.rcomplex equals the default method", {
  ed <- make_gcg_edges()
  rcx <- gene_clique_graph(make_gcg_rcx())
  expect_equal(rcx$gene_cliques, gene_clique_graph(ed))
})


test_that("gene_clique_graph errors when edges not yet computed", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(gene_clique_graph(rcx), "run find_coexpressologs")
})


test_that("classify_gene_cliques dispatches on rcomplex", {
  rcx <- gene_clique_graph(make_gcg_rcx())
  # The fixture's q-values sit well under the default alpha_graph = 0.9
  # ceiling on purpose; suppress the expected warning that triggers.
  rcx <- suppressWarnings(classify_gene_cliques(rcx))
  expect_s3_class(rcx, "rcomplex")
  expect_true(is.data.frame(rcx$gene_classification))
  expect_equal(rcx$gene_classification$clique_id, "HOG1_1")
  expect_equal(
    rcx$gene_classification$classification, "complete_conserved"
  )
})


test_that("classify_gene_cliques.rcomplex equals the default method", {
  ed <- make_gcg_edges()
  rcx <- suppressWarnings(
    classify_gene_cliques(gene_clique_graph(make_gcg_rcx()))
  )
  # The default call must be handed the same lineage the container
  # supplies from x$traits, or this compares the container against a
  # run with the lineage tiers switched off.
  expect_equal(
    rcx$gene_classification,
    classify_gene_cliques(
      gene_clique_graph(ed), ed, c("SP_A", "SP_B", "SP_C"),
      lineage = make_rcx_fixtures()$traits
    )
  )
})


test_that("classify_gene_cliques errors on empty slots", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_error(classify_gene_cliques(rcx), "run find_coexpressologs")
  rcx$edges <- make_gcg_edges()
  expect_error(classify_gene_cliques(rcx), "run gene_clique_graph")
})


test_that("the taxonomy slots surface in print and summary", {
  rcx <- suppressWarnings(
    classify_gene_cliques(gene_clique_graph(make_gcg_rcx()))
  )
  out <- capture.output(print(rcx))
  expect_true(any(grepl("Gene cliques:   1 cliques, 3 members", out)))
  expect_true(any(grepl("Gene taxonomy:  1 cliques", out)))

  s <- summary(rcx)
  expect_equal(s$n_gene_cliques, 1L)
  expect_equal(
    as.integer(s$gene_taxonomy_table[["complete_conserved"]]), 1L
  )
  sout <- capture.output(print(s))
  expect_true(any(grepl("Gene cliques: 1", sout)))
  expect_true(any(grepl("Gene clique taxonomy", sout)))
})


test_that("a fresh container has empty taxonomy slots", {
  fix <- make_rcx_fixtures()
  rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
  expect_null(rcx$gene_cliques)
  expect_null(rcx$gene_classification)
  expect_true(all(c("gene_cliques", "gene_classification") %in% names(rcx)))
  out <- capture.output(print(rcx))
  expect_true(any(grepl("Gene cliques:   \\[pending\\]", out)))
  expect_true(any(grepl("Gene taxonomy:  \\[pending\\]", out)))
})


test_that("the taxonomy chains through the container pipeline", {
  rcx <- suppressWarnings(
    make_gcg_rcx() |>
      gene_clique_graph() |>
      classify_gene_cliques()
  )
  expect_s3_class(rcx, "rcomplex")
  expect_true(!is.null(rcx$gene_cliques))
  expect_true(!is.null(rcx$gene_classification))
})


test_that("classify_gene_cliques.rcomplex supplies the container traits", {
  # classify_cliques.rcomplex passes x$traits; this method must too, or
  # the container path runs at lineage = NULL and the lineage_specific
  # and differentiated tiers are silently unreachable on an object that
  # is already holding the trait vector they need.
  rcx <- gene_clique_graph(make_gcg_rcx())
  from_container <- suppressWarnings(
    classify_gene_cliques(rcx)
  )$gene_classification
  explicit <- classify_gene_cliques.default(
    rcx$gene_cliques, rcx$edges, rcx$species,
    lineage = rcx$traits
  )
  expect_equal(from_container, explicit)

  # ...and an explicit lineage in ... still wins, without colliding with
  # the supplied default.
  flat <- stats::setNames(
    rep("one", length(rcx$species)), rcx$species
  )
  overridden <- suppressWarnings(
    classify_gene_cliques(rcx, lineage = flat)
  )
  expect_equal(
    overridden$gene_classification,
    classify_gene_cliques.default(
      rcx$gene_cliques, rcx$edges, rcx$species,
      lineage = flat
    )
  )
})


test_that("classify_gene_cliques.rcomplex actually reaches lineage tiers", {
  # The previous test only tells lineage = NULL apart from lineage =
  # x$traits through count columns; it never shows the container path
  # reach a tier other than complete_conserved. make_rcx_fixtures(3)
  # gives SP_A/SP_C "annual" and SP_B "perennial" -- split a
  # within-lineage-complete clique (HOG4: SP_A-SP_C significant, SP_B
  # absent) from a cross-lineage clique that is present at every pair
  # but significant only within the annual lineage (HOG3).
  rcx <- make_gcg_rcx()
  rcx$edges <- rbind(
    rcx$edges,
    data.frame(
      gene1 = c("SP_A_G3", "SP_A_G3", "SP_B_G3", "SP_A_G4"),
      gene2 = c("SP_B_G3", "SP_C_G3", "SP_C_G3", "SP_C_G4"),
      species1 = c("SP_A", "SP_A", "SP_B", "SP_A"),
      species2 = c("SP_B", "SP_C", "SP_C", "SP_C"),
      hog = c("HOG3", "HOG3", "HOG3", "HOG4"),
      q.value = c(0.5, 0.01, 0.5, 0.01),
      effect_size = c(1, 3, 1, 3),
      stringsAsFactors = FALSE
    )
  )
  # min_size = 2 so the 2-member HOG4 clique survives; alpha_graph = 0.9
  # so the two non-significant HOG3 edges still build the triangle.
  rcx <- suppressWarnings(
    gene_clique_graph(rcx, min_size = 2L, alpha_graph = 0.9)
  )
  from_container <- suppressWarnings(
    classify_gene_cliques(rcx)$gene_classification
  )
  no_lineage <- suppressWarnings(
    classify_gene_cliques(rcx, lineage = NULL)$gene_classification
  )

  expect_true("differentiated" %in% from_container$classification)
  expect_true("lineage_specific" %in% from_container$classification)
  # Without the container's traits, neither lineage-aware tier is
  # reachable: the same cliques fall back to a lineage-blind tier.
  expect_false(any(no_lineage$classification %in%
                     c("differentiated", "lineage_specific")))
})


test_that(
  "supplying the container's traits reaches a tier lineage = NULL cannot",
  {
    # The previous test only ever produces complete_conserved, which is
    # reachable with or without traits and so cannot show the container
    # path actually changes the tier. Five species split 3 annual (A, C, E)
    # / 2 perennial (B, D) lets an all-annual triangle clear min_size = 3
    # while B and D never appear in this HOG's edges at all -- an absence,
    # not a tested-and-rejected pair. With traits supplied that clique is a
    # complete lineage (lineage_specific); with lineage = NULL the same
    # clique is 2 species short of complete_conserved, past max_gap for
    # partial_present, and cannot qualify for any tier.
    fix <- make_rcx_fixtures(n_sp = 5)
    rcx <- rcomplex(fix$species, fix$traits, fix$networks, fix$orthologs)
    rcx$edges <- data.frame(
      gene1 = c("SP_A_G1", "SP_A_G1", "SP_C_G1"),
      gene2 = c("SP_C_G1", "SP_E_G1", "SP_E_G1"),
      species1 = c("SP_A", "SP_A", "SP_C"),
      species2 = c("SP_C", "SP_E", "SP_E"),
      hog = "HOG_LIN", q.value = c(0.01, 0.01, 0.01),
      effect_size = c(3, 3, 3), stringsAsFactors = FALSE
    )
    rcx <- suppressWarnings(gene_clique_graph(rcx))
    expect_equal(nrow(rcx$gene_cliques), 3L)

    with_traits <- suppressWarnings(classify_gene_cliques(rcx))
    expect_identical(
      with_traits$gene_classification$classification, "lineage_specific"
    )

    without_traits <- suppressWarnings(
      classify_gene_cliques(rcx, lineage = NULL)
    )
    expect_false(identical(
      without_traits$gene_classification$classification, "lineage_specific"
    ))
  }
)
