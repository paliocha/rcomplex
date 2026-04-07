# Shared test fixtures for clique pipeline tests
# (perturbation test, intensity test, etc.)

# Build a 2-species clique test fixture.
# 20 genes per species. A1 strongly connected to A2..A10,
# moderate to A11, A12 connected to A13..A15. SP_B mirrors SP_A.
# 1:1 orthologs.
make_clique_fixture <- function(n_genes = 20L) {
  ga <- paste0("A", seq_len(n_genes))
  gb <- paste0("B", seq_len(n_genes))

  make_net <- function(genes) {
    n <- length(genes)
    mat <- matrix(0, n, n, dimnames = list(genes, genes))
    for (i in 2:min(10, n)) mat[1, i] <- mat[i, 1] <- 10
    if (n >= 11) mat[1, 11] <- mat[11, 1] <- 3
    if (n >= 15) for (i in 13:15) mat[12, i] <- mat[i, 12] <- 4
    mat
  }

  networks <- list(
    SP_A = list(network = make_net(ga), threshold = 2.0),
    SP_B = list(network = make_net(gb), threshold = 2.0)
  )

  orthologs <- data.frame(
    Species1 = ga, Species2 = gb,
    hog = paste0("HOG", seq_len(n_genes)),
    stringsAsFactors = FALSE
  )

  target <- c("SP_A", "SP_B")
  edges <- find_coexpressologs(networks, orthologs, method = "analytical")
  cliques <- find_cliques(edges, target, min_species = 2L)

  list(networks = networks, orthologs = orthologs, edges = edges,
       cliques = cliques, target_species = target)
}


# Build a 3-species clique test fixture.
make_clique_fixture_3sp <- function(n_genes = 15L) {
  ga <- paste0("A", seq_len(n_genes))
  gb <- paste0("B", seq_len(n_genes))
  gc <- paste0("C", seq_len(n_genes))

  make_net <- function(genes) {
    n <- length(genes)
    mat <- matrix(0, n, n, dimnames = list(genes, genes))
    for (i in 2:min(8, n)) mat[1, i] <- mat[i, 1] <- 10
    if (n >= 9) mat[1, 9] <- mat[9, 1] <- 3
    mat
  }

  networks <- list(
    SP_A = list(network = make_net(ga), threshold = 2.0),
    SP_B = list(network = make_net(gb), threshold = 2.0),
    SP_C = list(network = make_net(gc), threshold = 2.0)
  )

  orthologs <- rbind(
    data.frame(Species1 = ga, Species2 = gb,
               hog = paste0("HOG", seq_len(n_genes)),
               stringsAsFactors = FALSE),
    data.frame(Species1 = ga, Species2 = gc,
               hog = paste0("HOG", seq_len(n_genes)),
               stringsAsFactors = FALSE),
    data.frame(Species1 = gb, Species2 = gc,
               hog = paste0("HOG", seq_len(n_genes)),
               stringsAsFactors = FALSE)
  )

  target <- c("SP_A", "SP_B", "SP_C")
  edges <- find_coexpressologs(networks, orthologs, method = "analytical")
  cliques <- find_cliques(edges, target, min_species = 2L)

  list(networks = networks, orthologs = orthologs, edges = edges,
       cliques = cliques, target_species = target)
}
