# make_fixture.R - synthetic 8-module data, seed 20260901
set.seed(20260901)
n_mod <- 8; n_samp <- 30; n_extra <- 60
make_sp <- function(prefix, per_mod) {
  latent <- matrix(rnorm(n_mod * n_samp), n_mod, n_samp)
  rows <- lapply(seq_len(n_mod), function(m)
    t(replicate(per_mod, 0.9 * latent[m, ] + rnorm(n_samp, sd = 0.7))))
  x <- rbind(do.call(rbind, rows), matrix(rnorm(n_extra * n_samp), n_extra)) + 8
  dimnames(x) <- list(paste0(prefix, "_", seq_len(nrow(x))), paste0("s", seq_len(n_samp)))
  list(x = x, module = c(rep(seq_len(n_mod), each = per_mod), rep(NA, n_extra)))
}
sp1 <- make_sp("sp1", 40); sp2 <- make_sp("sp2", 44)
g1 <- rownames(sp1$x)[1:320]
pos2 <- unlist(lapply(seq_len(n_mod), function(m) (m - 1) * 44 + 1:40))
g2 <- rownames(sp2$x)[pos2]
og <- data.frame(gene = c(g1, g2), species = rep(c("sp1", "sp2"), each = 320),
                 Ortholog_Group = rep(paste0("OG", 1:320), 2))
para2 <- rownames(sp2$x)[unlist(lapply(seq_len(n_mod), function(m) (m - 1) * 44 + 41:44))]
para_og <- paste0("OG", unlist(lapply(seq_len(n_mod), function(m) (m - 1) * 40 + 1:4)))
og <- rbind(og, data.frame(gene = para2, species = "sp2", Ortholog_Group = para_og))
sh <- sample(320, 64); og$gene[320 + sh] <- sample(og$gene[320 + sh])
write.table(og, "ortho_long.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
s1 <- og[og$species == "sp1", c("gene", "Ortholog_Group")]; names(s1) <- c("Species1", "hog")
s2 <- og[og$species == "sp2", c("gene", "Ortholog_Group")]; names(s2) <- c("Species2", "hog")
pairs <- merge(s1, s2, by = "hog")[, c("Species1", "Species2", "hog")]
write.table(pairs, "ortho_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
w <- function(x, f) write.table(data.frame(Genes = rownames(x), x, check.names = FALSE), f,
                                sep = "\t", quote = FALSE, row.names = FALSE)
w(sp1$x, "sp1_expr.tsv"); w(sp2$x, "sp2_expr.tsv")
