# Acceptance checks for dev/probes/structify.R, run before the grid.
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-stf Rscript dev/probes/structify_check.R
source("dev/probes/structify.R")
tree <- ape::read.tree(
  text = "(((A,B)ab,(C,D)cd)abcd,((E,F)ef,(G,H)gh)efgh)root;"
)
ok <- function(what, pass, extra = "") {
  cat(sprintf("%-4s %s %s\n", if (pass) "PASS" else "FAIL", what, extra))
}
key <- function(a, b) paste(pmin(a, b), pmax(a, b))

# 1. epsilon 0: the edges are the m lowest-ranked pairs
set.seed(1)
a <- structify_root(400, "spatial", p_none = 0)
d <- structify_ranks(a, "spatial")
m <- round(0.03 * length(d))
el <- structify_edges(d, m, 0, a$copy)
ij <- which(lower.tri(matrix(0, 400, 400)), arr.ind = TRUE)[order(d)[1:m], ]
ok(
  "eps 0 = m lowest ranks",
  setequal(key(el$gene1, el$gene2), key(a$copy[ij[, 2]], a$copy[ij[, 1]]))
)

# 2. epsilon 0.5: edge count and structured share over 20 draws
L <- length(d)
st <- replicate(20, {
  e <- structify_edges(d, m, 0.5, a$copy)
  c(nrow(e), mean(e$weight >= (L - m + 1) / L))
})
p1 <- 0.5 + 0.5 * m / L
ok(
  "eps 0.5 edges within 3% of m", abs(mean(st[1, ]) / m - 1) < 0.03,
  sprintf(
    "(mean %.0f, m %d; top-rank share %.3f, expected %.3f)",
    mean(st[1, ]), m, mean(st[2, ]), p1
  )
)

# 3. blocks at epsilon 0: 20 x 75 genes, density 0.03
b <- structify_root(1500, "blocks", p_none = 0)
d <- structify_ranks(b, "blocks")
m <- round(0.03 * length(d))
el <- structify_edges(d, m, 0, b$copy)
blk <- stats::setNames(b$block, b$copy)
within <- blk[el$gene1] == blk[el$gene2]
n_in <- N_BLOCKS * choose(75, 2)
ok(
  "blocks eps 0: all edges within blocks", all(within),
  sprintf(
    "(m %d, within pairs %d, %.3f of within pairs are edges)",
    m, n_in, sum(within) / n_in
  )
)

# 4. a tip with no events reproduces the root; 5. loss on abcd
s <- structify_tree(tree, 300, "blocks", 0.03, 0,
  data.frame(
    gene = "g0001", kind = "loss",
    branch = "abcd"
  ),
  seed = 2
)
s$attr$root$none[1] <- FALSE # make g0001 structured, rebuild
root <- s$attr$root
s <- structify_tree(tree, root, "blocks", 0.03, 0,
  data.frame(
    gene = "g0001", kind = "loss",
    branch = "abcd"
  ),
  seed = 2
)
ok(
  "tip E without events = root attributes",
  identical(s$attr$E[, -1], root[, -1])
)
none1 <- vapply(LETTERS[1:8], function(sp) s$attr[[sp]]$none[1], NA)
ok(
  "loss on abcd: none in A-D, structured in E-H",
  all(none1 == rep(c(TRUE, FALSE), each = 4)),
  paste(names(none1)[none1], collapse = "")
)
deg <- vapply(s$networks, function(x) {
  sum(x$network["g0001", ] >= x$threshold)
}, 1)
cat("     g0001 degree per tip:", deg, "\n")

# 6. dup_switch on ab: two copies in A, B; .2 moved; one copy elsewhere
s <- structify_tree(tree, root, "blocks", 0.03, 0,
  data.frame(
    gene = "g0001", kind = "dup_switch",
    branch = "ab"
  ),
  seed = 2
)
A <- s$attr$A
ok(
  "dup_switch on ab: copies, .1 kept, .2 moved",
  s$truth$copies[1] == "2,2,1,1,1,1,1,1" &&
    A$block[A$copy == "g0001.1"] == root$block[1] &&
    A$block[A$copy == "g0001.2"] != root$block[1]
)
