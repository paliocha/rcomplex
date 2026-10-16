# Run rcomplex() and edge_history() on structify toy data; score the
# recovery of planted gains, losses, switches and duplications.
# Usage, from the worktree root:
#   R_LIBS=~/.cache/rcomplex/lib-stf Rscript dev/probes/structify_run.R \
#     [structure] [epsilon] [seed]       (defaults blocks 0.2 1)
#   ... structify_run.R --all            (grid, writes summary.tsv)
# Per run: <key>.rds (cache, skipped if present), <key>_hogs.tsv,
# <key>_confusion.tsv in OUT. Design: 8 species, 1,500 genes,
# density 0.03; 60 each of gain, loss, switch, dup, dup_switch, plus 30
# dups at the root.
source("dev/probes/structify.R")
OUT <- "/Users/martin/Documents/R/rcomplex/prepare_data/cache/structify"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TREE <- ape::read.tree(
  text = "(((A,B)ab,(C,D)cd)abcd,((E,F)ef,(G,H)gh)efgh)root;"
)
SP <- TREE$tip.label
CLADES <- clades_from_tree(TREE)
N <- 1500L
DENSITY <- 0.03
MAX_WALL <- 300

simulate <- function(structure, epsilon, seed) {
  root <- withr::with_seed(seed, structify_root(N, structure))
  ev <- structify_events(TREE, root, 60, 60, 60, 60, 60, 30, seed = seed)
  structify_tree(TREE, root, structure, DENSITY, epsilon, ev, seed = seed)
}

run_one <- function(structure, epsilon, seed) {
  key <- sprintf("%s_e%.1f_s%d", structure, epsilon, seed)
  f <- file.path(OUT, paste0(key, ".rds"))
  if (file.exists(f)) {
    return(readRDS(f))
  }
  t0 <- Sys.time()
  sim <- simulate(structure, epsilon, seed)
  rc <- function(o) {
    rcomplex(
      networks = sim$networks, orthologs = o, clades = CLADES,
      density = DENSITY, alpha = 0.1, null = FALSE, n_cores = 4,
      seed = 1
    )
  }
  res <- rc(sim$orthologs)
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  eh <- edge_history(res$classification, TREE, SP)
  # (g) null: gene labels permuted within species
  o <- sim$orthologs
  o$gene <- withr::with_seed(seed, ave(o$gene, o$species, FUN = sample))
  nul <- rc(o)
  res$networks <- NULL
  out <- list(
    key = key, structure = structure, epsilon = epsilon, seed = seed,
    truth = sim$truth, events = sim$events, res = res, eh = eh,
    null_cliques = sum(nul$cliques$run == 1L),
    null_hogs = length(unique(nul$cliques$hog[nul$cliques$run == 1L])),
    wall = wall,
    wall_total = as.numeric(Sys.time() - t0, units = "secs")
  )
  saveRDS(out, f)
  out
}

# Branches that read as the planted one: b, plus its sister when b hangs
# from the root (a root tie resolves to absent, so a loss on abcd reads
# as a gain on efgh).
.parent <- function(tree) {
  nm <- rcomplex:::.node_names(tree)
  stats::setNames(nm[tree$edge[, 1]], nm[tree$edge[, 2]])
}
PAR <- .parent(TREE)
sister <- function(b) {
  vapply(b, function(x) {
    if (is.na(x) || x == "root" || PAR[[x]] != "root") {
      return(NA_character_)
    }
    setdiff(names(PAR)[PAR == "root"], x)
  }, "")
}
has <- function(lst, b) {
  !is.na(b) & mapply(function(l, x) x %in% strsplit(l, ",")[[1]], lst, b)
}

score <- function(r) {
  cl <- cbind(
    r$res$classification,
    r$eh$cliques[c("root_tied", "n_gain", "n_loss", "gains", "losses")]
  )
  c1 <- cl[cl$run == 1L, ]
  n_cl <- table(c1$hog)
  c1 <- c1[order(c1$hog, -c1$n_members, c1$mean_q), ]
  c1 <- c1[!duplicated(c1$hog), ]
  tr <- r$truth
  h <- merge(tr, c1[c(
    "hog", "clique_id", "hog_class", "n_members",
    "multiplicity", "capped", "root_tied", "n_gain",
    "n_loss", "gains", "losses"
  )],
  by = "hog", all.x = TRUE
  )
  h$n_cliques <- as.integer(n_cl[h$hog])
  h$n_cliques[is.na(h$n_cliques)] <- 0L
  h$hog_class[is.na(h$hog_class)] <- "no_clique"
  sis <- sister(h$branch)
  h$hit_loss <- has(h$losses, h$branch) | has(h$gains, sis)
  h$hit_gain <- has(h$gains, h$branch) | has(h$losses, sis)
  # members of the chosen clique in the species that carry two copies
  m1 <- r$res$members[r$res$members$run == 1L, ]
  mem <- merge(h[c("hog", "clique_id", "copies", "kind")], m1,
    by = "clique_id"
  )
  two <- mapply(
    function(cp, s) strsplit(cp, ",")[[1]][match(s, SP)] == "2",
    mem$copies, as.character(mem$species)
  )
  mem <- mem[two, ]
  dsw <- mem[mem$kind == "dup_switch", ]
  dup <- mem[mem$kind == "dup", ]
  twin <- tapply(dup$n_copies == 2L, dup$hog, all)

  kinds <- c("none", "gain", "loss", "switch", "dup", "dup_switch")
  k <- function(x) h$kind == x
  rec <- function(x, hit) mean(hit[k(x)])
  # (b) precision over every reported change
  rep_br <- function(col, dir) {
    s <- strsplit(ifelse(is.na(h[[col]]), "", h[[col]]), ",")
    data.frame(hog = rep(h$hog, lengths(s)), br = unlist(s),
               dir = rep(dir, sum(lengths(s))))
  }
  ch <- rbind(rep_br("gains", "gain"), rep_br("losses", "loss"))
  ch <- merge(ch, h[c("hog", "kind", "branch")], by = "hog")
  ch$sis <- sister(ch$branch)
  ok_br <- ch$kind %in% c("gain", "loss", "switch") &
    (ch$br == ch$branch | ch$br %in% ch$sis)
  stem <- h$branch %in% names(CLADES)
  # root children: the only stems whose complement is also a clade
  rc <- !is.na(sister(h$branch))
  spec <- c("trait_specific", "lineage_specific")
  tier <- function(x, cls) mean(h$hog_class[x] %in% cls)
  fe <- function(x) mean((h$n_gain + h$n_loss)[x & !is.na(h$n_gain)])
  reach <- vapply(kinds, function(x) sum(k(x) & h$n_cliques > 0), 1L)
  row <- data.frame(
    structure = r$structure, epsilon = r$epsilon, seed = r$seed,
    wall = round(r$wall), wall_null = round(r$wall_total - r$wall),
    cliques = nrow(cl[cl$run == 1L, ]),
    rec_loss = rec("loss", h$hit_loss), rec_gain = rec("gain", h$hit_gain),
    rec_switch = rec("switch", h$hit_loss),
    rec_gain_rc = mean(h$hit_gain[k("gain") & rc]),
    precision = mean(ok_br), n_changes = nrow(ch),
    fe_none = fe(k("none")), fe_dup = fe(k("dup") | k("dup_switch")),
    tier_none = tier(k("none"), "complete_conserved"),
    tier_loss_stem = tier(
      k("loss") & stem,
      c("trait_specific", "lineage_specific")
    ),
    tier_switch_stem = tier(k("switch") & stem, "differentiated"),
    tier_gain_rc = tier(k("gain") & rc, spec),
    tier_loss_rc = tier(k("loss") & rc, spec),
    tier_switch_rc = tier(k("switch") & rc, "differentiated"),
    tier_dup = tier(k("dup"), "complete_conserved"),
    tier_dupsw = tier(k("dup_switch"), "complete_conserved"),
    root_tied = mean(h$root_tied, na.rm = TRUE),
    copy_recall = mean(grepl("\\.1$", dsw$gene)),
    twin_frac = mean(twin),
    dup_split = sum(k("dup") & h$n_cliques > 1L),
    capped = length(unique(r$res$capped$hog[r$res$capped$run == 1L])),
    null_cliques = r$null_cliques
  )
  names(reach) <- paste0("reach_", kinds)
  row <- cbind(row, t(reach))
  conf <- as.data.frame(table(
    kind = factor(h$kind, kinds),
    hog_class = h$hog_class
  ))
  list(row = row, hogs = h, conf = conf[conf$Freq > 0, ])
}

do_run <- function(structure, epsilon, seed) {
  r <- run_one(structure, epsilon, seed)
  s <- score(r)
  w <- function(x, suf) {
    utils::write.table(x, file.path(OUT, paste0(r$key, suf)),
      sep = "\t",
      quote = FALSE, row.names = FALSE
    )
  }
  w(s$hogs, "_hogs.tsv")
  w(s$conf, "_confusion.tsv")
  cat(sprintf(
    "%s: %.0f s (+%.0f s null)\n", r$key, r$wall,
    r$wall_total - r$wall
  ))
  s
}

args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--all")) {
  grid <- expand.grid(seed = 1:3, epsilon = c(0, 0.2, 0.5, 0.8),
                      structure = STRUCTURES, stringsAsFactors = FALSE)
  rows <- list()
  for (g in seq_len(nrow(grid))) {
    s <- do_run(grid$structure[g], grid$epsilon[g], grid$seed[g])
    rows[[g]] <- s$row
    if (s$row$wall > MAX_WALL) {
      cat("run over", MAX_WALL, "s; grid stopped\n")
      break
    }
  }
  sm <- do.call(rbind, rows)
  utils::write.table(sm, file.path(OUT, "summary.tsv"),
    sep = "\t",
    quote = FALSE, row.names = FALSE
  )
  print(format(sm, digits = 2), row.names = FALSE)
} else {
  a <- c(args, c("blocks", "0.2", "1")[-seq_along(args)])
  s <- do_run(a[1], as.numeric(a[2]), as.integer(a[3]))
  print(t(s$row))
  print(s$conf, row.names = FALSE)
}
