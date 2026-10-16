# Score the finished runs of one PEO experiment.
# Usage: Rscript score.R <exp>
# Reads runs/<run>/ from run.R; writes score/<run>_{species,tiers,pairs,
# branches,hogs}.tsv, the experiment tables and score/summary.tsv.
# score/sacct.tsv (written by jobs.sh) adds wall time and MaxRSS.
suppressPackageStartupMessages(library(rcomplex))
args <- commandArgs(TRUE)
EXP <- args[1]
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
ROOTHOGS <- Sys.getenv("PEO_ROOTHOGS", paste0(
  "/mnt/project/FjellheimLab/martpali/AnnualPerennial/peo-fastoma/",
  "output/RootHOGs.tsv"
))
ALPHA <- 0.1
CONS <- c("complete_conserved", "partial_significant", "partial_present")
TIERS <- rcomplex:::.gcg_tiers
FOCAL <- c(C3vC4 = "C4", nodulation = "nodulating")
E <- file.path(ROOT, paste0("complex_", EXP))
D <- file.path(E, "inputs")
S <- file.path(E, "score")
dir.create(S, showWarnings = FALSE)
info <- read.delim(file.path(D, "species.tsv"))
tree <- ape::read.tree(file.path(D, "tree.nwk"))
ortho <- read_orthologs(file.path(D, "orthologs.tsv"))
runs <- intersect(
  c("leaf", "root", "organ_block"), list.files(file.path(E, "runs"))
)
tsv <- function(x, f) {
  data.table::fwrite(x, file.path(S, f),
    sep = "\t", na = "NA", quote = FALSE
  )
}
rd <- function(run, f) {
  p <- file.path(E, "runs", run, f)
  if (file.exists(p)) data.table::fread(p, data.table = FALSE)
}

# Species present in each clique: all species minus missing_species.
present <- function(cls, sp) {
  m <- cls$missing_species
  m[is.na(m)] <- ""
  lapply(strsplit(m, ",", fixed = TRUE), function(x) setdiff(sp, x))
}

# Tip pairs whose parent node has exactly these two tips as children.
sisters <- function(tree) {
  e <- tree$edge
  n <- ape::Ntip(tree)
  out <- character(0)
  for (p in unique(e[, 1L])) {
    ch <- e[e[, 1L] == p, 2L]
    if (length(ch) == 2L && all(ch <= n)) {
      out <- c(out, paste(sort(tree$tip.label[ch]), collapse = "|"))
    }
  }
  out
}

score_run <- function(run) {
  res <- readRDS(file.path(E, "runs", run, "res.rds"))
  nets <- rd(run, "networks.tsv")
  eh <- rd(run, "edge_history_cliques.tsv")
  tp <- rd(run, "trait_pass.tsv")
  sp <- info$species[info$run == run]
  cls <- res$classification
  e <- res$edges
  called <- !is.na(e$q_value) & e$q_value < ALPHA

  # species
  cl1 <- res$cliques[res$cliques$run == 1L, c("clique_id", "hog")]
  mem <- res$members[res$members$run == 1L, c("clique_id", "species")]
  mem$hog <- cl1$hog[match(mem$clique_id, cl1$clique_id)]
  spt <- do.call(rbind, lapply(sp, function(s) {
    g <- rownames(readRDS(file.path(D, run, paste0(s, ".rds"))))
    h <- table(ortho$hog[ortho$species == s & ortho$gene %in% g])
    inc <- length(unique(mem$hog[mem$species == s]))
    data.frame(
      species = s, n_runs = ncol(readRDS(file.path(
        D, run,
        paste0(s, ".rds")
      ))), genes = length(g), hogs = length(h),
      share_multi_copy = round(mean(h > 1), 3),
      r_threshold = nets$r_threshold[nets$species == s],
      edges_tested = sum(e$species1 == s | e$species2 == s),
      edges_called = sum(called & (e$species1 == s | e$species2 == s)),
      hogs_in_clique = inc, share_in_clique = round(inc / length(h), 3)
    )
  }))
  tsv(spt, paste0(run, "_species.tsv"))

  # tiers
  hc <- function(x) table(factor(x$hog_class[!duplicated(x$hog)], TIERS))
  cc <- function(x) table(factor(x$classification, TIERS))
  c0 <- res$classification_null
  tiers <- data.frame(
    tier = TIERS, hogs = as.vector(hc(cls)), hogs_null = as.vector(hc(c0)),
    cliques = as.vector(cc(cls)), cliques_null = as.vector(cc(c0))
  )
  tsv(tiers, paste0(run, "_tiers.tsv"))

  # pairs, on the strict graph
  pr <- present(cls[cls$run == 1L, ], sp)
  m <- Matrix::sparseMatrix(
    i = rep(seq_along(pr), lengths(pr)), j = match(unlist(pr), sp),
    x = 1, dims = c(length(pr), length(sp)), dimnames = list(NULL, sp)
  )
  obs <- as.matrix(Matrix::crossprod(m))
  ij <- which(upper.tri(obs), arr.ind = TRUE)
  n_sp <- diag(obs)
  pairs <- data.frame(
    species1 = sp[ij[, 1L]], species2 = sp[ij[, 2L]], cliques = obs[ij],
    expected = n_sp[ij[, 1L]] * n_sp[ij[, 2L]] / max(1, length(pr))
  )
  pairs$ratio <- pairs$cliques / pairs$expected
  pairs$rank <- rank(-pairs$ratio, ties.method = "min")
  key <- apply(pairs[, 1:2], 1L, function(x) paste(sort(x), collapse = "|"))
  pairs$sister <- key %in% sisters(tree)
  pairs$top_quartile <- pairs$rank <= ceiling(nrow(pairs) / 4)
  pairs <- pairs[order(pairs$rank), ]
  tsv(pairs, paste0(run, "_pairs.tsv"))

  # branches
  br <- rd(run, "edge_history_branches.tsv")
  cl <- clades_from_tree(tree, 1L)
  br$tips <- vapply(br$branch, function(b) paste(cl[[b]], collapse = ","), "")
  tsv(br, paste0(run, "_branches.tsv"))

  # hogs: best clique per hog
  cls$n_gain <- eh$n_gain
  cls$n_loss <- eh$n_loss
  cls$species <- vapply(present(cls, sp), paste, "", collapse = ",")
  o <- order(match(cls$hog, unique(cls$hog)), match(
    cls$classification,
    TIERS
  ), -cls$n_species, cls$mean_q)
  best <- cls[o, ][!duplicated(cls$hog[o]), ]
  hogs <- data.frame(
    hog = best$hog, graph = best$run, clique_id = best$clique_id,
    tier = best$hog_class, best_class = best$classification,
    clade = best$clade, n_species = best$n_species, species = best$species,
    n_gain = best$n_gain, n_loss = best$n_loss,
    trait_class = if (is.null(tp)) {
      NA
    } else {
      tp$hog_class[match(best$hog, tp$hog)]
    }
  )
  tsv(hogs, paste0(run, "_hogs.tsv"))
  tm <- rd(run, "timing.tsv")
  summ <- data.frame(
    run = run, species = length(sp), genes = sum(spt$genes),
    hogs_tested = length(unique(e$hog)), edges_tested = nrow(e),
    edges_called = sum(called),
    edges_called_null = sum(res$edges_null$q_value < ALPHA, na.rm = TRUE),
    cliques_strict = sum(res$cliques$run == 1L),
    hogs_conserved = sum(hogs$tier %in% CONS),
    hogs_complete = sum(hogs$tier == "complete_conserved"),
    hogs_lineage = sum(hogs$tier == "lineage_specific"),
    hogs_conserved_null = sum(hc(c0)[CONS]),
    sister_top_quartile = sum(pairs$sister & pairs$top_quartile),
    sisters = sum(pairs$sister),
    rcomplex_s = tm$seconds[tm$step == "rcomplex"],
    max_mb = max(tm$max_mb), job_id = tm$job_id[1L]
  )
  list(hogs = hogs, cls = cls, sp = sp, summ = summ, tested = unique(e$hog))
}

out <- sapply(runs, score_run, simplify = FALSE)

# LeafvRoot: conserved in leaf, root, both, neither; organ_block alongside.
if (all(c("leaf", "root") %in% runs)) {
  u <- union(out$leaf$tested, out$root$tested)
  tier <- function(run) out[[run]]$hogs$tier[match(u, out[[run]]$hogs$hog)]
  lf <- tier("leaf") %in% CONS
  rt <- tier("root") %in% CONS
  cat_ <- ifelse(lf & rt, "both", ifelse(lf, "leaf", ifelse(rt, "root",
    "neither"
  )))
  ob <- if ("organ_block" %in% runs) tier("organ_block") %in% CONS else NA
  tab <- as.data.frame(table(category = cat_, organ_block_conserved = ob))
  tsv(tab, "leaf_root.tsv")
  print(tab)
  top <- function(a, b, k) {
    h <- out[[a]]$hogs
    h <- h[h$tier %in% CONS & !h$hog %in% out[[b]]$hogs$hog[
      out[[b]]$hogs$tier %in% CONS
    ], ]
    h <- h[order(match(h$tier, TIERS), -h$n_species), ][seq_len(min(
      k,
      nrow(h)
    )), ]
    data.frame(
      hog = h$hog, n_species = h$n_species, species = h$species,
      tier_leaf = tier1(h$hog, "leaf"), tier_root = tier1(h$hog, "root"),
      tier_organ_block = tier1(h$hog, "organ_block")
    )
  }
  tier1 <- function(h, run) {
    if (is.null(out[[run]])) {
      NA
    } else {
      out[[run]]$hogs$tier[match(h, out[[run]]$hogs$hog)]
    }
  }
  tsv(top("leaf", "root", 30L), "leaf_not_root.tsv")
  tsv(top("root", "leaf", 30L), "root_not_leaf.tsv")
}

# Trait experiments: tier and trait class of anchor genes.
ANCHORS <- list(
  C3vC4 = data.frame(
    symbol = c(
      "PPC1", "PPC2", "PPC3", "PPDK", "PPDK-RP", "NADP-ME1", "NADP-ME2",
      "NADP-ME3", "NADP-ME4", "NAD-ME1", "NAD-ME2", "PCK1", "BCA1", "BCA2",
      "ASP5", "ALAAT1", "DiT1", "DiT2.1", "PPT1", "BASS2", "NHD1",
      "RBCS1A", "RCA", "PRK", "SBPase", "GAPA", "FBA1"
    ),
    id = c(
      "AT1G53310", "AT2G42600", "AT3G14940", "AT4G15530", "AT4G21210",
      "AT2G19900", "AT5G11670", "AT5G25880", "AT1G79750", "AT2G13560",
      "AT4G00570", "AT4G37870", "AT3G01500", "AT5G14740", "AT4G31990",
      "AT1G17290", "AT5G12860", "AT5G64290", "AT5G33320", "AT2G26900",
      "AT3G19490", "AT1G67090", "AT2G39730", "AT1G32060", "AT3G55800",
      "AT3G26650", "AT2G21330"
    ),
    role = rep(c("C4_cycle", "Calvin_control"), c(21L, 6L)),
    source = "AGI"
  ),
  # MTR_ ids from UniProt (organism 3880, ordered locus names); NIN, DMI1,
  # DMI2 and RPG have no locus there and take their Mt4.0 literature ids.
  # ENOD40 has no protein locus and is not placed.
  nodulation = data.frame(
    symbol = c(
      "NIN", "NSP1", "NSP2", "ERN1", "ENOD11", "NFP", "LYK3", "DMI1",
      "DMI2", "DMI3", "IPD3", "VPY", "RPG", "LB1_a", "LB1_b", "SUNN",
      "RAM1", "PT4"
    ),
    id = c(
      "MTR_5g099060", "MTR_8g020840", "MTR_3g072710", "MTR_7g085810",
      "MTR_3g415670", "MTR_5g019040", "MTR_5g086130", "MTR_2g005870",
      "MTR_5g030920", "MTR_8g043970", "MTR_5g026850", "MTR_6g027840",
      "MTR_1g090807", "MTR_1g011540", "MTR_5g066070", "MTR_4g070970",
      "MTR_7g027190", "MTR_1g028600"
    ),
    role = c(
      "nodulation", "other", "other", "other", "nodulation", "nodulation",
      "other", "common", "common", "common", "common", "common",
      "nodulation", "nodulation", "nodulation", "other", "common",
      "common"
    ),
    source = c(
      "literature", rep("uniprot", 6L), "literature", "literature",
      rep("uniprot", 3L), "literature", rep("uniprot", 5L)
    )
  )
)

if (EXP %in% names(ANCHORS)) {
  a <- ANCHORS[[EXP]]
  r <- out[[1L]]
  strip <- function(x) sub("\\.[0-9]+$", "", x)
  rh <- data.table::fread(ROOTHOGS, data.table = FALSE)
  rh$id <- strip(rh$Protein)
  focal <- info$species[info$run == runs[1L] & info$trait == FOCAL[[EXP]]]
  focal <- unique(focal)
  sp <- r$sp
  rows <- lapply(seq_len(nrow(a)), function(k) {
    root <- unique(rh$RootHOG[rh$id == a$id[k]])
    own <- ortho$hog[strip(ortho$gene) == a$id[k]]
    hg <- if (length(own) > 0L) {
      own
    } else {
      unique(ortho$hog[ortho$gene %in% rh$Protein[rh$RootHOG %in% root]])
    }
    if (length(hg) == 0L) {
      return(data.frame(a[k, ],
        root_hog = paste(root, collapse = ","),
        hog = NA, n_genes = 0L, tier = NA, trait_class = NA, n_gain = NA,
        n_loss = NA, focal_in_best = NA, other_in_best = NA,
        focal_own_clique = NA, focal_all_own = NA, absent_species = NA
      ))
    }
    do.call(rbind, lapply(hg, function(h) {
      b <- r$hogs[match(h, r$hogs$hog), ]
      bs <- strsplit(b$species, ",", fixed = TRUE)[[1L]]
      pc <- present(r$cls[r$cls$hog == h, ], sp)
      big <- pc[lengths(pc) >= min(3L, length(focal))]
      data.frame(a[k, ],
        root_hog = paste(root, collapse = ","), hog = h,
        n_genes = sum(ortho$hog == h), tier = b$tier,
        trait_class = b$trait_class, n_gain = b$n_gain, n_loss = b$n_loss,
        focal_in_best = sum(bs %in% focal),
        other_in_best = sum(!bs %in% focal & bs %in% sp),
        focal_own_clique = any(vapply(big, function(p) all(p %in% focal), NA)),
        focal_all_own = any(vapply(big, setequal, NA, focal)),
        absent_species = paste(setdiff(sp, ortho$species[ortho$hog == h]),
          collapse = ","
        )
      )
    }))
  })
  anc <- do.call(rbind, rows)
  rownames(anc) <- NULL
  tsv(anc, "anchors.tsv")
  print(anc[, c(
    "symbol", "role", "hog", "tier", "trait_class",
    "focal_in_best", "other_in_best", "focal_own_clique"
  )], row.names = FALSE)
  cat(
    "placed:", length(unique(anc$symbol[!is.na(anc$hog)])), "of",
    nrow(a), "anchors\n"
  )
}

summ <- do.call(rbind, lapply(out, `[[`, "summ"))
sa <- file.path(S, "sacct.tsv")
if (file.exists(sa) && file.size(sa) > 0L) {
  s <- read.delim(sa,
    sep = "|", header = FALSE,
    col.names = c("id", "state", "elapsed", "maxrss")
  )
  s$base <- sub("\\..*$", "", s$id)
  mb <- function(x) {
    u <- c(K = 1 / 1024, M = 1, G = 1024, T = 1024^2)
    v <- suppressWarnings(as.numeric(sub("[KMGT]$", "", x)))
    v * ifelse(grepl("[KMGT]$", x), u[substring(x, nchar(x))], 1 / 1024^2)
  }
  s$mb <- mb(s$maxrss)
  main <- s[s$id == s$base, ]
  summ$state <- main$state[match(summ$job_id, main$base)]
  summ$elapsed <- main$elapsed[match(summ$job_id, main$base)]
  summ$maxrss_gb <- round(vapply(as.character(summ$job_id), function(j) {
    max(c(0, s$mb[s$base == j]), na.rm = TRUE) / 1024
  }, 1), 1)
}
tsv(summ, "summary.tsv")
print(summ, row.names = FALSE)
