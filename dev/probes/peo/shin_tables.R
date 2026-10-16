# Tidy the Shin et al. 2021 supplementary tables (NAR 49:e3).
# Usage: Rscript shin_tables.R <supp2 dir> <out dir>
# Writes modules.tsv, genesets.tsv, muscari.tsv, nodulation_curated.tsv
# and geneset_patterns.tsv; prints genes per species against the sheets'
# gene-count rows.
suppressPackageStartupMessages(library(readxl))
args <- commandArgs(TRUE)
IN <- args[1]
OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SP <- c(
  P.patens = "Physcomitrella_patens", O.sativa = "Oryza_sativa",
  Z.mays = "Zea_mays", S.tuberosum = "Solanum_tuberosum",
  M.truncatula = "Medicago_truncatula", A.thaliana = "Arabidopsis_thaliana"
)
sheet <- function(f, s) {
  x <- read_excel(file.path(IN, f), s,
    col_names = FALSE, col_types = "text", .name_repair = "minimal"
  )
  as.matrix(x)
}
tsv <- function(x, f) {
  write.table(x, file.path(OUT, f),
    sep = "\t", quote = FALSE,
    row.names = FALSE
  )
}

# One block per "Species" cell: name above-right, species names to the
# right, the gene-count row below, gene ids down each species column.
blocks <- function(m, set_col) {
  hit <- which(m == "Species", arr.ind = TRUE)
  out <- list()
  chk <- list()
  for (k in seq_len(nrow(hit))) {
    r <- hit[k, 1L]
    c <- hit[k, 2L]
    set <- m[r - 1L, c + 1L]
    for (j in c + seq_len(6L)) {
      sp <- SP[[m[r, j]]]
      g <- m[(r + 2L):nrow(m), j]
      end <- which(is.na(g))[1L]
      g <- if (is.na(end)) g else g[seq_len(end - 1L)]
      out[[length(out) + 1L]] <- data.frame(
        species = rep(sp, length(g)), gene = g,
        set = rep(set, length(g))
      )
      chk[[length(chk) + 1L]] <- data.frame(
        set = set, species = sp,
        n = length(g), n_sheet = as.integer(m[r + 1L, j])
      )
    }
  }
  chk <- do.call(rbind, chk)
  bad <- chk[chk$n != chk$n_sheet, ]
  if (nrow(bad) > 0L) {
    print(head(bad))
    stop(nrow(bad), " blocks differ from their gene-count row")
  }
  names(out) <- NULL
  res <- do.call(rbind, out)
  names(res)[3L] <- set_col
  res
}
report <- function(x, what) {
  cat("==", what, nrow(x), "rows\n")
  print(table(x$species))
}

mod <- blocks(
  sheet("S4.Arboretum-proteome_result.xlsx", "S4B.Gene list"),
  "module"
)
mod$module <- as.integer(sub("^Module", "", mod$module))
mod <- mod[, c("species", "gene", "module")]
report(mod, "modules")
tsv(mod, "modules.tsv")

S5 <- "S5.Clade-specific_gene_sets_result.xlsx"
pat <- function(s, method) {
  m <- sheet(S5, s)
  h <- which(m == "ppa", arr.ind = TRUE)[1L, ]
  cols <- h[2L] + 0:10
  rows <- (h[1L] + 1L):nrow(m)
  rows <- rows[!is.na(m[rows, 1L])]
  p <- matrix(as.integer(m[rows, cols]), length(rows))
  colnames(p) <- m[h[1L], cols]
  type <- apply(p, 1L, function(v) {
    nz <- unique(v[v != 0L])
    if (length(nz) > 1L) {
      "change"
    } else if (any(v == 0L)) {
      "loss"
    } else {
      "conserved"
    }
  })
  data.frame(
    method = method, set = m[rows, 1L], p,
    pattern_type = type, check.names = FALSE
  )
}
pats <- rbind(
  pat("S5A.Summary(rule-based)", "rule"),
  pat("S5B.Summary(de novo)", "de_novo")
)
cat("== pattern types\n")
print(table(pats$method, pats$pattern_type))
tsv(pats, "geneset_patterns.tsv")
gs <- rbind(
  cbind(method = "rule", blocks(
    sheet(S5, "S5C.Gene list(rule-based)"),
    "set"
  )),
  cbind(method = "de_novo", blocks(
    sheet(S5, "S5D.Gene list(de novo)"),
    "set"
  ))
)
key <- paste(gs$method, gs$set)
miss <- setdiff(unique(key), paste(pats$method, pats$set))
if (length(miss) > 0L) stop("sets without a pattern: ", miss[1L])
gs$pattern_type <- pats$pattern_type[match(key, paste(
  pats$method,
  pats$set
))]
gs <- gs[, c("method", "set", "pattern_type", "species", "gene")]
report(gs, "genesets")
cat("sets:", length(unique(key)), "\n")
tsv(gs, "genesets.tsv")

mus <- blocks(sheet("S6.Muscari_result.xlsx", "S6B.Gene list"), "module")
mus <- mus[, c("species", "gene", "module")]
report(mus, "muscari")
tsv(mus, "muscari.tsv")

# Four lists side by side: id column, then name; header row 5.
m <- sheet(
  "S8.Curated_gene_list_for_Medicago_nodulation.xlsx",
  "Curated gene list"
)
lists <- c(
  nodulation1 = 1L, nodulation2 = 4L, mycorrhizal = 7L,
  lost_nonmycorrhizal = 10L
)
cur <- do.call(rbind, lapply(names(lists), function(l) {
  c <- lists[[l]]
  id <- m[6:nrow(m), c]
  nm <- m[6:nrow(m), c + 1L]
  ok <- !is.na(id)
  data.frame(
    list = l, gene_v4 = id[ok], name = nm[ok],
    gene = sub("^Medtr", "MTR_", id[ok])
  )
}))
cat("== curated\n")
print(table(cur$list))
tsv(cur, "nodulation_curated.tsv")
