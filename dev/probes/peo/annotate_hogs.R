# Annotate the trait_specific hogs of a trait pass, with the
# differentiated and lineage_specific hogs of run 1 for context.
# Usage: Rscript annotate_hogs.R <exp> <run>
# Members come from RootHOGs.tsv of the 103-species FastOMA run: the
# Arabidopsis and rice genes of the hog's root hog (Arabidopsis shares
# the root hog even when it is not in the species set). UniProt is read
# in two bulk stream calls (organism 3702 and 39947), matched on the
# ordered locus name; the reviewed entry wins, else the first. OMA
# descriptions come from omabrowser.org/api/hog/<OMAmerRootHOG>/; the
# call is skipped on error. PLAZA is not used: no documented REST call
# was reachable (403).
# Writes score/<run>_annotated.tsv and score/<run>_annotated_go.tsv.
args <- commandArgs(TRUE)
EXP <- args[1]
RUN <- args[2]
ROOT <- Sys.getenv(
  "PEO_ROOT", "/mnt/project/FjellheimLab/martpali/rcomplex-testrun"
)
ROOTHOGS <- Sys.getenv("PEO_ROOTHOGS", paste0(
  "/mnt/project/FjellheimLab/martpali/AnnualPerennial/peo-fastoma/",
  "output/RootHOGs.tsv"
))
FOCAL <- c(C3vC4 = "C4", nodulation = "nodulating")
CAP <- 3L
E <- file.path(ROOT, paste0("complex_", EXP))
D <- file.path(E, "inputs")
S <- file.path(E, "score")
rd <- function(...) data.table::fread(file.path(...), data.table = FALSE)

info <- read.delim(file.path(D, "species.tsv"))
info <- info[info$run == RUN, ]
focal <- info$species[info$trait == FOCAL[[EXP]]]
sp <- info$species
ortho <- rd(D, "orthologs.tsv")
age <- rd(D, "hog_age.tsv")
hogs <- rd(S, paste0(RUN, "_hogs.tsv"))
tp <- rd(E, "runs", RUN, "trait_pass.tsv")

# hogs to annotate
ts <- unique(tp$hog[tp$hog_class == "trait_specific"])
ctx <- hogs$hog[hogs$tier %in% c("differentiated", "lineage_specific")]
sel <- data.frame(hog = c(ts, setdiff(ctx, ts)))
sel$class <- ifelse(sel$hog %in% ts, "trait_specific",
  hogs$tier[match(sel$hog, hogs$hog)]
)
cat("hogs:", nrow(sel), "\n")
print(table(sel$class))

# species of the best clique: the trait-pass clique for trait_specific
present <- function(m) {
  m[is.na(m)] <- ""
  lapply(strsplit(m, ",", fixed = TRUE), function(x) setdiff(sp, x))
}
tpb <- tp[tp$classification == "trait_specific", ]
tpb <- tpb[order(-tpb$n_species), ]
tpb <- tpb[!duplicated(tpb$hog), ]
best <- strsplit(hogs$species[match(sel$hog, hogs$hog)], ",", fixed = TRUE)
i <- match(sel$hog, tpb$hog)
best[!is.na(i)] <- present(tpb$missing_species[i[!is.na(i)]])
sel$n_focal <- vapply(best, function(b) sum(b %in% focal), 1L)
sel$n_other <- vapply(best, function(b) sum(b %in% setdiff(sp, focal)), 1L)
sel$tier <- hogs$tier[match(sel$hog, hogs$hog)]
sel$trait_class <- tp$hog_class[match(sel$hog, tp$hog)]
ab <- function(s) {
  paste0(substr(s, 1L, 1L), substr(sub("^[^_]*_", "", s), 1L, 2L))
}
cnt <- table(ortho$hog, factor(ortho$species, sp))
sel$members <- vapply(sel$hog, function(h) {
  n <- cnt[h, ]
  paste0(ab(names(n))[n > 0], ":", n[n > 0], collapse = ",")
}, "")

# root hog members from the 103-species run
sel$root_hog <- age$root_hog[match(sel$hog, age$hog)]
rh <- rd(ROOTHOGS)
rh <- rh[rh$RootHOG %in% age$root_hog, ]
strip <- function(x) sub("\\.[0-9]+$", "", x)
rh$id <- strip(rh$Protein)
sel$omamer <- rh$OMAmerRootHOG[match(sel$root_hog, rh$RootHOG)]
by_root <- function(pat) {
  k <- grepl(pat, rh$id)
  lapply(split(rh$id[k], rh$RootHOG[k]), function(v) head(unique(v), CAP))
}
mem_at <- by_root("^AT[1-5MC]G")
mem_os <- by_root("^LOC_Os")

# PEO annotations (peo.ku.dk): one gene call per member gene of the
# selected hogs, 0.2 s apart, cached as json.gz under score/peo_api/
# (each response also holds the gene's 50 PEO co-expression
# neighbours, kept for a later edge check). Bulk pages of 5,000 genes
# drew HTTP 403 after about ten pages, so they are not used. After a
# 403 the script makes no further calls and records NA.
API <- "https://peo.ku.dk/api"
CACHE <- file.path(S, "peo_api")
dir.create(file.path(CACHE, "genes"), showWarnings = FALSE,
  recursive = TRUE)
blocked <- FALSE
get_json <- function(path, name) {
  f <- file.path(CACHE, paste0(name, ".json.gz"))
  if (!file.exists(f)) {
    if (blocked) return(NULL)
    Sys.sleep(0.2)
    txt <- tryCatch(
      paste(readLines(url(paste0(API, path)), warn = FALSE),
        collapse = ""),
      error = function(e) {
        if (grepl("403", conditionMessage(e))) blocked <<- TRUE
        NA_character_
      }
    )
    if (is.na(txt)) return(NULL)
    con <- gzfile(f, "w")
    writeLines(txt, con)
    close(con)
  }
  jsonlite::fromJSON(f, simplifyVector = FALSE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
dict <- function(kind, size, pages) {
  rows <- unlist(lapply(seq_len(pages) - 1L, function(i) {
    get_json(sprintf("/%s?pageSize=%d&pageIndex=%d", kind, size, i),
      sprintf("%s_%03d", kind, i))$geneAnnotations
  }), recursive = FALSE)
  data.frame(
    id = vapply(rows, `[[`, "", "_id"),
    label = vapply(rows, `[[`, "", "label"),
    name = vapply(rows, function(r) r$name %||% NA_character_, ""),
    desc = vapply(rows, function(r) r$details$desc %||% NA_character_, "")
  )
}
mapman <- dict("mapman", 50L, 125L)
pfam <- dict("interpro", 5000L, 2L)
cat("dictionaries: mapman", nrow(mapman), "interpro", nrow(pfam), "\n")
map <- read.delim(file.path(
  "/mnt/users/martpali/AnnualPerennial/coexpr/PEO",
  "taxid_species_mapping.tsv"
))
om <- ortho[ortho$species %in% sp & ortho$hog %in% sel$hog, ]
om <- om[order(!om$hog %in% sel$hog[sel$class == "trait_specific"]), ]
tx <- map$taxid[match(om$species, map$species)]
ga <- lapply(seq_len(nrow(om)), function(k) {
  g <- om$gene[k]
  j <- get_json(
    sprintf("/species/%s/genes?pageSize=1&queryFilter=%s", tx[k],
      utils::URLencode(g, reserved = TRUE)),
    file.path("genes", paste0(tx[k], "_", gsub("[^A-Za-z0-9._-]", "_", g)))
  )
  if (is.null(j) || length(j$genes) == 0L) return(NA)
  h <- j$genes[[1L]]
  ids <- toupper(c(h$label, vapply(h$alias, `[[`, "", "label")))
  if (!toupper(g) %in% c(ids, sub("\\.[0-9]+$", "", ids))) return(NA)
  unlist(h$ga_ids) %||% character(0)
})
cat("PEO gene calls: genes", nrow(om), "found",
  sum(!vapply(ga, function(x) identical(x, NA), NA)), "blocked", blocked,
  "\n")
ga <- lapply(ga, function(x) if (identical(x, NA)) NULL else x)
by_hog <- split(ga, om$hog)
fmt <- function(lab, nm, n) {
  if (length(lab) == 0L) return("")
  o <- order(-n)
  paste0(lab[o], " ", nm[o], " (", n[o], ")", collapse = "; ")
}
peo_hog <- function(h) {
  g <- by_hog[[h]]
  ids <- lapply(g, unique)
  mm <- mapman[mapman$id %in% unlist(ids), ]
  pf <- pfam[pfam$id %in% unlist(ids), ]
  n_of <- function(d) vapply(d$id, function(x) {
    sum(vapply(ids, function(v) x %in% v, NA))
  }, 1L)
  c(
    peo_annotated = round(mean(lengths(g) > 0L), 3),
    peo_mapman = fmt(mm$label, mm$desc, n_of(mm)),
    peo_pfam = fmt(pf$label, pf$name, n_of(pf))
  )
}
pa <- t(vapply(sel$hog, peo_hog, character(3)))
sel <- cbind(sel, pa)

# UniProt, two bulk calls
uni <- function(org) {
  f <- tempfile(fileext = ".tsv.gz")
  u <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?query=organism_id:", org,
    "&fields=accession,reviewed,gene_oln,gene_names,protein_name,go_p",
    "&format=tsv&compressed=true"
  )
  utils::download.file(u, f, quiet = TRUE, mode = "wb")
  x <- data.table::fread(f, data.table = FALSE, quote = "")
  names(x) <- c("acc", "reviewed", "oln", "genes", "protein", "go_p")
  x <- x[order(x$reviewed != "reviewed"), ]
  tok <- strsplit(toupper(paste(x$oln, x$genes)), "[ ;/]+")
  k <- data.frame(id = unlist(tok), row = rep(
    seq_len(nrow(x)),
    lengths(tok)
  ))
  k <- k[!duplicated(k$id), ]
  cat("UniProt", org, nrow(x), "entries\n")
  list(x = x, row = setNames(k$row, k$id))
}
u_at <- uni(3702)
u_os <- uni(39947)
look <- function(u, ids) {
  r <- u$row[toupper(ids)]
  u$x[r[!is.na(r)], , drop = FALSE]
}
go_terms <- function(go) {
  t <- unlist(strsplit(go[!is.na(go) & nzchar(go)], "; ", fixed = TRUE))
  sub(" \\[GO:[0-9]+\\]$", "", t)
}
annot <- function(root) {
  at <- mem_at[[root]]
  os <- mem_os[[root]]
  a <- look(u_at, at)
  o <- look(u_os, os)
  go <- go_terms(c(a$go_p, o$go_p))
  top <- names(sort(table(go), decreasing = TRUE))
  list(
    at = at, os = os,
    at_symbol = paste(unique(sub(" .*", "", a$genes)), collapse = ","),
    protein = paste(unique(c(a$protein, o$protein)), collapse = " | "),
    go = unique(go), go_top = paste(head(top, 3L), collapse = "; ")
  )
}
roots <- unique(sel$root_hog)
an <- lapply(roots, annot)
names(an) <- roots
pick <- function(f) {
  vapply(sel$root_hog, function(r) {
    v <- an[[r]][[f]]
    if (length(v) == 0L) "" else paste(v, collapse = ",")
  }, "")
}
sel$at_ids <- pick("at")
sel$at_symbols <- pick("at_symbol")
sel$rice_loci <- pick("os")
sel$protein <- pick("protein")
sel$go_bp <- pick("go_top")

# OMA, one call per hog
oma <- function(id) {
  if (is.na(id) || !nzchar(id)) {
    return(NA_character_)
  }
  Sys.sleep(0.2)
  j <- tryCatch(
    paste(readLines(url(paste0(
      "https://omabrowser.org/api/hog/", id,
      "/?format=json"
    )), warn = FALSE), collapse = ""),
    error = function(e) ""
  )
  d <- regmatches(j, regexpr('"description":"[^"]*"', j))
  if (length(d) == 0L) {
    NA_character_
  } else {
    sub(
      '^"description":"', "",
      sub('"$', "", d)
    )
  }
}
sel$oma_description <- vapply(sel$omamer, oma, "")
sel$plaza_family <- NA
sel$plaza_description <- NA

sel <- sel[order(sel$class != "trait_specific", -sel$n_focal), ]
data.table::fwrite(sel, file.path(S, paste0(RUN, "_annotated.tsv")),
  sep = "\t", na = "NA", quote = FALSE
)
cat(
  "trait_specific:", sum(sel$class == "trait_specific"),
  "with AT member", sum(sel$class == "trait_specific" & nzchar(sel$at_ids)),
  "with rice member",
  sum(sel$class == "trait_specific" & nzchar(sel$rice_loci)),
  "with OMA description",
  sum(sel$class == "trait_specific" & !is.na(sel$oma_description)), "\n"
)

# GO BP: hogs per term, the trait_specific set against all hogs with a
# clique in run 1
cl <- data.table::fread(file.path(E, "runs", RUN, "cliques.tsv"),
  select = c("run", "hog"), data.table = FALSE
)
bg <- unique(age$root_hog[match(unique(cl$hog[cl$run == 1L]), age$hog)])
cl_hogs <- unique(cl$hog[cl$run == 1L])
rm(cl)
bg_go <- lapply(bg, function(r) {
  if (!is.null(an[[r]])) {
    return(an[[r]]$go)
  }
  unique(go_terms(c(
    look(u_at, mem_at[[r]])$go_p,
    look(u_os, mem_os[[r]])$go_p
  )))
})
fg_go <- lapply(
  unique(sel$root_hog[sel$class == "trait_specific"]),
  function(r) an[[r]]$go
)
tab <- function(l) table(unlist(l))
fg <- tab(fg_go)
bgt <- tab(bg_go)
go <- data.frame(
  term = names(fg), n_trait_specific = as.vector(fg),
  n_background = as.vector(bgt[names(fg)])
)
go$share_ts <- go$n_trait_specific / length(fg_go)
go$share_bg <- go$n_background / length(bg_go)
go$ratio <- round(go$share_ts / go$share_bg, 3)
go <- go[order(-go$ratio, -go$n_trait_specific), ]
cat(
  "root hogs: trait_specific", length(fg_go), "background", length(bg_go),
  "\n"
)
data.table::fwrite(go, file.path(S, paste0(RUN, "_annotated_go.tsv")),
  sep = "\t", quote = FALSE
)
print(head(go, 10L), row.names = FALSE)

# Mapman level-1 and level-2 bins: hogs per bin, the trait_specific set
# against the annotated hogs (trait_specific, differentiated,
# lineage_specific) with a clique in run 1. PEO refused bulk calls, so
# there is no annotation for all cliqued hogs.
bin_sets <- function(h, lv) {
  ids <- unique(unlist(by_hog[[h]]))
  nm <- mapman$name[mapman$id %in% ids]
  unique(vapply(strsplit(nm, ".", fixed = TRUE), function(x) {
    paste(head(x, lv), collapse = ".")
  }, ""))
}
bg_h <- intersect(cl_hogs, names(by_hog))
fg_h <- intersect(sel$hog[sel$class == "trait_specific"], names(by_hog))
mm_out <- do.call(rbind, lapply(1:2, function(lv) {
  fg <- table(unlist(lapply(fg_h, bin_sets, lv)))
  bg <- table(unlist(lapply(bg_h, bin_sets, lv)))
  if (length(bg) == 0L) return(NULL)
  d <- data.frame(level = lv, bin = names(bg), n_background = as.vector(bg),
    n_trait_specific = as.vector(fg[names(bg)])
  )
  d$n_trait_specific[is.na(d$n_trait_specific)] <- 0L
  d$ratio <- round((d$n_trait_specific / length(fg_h)) /
    (d$n_background / length(bg_h)), 3)
  d[order(-d$ratio, -d$n_trait_specific), ]
}))
if (is.null(mm_out)) stop("no PEO Mapman bins: every PEO call failed")
data.table::fwrite(mm_out, file.path(S, paste0(RUN, "_annotated_mapman.tsv")),
  sep = "\t", quote = FALSE
)
cat("hogs: trait_specific", length(fg_h), "background", length(bg_h), "\n")
print(head(mm_out[mm_out$level == 2L, ], 10L), row.names = FALSE)
