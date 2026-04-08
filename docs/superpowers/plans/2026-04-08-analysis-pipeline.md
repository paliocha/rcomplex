# Pooideae Clique Analysis Pipeline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Nextflow DSL2 pipeline that runs the full rcomplex clique analysis on 8 Pooideae grass species on Orion HPC at NMBU.

**Architecture:** 15 R wrapper scripts call rcomplex exported functions; 12 Nextflow DSL2 modules wire them into a DAG with SLURM scheduling; main.nf orchestrates the workflow. GPU jobs (network construction, co-expressolog detection, density/threshold sweeps) use RTX PRO 6000 on gn-41; CPU jobs (stability, perturbation, intensity, assembly) run on the orion partition.

**Tech Stack:** Nextflow DSL2, R 4.4.2, rcomplex package, SLURM, CUDA/torch

---

## Key Constraints

- Do NOT modify any R or C++ source code in the rcomplex package (R/, src/ directories). Only create files under `nextflow/`.
- The pipeline is parameterized by `params.tissue` ("leaf", "root", or "both"). Tissue filtering happens in `reduce_orthogroups.R` before paralog correlation. Three separate `nextflow run` invocations produce independent results in `results/leaf/`, `results/root/`, `results/both/`.
- Correlation method is **pearson** (not spearman).

Read these rcomplex source files for API context:
- `R/cliques.R` -- find_cliques, clique_stability, clique_persistence, clique_threshold_sweep, clique_perturbation_test, clique_intensity_test, classify_cliques
- `R/comparison.R` -- find_coexpressologs, density_sweep, compare_neighborhoods, comparison_to_edges
- `R/network.R` -- compute_network (use_torch, mr_log_transform)
- `R/orthologs.R` -- reduce_orthogroups, prepare_orthologs
- `R/se_methods.R` -- extract_orthologs
- `CLAUDE.md` -- package architecture

---

## Species

| Abbr | Species | Trait | Genes | Samples | HOG coverage |
|------|---------|-------|-------|---------|--------------|
| BDIS | *Brachypodium distachyon* | annual | 24,669 | 40 | 85.0% |
| BSYL | *Brachypodium sylvaticum* | perennial | 26,031 | 40 | 83.9% |
| HVUL | *Hordeum vulgare* | annual | 21,874 | 40 | 90.3% |
| HJUB | *Hordeum jubatum* | perennial | 53,106 | 40 | 82.8% |
| BMAX | *Briza maxima* | annual | 30,293 | 40 | 66.4% |
| BMED | *Briza media* | perennial | 39,937 | 40 | 61.1% |
| VBRO | *Vulpia bromoides* | annual | 28,970 | 38 | 77.0% |
| FPRA | *Festuca pratensis* | perennial | 33,804 | 40 | 65.9% |

4 phylogenetic pairs: Brachypodium, Hordeum, Briza, Festuca/Vulpia.
Each species: 2 tissues (leaf, root) x 5 time points (T1-T5) x 4 bio reps.

Source data: `nf-rcomplex/vst_hog.RDS` (14.6M-row long-format tibble, 11 species).
Three additional species not used: PANN, PSUP (Poa pair), BMED2BMAX (BMED re-annotated on BMAX genome).

**Input SE files** are pre-generated in `prepare_data/data/` via `Rscript prepare_data/convert_to_se.R`.
Each SE has assay `"vst"` (gene x sample matrix), `rowData$hog` (HOG:NNNNNNN), and `colData` with sample_id, tissue, time_point, day, real_day.

---

## Parameters

| Parameter | Value |
|-----------|-------|
| `tissue` | "leaf", "root", or "both" (three separate runs) |
| `density` | 0.03 |
| `min_var` | 1e-3 |
| `cor_method` | "pearson" |
| `norm_method` | "MR" |
| `mr_log_transform` | TRUE |
| `use_torch` | TRUE (CUDA on Orion) |
| `min_species` | 2 |
| `max_genes_per_sp` | 10 |
| `cost_weights` | default (q=1, effect=0) -- do NOT pass explicitly |
| `coexpressolog method` | "permutation" (Besag-Clifford + Liang q-values) |
| `paralog cor_threshold` | 0.7 |
| `density_sweep multipliers` | seq(0.7, 1.3, by = 0.05) -- 13 levels, eff. density ~10% to ~0.9% |
| `density_sweep method` | "permutation" |
| `threshold_sweep multipliers` | c(1.2, 1.5, 2.0, 2.5, 3.0) -- eff. density ~1.4% to ~0.001% |
| `threshold_sweep method` | "permutation" |
| `stability max_k` | 6 (= 8 species - 2) |
| `perturbation n_boot` | 1000 (4 batches x 250) |
| `perturbation noise_sd` | 0.01 |
| `intensity n_perm` | 1000 (5 batches x 200) |
| `seed` | 42 |

---

## Output: `results/{tissue}/cliques_ranked.tsv`

Every clique gets a continuous trait-composition score and multi-dimensional robustness ranking. No discrete waterfall classification.

| Column group | Columns |
|-------------|---------|
| Identity | `hog`, per-species gene columns, `n_species` |
| Trait composition | `n_annual`, `n_perennial`, `annual_fraction`, `annual_enrichment` |
| Edge quality | `mean_q`, `max_q`, `mean_effect_size`, `intensity`, `coherence`, `min_effect_size` |
| Stability | `stability_class` (highest k where ALL subsets are stable) |
| Threshold sweep | `birth`, `death`, `threshold_persistence` |
| Density sweep | `n_densities_survived` (out of 13 levels) |
| Weakest link | `weakest_link_persistence`, `mean_persistence` |
| Noise robustness | `survival_rate`, `mean_jaccard`, `n_matched` |
| Null model | `observed_intensity`, `z_score`, `p_value`, `q_value` |

---

## R Script Conventions

Every R script in `nextflow/scripts/` must:
1. `#!/usr/bin/env Rscript` shebang
2. Parse args via `commandArgs(trailingOnly = TRUE)`
3. `library(rcomplex)` (plus Bioc packages if needed)
4. Read inputs from working directory (Nextflow stages them)
5. Write `.tsv` (sep="\t", row.names=FALSE, quote=FALSE), `.rds` for lists/matrices
6. Guard empty outputs: header-only file so Nextflow sees the output
7. Scripts staged as `file()` inputs, NOT `bin/`

## Orion Hard Rules (non-negotiable)

1. `stageInMode = 'copy'` -- symlinks break on `/mnt/project/`
2. `beforeScript` boilerplate on every process (CPU or GPU variant)
3. `unset R_HOME` after module load
4. Torch symlinks in every GPU beforeScript
5. `LD_LIBRARY_PATH` via `torchLd()` in every GPU script block
6. `optional: true` on outputs that might be empty
7. Never write-then-read empty data frames
8. `cache = 'lenient'` globally
9. R scripts as `file()` inputs, not `bin/`
10. GPU: `--nodelist=gn-41 --gres=gpu:rtxpro6000:1`, partition `TestGPU`
11. CPU: partition `orion`, no `--account` flag
12. Never install R packages from login nodes

---

## File Structure

All files live under `nextflow/` in the repo root:

```
nextflow/
├── main.nf                          # Workflow orchestration + boilerplate
├── nextflow.config                  # SLURM profiles, resource allocations
├── samplesheet.csv                  # Species, SE paths, traits
├── modules/
│   ├── reduce.nf                    # REDUCE process
│   ├── network.nf                   # NETWORK process
│   ├── orthologs.nf                 # ORTHOLOGS process
│   ├── coexpressologs.nf            # COEXPRESSOLOGS process
│   ├── density_sweep.nf             # DENSITY_SWEEP process
│   ├── cliques.nf                   # FIND_CLIQUES process
│   ├── stability.nf                 # STABILITY process
│   ├── threshold_edges.nf           # THRESHOLD_EDGES process (per multiplier)
│   ├── assemble_threshold.nf        # ASSEMBLE_THRESHOLD process
│   ├── persistence.nf               # PERSISTENCE process
│   ├── perturbation.nf              # PERTURBATION_BATCH + MERGE_PERTURBATION
│   ├── intensity.nf                 # INTENSITY_BATCH + MERGE_INTENSITY
│   └── assemble.nf                  # ASSEMBLE_RANKED process
├── scripts/
│   ├── reduce_orthogroups.R
│   ├── compute_network.R
│   ├── prepare_orthologs.R
│   ├── find_coexpressologs.R
│   ├── density_sweep.R
│   ├── find_cliques.R
│   ├── clique_stability.R
│   ├── threshold_edges.R
│   ├── assemble_threshold.R
│   ├── clique_persistence.R
│   ├── perturbation_batch.R
│   ├── merge_perturbation.R
│   ├── intensity_batch.R
│   ├── merge_intensity.R
│   └── assemble_ranked.R
```

---

### Task 1: Scaffold, samplesheet, and nextflow.config

**Files:**
- Create: `nextflow/samplesheet.csv`
- Create: `nextflow/nextflow.config`

- [ ] **Step 1: Create the `nextflow/` directory tree**

```bash
mkdir -p nextflow/modules nextflow/scripts
```

- [ ] **Step 2: Write `nextflow/samplesheet.csv`**

SE files are pre-generated in `prepare_data/data/` (see `prepare_data/convert_to_se.R`).
Each SE has assay `"vst"` (gene x sample), `rowData$hog` (HOG:NNNNNNN), and
`colData` with sample_id, tissue, time_point, day, real_day.

```csv
species,se_path,trait
BDIS,prepare_data/data/BDIS_se.rds,annual
BSYL,prepare_data/data/BSYL_se.rds,perennial
HVUL,prepare_data/data/HVUL_se.rds,annual
HJUB,prepare_data/data/HJUB_se.rds,perennial
BMAX,prepare_data/data/BMAX_se.rds,annual
BMED,prepare_data/data/BMED_se.rds,perennial
VBRO,prepare_data/data/VBRO_se.rds,annual
FPRA,prepare_data/data/FPRA_se.rds,perennial
```

- [ ] **Step 3: Write `nextflow/nextflow.config`**

```groovy
cache = 'lenient'

profiles {
  orion {
    process {
      executor     = 'slurm'
      stageInMode  = 'copy'
      beforeScript = cpuBeforeScript  // defined in main.nf

      // CPU defaults
      queue   = 'orion'
      cpus    = 8
      memory  = '64 GB'
      time    = '24 h'
    }
    withLabel: 'gpu' {
      queue          = 'TestGPU'
      clusterOptions = '--nodelist=gn-41 --gres=gpu:rtxpro6000:1'
      cpus           = 16
      memory         = '128 GB'
      time           = '8 h'
      beforeScript   = gpuBeforeScript
    }
    withName: 'REDUCE' {
      cpus   = 1
      memory = '8 GB'
      time   = '10 min'
    }
    withName: 'NETWORK' {
      label  = 'gpu'
      cpus   = 1
      memory = '16 GB'
      time   = '30 min'
    }
    withName: 'ORTHOLOGS' {
      cpus   = 1
      memory = '4 GB'
      time   = '10 min'
    }
    withName: 'COEXPRESSOLOGS' {
      label  = 'gpu'
      cpus   = 16
      memory = '32 GB'
      time   = '1 h'
    }
    withName: 'DENSITY_SWEEP' {
      label  = 'gpu'
      cpus   = 16
      memory = '32 GB'
      time   = '4 h'
    }
    withName: 'FIND_CLIQUES' {
      cpus   = 1
      memory = '4 GB'
      time   = '10 min'
    }
    withName: 'STABILITY' {
      cpus   = 16
      memory = '32 GB'
      time   = '2 h'
    }
    withName: 'THRESHOLD_EDGES' {
      label  = 'gpu'
      cpus   = 16
      memory = '32 GB'
      time   = '1 h'
    }
    withName: 'ASSEMBLE_THRESHOLD' {
      cpus   = 1
      memory = '4 GB'
      time   = '10 min'
    }
    withName: 'PERSISTENCE' {
      cpus   = 1
      memory = '4 GB'
      time   = '10 min'
    }
    withName: 'PERTURBATION_BATCH' {
      cpus   = 8
      memory = '16 GB'
      time   = '2 h'
    }
    withName: 'INTENSITY_BATCH' {
      cpus   = 8
      memory = '16 GB'
      time   = '2 h'
    }
    withName: 'MERGE_PERTURBATION' {
      cpus   = 1
      memory = '4 GB'
      time   = '10 min'
    }
    withName: 'MERGE_INTENSITY' {
      cpus   = 1
      memory = '4 GB'
      time   = '10 min'
    }
    withName: 'ASSEMBLE_RANKED' {
      cpus   = 1
      memory = '8 GB'
      time   = '10 min'
    }
  }
}
```

- [ ] **Step 4: Verify**

```bash
cd nextflow && ls modules/ scripts/ samplesheet.csv nextflow.config
```

- [ ] **Step 5: Commit**

```bash
git add nextflow/samplesheet.csv nextflow/nextflow.config
git commit -m "feat(pipeline): scaffold + samplesheet + nextflow.config"
```

---

### Task 2: R scripts — data preparation (reduce, network, orthologs)

**Files:**
- Create: `nextflow/scripts/reduce_orthogroups.R`
- Create: `nextflow/scripts/compute_network.R`
- Create: `nextflow/scripts/prepare_orthologs.R`

Read `R/orthologs.R` (reduce_orthogroups, prepare_orthologs) and `R/network.R` (compute_network) for API context. Follow the R Script Conventions listed above.

- [ ] **Step 1: Write `nextflow/scripts/reduce_orthogroups.R`**

```r
#!/usr/bin/env Rscript
# Paralog reduction for one species with tissue filtering.
# Args: <species_abbrev> <se_path> <cor_threshold> <tissue>
# Input: SE .rds file (contains both leaf and root samples)
# Output: {species}_reduction.rds
#
# The tissue argument ("leaf", "root", or "both") filters samples
# BEFORE computing paralog correlations. This means different tissues
# may produce different paralog merging decisions.

args <- commandArgs(trailingOnly = TRUE)
species   <- args[1]
se_path   <- args[2]
cor_thr   <- as.numeric(args[3])
tissue    <- args[4]             # "leaf", "root", or "both"

library(rcomplex)

se <- readRDS(se_path)

# Tissue filter: subset samples before extracting expression matrix
if (tissue != "both") {
  keep <- SummarizedExperiment::colData(se)$tissue == tissue
  se <- se[, keep]
  message(species, ": filtered to ", tissue, " (", sum(keep), " samples)")
}

expr <- as.matrix(SummarizedExperiment::assay(se))

# Build species-specific ortholog table from rowData
rd <- as.data.frame(SummarizedExperiment::rowData(se))
ortho_sp <- data.frame(
  Species1 = rownames(rd),
  hog      = rd$hog,
  stringsAsFactors = FALSE
)
ortho_sp <- ortho_sp[!is.na(ortho_sp$hog), ]

reduction <- reduce_orthogroups(expr, ortho_sp,
                                 gene_col = "Species1",
                                 cor_threshold = cor_thr)

saveRDS(reduction, paste0(species, "_reduction.rds"))
message("Reduced ", species, ": ", reduction$n_original, " -> ",
        reduction$n_reduced, " genes (", reduction$n_merged, " merged)")
```

- [ ] **Step 2: Write `nextflow/scripts/compute_network.R`**

```r
#!/usr/bin/env Rscript
# Network construction for one species (GPU-accelerated).
# Args: <species_abbrev> <reduction_rds> <cor_method> <norm_method>
#       <density> <min_var> <mr_log_transform> <use_torch>
# Output: {species}_network.rds

args <- commandArgs(trailingOnly = TRUE)
species          <- args[1]
reduction_path   <- args[2]
cor_method       <- args[3]
norm_method      <- args[4]
density          <- as.numeric(args[5])
min_var          <- as.numeric(args[6])
mr_log_transform <- as.logical(args[7])
use_torch        <- as.logical(args[8])

library(rcomplex)

reduction <- readRDS(reduction_path)

net <- compute_network(reduction$expr_matrix,
                       cor_method       = cor_method,
                       norm_method      = norm_method,
                       density          = density,
                       min_var          = min_var,
                       mr_log_transform = mr_log_transform,
                       use_torch        = use_torch)

saveRDS(net, paste0(species, "_network.rds"))
message("Network ", species, ": ", net$n_genes, " genes, threshold=",
        round(net$threshold, 4), ", removed=", net$n_removed)
```

- [ ] **Step 3: Write `nextflow/scripts/prepare_orthologs.R`**

```r
#!/usr/bin/env Rscript
# Ortholog preparation: extract pairwise orthologs from SEs, map through reductions.
# Args: <species_csv> <se_dir_or_paths> <reduction_dir_or_paths> <hog_col>
# Stdin: space-separated species names, SE paths, reduction paths
# Output: orthologs.tsv

# This script expects staged files:
#   {species}_se.rds       — per species
#   {species}_reduction.rds — per species
# Plus a species list file: species_list.txt (one species per line)

args <- commandArgs(trailingOnly = TRUE)
hog_col <- if (length(args) >= 1) args[1] else "hog"

library(rcomplex)

# Read species list
species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

# Load SEs and reductions
se_list <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_se.rds"))),
  species_names
)
reductions <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_reduction.rds"))),
  species_names
)

orthologs <- prepare_orthologs(se_list, reductions, hog_col = hog_col)

if (nrow(orthologs) > 0) {
  write.table(orthologs, "orthologs.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  # Write header-only file so Nextflow sees output
  writeLines("Species1\tSpecies2\thog", "orthologs.tsv")
}
message("Orthologs: ", nrow(orthologs), " rows across ",
        length(unique(orthologs$hog)), " HOGs")
```

- [ ] **Step 4: Verify R syntax**

```bash
for f in nextflow/scripts/reduce_orthogroups.R nextflow/scripts/compute_network.R nextflow/scripts/prepare_orthologs.R; do
  Rscript -e "parse('$f')" && echo "OK: $f" || echo "FAIL: $f"
done
```

Expected: `OK` for all three.

- [ ] **Step 5: Commit**

```bash
git add nextflow/scripts/reduce_orthogroups.R nextflow/scripts/compute_network.R nextflow/scripts/prepare_orthologs.R
git commit -m "feat(pipeline): R scripts for data preparation (reduce, network, orthologs)"
```

---

### Task 3: R scripts — comparison (coexpressologs, density sweep)

**Files:**
- Create: `nextflow/scripts/find_coexpressologs.R`
- Create: `nextflow/scripts/density_sweep.R`

Read `R/comparison.R` (find_coexpressologs, density_sweep) for API signatures.

- [ ] **Step 1: Write `nextflow/scripts/find_coexpressologs.R`**

```r
#!/usr/bin/env Rscript
# Co-expressolog detection across all species pairs (GPU-accelerated).
# Args: <method> <alternative> <alpha> <n_cores> <use_torch>
# Staged inputs: {species}_network.rds (x8), orthologs.tsv, species_list.txt
# Output: edges.tsv

args <- commandArgs(trailingOnly = TRUE)
method      <- args[1]          # "permutation" or "analytical"
alternative <- args[2]          # "greater"
alpha       <- as.numeric(args[3])
n_cores     <- as.integer(args[4])
use_torch   <- as.logical(args[5])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

networks <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_network.rds"))),
  species_names
)
orthologs <- read.table("orthologs.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

edges <- find_coexpressologs(networks, orthologs,
                              method      = method,
                              alternative = alternative,
                              alpha       = alpha,
                              n_cores     = n_cores,
                              use_torch   = use_torch)

if (nrow(edges) > 0) {
  write.table(edges, "edges.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  writeLines(paste(names(edges), collapse = "\t"), "edges.tsv")
}
message("Edges: ", nrow(edges), " (", sum(edges$type == "conserved"), " conserved)")
```

- [ ] **Step 2: Write `nextflow/scripts/density_sweep.R`**

```r
#!/usr/bin/env Rscript
# Density sweep: re-threshold networks at multiple multipliers.
# Args: <multipliers_csv> <method> <n_cores> <use_torch>
# Staged inputs: {species}_network.rds (x8), orthologs.tsv, species_list.txt
# Output: density_sweep.rds

args <- commandArgs(trailingOnly = TRUE)
multipliers <- as.numeric(strsplit(args[1], ",")[[1]])
method      <- args[2]
n_cores     <- as.integer(args[3])
use_torch   <- as.logical(args[4])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

networks <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_network.rds"))),
  species_names
)
orthologs <- read.table("orthologs.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

sweep <- density_sweep(networks, orthologs,
                        multipliers = multipliers,
                        method      = method,
                        use_torch   = use_torch,
                        n_cores     = n_cores)

saveRDS(sweep, "density_sweep.rds")
message("Density sweep: ", nrow(sweep), " multipliers, ",
        "range ", min(sweep$n_significant), "-", max(sweep$n_significant),
        " significant edges")
```

- [ ] **Step 3: Verify R syntax**

```bash
for f in nextflow/scripts/find_coexpressologs.R nextflow/scripts/density_sweep.R; do
  Rscript -e "parse('$f')" && echo "OK: $f" || echo "FAIL: $f"
done
```

- [ ] **Step 4: Commit**

```bash
git add nextflow/scripts/find_coexpressologs.R nextflow/scripts/density_sweep.R
git commit -m "feat(pipeline): R scripts for comparison (coexpressologs, density sweep)"
```

---

### Task 4: R scripts — clique core (find, stability, persistence)

**Files:**
- Create: `nextflow/scripts/find_cliques.R`
- Create: `nextflow/scripts/clique_stability.R`
- Create: `nextflow/scripts/clique_persistence.R`

Read `R/cliques.R` (find_cliques, clique_stability, clique_persistence) for API signatures.

**Critical:** `min_species = 2L` must be passed explicitly to both `find_cliques` and `clique_stability`. The default is `length(target_species)` = 8, which would miss all smaller cliques. Similarly, `max_genes_per_sp = 10L` must be consistent across all clique functions.

- [ ] **Step 1: Write `nextflow/scripts/find_cliques.R`**

```r
#!/usr/bin/env Rscript
# Find cliques from co-expressolog edges.
# Args: <min_species> <max_genes_per_sp>
# Staged inputs: edges.tsv, species_list.txt
# Output: cliques.tsv

args <- commandArgs(trailingOnly = TRUE)
min_species     <- as.integer(args[1])
max_genes_per_sp <- as.integer(args[2])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

edges <- read.table("edges.tsv", header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)

cliques <- find_cliques(edges, species_names,
                         min_species      = min_species,
                         max_genes_per_sp = max_genes_per_sp)

if (nrow(cliques) > 0) {
  write.table(cliques, "cliques.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  writeLines(paste(names(cliques), collapse = "\t"), "cliques.tsv")
}
message("Cliques: ", nrow(cliques), " across ",
        length(unique(cliques$hog)), " HOGs")
```

- [ ] **Step 2: Write `nextflow/scripts/clique_stability.R`**

```r
#!/usr/bin/env Rscript
# Leave-k-out jackknife stability for cliques.
# Args: <min_species> <max_genes_per_sp> <max_k> <n_cores>
# Staged inputs: edges.tsv, cliques.tsv, species_list.txt, traits.tsv
# Output: stability.rds

args <- commandArgs(trailingOnly = TRUE)
min_species      <- as.integer(args[1])
max_genes_per_sp <- as.integer(args[2])
max_k            <- as.integer(args[3])
n_cores          <- as.integer(args[4])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

edges   <- read.table("edges.tsv", header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
cliques <- read.table("cliques.tsv", header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)

# Load trait vector
traits_df <- read.table("traits.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
trait_vec <- setNames(traits_df$trait, traits_df$species)

stab <- clique_stability(edges, species_names,
                          species_trait    = trait_vec,
                          all_species      = species_names,
                          full_cliques     = cliques,
                          min_species      = min_species,
                          max_genes_per_sp = max_genes_per_sp,
                          max_k            = max_k,
                          n_cores          = n_cores)

saveRDS(stab, "stability.rds")
message("Stability: ", length(stab$stability_class), " cliques classified, ",
        "max class=", max(stab$stability_class, na.rm = TRUE))
```

- [ ] **Step 3: Write `nextflow/scripts/clique_persistence.R`**

```r
#!/usr/bin/env Rscript
# Weakest-link co-expressolog persistence scores.
# Args: (none)
# Staged inputs: cliques.tsv, edges.tsv, {species}_network.rds, species_list.txt
# Output: clique_persistence.tsv

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

cliques  <- read.table("cliques.tsv", header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
edges    <- read.table("edges.tsv", header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
networks <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_network.rds"))),
  species_names
)

# clique_persistence returns cliques df with persistence + mean_persistence appended
persist <- clique_persistence(cliques, species_names, networks, edges)

write.table(persist, "clique_persistence.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Persistence: ", sum(!is.na(persist$persistence)), "/",
        nrow(persist), " cliques scored")
```

- [ ] **Step 4: Verify R syntax**

```bash
for f in nextflow/scripts/find_cliques.R nextflow/scripts/clique_stability.R nextflow/scripts/clique_persistence.R; do
  Rscript -e "parse('$f')" && echo "OK: $f" || echo "FAIL: $f"
done
```

- [ ] **Step 5: Commit**

```bash
git add nextflow/scripts/find_cliques.R nextflow/scripts/clique_stability.R nextflow/scripts/clique_persistence.R
git commit -m "feat(pipeline): R scripts for clique core (find, stability, persistence)"
```

---

### Task 5: R scripts — threshold sweep (per-multiplier edges + assembly)

**Files:**
- Create: `nextflow/scripts/threshold_edges.R`
- Create: `nextflow/scripts/assemble_threshold.R`

The threshold sweep is restructured as parallel per-multiplier GPU jobs. Each multiplier re-thresholds all networks (threshold * multiplier) and runs `find_coexpressologs(method = "permutation")` + `find_cliques`. Assembly computes birth/death persistence by Jaccard-matching to baseline cliques. `jaccard_clique_match` is not exported from rcomplex, so it's inlined.

- [ ] **Step 1: Write `nextflow/scripts/threshold_edges.R`**

```r
#!/usr/bin/env Rscript
# Threshold sweep: re-threshold networks at one multiplier and run
# permutation-based co-expressolog detection + clique finding.
# Args: <multiplier> <min_species> <max_genes_per_sp> <n_cores> <use_torch>
# Staged inputs: {species}_network.rds (x8), orthologs.tsv, species_list.txt
# Output: threshold_{m}_edges.tsv, threshold_{m}_cliques.tsv

args <- commandArgs(trailingOnly = TRUE)
multiplier       <- as.numeric(args[1])
min_species      <- as.integer(args[2])
max_genes_per_sp <- as.integer(args[3])
n_cores          <- as.integer(args[4])
use_torch        <- as.logical(args[5])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

networks <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_network.rds"))),
  species_names
)
orthologs <- read.table("orthologs.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

# Re-threshold networks (same MR matrix, stricter threshold)
thr_networks <- lapply(networks, function(net) {
  list(network = net$network, threshold = net$threshold * multiplier)
})
names(thr_networks) <- names(networks)

# Permutation-based co-expressolog detection at this threshold
thr_edges <- tryCatch(
  find_coexpressologs(thr_networks, orthologs,
                       method      = "permutation",
                       alternative = "greater",
                       alpha       = 0.05,
                       n_cores     = n_cores,
                       use_torch   = use_torch),
  error = function(e) {
    warning("Threshold ", multiplier, "x: co-expressolog detection failed: ",
            conditionMessage(e))
    data.frame(gene1 = character(0), gene2 = character(0),
               species1 = character(0), species2 = character(0),
               hog = character(0), q.value = numeric(0),
               effect_size = numeric(0), jaccard = numeric(0),
               type = character(0))
  }
)

m_tag <- as.character(multiplier)

if (nrow(thr_edges) > 0) {
  # Find cliques at this threshold
  thr_cliques <- find_cliques(thr_edges, species_names,
                                min_species      = min_species,
                                max_genes_per_sp = max_genes_per_sp)

  write.table(thr_edges,
    paste0("threshold_", m_tag, "_edges.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  if (nrow(thr_cliques) > 0) {
    write.table(thr_cliques,
      paste0("threshold_", m_tag, "_cliques.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    writeLines(paste(names(thr_cliques), collapse = "\t"),
               paste0("threshold_", m_tag, "_cliques.tsv"))
  }
  message("Threshold ", m_tag, "x: ", nrow(thr_edges), " edges, ",
          nrow(thr_cliques), " cliques")
} else {
  writeLines(paste(names(thr_edges), collapse = "\t"),
             paste0("threshold_", m_tag, "_edges.tsv"))
  writeLines("hog", paste0("threshold_", m_tag, "_cliques.tsv"))
  message("Threshold ", m_tag, "x: 0 edges")
}
```

- [ ] **Step 2: Write `nextflow/scripts/assemble_threshold.R`**

This script includes the inline `jaccard_match` function (since `rcomplex:::jaccard_clique_match` is not exported), survival matrix construction, and birth/death/persistence computation.

```r
#!/usr/bin/env Rscript
# Assemble threshold sweep results: Jaccard-match baseline cliques to
# per-multiplier cliques and compute birth/death/persistence.
# Args: <multipliers_csv>
# Staged inputs: cliques.tsv, threshold_{m}_cliques.tsv (one per multiplier),
#                species_list.txt
# Output: threshold_persistence.tsv

args <- commandArgs(trailingOnly = TRUE)
multipliers <- as.numeric(strsplit(args[1], ",")[[1]])

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

baseline <- read.table("cliques.tsv", header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)

# Inline Jaccard match (rcomplex:::jaccard_clique_match is not exported)
jaccard_match <- function(row1, row2, species) {
  intersect_n <- 0L; union_n <- 0L
  for (sp in species) {
    g1 <- row1[[sp]]; g2 <- row2[[sp]]
    has1 <- !is.na(g1); has2 <- !is.na(g2)
    if (has1 || has2) {
      union_n <- union_n + 1L
      if (has1 && has2 && g1 == g2) intersect_n <- intersect_n + 1L
    }
  }
  if (union_n == 0L) 0 else intersect_n / union_n
}

n_cliques <- nrow(baseline)
survived <- matrix(FALSE, nrow = n_cliques, ncol = length(multipliers),
                    dimnames = list(NULL, as.character(multipliers)))

for (mi in seq_along(multipliers)) {
  m <- multipliers[mi]
  f <- paste0("threshold_", as.character(m), "_cliques.tsv")
  if (!file.exists(f)) next
  thr_cl <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(thr_cl) == 0L) next

  for (i in seq_len(n_cliques)) {
    candidates <- which(thr_cl$hog == baseline$hog[i])
    for (j in candidates) {
      if (jaccard_match(baseline[i, ], thr_cl[j, ], species_names) > 0) {
        survived[i, mi] <- TRUE
        break
      }
    }
  }
}

# Compute birth/death/persistence
all_mult <- c(1.0, multipliers)
persist_list <- vector("list", n_cliques)
for (i in seq_len(n_cliques)) {
  surv_at <- c(TRUE, survived[i, ])
  names(surv_at) <- as.character(all_mult)
  survived_mult <- all_mult[surv_at]
  birth <- if (length(survived_mult) > 0) min(survived_mult) else NA_real_

  death <- NA_real_
  if (!is.na(birth)) {
    above_birth <- all_mult[all_mult > birth]
    for (cand in above_birth) {
      if (!surv_at[as.character(cand)]) { death <- cand; break }
    }
  }

  persist_list[[i]] <- data.frame(
    clique_idx = i, hog = baseline$hog[i],
    birth = birth, death = death,
    threshold_persistence = if (!is.na(death)) death - birth else NA_real_,
    stringsAsFactors = FALSE)
}

threshold_persistence <- do.call(rbind, persist_list)
write.table(threshold_persistence, "threshold_persistence.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Threshold persistence: ", sum(!is.na(threshold_persistence$threshold_persistence)),
        "/", nrow(threshold_persistence), " cliques with finite persistence")
```

- [ ] **Step 3: Verify R syntax**

```bash
for f in nextflow/scripts/threshold_edges.R nextflow/scripts/assemble_threshold.R; do
  Rscript -e "parse('$f')" && echo "OK: $f" || echo "FAIL: $f"
done
```

- [ ] **Step 4: Commit**

```bash
git add nextflow/scripts/threshold_edges.R nextflow/scripts/assemble_threshold.R
git commit -m "feat(pipeline): R scripts for threshold sweep (per-multiplier + assembly)"
```

---

### Task 6: R scripts — perturbation and intensity (batch + merge)

**Files:**
- Create: `nextflow/scripts/perturbation_batch.R`
- Create: `nextflow/scripts/merge_perturbation.R`
- Create: `nextflow/scripts/intensity_batch.R`
- Create: `nextflow/scripts/merge_intensity.R`

Perturbation: `clique_perturbation_test()` returns `survival_rate` (= n_survived/n_boot) and `mean_jaccard` (= sum_jaccard/n_matched). Merge back-computes numerators to aggregate across batches.

Intensity: `clique_intensity_test()` returns summary stats (`null_mean`, `null_sd`, `n_matched`) not raw null distributions. Merge uses law of total variance: `Var(X) = E[Var(X|batch)] + Var(E[X|batch])` to pool batch summaries, then normal approximation for p-values.

- [ ] **Step 1: Write `nextflow/scripts/perturbation_batch.R`**

```r
#!/usr/bin/env Rscript
# Perturbation test batch: noise-injection bootstrap.
# Args: <n_boot> <noise_sd> <min_species> <max_genes_per_sp> <n_cores> <seed>
# Staged inputs: cliques.tsv, orthologs.tsv, {species}_network.rds, species_list.txt
# Output: perturbation_batch.tsv

args <- commandArgs(trailingOnly = TRUE)
n_boot           <- as.integer(args[1])
noise_sd         <- as.numeric(args[2])
min_species      <- as.integer(args[3])
max_genes_per_sp <- as.integer(args[4])
n_cores          <- as.integer(args[5])
seed             <- as.integer(args[6])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

cliques   <- read.table("cliques.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
orthologs <- read.table("orthologs.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
networks  <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_network.rds"))),
  species_names
)

result <- clique_perturbation_test(
  cliques, species_names, networks, orthologs,
  n_boot           = n_boot,
  noise_sd         = noise_sd,
  min_species      = min_species,
  max_genes_per_sp = max_genes_per_sp,
  jaccard_threshold = 0.5,
  n_cores          = n_cores,
  seed             = seed)

write.table(result, "perturbation_batch.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Perturbation batch (seed=", seed, "): ",
        sum(result$survival_rate > 0), "/", nrow(result), " survived")
```

- [ ] **Step 2: Write `nextflow/scripts/merge_perturbation.R`**

```r
#!/usr/bin/env Rscript
# Merge perturbation test batches.
# Args: (none, reads all perturbation_batch_*.tsv from working dir)
# Output: perturbation_merged.tsv

batch_files <- sort(Sys.glob("perturbation_batch_*.tsv"))
if (length(batch_files) == 0) stop("No perturbation batch files found")

batches <- lapply(batch_files, function(f)
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE))

n_cliques <- nrow(batches[[1]])
merged <- batches[[1]][, c("clique_idx", "hog")]
merged$survival_rate <- NA_real_
merged$mean_jaccard  <- NA_real_
merged$n_boot        <- 0L
merged$n_matched     <- 0L

for (i in seq_len(n_cliques)) {
  total_boot     <- 0L
  total_survived <- 0L
  total_matched  <- 0L
  sum_jaccard    <- 0

  for (b in batches) {
    nb <- b$n_boot[i]
    nm <- b$n_matched[i]
    total_boot     <- total_boot + nb
    total_survived <- total_survived + round(b$survival_rate[i] * nb)
    total_matched  <- total_matched + nm
    if (nm > 0 && !is.na(b$mean_jaccard[i]))
      sum_jaccard <- sum_jaccard + b$mean_jaccard[i] * nm
  }

  merged$n_boot[i]        <- total_boot
  merged$n_matched[i]     <- total_matched
  merged$survival_rate[i] <- total_survived / total_boot
  merged$mean_jaccard[i]  <- if (total_matched > 0) sum_jaccard / total_matched else NA_real_
}

write.table(merged, "perturbation_merged.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Perturbation merged: ", length(batch_files), " batches, ",
        merged$n_boot[1], " total bootstraps")
```

- [ ] **Step 3: Write `nextflow/scripts/intensity_batch.R`**

```r
#!/usr/bin/env Rscript
# Intensity test batch: ortholog-shuffle null model.
# Args: <n_perm> <min_species> <max_genes_per_sp> <n_cores> <seed>
# Staged inputs: cliques.tsv, edges.tsv, orthologs.tsv,
#                {species}_network.rds, species_list.txt
# Output: intensity_batch.tsv

args <- commandArgs(trailingOnly = TRUE)
n_perm           <- as.integer(args[1])
min_species      <- as.integer(args[2])
max_genes_per_sp <- as.integer(args[3])
n_cores          <- as.integer(args[4])
seed             <- as.integer(args[5])

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

cliques   <- read.table("cliques.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
edges     <- read.table("edges.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
orthologs <- read.table("orthologs.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
networks  <- setNames(
  lapply(species_names, function(sp) readRDS(paste0(sp, "_network.rds"))),
  species_names
)

result <- clique_intensity_test(
  cliques, species_names, networks, orthologs,
  n_perm           = n_perm,
  alternative      = "greater",
  edges            = edges,
  min_species      = min_species,
  max_genes_per_sp = max_genes_per_sp,
  n_cores          = n_cores,
  seed             = seed)

write.table(result, "intensity_batch.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Intensity batch (seed=", seed, "): ",
        sum(!is.na(result$z_score)), "/", nrow(result), " scored")
```

- [ ] **Step 4: Write `nextflow/scripts/merge_intensity.R`**

```r
#!/usr/bin/env Rscript
# Merge intensity test batches using law of total variance.
# Args: (none, reads all intensity_batch_*.tsv from working dir)
# Output: intensity_merged.tsv

batch_files <- sort(Sys.glob("intensity_batch_*.tsv"))
if (length(batch_files) == 0) stop("No intensity batch files found")

batches <- lapply(batch_files, function(f)
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE))

n_cliques <- nrow(batches[[1]])
merged <- batches[[1]][, c("clique_idx", "hog", "observed_intensity")]

pooled_mean <- pooled_sd <- z_score <- p_value <- rep(NA_real_, n_cliques)
total_perm <- total_matched <- rep(0L, n_cliques)

for (i in seq_len(n_cliques)) {
  bm <- sapply(batches, function(b) b$null_mean[i])
  bs <- sapply(batches, function(b) b$null_sd[i])
  bn <- sapply(batches, function(b) b$n_matched[i])
  bp <- sapply(batches, function(b) b$n_perm[i])

  valid <- !is.na(bm) & bn > 0
  total_perm[i] <- sum(bp)
  total_matched[i] <- sum(bn[valid])
  if (!any(valid)) next

  w <- bn[valid]
  N <- sum(w)

  # Pooled mean (weighted by n_matched)
  pm <- sum(bm[valid] * w) / N
  # Law of total variance: Var(X) = E[Var(X|batch)] + Var(E[X|batch])
  within_v  <- sum(bs[valid]^2 * w) / N
  between_v <- sum(w * (bm[valid] - pm)^2) / N
  ps <- sqrt(within_v + between_v)

  pooled_mean[i] <- pm
  pooled_sd[i]   <- ps

  obs <- merged$observed_intensity[i]
  z_score[i] <- if (!is.na(ps) && ps > 0) (obs - pm) / ps else NA_real_
  p_value[i] <- if (!is.na(z_score[i])) pnorm(z_score[i], lower.tail = FALSE) else NA_real_
}

merged$null_mean  <- pooled_mean
merged$null_sd    <- pooled_sd
merged$z_score    <- z_score
merged$p_value    <- p_value
merged$n_perm     <- total_perm
merged$n_matched  <- total_matched

# Storey q-values over cliques with sufficient matched permutations
reliable <- merged$n_matched >= 50 & !is.na(merged$p_value)
merged$q_value <- NA_real_
if (sum(reliable) >= 2) {
  merged$q_value[reliable] <- qvalue::qvalue(merged$p_value[reliable])$qvalues
}

write.table(merged, "intensity_merged.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Intensity merged: ", length(batch_files), " batches, ",
        total_perm[1], " total perms, ",
        sum(reliable), " cliques with q-values")
```

- [ ] **Step 5: Verify R syntax**

```bash
for f in nextflow/scripts/perturbation_batch.R nextflow/scripts/merge_perturbation.R nextflow/scripts/intensity_batch.R nextflow/scripts/merge_intensity.R; do
  Rscript -e "parse('$f')" && echo "OK: $f" || echo "FAIL: $f"
done
```

- [ ] **Step 6: Commit**

```bash
git add nextflow/scripts/perturbation_batch.R nextflow/scripts/merge_perturbation.R nextflow/scripts/intensity_batch.R nextflow/scripts/merge_intensity.R
git commit -m "feat(pipeline): R scripts for perturbation and intensity (batch + merge)"
```

---

### Task 7: R script — final assembly

**Files:**
- Create: `nextflow/scripts/assemble_ranked.R`

This script merges cliques with stability, threshold persistence, density sweep survival, weakest-link persistence, perturbation, and intensity test results into the final ranked table. See the "Output" section above for the column specification.

Key merge notes:
- `clique_persistence()` returns the full cliques df with `persistence` + `mean_persistence` appended (not a separate table). Rename `persistence` to `weakest_link_persistence` to avoid collision with threshold persistence.
- Density sweep survival: count HOGs present in conserved edges (`type == "conserved"`) across the 13 density levels.
- Deduplicate per HOG (keep best) for weakest-link, perturbation, and intensity merges.

- [ ] **Step 1: Write `nextflow/scripts/assemble_ranked.R`**

This script reads all staged result files and produces `cliques_ranked.tsv`. See the "Output" section at the top for the column specification.

```r
#!/usr/bin/env Rscript
# Assemble the final continuous clique ranking table.
# Args: (none)
# Staged inputs: cliques.tsv, stability.rds, threshold_persistence.tsv,
#   density_sweep.rds, clique_persistence.tsv, perturbation_merged.tsv,
#   intensity_merged.tsv, species_list.txt, traits.tsv
# Output: cliques_ranked.tsv

library(rcomplex)

species_names <- readLines("species_list.txt")
species_names <- trimws(species_names[nzchar(species_names)])

traits_df <- read.table("traits.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
annual_sp <- traits_df$species[traits_df$trait == "annual"]

# --- Load all results ---
ranked       <- read.table("cliques.tsv", header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
stab         <- readRDS("stability.rds")
thr_persist  <- read.table("threshold_persistence.tsv", header = TRUE,
                            sep = "\t", stringsAsFactors = FALSE)
ds           <- readRDS("density_sweep.rds")
wl_persist   <- read.table("clique_persistence.tsv", header = TRUE,
                            sep = "\t", stringsAsFactors = FALSE)
perturbation <- read.table("perturbation_merged.tsv", header = TRUE,
                            sep = "\t", stringsAsFactors = FALSE)
intensity    <- read.table("intensity_merged.tsv", header = TRUE,
                            sep = "\t", stringsAsFactors = FALSE)

# --- Trait composition ---
for (i in seq_len(nrow(ranked))) {
  present <- species_names[!is.na(ranked[i, species_names])]
  ranked$n_annual[i]    <- sum(present %in% annual_sp)
  ranked$n_perennial[i] <- sum(!present %in% annual_sp)
}
ranked$annual_fraction   <- ranked$n_annual / ranked$n_species
n_annual_total <- length(annual_sp)
n_perennial_total <- length(species_names) - n_annual_total
ranked$annual_enrichment <- ranked$n_annual / n_annual_total -
                            ranked$n_perennial / n_perennial_total

# --- Stability class ---
ranked$stability_class <- stab$stability_class

# --- Threshold persistence ---
ranked <- merge(ranked,
  thr_persist[, c("hog", "birth", "death", "threshold_persistence")],
  by = "hog", all.x = TRUE)

# --- Density sweep survival ---
hog_at_density <- lapply(ds$edges, function(e) {
  if (is.null(e) || nrow(e) == 0L) return(character(0))
  unique(e$hog[e$type == "conserved"])
})
all_hogs_ds <- unique(unlist(hog_at_density))
if (length(all_hogs_ds) > 0) {
  n_survived <- vapply(all_hogs_ds, function(h) {
    sum(vapply(hog_at_density, function(hogs) h %in% hogs, logical(1)))
  }, integer(1))
  hog_survival <- data.frame(hog = all_hogs_ds,
                              n_densities_survived = n_survived,
                              stringsAsFactors = FALSE)
} else {
  hog_survival <- data.frame(hog = character(0),
                              n_densities_survived = integer(0))
}
ranked <- merge(ranked, hog_survival, by = "hog", all.x = TRUE)
ranked$n_densities_survived[is.na(ranked$n_densities_survived)] <- 0L

# --- Weakest-link persistence ---
wl_best <- wl_persist[order(wl_persist$persistence, decreasing = TRUE), ]
wl_best <- wl_best[!duplicated(wl_best$hog), ]
ranked <- merge(ranked,
  wl_best[, c("hog", "persistence", "mean_persistence")],
  by = "hog", all.x = TRUE)
names(ranked)[names(ranked) == "persistence"] <- "weakest_link_persistence"

# --- Perturbation ---
pert_best <- perturbation[order(perturbation$survival_rate, decreasing = TRUE), ]
pert_best <- pert_best[!duplicated(pert_best$hog), ]
ranked <- merge(ranked,
  pert_best[, c("hog", "survival_rate", "mean_jaccard", "n_matched")],
  by = "hog", all.x = TRUE)

# --- Intensity test ---
int_best <- intensity[order(intensity$z_score, decreasing = TRUE, na.last = TRUE), ]
int_best <- int_best[!duplicated(int_best$hog), ]
ranked <- merge(ranked,
  int_best[, c("hog", "observed_intensity", "z_score", "p_value", "q_value")],
  by = "hog", all.x = TRUE)

# --- Write final table ---
write.table(ranked, "cliques_ranked.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
message("Final ranked table: ", nrow(ranked), " cliques, ",
        ncol(ranked), " columns")
```

- [ ] **Step 2: Verify R syntax**

```bash
Rscript -e "parse('nextflow/scripts/assemble_ranked.R')" && echo "OK" || echo "FAIL"
```

- [ ] **Step 3: Commit**

```bash
git add nextflow/scripts/assemble_ranked.R
git commit -m "feat(pipeline): R script for final ranked table assembly"
```

---

### Task 8: Nextflow modules — data preparation

**Files:**
- Create: `nextflow/modules/reduce.nf`
- Create: `nextflow/modules/network.nf`
- Create: `nextflow/modules/orthologs.nf`

Follow the Orion Hard Rules listed at the top of this plan. Every process must use `file()` for R scripts, `stageInMode = 'copy'` is in config. GPU processes need `torchLd()` in script block.

A shared `traits.tsv` and `species_list.txt` will be generated in main.nf from the samplesheet and passed as channel inputs. Each module should assume these are staged.

- [ ] **Step 1: Write `nextflow/modules/reduce.nf`**

```groovy
process REDUCE {
    tag "${species}"
    publishDir "${params.outdir}/reduce", mode: 'copy'

    input:
    tuple val(species), path(se_rds)
    path script

    output:
    tuple val(species), path("${species}_reduction.rds"), emit: reduction

    script:
    """
    Rscript ${script} ${species} ${se_rds} ${params.cor_threshold} ${params.tissue}
    """
}
```

- [ ] **Step 2: Write `nextflow/modules/network.nf`**

```groovy
process NETWORK {
    tag "${species}"
    label 'gpu'
    publishDir "${params.outdir}/network", mode: 'copy'

    input:
    tuple val(species), path(reduction_rds)
    path script

    output:
    tuple val(species), path("${species}_network.rds"), emit: network

    script:
    """
    ${torchLd()}
    Rscript ${script} ${species} ${reduction_rds} \
        pearson MR ${params.density} ${params.min_var} TRUE TRUE
    """
}
```

- [ ] **Step 3: Write `nextflow/modules/orthologs.nf`**

```groovy
process ORTHOLOGS {
    tag "all_species"
    publishDir "${params.outdir}/orthologs", mode: 'copy'

    input:
    path se_files            // all {species}_se.rds
    path reduction_files     // all {species}_reduction.rds
    path species_list        // species_list.txt
    path script

    output:
    path "orthologs.tsv", emit: orthologs

    script:
    """
    Rscript ${script} hog
    """
}
```

- [ ] **Step 4: Commit**

```bash
git add nextflow/modules/reduce.nf nextflow/modules/network.nf nextflow/modules/orthologs.nf
git commit -m "feat(pipeline): Nextflow modules for data preparation"
```

---

### Task 9: Nextflow modules — comparison

**Files:**
- Create: `nextflow/modules/coexpressologs.nf`
- Create: `nextflow/modules/density_sweep.nf`

- [ ] **Step 1: Write `nextflow/modules/coexpressologs.nf`**

```groovy
process COEXPRESSOLOGS {
    tag "permutation"
    label 'gpu'
    publishDir "${params.outdir}/coexpressologs", mode: 'copy'

    input:
    path network_files       // all {species}_network.rds
    path orthologs           // orthologs.tsv
    path species_list        // species_list.txt
    path script

    output:
    path "edges.tsv", emit: edges

    script:
    """
    ${torchLd()}
    Rscript ${script} permutation greater 0.05 ${task.cpus} TRUE
    """
}
```

- [ ] **Step 2: Write `nextflow/modules/density_sweep.nf`**

```groovy
process DENSITY_SWEEP {
    tag "sweep"
    label 'gpu'
    publishDir "${params.outdir}/density_sweep", mode: 'copy'

    input:
    path network_files
    path orthologs
    path species_list
    path script

    output:
    path "density_sweep.rds", emit: sweep

    script:
    """
    ${torchLd()}
    Rscript ${script} "0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3" \
        permutation ${task.cpus} TRUE
    """
}
```

- [ ] **Step 3: Commit**

```bash
git add nextflow/modules/coexpressologs.nf nextflow/modules/density_sweep.nf
git commit -m "feat(pipeline): Nextflow modules for comparison (coexpressologs, density sweep)"
```

---

### Task 10: Nextflow modules — clique analysis

**Files:**
- Create: `nextflow/modules/cliques.nf`
- Create: `nextflow/modules/stability.nf`
- Create: `nextflow/modules/persistence.nf`

- [ ] **Step 1: Write `nextflow/modules/cliques.nf`**

```groovy
process FIND_CLIQUES {
    tag "cliques"
    publishDir "${params.outdir}/cliques", mode: 'copy'

    input:
    path edges
    path species_list
    path script

    output:
    path "cliques.tsv", emit: cliques

    script:
    """
    Rscript ${script} 2 10
    """
}
```

- [ ] **Step 2: Write `nextflow/modules/stability.nf`**

```groovy
process STABILITY {
    tag "jackknife"
    publishDir "${params.outdir}/stability", mode: 'copy'

    input:
    path edges
    path cliques
    path species_list
    path traits
    path script

    output:
    path "stability.rds", emit: stability

    script:
    """
    Rscript ${script} 2 10 6 ${task.cpus}
    """
}
```

- [ ] **Step 3: Write `nextflow/modules/persistence.nf`**

```groovy
process PERSISTENCE {
    tag "weakest_link"
    publishDir "${params.outdir}/persistence", mode: 'copy'

    input:
    path cliques
    path edges
    path network_files
    path species_list
    path script

    output:
    path "clique_persistence.tsv", emit: persistence

    script:
    """
    Rscript ${script}
    """
}
```

- [ ] **Step 4: Commit**

```bash
git add nextflow/modules/cliques.nf nextflow/modules/stability.nf nextflow/modules/persistence.nf
git commit -m "feat(pipeline): Nextflow modules for clique analysis"
```

---

### Task 11: Nextflow modules — threshold, perturbation, intensity, assembly

**Files:**
- Create: `nextflow/modules/threshold_edges.nf`
- Create: `nextflow/modules/assemble_threshold.nf`
- Create: `nextflow/modules/perturbation.nf`
- Create: `nextflow/modules/intensity.nf`
- Create: `nextflow/modules/assemble.nf`

- [ ] **Step 1: Write `nextflow/modules/threshold_edges.nf`**

```groovy
process THRESHOLD_EDGES {
    tag "m=${multiplier}"
    label 'gpu'
    publishDir "${params.outdir}/threshold_sweep", mode: 'copy'

    input:
    val multiplier
    path network_files
    path orthologs
    path species_list
    path script

    output:
    path "threshold_${multiplier}_edges.tsv",   emit: edges,   optional: true
    path "threshold_${multiplier}_cliques.tsv", emit: cliques, optional: true

    script:
    """
    ${torchLd()}
    Rscript ${script} ${multiplier} 2 10 ${task.cpus} TRUE
    """
}
```

- [ ] **Step 2: Write `nextflow/modules/assemble_threshold.nf`**

```groovy
process ASSEMBLE_THRESHOLD {
    tag "persistence"
    publishDir "${params.outdir}/threshold_sweep", mode: 'copy'

    input:
    path cliques                   // baseline cliques.tsv
    path threshold_cliques_files   // all threshold_*_cliques.tsv
    path species_list
    path script

    output:
    path "threshold_persistence.tsv", emit: persistence

    script:
    """
    Rscript ${script} "${params.threshold_multipliers}"
    """
}
```

- [ ] **Step 3: Write `nextflow/modules/perturbation.nf`**

This module defines both the batch process and the merge process.

```groovy
process PERTURBATION_BATCH {
    tag "batch_${batch_id}"
    publishDir "${params.outdir}/perturbation", mode: 'copy'

    input:
    val batch_id
    path cliques
    path orthologs
    path network_files
    path species_list
    path script

    output:
    path "perturbation_batch.tsv", emit: batch

    script:
    def n_boot_per = (params.n_boot / params.n_boot_batches) as int
    def seed = params.seed + batch_id
    """
    Rscript ${script} ${n_boot_per} ${params.noise_sd} 2 10 ${task.cpus} ${seed}
    mv perturbation_batch.tsv perturbation_batch_${batch_id}.tsv
    ln -s perturbation_batch_${batch_id}.tsv perturbation_batch.tsv
    """
}

process MERGE_PERTURBATION {
    tag "merge"
    publishDir "${params.outdir}/perturbation", mode: 'copy'

    input:
    path batch_files    // all perturbation_batch_*.tsv
    path script

    output:
    path "perturbation_merged.tsv", emit: merged

    script:
    """
    Rscript ${script}
    """
}
```

- [ ] **Step 4: Write `nextflow/modules/intensity.nf`**

```groovy
process INTENSITY_BATCH {
    tag "batch_${batch_id}"
    publishDir "${params.outdir}/intensity", mode: 'copy'

    input:
    val batch_id
    path cliques
    path edges
    path orthologs
    path network_files
    path species_list
    path script

    output:
    path "intensity_batch.tsv", emit: batch

    script:
    def n_perm_per = (params.n_perm / params.n_perm_batches) as int
    def seed = params.seed + batch_id
    """
    Rscript ${script} ${n_perm_per} 2 10 ${task.cpus} ${seed}
    mv intensity_batch.tsv intensity_batch_${batch_id}.tsv
    ln -s intensity_batch_${batch_id}.tsv intensity_batch.tsv
    """
}

process MERGE_INTENSITY {
    tag "merge"
    publishDir "${params.outdir}/intensity", mode: 'copy'

    input:
    path batch_files    // all intensity_batch_*.tsv
    path script

    output:
    path "intensity_merged.tsv", emit: merged

    script:
    """
    Rscript ${script}
    """
}
```

- [ ] **Step 5: Write `nextflow/modules/assemble.nf`**

```groovy
process ASSEMBLE_RANKED {
    tag "final"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path cliques
    path stability
    path threshold_persistence
    path density_sweep
    path clique_persistence
    path perturbation_merged
    path intensity_merged
    path species_list
    path traits
    path script

    output:
    path "cliques_ranked.tsv", emit: ranked

    script:
    """
    Rscript ${script}
    """
}
```

- [ ] **Step 6: Commit**

```bash
git add nextflow/modules/threshold_edges.nf nextflow/modules/assemble_threshold.nf nextflow/modules/perturbation.nf nextflow/modules/intensity.nf nextflow/modules/assemble.nf
git commit -m "feat(pipeline): Nextflow modules for threshold, perturbation, intensity, assembly"
```

---

### Task 12: main.nf — workflow orchestration

**Files:**
- Create: `nextflow/main.nf`

This is the most complex file. It must:
1. Define `cpuBeforeScript`, `gpuBeforeScript`, `torchLd()`
2. Parse the samplesheet into channels
3. Generate `species_list.txt` and `traits.tsv` from samplesheet
4. Stage R scripts as `file()` inputs
5. Wire all processes in the correct DAG order (see Appendix A)
6. Fan-out threshold multipliers and batch IDs into channels

- [ ] **Step 1: Write `nextflow/main.nf`**

The file must contain:

**Header + boilerplate**: `cpuBeforeScript`, `gpuBeforeScript`, `torchLd()`, and all `params.*` declarations:

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def cpuBeforeScript = '''
    module load R/4.4.2
    unset R_HOME
    export RSTUDIO_APPTAINER_BINDPATH=$(
        echo "$RSTUDIO_APPTAINER_BINDPATH" | tr ',' '\\n' |
        grep -v 'rstudio-server' | paste -sd, -
    )
    export RSTUDIO_APPTAINER_BINDPATH="${RSTUDIO_APPTAINER_BINDPATH},/mnt/project:/mnt/project,/net/fs-2/scale:/net/fs-2/scale,/work:/work"
'''

def gpuBeforeScript = '''
    module load R/4.4.2
    unset R_HOME
    export RSTUDIO_APPTAINER_BINDPATH=$(
        echo "$RSTUDIO_APPTAINER_BINDPATH" | tr ',' '\\n' |
        grep -v 'rstudio-server' | paste -sd, -
    )
    export RSTUDIO_APPTAINER_BINDPATH="${RSTUDIO_APPTAINER_BINDPATH},/mnt/project:/mnt/project,/net/fs-2/scale:/net/fs-2/scale,/work:/work"
    export APPTAINER_NV=1
    export CUDA_MODULE_LOADING=LAZY
    TORCH_LIB="/mnt/users/martpali/R/library-4.4/torch/lib"
    if [[ -d "$TORCH_LIB" ]]; then
        for f in "$TORCH_LIB"/lib*-*.so.*; do
            [[ -f "$f" ]] || continue
            base=$(basename "$f")
            link=$(echo "$base" | sed 's/-[0-9a-f]\\{8\\}//')
            [[ -e "$TORCH_LIB/$link" ]] || ln -sf "$base" "$TORCH_LIB/$link"
        done
    fi
'''

def torchLd() {
    'TORCH_LIB="/mnt/users/martpali/R/library-4.4/torch/lib"; ' +
    'export LD_LIBRARY_PATH="${TORCH_LIB}:${LD_LIBRARY_PATH:-}"'
}

params.samplesheet    = 'samplesheet.csv'
params.tissue         = 'both'       // 'leaf', 'root', or 'both'
params.density        = 0.03
params.min_var        = 1e-3
params.cor_threshold  = 0.7
params.n_cores        = 8
params.seed           = 42
params.n_boot         = 1000
params.n_perm         = 1000
params.n_boot_batches = 4
params.n_perm_batches = 5
params.noise_sd       = 0.01
params.threshold_multipliers = '1.2,1.5,2.0,2.5,3.0'
params.outdir         = "results/${params.tissue}"
```

**Include statements**: one per module file.

```groovy
include { REDUCE }               from './modules/reduce'
include { NETWORK }              from './modules/network'
include { ORTHOLOGS }            from './modules/orthologs'
include { COEXPRESSOLOGS }       from './modules/coexpressologs'
include { DENSITY_SWEEP }        from './modules/density_sweep'
include { FIND_CLIQUES }         from './modules/cliques'
include { STABILITY }            from './modules/stability'
include { THRESHOLD_EDGES }      from './modules/threshold_edges'
include { ASSEMBLE_THRESHOLD }   from './modules/assemble_threshold'
include { PERSISTENCE }          from './modules/persistence'
include { PERTURBATION_BATCH }   from './modules/perturbation'
include { MERGE_PERTURBATION }   from './modules/perturbation'
include { INTENSITY_BATCH }      from './modules/intensity'
include { MERGE_INTENSITY }      from './modules/intensity'
include { ASSEMBLE_RANKED }      from './modules/assemble'
```

**Workflow block**: Parse samplesheet, create derived channels (`species_list.txt`, `traits.tsv`), then wire processes:

```groovy
workflow {
    // Parse samplesheet
    Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.species, file(row.se_path), row.trait) }
        .set { samples_ch }

    // Extract species + SE for reduction
    samples_ch
        .map { species, se, trait -> tuple(species, se) }
        .set { reduce_input_ch }

    // Generate species_list.txt
    samples_ch
        .map { species, se, trait -> species }
        .collect()
        .map { species_list ->
            def f = file("${workDir}/species_list.txt")
            f.text = species_list.join('\n') + '\n'
            f
        }
        .set { species_list_ch }

    // Generate traits.tsv
    samples_ch
        .map { species, se, trait -> "${species}\t${trait}" }
        .collect()
        .map { lines ->
            def f = file("${workDir}/traits.tsv")
            f.text = "species\ttrait\n" + lines.join('\n') + '\n'
            f
        }
        .set { traits_ch }

    // Collect SE files for ORTHOLOGS (renamed to {species}_se.rds)
    samples_ch
        .map { species, se, trait -> tuple(species, se) }
        .set { se_for_ortho_ch }

    // --- REDUCE + NETWORK per species ---
    REDUCE(reduce_input_ch, file("scripts/reduce_orthogroups.R"))
    NETWORK(REDUCE.out.reduction, file("scripts/compute_network.R"))

    // Collect all network files
    NETWORK.out.network
        .map { species, net_rds -> net_rds }
        .collect()
        .set { all_networks_ch }

    // Collect all reduction files
    REDUCE.out.reduction
        .map { species, red_rds -> red_rds }
        .collect()
        .set { all_reductions_ch }

    // Collect SE files (rename for ORTHOLOGS staging)
    se_for_ortho_ch
        .map { species, se -> se }
        .collect()
        .set { all_se_ch }

    // --- ORTHOLOGS ---
    ORTHOLOGS(all_se_ch, all_reductions_ch, species_list_ch,
              file("scripts/prepare_orthologs.R"))

    // --- COEXPRESSOLOGS ---
    COEXPRESSOLOGS(all_networks_ch, ORTHOLOGS.out.orthologs,
                   species_list_ch, file("scripts/find_coexpressologs.R"))

    // --- DENSITY_SWEEP (parallel with FIND_CLIQUES) ---
    DENSITY_SWEEP(all_networks_ch, ORTHOLOGS.out.orthologs,
                  species_list_ch, file("scripts/density_sweep.R"))

    // --- FIND_CLIQUES ---
    FIND_CLIQUES(COEXPRESSOLOGS.out.edges, species_list_ch,
                 file("scripts/find_cliques.R"))

    // --- STABILITY ---
    STABILITY(COEXPRESSOLOGS.out.edges, FIND_CLIQUES.out.cliques,
              species_list_ch, traits_ch, file("scripts/clique_stability.R"))

    // --- PERSISTENCE ---
    PERSISTENCE(FIND_CLIQUES.out.cliques, COEXPRESSOLOGS.out.edges,
                all_networks_ch, species_list_ch,
                file("scripts/clique_persistence.R"))

    // --- THRESHOLD EDGES (fan-out over multipliers) ---
    threshold_mult_ch = Channel.of(
        params.threshold_multipliers.tokenize(',').collect { it.trim() as float }
    ).flatten()

    THRESHOLD_EDGES(threshold_mult_ch, all_networks_ch,
                    ORTHOLOGS.out.orthologs, species_list_ch,
                    file("scripts/threshold_edges.R"))

    // Collect threshold clique files for assembly
    THRESHOLD_EDGES.out.cliques
        .collect()
        .set { all_threshold_cliques_ch }

    ASSEMBLE_THRESHOLD(FIND_CLIQUES.out.cliques, all_threshold_cliques_ch,
                       species_list_ch, file("scripts/assemble_threshold.R"))

    // --- PERTURBATION BATCHES ---
    pert_batch_ch = Channel.of(1..params.n_boot_batches)

    PERTURBATION_BATCH(pert_batch_ch, FIND_CLIQUES.out.cliques,
                       ORTHOLOGS.out.orthologs, all_networks_ch,
                       species_list_ch, file("scripts/perturbation_batch.R"))

    PERTURBATION_BATCH.out.batch
        .collect()
        .set { all_pert_batches_ch }

    MERGE_PERTURBATION(all_pert_batches_ch,
                       file("scripts/merge_perturbation.R"))

    // --- INTENSITY BATCHES ---
    int_batch_ch = Channel.of(1..params.n_perm_batches)

    INTENSITY_BATCH(int_batch_ch, FIND_CLIQUES.out.cliques,
                    COEXPRESSOLOGS.out.edges, ORTHOLOGS.out.orthologs,
                    all_networks_ch, species_list_ch,
                    file("scripts/intensity_batch.R"))

    INTENSITY_BATCH.out.batch
        .collect()
        .set { all_int_batches_ch }

    MERGE_INTENSITY(all_int_batches_ch,
                    file("scripts/merge_intensity.R"))

    // --- ASSEMBLE RANKED TABLE ---
    ASSEMBLE_RANKED(
        FIND_CLIQUES.out.cliques,
        STABILITY.out.stability,
        ASSEMBLE_THRESHOLD.out.persistence,
        DENSITY_SWEEP.out.sweep,
        PERSISTENCE.out.persistence,
        MERGE_PERTURBATION.out.merged,
        MERGE_INTENSITY.out.merged,
        species_list_ch,
        traits_ch,
        file("scripts/assemble_ranked.R")
    )
}
```

**Important:** The `torchLd()` function and `beforeScript` definitions must be at the top of main.nf (before the `workflow` block), along with all `params.*` defaults. Copy them exactly from the plan.

- [ ] **Step 2: Verify Nextflow syntax**

```bash
cd nextflow && nextflow -version 2>/dev/null || echo "nextflow not installed locally (OK, will test on Orion)"
```

If nextflow is available locally:
```bash
cd nextflow && nextflow run main.nf -profile orion -preview
```

- [ ] **Step 3: Commit**

```bash
git add nextflow/main.nf
git commit -m "feat(pipeline): main.nf workflow orchestration"
```

---

### Task 13: Validation and cleanup

- [ ] **Step 1: Verify all R scripts parse**

```bash
for f in nextflow/scripts/*.R; do
  Rscript -e "parse('$f')" 2>&1 | grep -q "Error" && echo "FAIL: $f" || echo "OK: $f"
done
```

Expected: `OK` for all 15 scripts.

- [ ] **Step 2: Verify file structure matches the plan**

```bash
echo "=== Expected files ==="
cat <<'EOF'
nextflow/main.nf
nextflow/nextflow.config
nextflow/samplesheet.csv
nextflow/modules/reduce.nf
nextflow/modules/network.nf
nextflow/modules/orthologs.nf
nextflow/modules/coexpressologs.nf
nextflow/modules/density_sweep.nf
nextflow/modules/cliques.nf
nextflow/modules/stability.nf
nextflow/modules/threshold_edges.nf
nextflow/modules/assemble_threshold.nf
nextflow/modules/persistence.nf
nextflow/modules/perturbation.nf
nextflow/modules/intensity.nf
nextflow/modules/assemble.nf
nextflow/scripts/reduce_orthogroups.R
nextflow/scripts/compute_network.R
nextflow/scripts/prepare_orthologs.R
nextflow/scripts/find_coexpressologs.R
nextflow/scripts/density_sweep.R
nextflow/scripts/find_cliques.R
nextflow/scripts/clique_stability.R
nextflow/scripts/threshold_edges.R
nextflow/scripts/assemble_threshold.R
nextflow/scripts/clique_persistence.R
nextflow/scripts/perturbation_batch.R
nextflow/scripts/merge_perturbation.R
nextflow/scripts/intensity_batch.R
nextflow/scripts/merge_intensity.R
nextflow/scripts/assemble_ranked.R
EOF

echo "=== Actual files ==="
find nextflow -type f | sort
```

Verify all expected files exist and no extras.

- [ ] **Step 3: Review the PERTURBATION_BATCH mv/ln-s pattern**

In `nextflow/modules/perturbation.nf`, the batch process renames the output to include the batch ID. Verify the glob pattern in `merge_perturbation.R` (`Sys.glob("perturbation_batch_*.tsv")`) will match the staged filenames. Same for `intensity.nf` / `merge_intensity.R`.

If the Nextflow staging renames files, adjust the merge scripts' glob patterns accordingly.

- [ ] **Step 4: Final commit if any fixes were needed**

```bash
git add -A nextflow/
git status
# Only commit if there are changes:
git diff --cached --quiet || git commit -m "fix(pipeline): validation fixes"
```

---

## Appendix A: Pipeline DAG

```
                    +--- REDUCE(BDIS) -> NETWORK(BDIS) ---+
                    +--- REDUCE(BSYL) -> NETWORK(BSYL) ---+
samplesheet --+---  +--- REDUCE(HVUL) -> NETWORK(HVUL) ---+
              |     +--- REDUCE(HJUB) -> NETWORK(HJUB) ---+
              |     +--- REDUCE(BMAX) -> NETWORK(BMAX) ---+
              |     +--- REDUCE(BMED) -> NETWORK(BMED) ---+
              |     +--- REDUCE(VBRO) -> NETWORK(VBRO) ---+
              |     +--- REDUCE(FPRA) -> NETWORK(FPRA) ---+
              |                                            |
              +---- ORTHOLOGS -------------------> COEXPRESSOLOGS --+-- DENSITY_SWEEP
                                                        |           |
                                                        v           |
                                                   FIND_CLIQUES    |
                                                        |           |
                                          +-------------+-------+   |
                                          v             v       v   |
                                    STABILITY    PERSISTENCE    |   |
                                          |             |       |   |
         +-- THRESHOLD_EDGES(1.2) -+      |             |       |   |
         +-- THRESHOLD_EDGES(1.5) -+      |             |       |   |
         +-- THRESHOLD_EDGES(2.0) -+--> ASSEMBLE_THR    |       |   |
         +-- THRESHOLD_EDGES(2.5) -+      |             |       |   |
         +-- THRESHOLD_EDGES(3.0) -+      |             |       |   |
                                          |             |       |   |
         +-- PERTURBATION_BATCH(1..4) --> MERGE_PERT    |       |   |
                                          |             |       |   |
         +-- INTENSITY_BATCH(1..5) -----> MERGE_INT     |       |   |
                                          |             |       |   |
                                          v             v       v   v
                                    ASSEMBLE_RANKED <-----------+---+
```

PERTURBATION_BATCH and INTENSITY_BATCH run in parallel (both depend only on cliques + networks + orthologs + edges). THRESHOLD_EDGES jobs also run in parallel (depend on networks + orthologs only).

---

## Appendix B: Resource Estimates (8 species, ~20K genes, per tissue run)

| Process | Wall time | CPUs | RAM | Partition | GPU |
|---------|-----------|------|-----|-----------|-----|
| REDUCE (x8) | 1 min each | 1 | 8 GB | orion | no |
| NETWORK (x8, 6 concurrent) | 2 min each | 1 | 16 GB | TestGPU | 1x RTX PRO 6000 |
| ORTHOLOGS | <10s | 1 | 4 GB | orion | no |
| COEXPRESSOLOGS | 30 min | 16 | 32 GB | TestGPU | 1x RTX PRO 6000 |
| DENSITY_SWEEP | 2 h | 16 | 32 GB | TestGPU | 1x RTX PRO 6000 |
| FIND_CLIQUES | <30s | 1 | 4 GB | orion | no |
| STABILITY | 30 min | 16 | 32 GB | orion | no |
| THRESHOLD_EDGES (x5) | 30 min each | 16 | 32 GB | TestGPU | 1x RTX PRO 6000 |
| ASSEMBLE_THRESHOLD | <1 min | 1 | 4 GB | orion | no |
| PERSISTENCE | <1 min | 1 | 4 GB | orion | no |
| PERTURBATION_BATCH (x4) | 30 min each | 8 | 16 GB | orion | no |
| INTENSITY_BATCH (x5) | 20 min each | 8 | 16 GB | orion | no |
| MERGE_PERTURBATION | <1 min | 1 | 4 GB | orion | no |
| MERGE_INTENSITY | <1 min | 1 | 4 GB | orion | no |
| ASSEMBLE_RANKED | <1 min | 1 | 8 GB | orion | no |
| **Total (parallel)** | **~4 h** | | | | |
| **Total (serial)** | **~9 h** | | | | |

GPU time budget: 8 NETWORK + 1 COEXPRESSOLOGS + 1 DENSITY_SWEEP + 5 THRESHOLD_EDGES = 15 GPU jobs total. gn-41 has 6 GPUs; SLURM packs concurrently.

Three tissue runs (leaf, root, both) = ~12 h parallel / ~27 h serial total.

---

## Appendix C: Running the Pipeline

```bash
# One-time data preparation (run locally or on login node)
Rscript prepare_data/convert_to_se.R

# Three tissue runs on Orion
cd nextflow
nextflow run main.nf -profile orion --tissue leaf  --outdir results/leaf
nextflow run main.nf -profile orion --tissue root  --outdir results/root
nextflow run main.nf -profile orion --tissue both  --outdir results/both
```

Nextflow caching (`cache = 'lenient'`) means the ORTHOLOGS step (HOG extraction from SEs) will be cached across runs since SE rowData is tissue-independent. REDUCE and everything downstream re-runs because paralog correlations differ by tissue.
