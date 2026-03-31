# rcomplex Audit: Ranked Improvement Plan

Consolidated from 5-agent audit (R idioms, C++ idioms, API consistency, data structures, parallelisation). Findings that conflict with CLAUDE.md intentional decisions have been discarded.

## Tier 1: Worth Fixing

### 1. ~~Fork + OpenMP oversubscription~~ DONE (d26cd5a)
`num_threads()` clause + `OMP_NUM_THREADS=1` in mclapply forks.

### 2. ~~Global `omp_set_num_threads()` leaks between functions~~ DONE (d26cd5a)
Replaced all 7 calls with `num_threads(n_cores)` pragma clause.

### 3. ~~Besag-Clifford adaptive stopping causes thread starvation~~ DONE (3638e3b)
HOGs sorted by descending size, `schedule(guided)` for permutation loops.

### 4. ~~Missing `comparison_to_edges()` helper~~ DONE (020de30)
Added `comparison_to_edges(comparison, sp1, sp2)` — bridges comparison→cliques.

### 5. ~~Inconsistent ortholog group naming (`OrthoGroup` vs `hog`)~~ DONE (020de30)
Renamed to `hog` throughout all R code, tests, and fixture files.

### 6. ~~Implicit species identity in comparison output~~ RESOLVED via item 4
Species identity is injected by `comparison_to_edges(sp1, sp2)` — no need to change `compare_neighborhoods()` output.

### 7. ~~Missing `@examples` on most exported functions~~ DONE
Added `\dontrun{}` examples to all 13 exported functions.

### 8. ~~CLR inner loop could use Armadillo vectorisation~~ DONE (f8de0ff)
Z-scoring vectorised with `each_col`/`transform`, in-place hypot.

## Tier 2: Minor

### 9. `schedule(dynamic)` on uniform-cost loops
**Files**: `src/neighborhood_comparison.cpp:146,158`, `src/clr.cpp:61`, `src/mutual_rank.cpp:105,112`
**Problem**: Neighbour extraction and rank computation have near-uniform cost per gene. `dynamic` incurs unnecessary synchronisation overhead vs `static`.
**Fix**: Switch to `schedule(static)` for uniform-cost loops; keep `dynamic` only for variable-cost work (permutations, clustering).

### 10. GPU-CPU sync bottleneck in FE tiling
**File**: `R/summary.R:229-278`
**Problem**: `$cpu()` is synchronous inside each tile iteration. Multiple GPU→CPU roundtrips where one would suffice.
**Fix**: Accumulate FE tiles on GPU, single `$cpu()` transfer after loop completes.

### 11. Function signature inconsistency across clique functions
**Files**: `R/cliques.R:90, 285, 482, 643`
**Problem**: `target_species` is position 2 in `find_cliques()` / `clique_stability()` but position 3 in `clique_persistence()`. `edges` is position 1 in `find_cliques()` but position 4 in `clique_persistence()`.
**Fix**: Harmonise parameter ordering. Consider: data objects first (`cliques`, `edges`, `networks`), then `target_species`, then tuning parameters.

### 12. GPU memory not freed on error paths
**File**: `R/summary.R:463-481`
**Problem**: If `fe_hog_permutation_test_cpp()` throws, the `combined` tensor stays on GPU.
**Fix**: `on.exit(rm(combined); gc())` before the C++ call.

### 13. MPS float32 precision not enforced at runtime
**Files**: `R/network.R:120-127`, `R/summary.R:318-333`
**Problem**: No warning if user runs `permutation_hog_test(use_torch=TRUE)` on MPS after computing Spearman + MR correlations in float32 (which has rank-swap artifacts).
**Fix**: Check `net$params$cor_method` and device; warn if the combination is lossy.

### 14. Bit-vector mode memory scaling
**File**: `src/hog_permutation.cpp:314-335`
**Problem**: Bit-vector mode allocates O(n^2/64) memory. At n=100K, this is ~6GB. No diagnostic when falling back to the slower flag-vector path.
**Fix**: `Rcpp::message()` when switching modes, or auto-tune threshold based on available memory.

### 15. `std::replace()` instead of manual cluster ID loop
**File**: `src/reduce_orthogroups.cpp:53-56`
**Problem**: Manual `for` loop to replace cluster IDs where `std::replace(cid.begin(), cid.end(), old_id, new_id)` is the STL idiom.
**Fix**: One-line change to `std::replace()`.

### 16. Error messages could include received type
**Files**: `R/cliques.R:486,645`, `R/summary.R:80,408`
**Problem**: `stop("cliques must be a data frame from find_cliques()")` doesn't say what was received.
**Fix**: Include `class(cliques)` in the error message.

## Tier 3: Cosmetic

### 17. Remove `stringsAsFactors = FALSE` throughout
**Files**: All R files with `data.frame()` calls
**Problem**: Redundant since `Depends: R (>= 4.1.0)` and R 4.0+ defaults to FALSE.

### 18. Inconsistent internal function documentation (`@noRd` vs `@keywords internal`)
**Files**: `R/network.R:13`, `R/summary.R:9`, `R/cliques.R:14`
**Fix**: Standardise on `@noRd`.

### 19. Redundant `as.integer()` on already-integer parameters
**Files**: Multiple; e.g., `n_cores = 1L` followed by `n_cores <- as.integer(n_cores)`
**Fix**: Accept numeric and coerce once, or remove the coercion if default is already `1L`.

### 20. `permutation_hog_test()` breaks verb-noun naming pattern
**Problem**: Other functions use verb-noun (`compute_network`, `detect_modules`, `find_cliques`). This one uses noun-noun-noun.
**Fix**: Consider `test_hog_permutation()` for consistency, with a deprecation alias.

### 21. Empty data frame construction patterns inconsistent
**Files**: `R/cliques.R:111-117` vs `R/cliques.R:335-341`
**Fix**: Standardise on `data.frame(col = type(0), ...)`.

### 22. No early stopping for zero T_obs in permutation tests
**Files**: `src/hog_permutation.cpp:388-396`
**Problem**: If T_obs = 0, every permutation exceeds it. Running max_permutations is wasted.
**Fix**: `if (T_obs <= 0.0) { p_value = 1.0; continue; }` after computing T_obs.

### 23. Co-classification could parallelise outer K loop
**File**: `src/coclassification.cpp:143-154`
**Problem**: K sequential edge-loops (one per resolution) where the outer K loop could be parallelised instead.
**Fix**: Restructure to `#pragma omp parallel for` over K with per-thread accumulation.
