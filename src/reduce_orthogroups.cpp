// reduce_orthogroups.cpp
// Ward.D2 agglomerative clustering to merge correlated paralogs within HOGs.
// Matches R's hclust(method="ward.D2") + cutree(h = 1 - cor_threshold).

#include <RcppArmadillo.h>
#include <numeric>
#include <vector>

using namespace Rcpp;

// Ward.D2 clustering with early stopping at cut height h.
// Input: unsquared Pearson distance matrix (1 - cor). Squares internally.
// Returns: cluster assignments (contiguous 0..k-1) and number of clusters.
static int ward_d2_cut(const arma::mat& dist_unsq, double h, std::vector<int>& cid) {
    int n = dist_unsq.n_rows;
    cid.resize(n);
    std::iota(cid.begin(), cid.end(), 0);
    if (n <= 1) return n;

    arma::mat D = dist_unsq % dist_unsq;
    double h_sq = h * h;

    std::vector<bool> active(n, true);
    std::vector<int> sizes(n, 1);

    for (int step = 0; step < n - 1; step++) {
        double min_d = arma::datum::inf;
        int mi = -1, mj = -1;
        for (int i = 0; i < n; i++) {
            if (!active[i]) continue;
            for (int j = i + 1; j < n; j++) {
                if (!active[j]) continue;
                if (D(i, j) < min_d) {
                    min_d = D(i, j);
                    mi = i; mj = j;
                }
            }
        }

        if (mi < 0 || min_d > h_sq) break;

        // Lance-Williams update for Ward on squared distances
        int ni = sizes[mi], nj = sizes[mj];
        for (int k = 0; k < n; k++) {
            if (!active[k] || k == mi || k == mj) continue;
            int nk = sizes[k];
            double new_d = ((nk + ni) * D(k, mi) +
                            (nk + nj) * D(k, mj) -
                            nk * D(mi, mj)) / (nk + ni + nj);
            D(k, mi) = D(mi, k) = new_d;
        }

        int old_id = cid[mj], new_id = cid[mi];
        for (int i = 0; i < n; i++) {
            if (cid[i] == old_id) cid[i] = new_id;
        }
        sizes[mi] += sizes[mj];
        active[mj] = false;
    }

    // Renumber to contiguous 0..k-1 using flat vector (IDs are in [0, n))
    std::vector<int> remap(n, -1);
    int n_clusters = 0;
    for (int i = 0; i < n; i++) {
        if (remap[cid[i]] < 0) remap[cid[i]] = n_clusters++;
        cid[i] = remap[cid[i]];
    }
    return n_clusters;
}


//' Reduce orthogroups by merging correlated paralogs
//'
//' For each HOG with multiple paralogs, computes pairwise Pearson correlation
//' on the expression matrix, clusters with Ward.D2, and merges paralogs
//' whose correlation exceeds \code{cor_threshold} into a single averaged
//' representative.
//'
//' @param expr Numeric matrix (genes x samples). Row order matches gene_ids.
//' @param hog_members List of integer vectors. Each element contains 1-indexed
//'   row positions in \code{expr} for genes belonging to one HOG.
//' @param non_hog_idx Integer vector of 1-indexed row positions for genes
//'   not assigned to any HOG (kept as-is).
//' @param cor_threshold Pearson correlation threshold for merging (default 0.7).
//' @return List with: expr_matrix (reduced matrix), out_row_source (1-indexed
//'   original row for each output row), map_from / map_to (1-indexed gene
//'   mapping), n_original, n_reduced, n_merged.
//' @keywords internal
// [[Rcpp::export]]
List reduce_orthogroups_cpp(const arma::mat& expr,
                            const List& hog_members,
                            const IntegerVector& non_hog_idx,
                            double cor_threshold) {
    int n_samples = expr.n_cols;
    int n_original = expr.n_rows;
    double h = 1.0 - cor_threshold;
    int n_merged = 0;

    // Pre-allocate output matrix at upper bound, fill rows directly
    arma::mat out_mat(n_original, n_samples);
    int out_row = 0;
    std::vector<int> out_source;
    std::vector<int> map_from, map_to;
    out_source.reserve(n_original);
    map_from.reserve(n_original);
    map_to.reserve(n_original);

    // Emit a single gene unchanged
    auto emit = [&](int gi) {
        out_mat.row(out_row++) = expr.row(gi);
        out_source.push_back(gi);
        map_from.push_back(gi);
        map_to.push_back(gi);
    };

    int n_hogs = hog_members.size();
    std::vector<int> clusters;  // reused across HOGs

    for (int hi = 0; hi < n_hogs; hi++) {
        if (hi % 500 == 0) Rcpp::checkUserInterrupt();

        IntegerVector members = hog_members[hi];
        int ng = members.size();
        if (ng == 0) continue;
        if (ng == 1) { emit(members[0] - 1); continue; }

        // Separate zero-variance genes (cannot correlate)
        std::vector<int> var_gi, zv_gi;
        for (int i = 0; i < ng; i++) {
            int gi = members[i] - 1;
            if (arma::as_scalar(arma::var(expr.row(gi))) > 0) var_gi.push_back(gi);
            else zv_gi.push_back(gi);
        }
        for (int gi : zv_gi) emit(gi);

        if (var_gi.size() <= 1) {
            if (var_gi.size() == 1) emit(var_gi[0]);
            continue;
        }

        // Pearson correlation → distance → Ward.D2 clustering
        int nv = static_cast<int>(var_gi.size());
        arma::mat sub_expr(nv, n_samples);
        for (int i = 0; i < nv; i++) sub_expr.row(i) = expr.row(var_gi[i]);

        arma::mat dist_mat = 1.0 - arma::cor(sub_expr.t());
        dist_mat.clamp(0.0, 2.0);

        int n_clusters = ward_d2_cut(dist_mat, h, clusters);

        // Group by cluster and emit (flat vector, no std::map)
        std::vector<std::vector<int>> groups(n_clusters);
        for (int i = 0; i < nv; i++) groups[clusters[i]].push_back(i);

        for (int c = 0; c < n_clusters; c++) {
            const auto& grp = groups[c];
            int rep_gi = var_gi[grp[0]];

            if (grp.size() == 1) {
                emit(rep_gi);
            } else {
                arma::rowvec avg = arma::zeros<arma::rowvec>(n_samples);
                for (int m : grp) avg += expr.row(var_gi[m]);
                avg /= static_cast<double>(grp.size());

                out_mat.row(out_row++) = avg;
                out_source.push_back(rep_gi);
                for (int m : grp) {
                    map_from.push_back(var_gi[m]);
                    map_to.push_back(rep_gi);
                }
                n_merged += static_cast<int>(grp.size()) - 1;
            }
        }
    }

    // Non-HOG genes: keep as-is
    for (int i = 0; i < static_cast<int>(non_hog_idx.size()); i++) {
        emit(non_hog_idx[i] - 1);
    }

    // Trim output matrix to actual size
    out_mat = out_mat.head_rows(out_row);

    // Convert to 1-indexed for R
    IntegerVector out_source_r(out_source.begin(), out_source.end());
    IntegerVector map_from_r(map_from.begin(), map_from.end());
    IntegerVector map_to_r(map_to.begin(), map_to.end());
    out_source_r = out_source_r + 1;
    map_from_r = map_from_r + 1;
    map_to_r = map_to_r + 1;

    return List::create(
        Named("expr_matrix") = out_mat,
        Named("out_row_source") = out_source_r,
        Named("map_from") = map_from_r,
        Named("map_to") = map_to_r,
        Named("n_original") = n_original,
        Named("n_reduced") = out_row,
        Named("n_merged") = n_merged
    );
}
