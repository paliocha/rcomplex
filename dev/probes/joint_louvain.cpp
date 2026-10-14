// Probe E2: serial Louvain on a species-blocked joint objective.
// Written from the formula in the probe spec (joint_louvain.R header).
//
//   Q = sum_s (1 / 2m_s) sum_{i,j in s} (A_ij - g k_i k_j / 2m_s) d(c_i, c_j)
//     + (l / 2m_bar) sum_{i in s, j in t, s != t}
//         (w_ij - h_i^t h_j^s / H_st) d(c_i, c_j)
//
// Both terms fold into one graph. An edge carries A_ij / 2m_s (same
// species) or l w_ij / 2m_bar (hog pair). Each node carries a profile u
// of length S + S^2: vol_s, then G[s, t] at S + s S + t. The null
// between node P and community C is the bilinear form
//   sum_s g vol_s(P) vol_s(C) / (2m_s)^2
//   + (l / 2m_bar) sum_{s != t} G_P[s, t] G_C[t, s] / H_st.
// The gain of P into C is e(P, C) - null(P, C); the change in Q is twice
// that. Supernodes keep their internal weight as a self-loop.

#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

namespace {

struct Graph {
  int n;
  std::vector<int> p, j;      // CSR, no self entries
  std::vector<double> w, self;
  std::vector<double> u;      // n x D profiles, row-major
};

struct Null {
  int S, D;
  // coefficient of u_P[a] * U_C[b]: idx a pairs with partner b
  std::vector<int> partner;
  std::vector<double> coef;
};

// Sparse view of node P's null row: (index into U_C, weight)
void null_row(const Null& nl, const double* uP, std::vector<int>& idx,
              std::vector<double>& val) {
  idx.clear();
  val.clear();
  for (int a = 0; a < nl.D; ++a) {
    if (uP[a] != 0.0 && nl.coef[a] != 0.0) {
      idx.push_back(nl.partner[a]);
      val.push_back(uP[a] * nl.coef[a]);
    }
  }
}

// One level of local moves from the partition comm (ids < n). Returns
// the number of moves.
int local_moves(const Graph& g, const Null& nl, std::vector<int>& comm,
                int max_passes, double tol) {
  const int n = g.n, D = nl.D;
  std::vector<double> U(static_cast<size_t>(n) * D, 0.0);
  std::vector<int> size(n, 0), empty;
  for (int i = 0; i < n; ++i) {
    size[comm[i]]++;
    for (int a = 0; a < D; ++a) {
      U[static_cast<size_t>(comm[i]) * D + a] +=
        g.u[static_cast<size_t>(i) * D + a];
    }
  }
  for (int c = n - 1; c >= 0; --c) {
    if (size[c] == 0) empty.push_back(c);
  }
  std::vector<double> nw(n, 0.0);
  std::vector<int> touched;
  std::vector<char> seen(n, 0);
  std::vector<int> order(n), vidx;
  std::vector<double> vval;
  for (int i = 0; i < n; ++i) order[i] = i;
  int total_moves = 0;
  for (int pass = 0; pass < max_passes; ++pass) {
    for (int i = n - 1; i > 0; --i) {
      int k = static_cast<int>(R::unif_rand() * (i + 1));
      if (k > i) k = i;
      std::swap(order[i], order[k]);
    }
    int moves = 0;
    double gain_sum = 0.0;
    for (int oi = 0; oi < n; ++oi) {
      const int i = order[oi], d = comm[i];
      const double* ui = &g.u[static_cast<size_t>(i) * D];
      double* Ud = &U[static_cast<size_t>(d) * D];
      for (int a = 0; a < D; ++a) Ud[a] -= ui[a];
      size[d]--;
      null_row(nl, ui, vidx, vval);
      touched.clear();
      for (int e = g.p[i]; e < g.p[i + 1]; ++e) {
        int c = comm[g.j[e]];
        if (!seen[c]) {
          seen[c] = 1;
          touched.push_back(c);
        }
        nw[c] += g.w[e];
      }
      auto gain = [&](int c) {
        const double* Uc = &U[static_cast<size_t>(c) * D];
        double nu = 0.0;
        for (size_t a = 0; a < vidx.size(); ++a) nu += vval[a] * Uc[vidx[a]];
        return nw[c] - nu;
      };
      const double g_d = size[d] > 0 ? gain(d) : 0.0;
      int best = d;
      double best_g = g_d;
      for (int c : touched) {
        if (c == d) continue;
        double gc = gain(c);
        if (gc > best_g) {
          best_g = gc;
          best = c;
        }
      }
      if (size[d] > 0 && best_g < 0.0 && !empty.empty()) {
        best = empty.back();
        best_g = 0.0;
      }
      for (int c : touched) {
        nw[c] = 0.0;
        seen[c] = 0;
      }
      if (best != d) {
        if (size[best] == 0) empty.pop_back();
        if (size[d] == 0) empty.push_back(d);
        ++moves;
        gain_sum += best_g - g_d;
      }
      comm[i] = best;
      size[best]++;
      double* Ub = &U[static_cast<size_t>(best) * D];
      for (int a = 0; a < D; ++a) Ub[a] += ui[a];
    }
    total_moves += moves;
    if (moves == 0 || 2.0 * gain_sum < tol) break;
  }
  return total_moves;
}

// Relabel comm to 0..K-1 in order of first appearance; return K
int compact(std::vector<int>& comm) {
  std::vector<int> map(comm.size(), -1);
  int k = 0;
  for (int& c : comm) {
    if (map[c] < 0) map[c] = k++;
    c = map[c];
  }
  return k;
}

Graph aggregate(const Graph& g, const std::vector<int>& comm, int K,
                int D) {
  Graph h;
  h.n = K;
  h.self.assign(K, 0.0);
  h.u.assign(static_cast<size_t>(K) * D, 0.0);
  h.p.assign(K + 1, 0);
  std::vector<std::vector<int>> mem(K);
  for (int i = 0; i < g.n; ++i) mem[comm[i]].push_back(i);
  std::vector<double> nw(K, 0.0);
  std::vector<char> seen(K, 0);
  std::vector<int> touched;
  for (int a = 0; a < K; ++a) {
    touched.clear();
    for (int i : mem[a]) {
      h.self[a] += g.self[i];
      const double* ui = &g.u[static_cast<size_t>(i) * D];
      double* ua = &h.u[static_cast<size_t>(a) * D];
      for (int x = 0; x < D; ++x) ua[x] += ui[x];
      for (int e = g.p[i]; e < g.p[i + 1]; ++e) {
        int b = comm[g.j[e]];
        if (b == a) {
          h.self[a] += g.w[e];
          continue;
        }
        if (!seen[b]) {
          seen[b] = 1;
          touched.push_back(b);
        }
        nw[b] += g.w[e];
      }
    }
    std::sort(touched.begin(), touched.end());
    for (int b : touched) {
      h.j.push_back(b);
      h.w.push_back(nw[b]);
      nw[b] = 0.0;
      seen[b] = 0;
    }
    h.p[a + 1] = static_cast<int>(h.j.size());
  }
  return h;
}

}  // namespace

// Quality of a label vector at the gene level, split into parts.
// co_* and hog_*: CSC slots of the symmetric gene x gene matrices.
// [[Rcpp::export]]
List jl_quality_cpp(IntegerVector lab, IntegerVector sp, int S,
                    IntegerVector co_p, IntegerVector co_i,
                    NumericVector co_x, IntegerVector hog_p,
                    IntegerVector hog_i, NumericVector hog_x,
                    NumericVector h, NumericVector H, double lambda,
                    double gamma) {
  const int n = lab.size();
  int K = 0;
  for (int i = 0; i < n; ++i) K = std::max(K, lab[i] + 1);
  std::vector<double> two_m(S, 0.0), in_co(S, 0.0);
  std::vector<double> vol(static_cast<size_t>(K) * S, 0.0);
  std::vector<double> G(static_cast<size_t>(K) * S * S, 0.0);
  double in_hog = 0.0;
  for (int i = 0; i < n; ++i) {
    const int s = sp[i], c = lab[i];
    double k = 0.0;
    for (int e = co_p[i]; e < co_p[i + 1]; ++e) {
      k += co_x[e];
      if (lab[co_i[e]] == c) in_co[s] += co_x[e];
    }
    two_m[s] += k;
    vol[static_cast<size_t>(c) * S + s] += k;
    for (int e = hog_p[i]; e < hog_p[i + 1]; ++e) {
      if (lab[hog_i[e]] == c) in_hog += hog_x[e];
    }
    for (int t = 0; t < S; ++t) {
      G[(static_cast<size_t>(c) * S + s) * S + t] +=
        h[static_cast<size_t>(t) * n + i];
    }
  }
  NumericVector parts(S);
  double bar = 0.0;
  for (int s = 0; s < S; ++s) bar += two_m[s] / S;
  for (int s = 0; s < S; ++s) {
    double nul = 0.0;
    for (int c = 0; c < K; ++c) {
      double v = vol[static_cast<size_t>(c) * S + s];
      nul += v * v;
    }
    parts[s] = in_co[s] / two_m[s] - gamma * nul / (two_m[s] * two_m[s]);
  }
  double nul = 0.0;
  for (int c = 0; c < K; ++c) {
    const double* Gc = &G[static_cast<size_t>(c) * S * S];
    for (int s = 0; s < S; ++s) {
      for (int t = 0; t < S; ++t) {
        double Hst = H[t * S + s];
        if (s != t && Hst > 0) nul += Gc[s * S + t] * Gc[t * S + s] / Hst;
      }
    }
  }
  double hog_raw = (in_hog - nul) / bar;
  double q = lambda * hog_raw;
  for (int s = 0; s < S; ++s) q += parts[s];
  return List::create(_["quality"] = q, _["species"] = parts,
                      _["hog_raw"] = hog_raw, _["hog"] = lambda * hog_raw);
}

// Louvain on the joint objective. With freeze, the first level is the
// aggregate of init (0-based labels), so single genes never move.
// Returns 0-based compact labels and the number of levels.
// [[Rcpp::export]]
List jl_louvain_cpp(IntegerVector sp, int S, IntegerVector co_p,
                    IntegerVector co_i, NumericVector co_x,
                    IntegerVector hog_p, IntegerVector hog_i,
                    NumericVector hog_x, NumericVector h, NumericVector H,
                    double lambda, double gamma, IntegerVector init,
                    bool freeze, int max_passes, double tol) {
  RNGScope scope;
  const int n = sp.size(), D = S + S * S;
  std::vector<double> two_m(S, 0.0), k(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int e = co_p[i]; e < co_p[i + 1]; ++e) k[i] += co_x[e];
    two_m[sp[i]] += k[i];
  }
  double bar = 0.0;
  for (int s = 0; s < S; ++s) bar += two_m[s] / S;
  const double lb = lambda / bar;

  Null nl;
  nl.S = S;
  nl.D = D;
  nl.partner.assign(D, 0);
  nl.coef.assign(D, 0.0);
  for (int s = 0; s < S; ++s) {
    nl.partner[s] = s;
    nl.coef[s] = gamma / (two_m[s] * two_m[s]);
    for (int t = 0; t < S; ++t) {
      int a = S + s * S + t;
      nl.partner[a] = S + t * S + s;
      double Hst = H[t * S + s];
      if (s != t && Hst > 0) nl.coef[a] = lb / Hst;
    }
  }

  Graph g;
  g.n = n;
  g.self.assign(n, 0.0);
  g.u.assign(static_cast<size_t>(n) * D, 0.0);
  g.p.assign(n + 1, 0);
  for (int i = 0; i < n; ++i) {
    const int s = sp[i];
    for (int e = co_p[i]; e < co_p[i + 1]; ++e) {
      if (co_i[e] == i) {
        g.self[i] += co_x[e] / two_m[s];
        continue;
      }
      g.j.push_back(co_i[e]);
      g.w.push_back(co_x[e] / two_m[s]);
    }
    if (lambda != 0.0) {
      for (int e = hog_p[i]; e < hog_p[i + 1]; ++e) {
        g.j.push_back(hog_i[e]);
        g.w.push_back(lb * hog_x[e]);
      }
    }
    g.p[i + 1] = static_cast<int>(g.j.size());
    double* ui = &g.u[static_cast<size_t>(i) * D];
    ui[s] = k[i];
    for (int t = 0; t < S; ++t) {
      ui[S + s * S + t] = h[static_cast<size_t>(t) * n + i];
    }
  }

  // orig: gene -> node of the base graph (genes, or init supernodes)
  std::vector<int> orig(n);
  for (int i = 0; i < n; ++i) orig[i] = freeze ? init[i] : i;
  if (freeze) g = aggregate(g, orig, compact(orig), D);

  // Coarsen: local moves, aggregate, repeat until a level makes no move.
  // part: base node -> community.
  std::vector<int> part(g.n);
  for (int i = 0; i < g.n; ++i) part[i] = i;
  Graph cur;
  int levels = 0;
  for (; levels < 1000; ++levels) {
    const Graph& gl = levels == 0 ? g : cur;
    std::vector<int> comm(gl.n);
    for (int i = 0; i < gl.n; ++i) comm[i] = i;
    if (local_moves(gl, nl, comm, max_passes, tol) == 0) break;
    int K = compact(comm);
    for (int& c : part) c = comm[c];
    Graph next = aggregate(gl, comm, K, D);
    cur = std::move(next);
  }
  // Refine: one round of local moves on the base graph from the result.
  // Plain Louvain cannot take back a node once its supernode merged.
  int refined = local_moves(g, nl, part, max_passes, tol);
  compact(part);
  for (int& c : orig) c = part[c];
  return List::create(_["labels"] = wrap(orig), _["levels"] = levels,
                      _["refine_moves"] = refined);
}
