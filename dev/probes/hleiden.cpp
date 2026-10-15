// Probe WP15b: serial Leiden on a species-blocked objective with a strict
// hyperedge term. Written from the papers only: the Leiden phases from
// Traag, Waltman & van Eck 2019 (fast local move with a queue,
// refinement of each community from singletons with randomised merges
// among well-connected subsets, aggregation on the refined partition);
// the alpha blend and its schedule from Kaminski, Pralat & Theberge 2024
// (h-Louvain). No library source was consulted.
//
//   q(a) = sum_s Q_s + l ((1 - a) P + a T)
//
// Q_s: per-species modularity, as E2 (joint_louvain.cpp).
// P:   E2's hog-pair term, (1 / 2m_bar) sum (w_ij - h_i^t h_j^s / H_st).
// T:   strict hyperedge term. For a hog h with species set S_h (two or
//      more species), c_hs copies in species s and n_hs(C) of them in
//      module C:
//        T = (1 / 2m_bar) sum_C sum_h w_h
//              (prod_{s in S_h} n_hs(C) / c_hs - prod_{s in S_h} p_s(C))
//      with w_h = |S_h| (|S_h| - 1), the credit E2 gives a hog whose
//      copies all share a module. Paralogs are collapsed: species s
//      counts as the fraction n_hs(C) / c_hs, so T credits w_h in full
//      only when every copy sits in C. p_s(C) = H_s(C) / H_s, where a
//      copy carries hog weight 1 / c_hs, H_s(C) sums it over species-s
//      genes in C and H_s over all of species s (species-blocked null).
//
// Node gains: pairwise part 2 (e(P, C) - null(P, C)) as E2; strict part
// a l / 2m_bar (change in observed - change in null), both exact.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <deque>
#include <numeric>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

typedef std::vector<int> ivec;
typedef std::vector<double> dvec;
const double EPS = 1e-14;

struct Csr {
  ivec p, j;
  dvec w;
};

struct Graph {
  int n;
  Csr co, hg;   // co-expression (A / 2m_s); hog pairs (l / 2m_bar) w_ij
  dvec u;       // n x D pairwise null profile (E2 layout)
  dvec hc;      // n x S strict profile: H_s(node) / H_s
  dvec self;    // co-expression self weight, base graph only
};

struct Model {
  int S, D, n0;
  double lb;               // lambda / 2m_bar
  ivec partner;
  dvec coef_co, coef_hog;  // pairwise null coefficients; hog part at a = 0
  ivec mask;               // species sets present among hogs
  dvec wmask;              // total w_h per species set
  ivec hsp, hss;           // hog -> its species (CSR)
  dvec hsc;                // copies of that species
  ivec hgp, hgg, hgl;      // hog -> genes and their local species index
  dvec hw;                 // w_h
  ivec ghp, ghh;           // gene -> hogs (CSR)
};

// Community state: profile sums and a cache of d T_null / d p_s
struct Comm {
  int S, D;
  dvec U, H, dn;
  ivec size, ver, dver;
  void reset(int n, int S_, int D_) {
    S = S_;
    D = D_;
    U.assign(static_cast<size_t>(n) * D, 0.0);
    H.assign(static_cast<size_t>(n) * S, 0.0);
    dn.assign(static_cast<size_t>(n) * S, 0.0);
    size.assign(n, 0);
    ver.assign(n, 0);
    dver.assign(static_cast<size_t>(n) * S, -1);
  }
  void add(int c, const double* u, const double* h, int sgn) {
    double* Uc = &U[static_cast<size_t>(c) * D];
    for (int a = 0; a < D; ++a) Uc[a] += sgn * u[a];
    double* Hc = &H[static_cast<size_t>(c) * S];
    for (int s = 0; s < S; ++s) Hc[s] += sgn * h[s];
    size[c] += sgn;
    ver[c]++;
  }
};

struct Work {
  dvec nw, nwc, sobs;
  std::vector<char> seen;
  ivec touched, hstamp;
  int stamp;
  std::vector<std::pair<int, int>> oth;
  Work(int n, int nh)
      : nw(n, 0.0), nwc(n, 0.0), sobs(n, 0.0), seen(n, 0),
        hstamp(nh, 0), stamp(0) {}
  void touch(int c) {
    if (!seen[c]) {
      seen[c] = 1;
      touched.push_back(c);
    }
  }
  void clear() {
    for (int c : touched) {
      nw[c] = nwc[c] = sobs[c] = 0.0;
      seen[c] = 0;
    }
    touched.clear();
  }
};

void shuffle(ivec& x) {
  for (int i = static_cast<int>(x.size()) - 1; i > 0; --i) {
    int k = static_cast<int>(R::unif_rand() * (i + 1));
    if (k > i) k = i;
    std::swap(x[i], x[k]);
  }
}

// Relabel to 0..K-1 by first appearance; labels must be < x.size()
int compact(ivec& x) {
  ivec map(x.size(), -1);
  int k = 0;
  for (int& c : x) {
    if (map[c] < 0) map[c] = k++;
    c = map[c];
  }
  return k;
}

// CSR of members: group g -> items with lab[item] == g
void members(const ivec& lab, int K, ivec& mp, ivec& mg) {
  mp.assign(K + 1, 0);
  for (int c : lab) mp[c + 1]++;
  for (int k = 0; k < K; ++k) mp[k + 1] += mp[k];
  mg.assign(lab.size(), 0);
  ivec pos(mp.begin(), mp.end() - 1);
  for (int i = 0; i < static_cast<int>(lab.size()); ++i) {
    mg[pos[lab[i]]++] = i;
  }
}

void build(const IntegerVector& sp, int S, const IntegerVector& co_p,
           const IntegerVector& co_i, const NumericVector& co_x,
           const IntegerVector& hog_p, const IntegerVector& hog_i,
           const NumericVector& hog_x, const NumericVector& h,
           const NumericVector& H, const IntegerVector& ln,
           const IntegerVector& lh, double lambda, double gamma,
           Model& M, Graph& g) {
  if (S > 20) stop("hleiden: at most 20 species");
  const int n = sp.size(), D = S + S * S;
  M.S = S;
  M.D = D;
  M.n0 = n;
  dvec two_m(S, 0.0), k(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int e = co_p[i]; e < co_p[i + 1]; ++e) k[i] += co_x[e];
    two_m[sp[i]] += k[i];
  }
  double bar = 0.0;
  for (int s = 0; s < S; ++s) bar += two_m[s] / S;
  M.lb = lambda / bar;
  M.partner.assign(D, 0);
  M.coef_co.assign(D, 0.0);
  M.coef_hog.assign(D, 0.0);
  for (int s = 0; s < S; ++s) {
    M.partner[s] = s;
    M.coef_co[s] = gamma / (two_m[s] * two_m[s]);
    for (int t = 0; t < S; ++t) {
      int a = S + s * S + t;
      M.partner[a] = S + t * S + s;
      double Hst = H[t * S + s];
      if (s != t && Hst > 0) M.coef_hog[a] = M.lb / Hst;
    }
  }
  g.n = n;
  g.self.assign(n, 0.0);
  g.u.assign(static_cast<size_t>(n) * D, 0.0);
  g.co.p.assign(n + 1, 0);
  g.hg.p.assign(n + 1, 0);
  for (int i = 0; i < n; ++i) {
    const int s = sp[i];
    for (int e = co_p[i]; e < co_p[i + 1]; ++e) {
      if (co_i[e] == i) {
        g.self[i] += co_x[e] / two_m[s];
        continue;
      }
      g.co.j.push_back(co_i[e]);
      g.co.w.push_back(co_x[e] / two_m[s]);
    }
    g.co.p[i + 1] = static_cast<int>(g.co.j.size());
    for (int e = hog_p[i]; e < hog_p[i + 1]; ++e) {
      g.hg.j.push_back(hog_i[e]);
      g.hg.w.push_back(M.lb * hog_x[e]);
    }
    g.hg.p[i + 1] = static_cast<int>(g.hg.j.size());
    double* ui = &g.u[static_cast<size_t>(i) * D];
    ui[s] = k[i];
    for (int t = 0; t < S; ++t) {
      ui[S + s * S + t] = h[static_cast<size_t>(t) * n + i];
    }
  }

  // Strict term: multi-species hogs from the (node, hog) link table
  const int nl = ln.size();
  ivec ord(nl);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b) {
    return lh[a] != lh[b] ? lh[a] < lh[b] : ln[a] < ln[b];
  });
  dvec dg(n, 0.0), Hs(S, 0.0), wm(static_cast<size_t>(1) << S, 0.0);
  ivec cnt(S), loc(S);
  M.hsp.assign(1, 0);
  M.hgp.assign(1, 0);
  for (int a = 0; a < nl;) {
    int b = a;
    while (b < nl && lh[ord[b]] == lh[ord[a]]) ++b;
    std::fill(cnt.begin(), cnt.end(), 0);
    for (int x = a; x < b; ++x) cnt[sp[ln[ord[x]]]]++;
    int ks = 0, mask = 0;
    for (int s = 0; s < S; ++s) {
      if (cnt[s]) {
        ++ks;
        mask |= 1 << s;
      }
    }
    if (ks >= 2) {
      const double w = ks * (ks - 1.0);
      M.hw.push_back(w);
      wm[mask] += w;
      for (int s = 0; s < S; ++s) {
        if (!cnt[s]) continue;
        loc[s] = static_cast<int>(M.hss.size()) - M.hsp.back();
        M.hss.push_back(s);
        M.hsc.push_back(cnt[s]);
      }
      M.hsp.push_back(static_cast<int>(M.hss.size()));
      for (int x = a; x < b; ++x) {
        const int gene = ln[ord[x]], s = sp[gene];
        M.hgg.push_back(gene);
        M.hgl.push_back(loc[s]);
        dg[gene] += 1.0 / cnt[s];
        Hs[s] += 1.0 / cnt[s];
      }
      M.hgp.push_back(static_cast<int>(M.hgg.size()));
    }
    a = b;
  }
  for (int m = 0; m < (1 << S); ++m) {
    if (wm[m] > 0) {
      M.mask.push_back(m);
      M.wmask.push_back(wm[m]);
    }
  }
  const int nh = static_cast<int>(M.hw.size());
  M.ghp.assign(n + 1, 0);
  for (int x = 0; x < static_cast<int>(M.hgg.size()); ++x) {
    M.ghp[M.hgg[x] + 1]++;
  }
  for (int i = 0; i < n; ++i) M.ghp[i + 1] += M.ghp[i];
  M.ghh.assign(M.hgg.size(), 0);
  ivec pos(M.ghp.begin(), M.ghp.end() - 1);
  for (int hh = 0; hh < nh; ++hh) {
    for (int x = M.hgp[hh]; x < M.hgp[hh + 1]; ++x) {
      M.ghh[pos[M.hgg[x]]++] = hh;
    }
  }
  g.hc.assign(static_cast<size_t>(n) * S, 0.0);
  for (int i = 0; i < n; ++i) {
    const int s = sp[i];
    if (Hs[s] > 0) g.hc[static_cast<size_t>(i) * S + s] = dg[i] / Hs[s];
  }
}

// d T_null / d p_s for community c (cached until c changes)
double dn_of(const Model& M, Comm& C, int c, int s) {
  const size_t idx = static_cast<size_t>(c) * M.S + s;
  if (C.dver[idx] == C.ver[c]) return C.dn[idx];
  const double* Hc = &C.H[static_cast<size_t>(c) * M.S];
  double sum = 0.0;
  for (size_t t = 0; t < M.mask.size(); ++t) {
    const int m = M.mask[t];
    if (!(m >> s & 1)) continue;
    double pr = M.wmask[t];
    for (int v = 0; v < M.S && pr != 0.0; ++v) {
      if (v != s && (m >> v & 1)) pr *= Hc[v];
    }
    sum += pr;
  }
  C.dn[idx] = sum;
  C.dver[idx] = C.ver[c];
  return sum;
}

// -1: no strict weight; s: one species; -2: several
int species_of(const double* p, int S) {
  int s = -1;
  for (int v = 0; v < S; ++v) {
    if (p[v] != 0.0) {
      if (s >= 0) return -2;
      s = v;
    }
  }
  return s;
}

// Increase of the strict null when p joins community c (c < 0: empty)
double snull(const Model& M, Comm& C, int c, const double* p, int single) {
  if (single == -1) return 0.0;
  if (single >= 0) return c < 0 ? 0.0 : p[single] * dn_of(M, C, c, single);
  double sum = 0.0;
  for (size_t t = 0; t < M.mask.size(); ++t) {
    const int m = M.mask[t];
    double a = M.wmask[t], b = M.wmask[t];
    for (int v = 0; v < M.S; ++v) {
      if (!(m >> v & 1)) continue;
      const double x = c < 0 ? 0.0 : C.H[static_cast<size_t>(c) * M.S + v];
      a *= x + p[v];
      b *= x;
    }
    sum += a - b;
  }
  return sum;
}

// Change in the strict observed term when node v joins each community
// holding other copies of its hogs (lab: node -> community). Adds
// w_h (val_X - base_h) to sobs[X]; returns sum_h w_h base_h, the change
// for a community with no other copy.
double strict_obs(const Model& M, const ivec& mp, const ivec& mg, int v,
                  const ivec& node_of, const ivec& lab, Work& W) {
  double cst = 0.0, p[32], cn[32];
  ++W.stamp;
  for (int e = mp[v]; e < mp[v + 1]; ++e) {
    const int gene = mg[e];
    for (int f = M.ghp[gene]; f < M.ghp[gene + 1]; ++f) {
      const int h = M.ghh[f];
      if (W.hstamp[h] == W.stamp) continue;
      W.hstamp[h] = W.stamp;
      const int k = M.hsp[h + 1] - M.hsp[h];
      const double* c = &M.hsc[M.hsp[h]];
      std::fill(p, p + k, 0.0);
      W.oth.clear();
      for (int x = M.hgp[h]; x < M.hgp[h + 1]; ++x) {
        const int nd = node_of[M.hgg[x]];
        if (nd == v) {
          p[M.hgl[x]] += 1.0;
        } else {
          W.oth.emplace_back(lab[nd], M.hgl[x]);
        }
      }
      double base = 1.0;
      for (int l = 0; l < k; ++l) base *= p[l] / c[l];
      cst += M.hw[h] * base;
      std::sort(W.oth.begin(), W.oth.end());
      for (size_t a = 0; a < W.oth.size();) {
        size_t b = a;
        std::fill(cn, cn + k, 0.0);
        while (b < W.oth.size() && W.oth[b].first == W.oth[a].first) {
          cn[W.oth[b].second] += 1.0;
          ++b;
        }
        double with = 1.0, without = 1.0;
        for (int l = 0; l < k; ++l) {
          with *= (cn[l] + p[l]) / c[l];
          without *= cn[l] / c[l];
        }
        const double val = with - without - base;
        if (val != 0.0) {
          const int X = W.oth[a].first;
          W.touch(X);
          W.sobs[X] += M.hw[h] * val;
        }
        a = b;
      }
    }
  }
  return cst;
}

// Everything the gain of node v needs, for one partition lab
struct Eval {
  const Model& M;
  Comm& C;
  Work& W;
  double a, cst;
  const double* p;
  int single;
  ivec vidx;
  dvec vval;
  Eval(const Model& M_, Comm& C_, Work& W_, double a_)
      : M(M_), C(C_), W(W_), a(a_), cst(0.0), p(nullptr), single(-1) {}

  // Fill W for node v; only neighbours with keep_lab[u] == keep count
  // (keep < 0: all)
  void scan(const Graph& g, const ivec& mp, const ivec& mg,
            const ivec& node_of, const ivec& lab, const dvec& coef, int v,
            const ivec* keep_lab, int keep) {
    const double* uv = &g.u[static_cast<size_t>(v) * M.D];
    vidx.clear();
    vval.clear();
    for (int x = 0; x < M.D; ++x) {
      if (uv[x] != 0.0 && coef[x] != 0.0) {
        vidx.push_back(M.partner[x]);
        vval.push_back(uv[x] * coef[x]);
      }
    }
    for (int e = g.co.p[v]; e < g.co.p[v + 1]; ++e) {
      const int u = g.co.j[e];
      if (keep >= 0 && (*keep_lab)[u] != keep) continue;
      const int c = lab[u];
      W.touch(c);
      W.nw[c] += g.co.w[e];
      W.nwc[c] += g.co.w[e];
    }
    for (int e = g.hg.p[v]; e < g.hg.p[v + 1]; ++e) {
      const int u = g.hg.j[e];
      if (keep >= 0 && (*keep_lab)[u] != keep) continue;
      const int c = lab[u];
      W.touch(c);
      W.nw[c] += (1.0 - a) * g.hg.w[e];
      W.nwc[c] += g.hg.w[e];
    }
    p = &g.hc[static_cast<size_t>(v) * M.S];
    single = species_of(p, M.S);
    cst = a > 0.0 ? strict_obs(M, mp, mg, v, node_of, lab, W) : 0.0;
  }

  double gain(int c) {
    double e = 0.0, nu = 0.0, st = 0.0;
    if (c >= 0) {
      e = W.nw[c];
      const double* Uc = &C.U[static_cast<size_t>(c) * M.D];
      for (size_t x = 0; x < vidx.size(); ++x) nu += vval[x] * Uc[vidx[x]];
    }
    if (a > 0.0) {
      st = a * M.lb *
           (cst + (c >= 0 ? W.sobs[c] : 0.0) - snull(M, C, c, p, single));
    }
    return 2.0 * (e - nu) + st;
  }
};

dvec coef_at(const Model& M, double a) {
  dvec c(M.D);
  for (int x = 0; x < M.D; ++x) {
    c[x] = x < M.S ? M.coef_co[x] : (1.0 - a) * M.coef_hog[x];
  }
  return c;
}

void enqueue(int u, std::deque<int>& q, std::vector<char>& inq) {
  if (!inq[u]) {
    inq[u] = 1;
    q.push_back(u);
  }
}

// Fast local move (Traag 2019, MoveNodesFast). Returns moves made.
int move_nodes(const Model& M, const Graph& g, const ivec& mp,
               const ivec& mg, const ivec& node_of, ivec& comm, double a,
               Work& W, Comm& C) {
  const int n = g.n, D = M.D, S = M.S;
  C.reset(n, S, D);
  for (int i = 0; i < n; ++i) {
    C.add(comm[i], &g.u[static_cast<size_t>(i) * D],
          &g.hc[static_cast<size_t>(i) * S], 1);
  }
  ivec empty;
  for (int c = n - 1; c >= 0; --c) {
    if (C.size[c] == 0) empty.push_back(c);
  }
  const dvec coef = coef_at(M, a);
  ivec order(n);
  std::iota(order.begin(), order.end(), 0);
  shuffle(order);
  std::deque<int> q(order.begin(), order.end());
  std::vector<char> inq(n, 1);
  Eval ev(M, C, W, a);
  int moves = 0;
  while (!q.empty()) {
    const int v = q.front();
    q.pop_front();
    inq[v] = 0;
    const int d = comm[v];
    const double* uv = &g.u[static_cast<size_t>(v) * D];
    const double* hv = &g.hc[static_cast<size_t>(v) * S];
    C.add(d, uv, hv, -1);
    ev.scan(g, mp, mg, node_of, comm, coef, v, nullptr, -1);
    const double ge = ev.gain(-1);
    const bool alone = C.size[d] == 0;
    int best = d;
    double bg = alone ? ge : ev.gain(d);
    for (int c : W.touched) {
      if (c == d) continue;
      const double gc = ev.gain(c);
      if (gc > bg + EPS) {
        bg = gc;
        best = c;
      }
    }
    // A node with no 2-section edge into the rest of its community leaves
    // for an empty one: the gain there is at least that of staying.
    if (!alone && !empty.empty() &&
        (ge > bg + EPS || (best == d && W.nwc[d] == 0.0))) {
      best = empty.back();
      bg = ge;
    }
    W.clear();
    if (best != d) {
      if (C.size[best] == 0) empty.pop_back();
      if (alone) empty.push_back(d);
      ++moves;
      for (int e = g.co.p[v]; e < g.co.p[v + 1]; ++e) {
        if (comm[g.co.j[e]] != best) enqueue(g.co.j[e], q, inq);
      }
      for (int e = g.hg.p[v]; e < g.hg.p[v + 1]; ++e) {
        if (comm[g.hg.j[e]] != best) enqueue(g.hg.j[e], q, inq);
      }
      if (a > 0.0) {  // paralogs share no edge but share the strict term
        for (int e = mp[v]; e < mp[v + 1]; ++e) {
          const int gene = mg[e];
          for (int f = M.ghp[gene]; f < M.ghp[gene + 1]; ++f) {
            const int h = M.ghh[f];
            for (int x = M.hgp[h]; x < M.hgp[h + 1]; ++x) {
              const int u = node_of[M.hgg[x]];
              if (u != v && comm[u] != best) enqueue(u, q, inq);
            }
          }
        }
      }
    }
    comm[v] = best;
    C.add(best, uv, hv, 1);
  }
  return moves;
}

// Refinement (Traag 2019, RefinePartition / MergeNodesSubset). ref gets
// refined labels (node ids); returns the number of merges. Well-connected
// is judged on the 2-section at a = 0: e(X, S - X) >= null(X, S - X).
// Merge gains are judged on q(0), not q(a): the strict term gives a lone
// copy nothing until its whole hog shares a module, so refinement on q(a)
// cannot grow a hog one node at a time; it stalls, and a complete hog
// with no edge to its module stays there (check 2 failed on the planted
// schedule). Local moves and aggregation still use q(a). The parameter a
// is kept for the call site; refinement ignores it.
int refine(const Model& M, const Graph& g, const ivec& mp, const ivec& mg,
           const ivec& node_of, const ivec& comm, ivec& ref, double a,
           double theta, Work& W, Comm& R) {
  const int n = g.n, D = M.D, S = M.S;
  R.reset(n, S, D);
  for (int i = 0; i < n; ++i) {
    ref[i] = i;
    R.add(i, &g.u[static_cast<size_t>(i) * D],
          &g.hc[static_cast<size_t>(i) * S], 1);
  }
  int K = 0;
  for (int c : comm) K = std::max(K, c + 1);
  ivec cp, cm;
  members(comm, K, cp, cm);
  dvec US(static_cast<size_t>(K) * D, 0.0), ext(n, 0.0);
  for (int i = 0; i < n; ++i) {
    const double* ui = &g.u[static_cast<size_t>(i) * D];
    double* Uc = &US[static_cast<size_t>(comm[i]) * D];
    for (int x = 0; x < D; ++x) Uc[x] += ui[x];
    for (int e = g.co.p[i]; e < g.co.p[i + 1]; ++e) {
      if (comm[g.co.j[e]] == comm[i]) ext[i] += g.co.w[e];
    }
    for (int e = g.hg.p[i]; e < g.hg.p[i + 1]; ++e) {
      if (comm[g.hg.j[e]] == comm[i]) ext[i] += g.hg.w[e];
    }
  }
  dvec extR(ext);
  const dvec coef0 = coef_at(M, 0.0);
  (void)a;
  auto nullconn = [&](const double* uX, const double* uS) {
    double s = 0.0;
    for (int x = 0; x < D; ++x) {
      if (uX[x] == 0.0 || coef0[x] == 0.0) continue;
      const int y = M.partner[x];
      s += coef0[x] * uX[x] * (uS[y] - uX[y]);
    }
    return s;
  };
  Eval ev(M, R, W, 0.0);
  struct Cand {
    int r;
    double dh, ec;
  };
  std::vector<Cand> cand;
  int merges = 0;
  for (int c = 0; c < K; ++c) {
    if (cp[c + 1] - cp[c] < 2) continue;
    ivec mem(cm.begin() + cp[c], cm.begin() + cp[c + 1]);
    shuffle(mem);
    const double* uS = &US[static_cast<size_t>(c) * D];
    for (int v : mem) {
      if (R.size[ref[v]] != 1) continue;
      const double* uv = &g.u[static_cast<size_t>(v) * D];
      if (ext[v] < nullconn(uv, uS)) continue;
      ev.scan(g, mp, mg, node_of, ref, coef0, v, &comm, c);
      const double ge = ev.gain(-1);
      cand.clear();
      cand.push_back({v, 0.0, 0.0});
      double top = 0.0;
      for (int r : W.touched) {
        if (r == v || comm[r] != c || W.nwc[r] <= 0.0) continue;
        if (extR[r] < nullconn(&R.U[static_cast<size_t>(r) * D], uS)) {
          continue;
        }
        const double dh = ev.gain(r) - ge;
        if (dh < 0.0) continue;
        cand.push_back({r, dh, W.nwc[r]});
        top = std::max(top, dh);
      }
      W.clear();
      if (cand.size() == 1) continue;
      double tot = 0.0;
      for (Cand& x : cand) tot += std::exp((x.dh - top) / theta);
      double draw = R::unif_rand() * tot;
      size_t pick = 0;
      for (; pick + 1 < cand.size(); ++pick) {
        draw -= std::exp((cand[pick].dh - top) / theta);
        if (draw < 0.0) break;
      }
      const Cand& ch = cand[pick];
      if (ch.r == v) continue;
      const double* hv = &g.hc[static_cast<size_t>(v) * S];
      extR[ch.r] += ext[v] - 2.0 * ch.ec;
      R.add(v, uv, hv, -1);
      R.add(ch.r, uv, hv, 1);
      ref[v] = ch.r;
      ++merges;
    }
  }
  return merges;
}

Csr agg_csr(const Csr& in, const ivec& map, const ivec& ap,
            const ivec& am, int K, Work& W) {
  Csr out;
  out.p.assign(K + 1, 0);
  for (int a = 0; a < K; ++a) {
    for (int x = ap[a]; x < ap[a + 1]; ++x) {
      const int i = am[x];
      for (int e = in.p[i]; e < in.p[i + 1]; ++e) {
        const int b = map[in.j[e]];
        if (b == a) continue;
        W.touch(b);
        W.nw[b] += in.w[e];
      }
    }
    for (int b : W.touched) {
      out.j.push_back(b);
      out.w.push_back(W.nw[b]);
    }
    W.clear();
    out.p[a + 1] = static_cast<int>(out.j.size());
  }
  return out;
}

Graph aggregate(const Graph& g, const ivec& map, int K, const Model& M,
                Work& W) {
  ivec ap, am;
  members(map, K, ap, am);
  Graph h;
  h.n = K;
  h.u.assign(static_cast<size_t>(K) * M.D, 0.0);
  h.hc.assign(static_cast<size_t>(K) * M.S, 0.0);
  for (int i = 0; i < g.n; ++i) {
    for (int x = 0; x < M.D; ++x) {
      h.u[static_cast<size_t>(map[i]) * M.D + x] +=
        g.u[static_cast<size_t>(i) * M.D + x];
    }
    for (int s = 0; s < M.S; ++s) {
      h.hc[static_cast<size_t>(map[i]) * M.S + s] +=
        g.hc[static_cast<size_t>(i) * M.S + s];
    }
  }
  h.co = agg_csr(g.co, map, ap, am, K, W);
  h.hg = agg_csr(g.hg, map, ap, am, K, W);
  return h;
}

}  // namespace

// Quality parts per unit lambda at the gene level, from the kernel's own
// structures: co (sum of Q_s), hog_raw (P) and strict_raw (T).
// [[Rcpp::export]]
List hl_quality_cpp(IntegerVector lab, IntegerVector sp, int S,
                    IntegerVector co_p, IntegerVector co_i,
                    NumericVector co_x, IntegerVector hog_p,
                    IntegerVector hog_i, NumericVector hog_x,
                    NumericVector h, NumericVector H, IntegerVector ln,
                    IntegerVector lh, double gamma) {
  Model M;
  Graph g;
  build(sp, S, co_p, co_i, co_x, hog_p, hog_i, hog_x, h, H, ln, lh, 1.0,
        gamma, M, g);
  const int n = g.n, D = M.D;
  int K = 0;
  for (int i = 0; i < n; ++i) K = std::max(K, lab[i] + 1);
  dvec U(static_cast<size_t>(K) * D, 0.0), P(static_cast<size_t>(K) * S, 0.0);
  double co_in = 0.0, hog_in = 0.0;
  for (int i = 0; i < n; ++i) {
    const int c = lab[i];
    co_in += g.self[i];
    for (int e = g.co.p[i]; e < g.co.p[i + 1]; ++e) {
      if (lab[g.co.j[e]] == c) co_in += g.co.w[e];
    }
    for (int e = g.hg.p[i]; e < g.hg.p[i + 1]; ++e) {
      if (lab[g.hg.j[e]] == c) hog_in += g.hg.w[e];
    }
    for (int x = 0; x < D; ++x) {
      U[static_cast<size_t>(c) * D + x] += g.u[static_cast<size_t>(i) * D + x];
    }
    for (int s = 0; s < S; ++s) {
      P[static_cast<size_t>(c) * S + s] += g.hc[static_cast<size_t>(i) * S + s];
    }
  }
  double co_nul = 0.0, hog_nul = 0.0, s_nul = 0.0, s_obs = 0.0;
  for (int c = 0; c < K; ++c) {
    const double* Uc = &U[static_cast<size_t>(c) * D];
    for (int x = 0; x < D; ++x) {
      const double v = Uc[x] * Uc[M.partner[x]];
      if (x < S) {
        co_nul += M.coef_co[x] * v;
      } else {
        hog_nul += M.coef_hog[x] * v;
      }
    }
    const double* Pc = &P[static_cast<size_t>(c) * S];
    for (size_t t = 0; t < M.mask.size(); ++t) {
      double pr = M.wmask[t];
      for (int s = 0; s < S; ++s) {
        if (M.mask[t] >> s & 1) pr *= Pc[s];
      }
      s_nul += pr;
    }
  }
  std::vector<std::pair<int, int>> oth;
  double cn[32];
  for (int hh = 0; hh < static_cast<int>(M.hw.size()); ++hh) {
    const int k = M.hsp[hh + 1] - M.hsp[hh];
    const double* c = &M.hsc[M.hsp[hh]];
    oth.clear();
    for (int x = M.hgp[hh]; x < M.hgp[hh + 1]; ++x) {
      oth.emplace_back(lab[M.hgg[x]], M.hgl[x]);
    }
    std::sort(oth.begin(), oth.end());
    for (size_t a = 0; a < oth.size();) {
      size_t b = a;
      std::fill(cn, cn + k, 0.0);
      while (b < oth.size() && oth[b].first == oth[a].first) {
        cn[oth[b].second] += 1.0;
        ++b;
      }
      double pr = M.hw[hh];
      for (int l = 0; l < k; ++l) pr *= cn[l] / c[l];
      s_obs += pr;
      a = b;
    }
  }
  return List::create(_["co"] = co_in - co_nul,
                      _["hog_raw"] = hog_in - hog_nul,
                      _["strict_raw"] = M.lb * (s_obs - s_nul));
}

// Leiden on q(alpha). alpha[l] is used at aggregation level l of the first
// iteration (last value held); later iterations use the last value.
// Iterates until an iteration leaves the partition unchanged.
// [[Rcpp::export]]
List hl_leiden_cpp(IntegerVector sp, int S, IntegerVector co_p,
                   IntegerVector co_i, NumericVector co_x,
                   IntegerVector hog_p, IntegerVector hog_i,
                   NumericVector hog_x, NumericVector h, NumericVector H,
                   IntegerVector ln, IntegerVector lh, double lambda,
                   double gamma, NumericVector alpha, double theta,
                   int max_iter) {
  RNGScope scope;
  Model M;
  Graph g0;
  build(sp, S, co_p, co_i, co_x, hog_p, hog_i, hog_x, h, H, ln, lh, lambda,
        gamma, M, g0);
  const int n0 = g0.n, L = alpha.size();
  Work W(n0, static_cast<int>(M.hw.size()));
  Comm C, R;
  ivec part(n0);
  std::iota(part.begin(), part.end(), 0);
  ivec levels;
  int stalls = 0, it = 0, moves0 = 0;
  bool converged = false;
  for (; it < max_iter && !converged; ++it) {
    ivec node_of(n0), mp, mg;
    std::iota(node_of.begin(), node_of.end(), 0);
    members(node_of, n0, mp, mg);
    ivec comm(part);
    Graph cur;
    int level = 0;
    for (;; ++level) {
      const Graph& gl = level == 0 ? g0 : cur;
      const double a = alpha[it == 0 ? std::min(level, L - 1) : L - 1];
      const int mv = move_nodes(M, gl, mp, mg, node_of, comm, a, W, C);
      if (level == 0) moves0 = mv;
      if (compact(comm) == gl.n) break;
      ivec ref(gl.n);
      if (refine(M, gl, mp, mg, node_of, comm, ref, a, theta, W, R) == 0) {
        ++stalls;
        break;
      }
      const int Kr = compact(ref);
      ivec next(Kr);
      for (int i = 0; i < gl.n; ++i) next[ref[i]] = comm[i];
      Graph ng = aggregate(gl, ref, Kr, M, W);
      for (int& x : node_of) x = ref[x];
      members(node_of, Kr, mp, mg);
      cur = std::move(ng);
      comm.swap(next);
    }
    ivec lab(n0);
    for (int i = 0; i < n0; ++i) lab[i] = comm[node_of[i]];
    compact(lab);
    converged = lab == part;
    part.swap(lab);
    levels.push_back(level + 1);
  }
  return List::create(_["labels"] = wrap(part), _["iterations"] = it,
                      _["converged"] = converged, _["levels"] = wrap(levels),
                      _["stalls"] = stalls, _["last_moves"] = moves0);
}
