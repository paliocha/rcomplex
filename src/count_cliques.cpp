// count_cliques.cpp
// Count the maximal cliques of a graph, with a ceiling. Bron-Kerbosch
// with the Tomita pivot; the count stops as soon as it passes the
// ceiling, so a dense hog costs no more than the budget it breaks.
//
// Integer indices and sorted vectors only. Single-threaded.

#include <Rcpp.h>
#include <algorithm>
#include <iterator>
#include <vector>

namespace {

struct Counter {
  std::vector<std::vector<int>> adj;
  int min_size;
  double ceiling;
  double count = 0;

  std::vector<int> meet(const std::vector<int>& s, int v) const {
    std::vector<int> out;
    const std::vector<int>& nb = adj[v];
    std::set_intersection(s.begin(), s.end(), nb.begin(), nb.end(),
                          std::back_inserter(out));
    return out;
  }

  std::size_t n_meet(const std::vector<int>& s, int v) const {
    const std::vector<int>& nb = adj[v];
    std::size_t n = 0;
    auto i = s.begin();
    auto j = nb.begin();
    while (i != s.end() && j != nb.end()) {
      if (*i < *j) {
        ++i;
      } else if (*j < *i) {
        ++j;
      } else {
        ++n;
        ++i;
        ++j;
      }
    }
    return n;
  }

  // Returns false once the count passes the ceiling.
  bool expand(int r, std::vector<int> p, std::vector<int> x) {
    if (p.empty()) {
      if (x.empty() && r >= min_size) count += 1;
      return count <= ceiling;
    }
    // No clique through r can reach min_size.
    if (r + static_cast<int>(p.size()) < min_size) return true;
    // Pivot: the vertex of p or x with the most neighbours in p.
    int pivot = -1;
    std::size_t best = 0;
    for (const std::vector<int>* s : {&p, &x}) {
      for (int u : *s) {
        std::size_t n = n_meet(p, u);
        if (pivot < 0 || n > best) {
          pivot = u;
          best = n;
        }
      }
    }
    std::vector<int> cand;
    const std::vector<int>& np = adj[pivot];
    std::set_difference(p.begin(), p.end(), np.begin(), np.end(),
                        std::back_inserter(cand));
    for (int v : cand) {
      if (!expand(r + 1, meet(p, v), meet(x, v))) return false;
      p.erase(std::lower_bound(p.begin(), p.end(), v));
      x.insert(std::lower_bound(x.begin(), x.end(), v), v);
    }
    return true;
  }
};

}  // namespace

//' Count maximal cliques of at least min_size vertices, up to a ceiling
//'
//' @param i1,i2 Edge end points, 1-based vertex indices.
//' @param n_v Number of vertices.
//' @param min_size Smallest clique counted.
//' @param ceiling The count stops once it exceeds this value.
//' @return The count; above `ceiling` when counting stopped early.
//' @noRd
// [[Rcpp::export]]
double count_cliques_cpp(Rcpp::IntegerVector i1, Rcpp::IntegerVector i2,
                         int n_v, int min_size, double ceiling) {
  Counter c;
  c.adj.assign(n_v, std::vector<int>());
  for (R_xlen_t k = 0; k < i1.size(); ++k) {
    int a = i1[k] - 1;
    int b = i2[k] - 1;
    if (a == b) continue;
    c.adj[a].push_back(b);
    c.adj[b].push_back(a);
  }
  for (std::vector<int>& nb : c.adj) {
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
  }
  c.min_size = min_size;
  c.ceiling = ceiling;
  std::vector<int> p(n_v);
  for (int v = 0; v < n_v; ++v) p[v] = v;
  c.expand(0, p, std::vector<int>());
  return c.count;
}
