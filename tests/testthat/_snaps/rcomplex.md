# print() shows species, edges and tiers

    Code
      print(res)
    Output
      rcomplex: 3 species, 120 hogs, 3 tiers
        SpA  120 genes  20 samples  density 0.03  r >= 0.26
        SpB  120 genes  20 samples  density 0.03  r >= 0.31
        SpC  120 genes  20 samples  density 0.03  r >= 0.30
        edges 360 tested, 33 called at q < 0.1
        graphs: 3 cliques at q < 0.1, 120 at q < 0.9, 109 on all tested pairs
        cliques 232: complete_conserved   6, partial_significant   8, unclassified 218

---

    Code
      print(rcomplex(networks = nets, orthologs = syn$ortho))
    Output
      rcomplex: 3 species, 120 hogs, 3 tiers
        SpA  120 genes  ? samples  density 0.03  r ?
        SpB  120 genes  ? samples  density 0.03  r ?
        SpC  120 genes  ? samples  density 0.03  r ?
        edges 360 tested, 31 called at q < 0.1
        graphs: 3 cliques at q < 0.1, 120 at q < 0.9, 111 on all tested pairs
        cliques 234: complete_conserved   6, partial_significant   6, unclassified 222

