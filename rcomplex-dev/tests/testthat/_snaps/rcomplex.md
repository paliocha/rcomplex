# print() shows species, edges and tiers

    Code
      print(res)
    Output
      rcomplex: 3 species, 120 hogs, 1 tier
        SpA  120 genes  20 samples  density 0.03  r >= 0.26
        SpB  120 genes  20 samples  density 0.03  r >= 0.31
        SpC  120 genes  20 samples  density 0.03  r >= 0.30
        edges 360 tested, 33 called at q < 0.1
        cliques 3: complete_conserved 3

---

    Code
      print(rcomplex(networks = nets, orthologs = syn$ortho))
    Output
      rcomplex: 3 species, 120 hogs, 1 tier
        SpA  120 genes  ? samples  density 0.03  r ?
        SpB  120 genes  ? samples  density 0.03  r ?
        SpC  120 genes  ? samples  density 0.03  r ?
        edges 360 tested, 31 called at q < 0.1
        cliques 3: complete_conserved 3

