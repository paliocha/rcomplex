## Signed networks: anticorrelation

`compute_network(sign = "negative")` builds the network of
anticorrelation. It negates every correlation between two genes and then
ranks, normalises and thresholds exactly as for the positive network, so
the strongest anticorrelations become the top edges; each gene keeps its
own correlation of 1 and ranks itself first in both signs. With a
`partition`, each level's correlation is negated before negative values
are set to zero, so the rectified average keeps anticorrelations only.
Every consumer reads membership, so `find_coexpressologs()` on a positive
network of one species and a negative network of another already tests
whether the first species' partners are the second's antipartners. The
default stays `sign = "positive"`. Read a negative network with care:
anticorrelations are rarer and weaker than correlations in RNA-seq data,
yet a density threshold keeps the top fraction of pairs at any sample
size, so a top-3 % negative network exists even when it holds only noise
(at n = 20 samples, r >= -0.3 is noise). The weakest correlation that
passed the threshold says whether the network holds anything.
