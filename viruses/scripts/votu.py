"""Cluster confident viral sequences into vOTUs: greedy, longest first, members at >= 95% ANI and >= 85% aligned fraction of the
shorter sequence (MIUViG-style), from skani triangle -i --sparse. Out: out/votu.tsv (seq, region, votu, rep)."""
import pandas as pd, numpy as np
V = '/data/dev/rnaseq-exploration/viruses/'
vc = pd.read_csv(V + 'out/viral_confident.tsv', sep='\t')
sk = pd.read_csv(V + 'grex/skani_sparse.tsv.gz', sep='\t')
name = lambda s: s.split()[0]
sk['a'] = sk.Ref_name.map(name); sk['b'] = sk.Query_name.map(name)
L = dict(zip(vc.region, vc.length))
# aligned fraction of the shorter sequence
sk['af_short'] = np.where(sk.a.map(L) <= sk.b.map(L), sk.Align_fraction_ref, sk.Align_fraction_query)
e = sk[(sk.ANI >= 95) & (sk.af_short >= 85)]
nb = {}
for a, b in zip(e.a, e.b): nb.setdefault(a, set()).add(b); nb.setdefault(b, set()).add(a)
rep = {}
for r in vc.sort_values('length', ascending=False).region:
    if r in rep: continue
    rep[r] = r
    for x in nb.get(r, ()):
        if x not in rep and L.get(x, 0) <= L[r]: rep[x] = r
vc['votu_rep'] = vc.region.map(rep)
vc.to_csv(V + 'out/viral_confident.tsv', sep='\t', index=False)
print(len(vc), 'sequences ->', vc.votu_rep.nunique(), 'vOTUs;', (vc.groupby('votu_rep').size() > 1).sum(), 'multi-member')
