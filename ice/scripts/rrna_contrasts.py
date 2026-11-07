"""Ice-vs-water and top-vs-bottom paired contrasts of rRNA-read taxonomy (SILVA genus), separately for eukaryotes
(share of eukaryotic rRNA) and prokaryotes (share of bacterial+archaeal rRNA, plastid/mito excluded). Out: ice/out/rrna_{euk,prok}_{icewater,topbot}.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
import contrasts as C   # reuses pairs + contrast()
T = pd.read_csv(D + 'marine_coassembly.taxonomy_rrna.tsv', sep='\t')
g = T[T['rank'] == 'genus'].drop(columns='rank').set_index('lineage').apply(pd.to_numeric, errors='coerce').fillna(0)
for dom, sel in (('euk', g.index.str.startswith('Eukaryota')),
                 ('prok', ~g.index.str.startswith('Eukaryota') & ~g.index.str.contains('Chloroplast|Mitochondria'))):
    X = g[sel]; X = X / X.sum(axis=0)
    X = X[(X > 0.005).sum(axis=1) >= 5]  # genera >=0.5% of the domain in >=5 libraries
    L = np.log2(X + 1e-4)
    for name, P in (('icewater', C.PI), ('topbot', C.PT)):
        r, _ = C.contrast(L, P)
        r['mean_share'] = X.mean(axis=1)
        r.sort_values('delta').to_csv(C.O + f'rrna_{dom}_{name}.tsv', sep='\t')
