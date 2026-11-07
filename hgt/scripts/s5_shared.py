"""Cross-taxon near-identical HC proteins (>=95% aa id, >=80% cov both) -> candidate HGT pairs + clusters."""
import sys; sys.path.insert(0, '.')
from hcdef import *
H = '/data/dev/rnaseq-exploration/hgt/'
h = pd.read_csv(H + 'hc_genes.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
tx = pd.read_csv(H + 'hc_host.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
h = h.join(tx)
a = pd.read_csv(H + 'grex/hc_ava.m8', sep='\t', header=None, names=['q', 't', 'fident', 'alnlen', 'qcov', 'tcov', 'evalue', 'bits'])
a = a[(a.q < a.t) & (a.fident >= 0.95) & (a.qcov >= 0.8) & (a.tcov >= 0.8)]
a = a[a.q.isin(h.index) & a.t.isin(h.index)]
Q, T = h.loc[a.q], h.loc[a.t]
a['q_contig'], a['t_contig'] = Q.contig.values, T.contig.values
a = a[a.q_contig != a.t_contig]
for side, X in (('q', Q), ('t', T)):
    for c in ['bin', 'category', 'tier', 'phylum', 'class', 'order', 'family', 'genus', 'tax_source']:
        a[f'{side}_{c}'] = X.loc[a[side]][c].values if False else h.loc[a[side], c].values
def diff_level(r):
    for lvl in ['phylum', 'class', 'order', 'family', 'genus']:
        x, y = r[f'q_{lvl}'], r[f't_{lvl}']
        if x and y and x != y: return lvl
    return ''
a['differ_at'] = a.apply(diff_level, axis=1)
a['both_mags'] = (a.q_bin != '') & (a.t_bin != '') & (a.q_bin != a.t_bin)
a['both_gtdb'] = (a.q_tax_source == 'GTDB(MAG)') & (a.t_tax_source == 'GTDB(MAG)')
a.to_csv(H + 'shared_pairs.tsv', sep='\t', index=False)
x = a[a.differ_at.isin(['phylum', 'class', 'order', 'family'])]
print('pairs >=95% id, different contigs:', len(a), '| hosts differ at family or above:', len(x), '| of which both GTDB MAGs:', x.both_gtdb.sum())
print(x.groupby('q_category').size().sort_values(ascending=False).to_string())
# union-find clusters over all >=95% pairs
par = {}
def f(u):
    while par.setdefault(u, u) != u: par[u] = par[par[u]]; u = par[u]
    return u
for q, t in zip(a.q, a.t): par[f(q)] = f(t)
cl = pd.Series({g: f(g) for g in par})
C = h.loc[cl.index].assign(cluster=cl.values)
summ = C.groupby('cluster').agg(n=('contig', 'size'), n_contigs=('contig', 'nunique'), n_families=('family', lambda s: s[s != ''].nunique()),
                                n_phyla=('phylum', lambda s: s[s != ''].nunique()), category=('category', 'first'),
                                families=('family', lambda s: ';'.join(sorted(set(s) - {''}))))
summ = summ.sort_values(['n_families', 'n'], ascending=False)
summ.to_csv(H + 'shared_clusters.tsv', sep='\t')
C.to_csv(H + 'shared_cluster_members.tsv', sep='\t')
print(summ[summ.n_families >= 2].head(30).to_string())
