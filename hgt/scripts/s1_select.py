"""Select HC genes -> ../hc_genes.tsv. Evidence: Bakta product, coordinate-remapped KO (if present), CANT-HYD (if present)."""
import sys, os; sys.path.insert(0, '.')
from hcdef import *
D = '/data/dev/rnaseq-exploration/data/'
g = pd.read_pickle(D + 'genemap.pkl').drop(columns=['KO'])
ko = None
if os.path.exists(D + 'ko_remap.tsv.gz'):
    k = pd.read_csv(D + 'ko_remap.tsv.gz', sep='\t', dtype=str, keep_default_na=False)
    ko = k.set_index('gene_id')['KO']; ko = ko[ko != '']
c, ev = classify(g, ko)
h = g[c != ''].copy(); h['category'] = c[c != '']; h['evidence'] = ev[c != '']
ch = '/data/dev/rnaseq-exploration/canthyd/canthyd_hits.tsv'
t = pd.read_csv(ch, sep='\t', dtype=str).drop_duplicates('gene_id').set_index('gene_id')
t = t[t.hmm.isin(CANTHYD)]
gi = g.reset_index().set_index('gene_id')
for gid, hmm in t.hmm.items():
    i = gi.loc[gid, 'index']; cc = CANTHYD[hmm]
    if c.iat[i] == '': c.iat[i] = cc; ev.iat[i] = 'CANT-HYD'
    else:
        if c.iat[i] != cc and hmm in ('AlkB', 'CYP153'): c.iat[i] = cc   # HMM wins for the specific alkane hydroxylases
        ev.iat[i] = ev.iat[i] + '+CANT-HYD'
h = g[c != ''].copy(); h['category'] = c[c != '']; h['evidence'] = ev[c != '']
h['canthyd'] = h.gene_id.map(t.hmm).fillna('')
h['tier'] = np.where(h.evidence.str.contains('product|CANT-HYD'), 'confident', 'putative')
h['pathway'] = h.category.map(PATHWAY); h['step'] = h.category.map(STEP)
h['row'] = g.index[c != ''] if 'row' not in h else h['row']
h.to_csv('../hc_genes.tsv', sep='\t', index=False)
h.contig.drop_duplicates().to_csv('../hc_contigs.txt', index=False, header=False)
h.gene_id.to_csv('../hc_gene_ids.txt', index=False, header=False)
print(len(h), 'genes on', h.contig.nunique(), 'contigs; binned', (h.bin != '').mean().round(3))
print(h.groupby(['category', 'evidence']).size().unstack(fill_value=0).to_string())
