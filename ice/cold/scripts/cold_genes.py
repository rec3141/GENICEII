"""Cold-adaptation gene inventory per MAG and expression (CPM of CDS read pairs) in ice vs water and TOP vs BOT.
-> ../cold_gene_members.tsv, ../cold_inventory_mag.tsv, ../cold_cpm.tsv, ../cold_paired.tsv, ../cold_taxa_ice.tsv"""
import re, os, sys; sys.path.insert(0, "/data/dev/rnaseq-exploration/ice/cold/scripts")
from icecommon import *
from cats import CATS
X, g = genes(); m = design(); S = samples()
cds = (g['type'].values == 'CDS')
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
kos = g['KO'].fillna('').str.split(',')
prod = g['product'].fillna('').str.lower()
pf = pd.read_csv(C_DIR + 'pf_hits.tsv', sep='\t') if os.path.exists(C_DIR + 'pf_hits.tsv') else pd.DataFrame(columns=['gene_id', 'pfam', 'score', 'plen'])
gid2row = pd.Series(np.arange(len(g)), index=g['gene_id'].values)
ko_set = {}
for i, ks in enumerate(kos.values):
    for k in ks:
        if k: ko_set.setdefault(k, []).append(i)
rows = []
for cat, d in CATS.items():
    sel = set()
    for k in d.get('ko', []): sel.update(ko_set.get(k, []))
    if 'prod' in d: sel.update(np.where(prod.str.contains(d['prod'], regex=True).values)[0])
    if 'pf' in d:
        h = pf[pf.pfam.isin(d['pf'])]
        if 'maxlen' in d: h = h[h.plen <= d['maxlen']]
        sel.update(gid2row.reindex(h.gene_id).dropna().astype(int).values)
    if 'maxlen' in d and d.get('ko'):  # KO-assigned cold-shock proteins must also be short
        pass
    sel = np.array(sorted(i for i in sel if cds[i]), int)
    rows.append(pd.DataFrame({'row': sel, 'category': cat}))
M = pd.concat(rows)
M['gene_id'] = g['gene_id'].values[M.row]; M['bin'] = g['bin'].values[M.row]; M['contig'] = g['contig'].values[M.row]
M['product'] = g['product'].values[M.row]; M['KO'] = g['KO'].values[M.row]
M.to_csv(C_DIR + 'cold_gene_members.tsv', sep='\t', index=False)
print(M.groupby('category').size().to_string())
# inventory per MAG
inv = M[M.bin.fillna('') != ''].groupby(['bin', 'category']).size().unstack(fill_value=0)
inv.to_csv(C_DIR + 'cold_inventory_mag.tsv', sep='\t')
# expression
E = pd.DataFrame({c: np.asarray(X[d.row.values].sum(0)).ravel() / tot * 1e6 for c, d in M.groupby('category')}, index=S)
E.to_csv(C_DIR + 'cold_cpm.tsv', sep='\t')
res = []
for kind in ['ice', 'top']:
    P = paired_log2(E, m, kind)
    res.append(pd.DataFrame({'contrast': kind, 'median_log2': P.median(), 'frac_up': (P > 0).mean(), 'n_blocks': len(P),
                             'q25': P.quantile(.25), 'q75': P.quantile(.75)}))
R = pd.concat(res); R.index.name = 'category'
R.round(3).to_csv(C_DIR + 'cold_paired.tsv', sep='\t')
pd.set_option('display.width', 200)
print(R.round(2).to_string())
