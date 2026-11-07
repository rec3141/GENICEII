"""(a) Per-KO/Pfam breakdown of the cold categories in ice vs water and TOP vs BOT;
(b) within-MAG regulation: category share of each MAG's own CDS transcripts, ice vs water, paired by block (removes community turnover).
-> ../cold_ko_paired.tsv, ../cold_within_mag.tsv, ../cold_within_summary.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
X, g = genes(); m = design(); S = samples()
cds = g['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
M = pd.read_csv(C_DIR + 'cold_gene_members.tsv', sep='\t', keep_default_na=False)
# (a) per-KO (first KO of the gene, else 'Pfam/product')
M['unit'] = np.where(M.KO != '', M.KO.str.split(',').str[0], M['product'].str.slice(0, 40))
KD = pd.read_csv('/data/dev/rnaseq-exploration/data/ko_list.tsv', sep='\t', header=None, usecols=[0, 1], names=['ko', 'def'], skiprows=1).set_index('ko')['def'] if True else None
rows = []
for (cat, u), d in M.groupby(['category', 'unit']):
    if len(d) < 5: continue
    e = pd.Series(np.asarray(X[d.row.values].sum(0)).ravel() / tot * 1e6, index=S)
    E = e.to_frame('v')
    r = {'category': cat, 'unit': u, 'definition': KD.get(u, ''), 'n_genes': len(d), 'mean_cpm': e.mean()}
    for kind in ['ice', 'top']:
        P = paired_log2(E, m, kind, pseudo=0.5)['v']
        r[f'{kind}_log2'] = P.median(); r[f'{kind}_frac_up'] = (P > 0).mean()
    rows.append(r)
K = pd.DataFrame(rows).sort_values(['category', 'mean_cpm'], ascending=[True, False])
K.round(3).to_csv(C_DIR + 'cold_ko_paired.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(K[K.mean_cpm >= 2].round(2).to_string(index=False, max_colwidth=45))
# (b) within-MAG
binned = (g['bin'].fillna('') != '').values & cds
bins = g['bin'].values
codes, uniq = pd.factorize(pd.Series(bins[binned]))
import scipy.sparse as sp
B = sp.csr_matrix((np.ones(binned.sum()), (codes, np.where(binned)[0])), shape=(len(uniq), len(g)))
mag_tot = pd.DataFrame((B @ X).toarray(), index=uniq, columns=S)
keep = mag_tot.index[(mag_tot.sum(1) >= 20000)]
out = []
for cat, d in M[M.bin != ''].groupby('category'):
    cc, _ = pd.factorize(pd.Series(d.bin.values).astype(pd.CategoricalDtype(uniq)))
    Bc = sp.csr_matrix((np.ones(len(d)), (pd.Categorical(d.bin, categories=uniq).codes, d.row.values)), shape=(len(uniq), len(g)))
    cat_cnt = pd.DataFrame((Bc @ X).toarray(), index=uniq, columns=S).loc[keep]
    T = mag_tot.loc[keep]
    for lab, a, b in blocks(m, 'ice'):
        ca, ta, cb, tb = cat_cnt[a].sum(1), T[a].sum(1), cat_cnt[b].sum(1), T[b].sum(1)
        ok = (ta >= 2000) & (tb >= 2000) & ((ca + cb) >= 10)
        lr = np.log2(((ca + .5) / ta) / ((cb + .5) / tb))[ok]
        out.append(pd.DataFrame({'category': cat, 'block': lab, 'bin': lr.index, 'log2': lr.values}))
W = pd.concat(out)
W.round(3).to_csv(C_DIR + 'cold_within_mag.tsv', sep='\t', index=False)
# per category: median over (MAG medians across blocks), fraction of MAGs up
per_mag = W.groupby(['category', 'bin']).log2.median().reset_index()
Sm = per_mag.groupby('category').log2.agg(n_mags='size', median='median', frac_up=lambda v: (v > 0).mean())
Sm.round(3).to_csv(C_DIR + 'cold_within_summary.tsv', sep='\t')
print(Sm.round(2).to_string())
