"""Ice-binding protein (DUF3494 etc.) genes: host taxon (GTDB bin or Kaiju contig), expression in ice vs water and TOP vs BOT.
-> ../ibp_genes.tsv, ../ibp_taxa.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
X, g = genes(); m = design(); S = samples()
cds = g['type'].values == 'CDS'; tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
M = pd.read_csv(C_DIR + 'cold_gene_members.tsv', sep='\t', keep_default_na=False)
I = M[M.category == 'Ice-binding proteins'].copy()
pf = pd.read_csv(C_DIR + 'pf_dom.tsv', sep='\t'); pf = pf[pf.pfam.isin(['DUF3494', 'AFP', 'CfAFP', 'Ice_nucleation'])].drop_duplicates('gene_id').set_index('gene_id')
I['plen'] = I.gene_id.map(pf.tlen); I['pf_score'] = I.gene_id.map(pf.score)
tx = pd.read_csv('/data/dev/rnaseq-exploration/data/gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1])
tx['bin'] = np.where(tx.user_genome.str.startswith('dastool-'), tx.user_genome, 'dastool-' + tx.user_genome)
tx = tx.drop_duplicates('bin').set_index('bin').classification
kj = pd.read_csv('/data/dev/rnaseq-exploration/data/kaiju_contigs.tsv', sep='\t', index_col=0)
def lineage(r):
    if r.bin and r.bin in tx.index:
        t = dict(x.split('__') for x in tx[r.bin].split(';'))
        return t.get('p', ''), t.get('f', ''), t.get('g', '') or kaiju_bin_label(r.bin) or t.get('f', ''), 'GTDB'
    k = kj.lineage.get(r.contig)
    if not isinstance(k, str): return '', '', 'unassigned', ''
    p = (k.split('; ') + [''] * 7)[:7]
    return p[1], p[4], (p[5] or p[4] or p[1] or 'unassigned') + ' (Kaiju)', 'Kaiju'
L = I.apply(lineage, axis=1, result_type='expand'); L.columns = ['phylum', 'family', 'genus', 'tax_src']
I = I.join(L)
ti = pd.read_csv(C_DIR + 'grex/ibp_tiara.tsv', sep='\t').set_index('sequence_id').class_fst_stage
I['tiara'] = I.contig.map(ti).fillna('n/a')
I['domain'] = np.where(I.tax_src == 'GTDB', 'bacteria (MAG)', np.where(I.tiara == 'eukarya', 'eukaryote', np.where(I.tiara.isin(['bacteria', 'prokarya', 'archaea']), 'prokaryote (unbinned)', 'unknown')))
E = pd.DataFrame(X[I.row.values].toarray().T / tot[:, None] * 1e6, index=S, columns=I.gene_id.values)
I['reads'] = np.asarray(X[I.row.values].sum(1)).ravel()
si, sw = m.medium == 'SI', m.medium == 'SW'
I['cpm_SI'] = E.loc[si].mean().values; I['cpm_SW'] = E.loc[sw].mean().values
P = paired_log2(E, m, 'ice', pseudo=0.05)
I['ice_log2'] = P.median().values
I.drop(columns='row').round(4).to_csv(C_DIR + 'ibp_genes.tsv', sep='\t', index=False)
T = I.groupby(['domain', 'phylum', 'genus']).agg(n=('gene_id', 'size'), bins=('bin', lambda b: len(set(b) - {''})), reads=('reads', 'sum'),
                                       cpm_SI=('cpm_SI', 'sum'), cpm_SW=('cpm_SW', 'sum')).sort_values('cpm_SI', ascending=False)
T['share_SI'] = T.cpm_SI / T.cpm_SI.sum() * 100
T.round(3).to_csv(C_DIR + 'ibp_taxa.tsv', sep='\t')
pd.set_option('display.width', 200)
D = I.groupby('domain').agg(n=('gene_id', 'size'), cpm_SI=('cpm_SI', 'sum'), cpm_SW=('cpm_SW', 'sum'))
D['share_SI'] = D.cpm_SI / D.cpm_SI.sum() * 100; D['share_SW'] = D.cpm_SW / D.cpm_SW.sum() * 100
D.round(3).to_csv(C_DIR + 'ibp_domain.tsv', sep='\t'); print(D.round(2))
print(len(I), 'IBP genes;', (I.bin != '').sum(), 'binned in', I.bin.replace('', np.nan).nunique(), 'MAGs')
print(T.head(25).round(2).to_string())
print(I.sort_values('cpm_SI', ascending=False).head(15)[['gene_id', 'bin', 'genus', 'plen', 'product', 'cpm_SI', 'cpm_SW', 'ice_log2']].round(2).to_string())
