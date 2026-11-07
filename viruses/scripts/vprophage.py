"""Prophage induction: geNomad proviruses (hgt/mge/provirus.tsv) on host contigs. Genes inside each provirus come from Grex
annotation.tsv (reconciled tags, Nov 2) translated back to the local (old extra) tags via locus_tag_map. For proviruses on binned
contigs: provirus share of the host MAG's CDS transcripts, diesel vs control (SW, after dosing) and, in 2025-26, vs each pool's own
pre-dosing level (BACI); activity = share / provirus fraction of host genes. Out: out/prophage.tsv, out/prophage_cpm.pkl"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
import scipy.sparse as sp
V = '/data/dev/rnaseq-exploration/viruses/'
pv = pd.read_csv('/data/dev/rnaseq-exploration/hgt/mge/provirus.tsv', sep='\t')
F = pd.read_csv(V + 'grex/provirus_features.tsv', sep='\t', header=None, comment='#',
                names=['contig', 'type', 'start', 'stop', 'strand', 'tag', 'gene', 'product', 'xref'], keep_default_na=False)
tm = pd.read_csv(V + 'grex/provirus_tagmap.tsv', sep='\t', header=None, names=['extra_tag', 'new_tag', 'matched']).set_index('new_tag').extra_tag
F['gene_id'] = F.tag.map(tm)          # old (local) tag
F = F[F.type == 'cds']
X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'; tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
row = pd.Series(np.arange(len(G)), index=G.gene_id.values)
c2b = G.drop_duplicates('contig').set_index('contig').bin
bins = G['bin'].fillna('').values; ub, bi = np.unique(bins, return_inverse=True)
ind = sp.csr_matrix((cds.astype(float), (bi, np.arange(len(bins)))), shape=(len(ub), len(bins)))
BC = pd.DataFrame((ind @ X).toarray() / tot * 1e6, index=ub, columns=S); nb = pd.Series(np.asarray(ind.sum(1)).ravel(), index=ub)
rows, cp = [], {}
oil, ctl, pre = m.phase.eq('oiled'), m.phase.eq('control'), m.phase.eq('pre')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
def lab(b):
    if not isinstance(b, str) or not b: return 'unbinned'
    if special_label(b): return special_label(b)
    t = [x.split('__', 1)[1] for x in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in x and x.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b
for p in pv.itertuples():
    f = F[(F.contig == p.source_seq) & ((F.start + F.stop) / 2 >= p.start) & ((F.start + F.stop) / 2 <= p.end)]
    r_ = f.gene_id.map(row).dropna().astype(int).values
    if not len(r_): continue
    e = pd.Series(np.asarray(X[r_].sum(0)).ravel() / tot * 1e6, index=S); cp[p.seq_name] = e
    b = c2b.get(p.source_seq, ''); b = b if isinstance(b, str) else ''
    rec = dict(provirus=p.seq_name, contig=p.source_seq, start=p.start, end=p.end, length=p.length, n_genes=len(r_), bin=b, host=lab(b),
               integrase=bool(isinstance(p.integrases, str) and p.integrases), cpm_max=e.max())
    host = BC.loc[b] if b in BC.index else None
    for s in ('2024-25', '2025-26'):
        a = m.index[oil & (m.season == s) & (m.medium == 'SW')]; c = m.index[ctl & (m.season == s) & (m.medium == 'SW')]
        rec[f'cpm_oil_{s}'] = e[a].median(); rec[f'cpm_ctl_{s}'] = e[c].median()
        if host is not None:
            sh = e / host.where(host >= 20)
            rec[f'host_cpm_oil_{s}'] = host[a].median(); rec[f'host_cpm_ctl_{s}'] = host[c].median()
            rec[f'share_oil_{s}'] = sh[a].median(); rec[f'share_ctl_{s}'] = sh[c].median()
            rec[f'activity_oil_{s}'] = sh[a].median() / (len(r_) / nb[b]) if nb.get(b, 0) else np.nan
            if s == '2025-26':
                pe = m.index[pre & (m.season == s) & (m.medium == 'SW') & (m.group == 'EXP')]; pc = m.index[pre & (m.season == s) & (m.medium == 'SW') & (m.group == 'CTL')]
                q = lambda ix: sh[ix].median()
                rec['baci_share_2025-26'] = np.log2((q(a) + 1e-4) / (q(pe) + 1e-4)) - np.log2((q(c) + 1e-4) / (q(pc) + 1e-4))
    rows.append(rec)
P = pd.DataFrame(rows); P.to_csv(V + 'out/prophage.tsv', sep='\t', index=False)
pd.to_pickle(pd.DataFrame(cp).T, V + 'out/prophage_cpm.pkl')
pd.set_option('display.width', 260)
c = ['provirus', 'host', 'length', 'n_genes', 'cpm_oil_2024-25', 'cpm_ctl_2024-25', 'cpm_oil_2025-26', 'cpm_ctl_2025-26', 'share_oil_2025-26', 'share_ctl_2025-26', 'baci_share_2025-26', 'activity_oil_2025-26', 'activity_oil_2024-25']
print(len(P), 'proviruses with genes;', (P.bin != '').sum(), 'in MAGs')
print(P.sort_values('cpm_max', ascending=False)[c].head(25).round(3).to_string())
