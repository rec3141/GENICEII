"""Mobile context of metal and toxicity genes, with the same flags hgt/scripts/s3_context.py uses for hydrocarbon genes, plus co-location
with confident hydrocarbon genes (within 10 kb). Gene coordinates come from the reconciled extra Bakta run on Grex (metals/grex/), translated
back to the pre-Nov-2 tags used locally via locus_tag_map. Out: metals/out/mobility.tsv (per gene), mobility_summary.tsv, cassettes.tsv."""
import sys, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
R = '/data/dev/rnaseq-exploration/'; O = R + 'metals/out/'; H = R + 'hgt/'
M = pd.read_csv(O + 'members.tsv', sep='\t', keep_default_na=False)
tm = pd.read_csv(R + 'metals/grex/tagmap_subset.tsv', sep='\t', header=None, names=['old', 'new', 'matched'])
F = pd.read_csv(R + 'metals/grex/contig_features.tsv.gz', sep='\t', header=None, comment='#',
                names=['contig', 'type', 'start', 'stop', 'strand', 'tag', 'gene', 'product'], keep_default_na=False)
co = F.set_index('tag')[['contig', 'start', 'stop', 'strand']]
new = dict(zip(tm.old, tm.new))
M['new_tag'] = M.gene_id.map(new)
M = M.join(co[['start', 'stop', 'strand']], on='new_tag')
print('members with coordinates:', M.start.notna().sum(), 'of', len(M))
hc = pd.read_csv(H + 'hc_genes.tsv', sep='\t', keep_default_na=False); hc = hc[hc.tier == 'confident']
hc['new_tag'] = hc.gene_id.map(new); hc = hc.join(co[['contig', 'start', 'stop']].rename(columns={'contig': 'c2'}), on='new_tag')
p = pd.read_csv(H + 'mge/plasmid_summary.tsv', sep='\t'); v = pd.read_csv(H + 'mge/virus_summary.tsv', sep='\t'); v = v[~v.seq_name.str.contains(r'\|provirus')]
pv = pd.read_csv(H + 'mge/provirus.tsv', sep='\t'); gi = pd.read_csv(H + 'mge/genomic_islands.tsv', sep='\t'); it = pd.read_csv(H + 'mge/integrons.tsv', sep='\t')
FALSE_PLASMID = {'contig_0030642'}
M = M[M.start.notna()].copy(); M['start'] = M.start.astype(int); M['stop'] = M.stop.astype(int)
M['plasmid_any'] = M.contig.isin(p.seq_name) & ~M.contig.isin(FALSE_PLASMID)
M['plasmid_strict'] = M.contig.isin(p[(p.plasmid_score >= 0.9) & ((p.n_hallmarks >= 1) | p.conjugation_genes.notna())].seq_name)
M['phage_contig'] = M.contig.isin(v[(v.n_hallmarks >= 2) & (v.virus_score >= 0.95)].seq_name)
def within(df, c, s, e, pad=0):
    grp = {k: d[[s, e]].values for k, d in df.groupby(c)}
    return np.array([(lambda r: r is not None and bool(((r[:, 0] - pad <= b) & (r[:, 1] + pad >= a)).any()))(grp.get(ct)) for ct, a, b in zip(M.contig, M.start, M.stop)])
M['provirus'] = within(pv, 'source_seq', 'start', 'end')
M['genomic_island'] = within(gi, 'contig', 'start', 'end')
ib = it.groupby(['ID_integron', 'ID_replicon']).agg(start=('pos_beg', 'min'), end=('pos_end', 'max')).reset_index()
M['integron'] = within(ib, 'ID_replicon', 'start', 'end', pad=1000)
rx_is = re.compile(r'transposase|insertion element|\bIS\d|IS[0-9]+ family|integrase|recombinase|resolvase|invertase', re.I)
house = F['product'].str.contains(r'RecA|Xer[CD]|Holliday|Ruv[ABCX]|RecB|RecG|RecN|RecO|RecR|DNA repair', case=False, regex=True)
isf = F[(F['product'].str.contains(rx_is) | F['gene'].str.contains(r'^(tnp|ins[A-Z]|intI?$)', case=False, regex=True)) & ~house]
isd = isf.groupby('contig')[['start', 'stop']].apply(lambda d: d.values).to_dict()
hcd = hc.dropna(subset=['start']).groupby('c2')[['start', 'stop']].apply(lambda d: d.values.astype(int)).to_dict()
hcc = hc.dropna(subset=['start']).groupby('c2').apply(lambda d: list(zip(d.start.astype(int), d.stop.astype(int), d.category))).to_dict()
def mind(d, ct, a, b):
    r = d.get(ct)
    if r is None: return np.nan
    return float(np.maximum(0, np.maximum(r[:, 0] - b, a - r[:, 1])).min())
M['IS_dist'] = [mind(isd, c, a, b) for c, a, b in zip(M.contig, M.start, M.stop)]
M['IS_within_10kb'] = M.IS_dist <= 10000
M['hc_dist'] = [mind(hcd, c, a, b) for c, a, b in zip(M.contig, M.start, M.stop)]
M['hc_within_10kb'] = M.hc_dist <= 10000
M['hc_near_cats'] = ['; '.join(sorted({cat for s, e, cat in hcc.get(c, []) if max(0, s - b, a - e) <= 10000})) for c, a, b in zip(M.contig, M.start, M.stop)]
M['mge_borne'] = M[['plasmid_strict', 'phage_contig', 'provirus', 'integron', 'genomic_island']].any(axis=1)
M['mobile_assoc'] = M.mge_borne | M.IS_within_10kb
M.to_csv(O + 'mobility.tsv', sep='\t', index=False)
cols = ['plasmid_strict', 'phage_contig', 'provirus', 'genomic_island', 'integron', 'IS_within_10kb', 'mge_borne', 'mobile_assoc', 'hc_within_10kb']
Mc = M[M.confident.astype(str) == 'True']
summ = Mc.groupby(['group', 'category'])[cols].mean().mul(100).round(2); summ.insert(0, 'n', Mc.groupby(['group', 'category']).size())
ref = pd.read_csv(H + 'hc_context.tsv', sep='\t', keep_default_na=False)
ref = ref[ref.tier == 'confident']
rr = {c: ref[c].astype(str).eq('True').mean() * 100 for c in cols if c in ref.columns}
summ.loc[('reference', 'Hydrocarbon genes (confident, hgt/)'), :] = [len(ref)] + [round(rr.get(c, np.nan), 2) for c in cols]
summ.to_csv(O + 'mobility_summary.tsv', sep='\t')
# cassettes: mobile toxic-metal genes grouped by contig, with nearby hydrocarbon genes
tm_ = Mc[(Mc.group == 'toxic metal') & Mc.mobile_assoc]
cas = tm_.groupby('contig').agg(bin=('bin', 'first'), n=('gene_id', 'size'), categories=('category', lambda x: '; '.join(sorted(set(x)))),
                                start=('start', 'min'), end=('stop', 'max'), plasmid=('plasmid_strict', 'any'), island=('genomic_island', 'any'),
                                integron=('integron', 'any'), IS=('IS_within_10kb', 'any'), hc_near=('hc_near_cats', lambda x: '; '.join(sorted({c for s in x for c in s.split('; ') if c}))),
                                products=('product', lambda x: '; '.join(pd.Series(x).value_counts().index[:6]))).reset_index()
cas.sort_values('n', ascending=False).to_csv(O + 'cassettes.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80)
print(summ.to_string())
print(cas[cas.hc_near != ''].sort_values('n', ascending=False).head(15).to_string())
