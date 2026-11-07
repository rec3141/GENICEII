"""Mobile-genetic context of each HC gene -> ../hc_context.tsv."""
import sys, gzip, re; sys.path.insert(0, '.')
from hcdef import *
H = '/data/dev/rnaseq-exploration/hgt/'
h = pd.read_csv(H + 'hc_genes.tsv', sep='\t', keep_default_na=False)
cg = pd.read_csv(H + 'grex/hc_contig_genes.tsv.gz', sep='\t', keep_default_na=False, dtype={'start': int, 'end': int})
co = cg.set_index('gene_id')[['start', 'end', 'strand']]
h = h.join(co, on='gene_id')
# --- MGE calls
p = pd.read_csv(H + 'mge/plasmid_summary.tsv', sep='\t')
v = pd.read_csv(H + 'mge/virus_summary.tsv', sep='\t'); v = v[~v.seq_name.str.contains(r'\|provirus')]
pv = pd.read_csv(H + 'mge/provirus.tsv', sep='\t')
gi = pd.read_csv(H + 'mge/genomic_islands.tsv', sep='\t')
it = pd.read_csv(H + 'mge/integrons.tsv', sep='\t')
pd_ = p.set_index('seq_name')
FALSE_PLASMID = {'contig_0030642'}   # Cycloclasticus chromosome (colinear with C. zancles; cyclo fork)
h['plasmid_any'] = h.contig.isin(p.seq_name) & ~h.contig.isin(FALSE_PLASMID)
strict_p = p[(p.plasmid_score >= 0.9) & ((p.n_hallmarks >= 1) | p.conjugation_genes.notna())]
h['plasmid_strict'] = h.contig.isin(strict_p.seq_name)
h['plasmid_conj'] = h.contig.isin(p[p.conjugation_genes.notna()].seq_name)
h['phage_contig'] = h.contig.isin(v[(v.n_hallmarks >= 2) & (v.virus_score >= 0.95)].seq_name)
def within(df, c, s, e, pad=0):
    out = np.zeros(len(h), bool)
    grp = {k: d[[s, e]].values for k, d in df.groupby(c)}
    for i, (ct, a, b) in enumerate(zip(h.contig, h.start, h.end)):
        r = grp.get(ct)
        if r is not None and ((r[:, 0] - pad <= b) & (r[:, 1] + pad >= a)).any(): out[i] = True
    return out
h['provirus'] = within(pv, 'source_seq', 'start', 'end')
h['genomic_island'] = within(gi, 'contig', 'start', 'end')
ib = it.groupby(['ID_integron', 'ID_replicon']).agg(start=('pos_beg', 'min'), end=('pos_end', 'max')).reset_index()
h['integron'] = within(ib, 'ID_replicon', 'start', 'end', pad=1000)
# --- IS/transposase/integrase within +-10 kb, conjugation machinery on same contig
rx_is = re.compile(r'transposase|insertion element|\bIS\d|IS[0-9]+ family|integrase|recombinase|resolvase|invertase', re.I)
rx_conj = re.compile(r'conjugal|conjugative|type IV secretion|relaxase|\bTra[A-Z]\b|\bTrb[A-Z]\b|mobilization protein|\bMob[A-Z]\b', re.I)
housekeeping = cg['product'].str.contains(r'RecA|Xer[CD]|Holliday|Ruv[ABCX]|RecB|RecG|RecN|RecO|RecR|DNA repair', case=False, regex=True)
cg['is'] = (cg['product'].str.contains(rx_is) | cg['gene'].str.contains(r'^(tnp|ins[A-Z]|intI?$)', case=False, regex=True)) & ~housekeeping
cg['conj'] = cg['product'].str.contains(rx_conj)
isd = cg[cg['is']].groupby('contig')[['start', 'end']].apply(lambda d: d.values).to_dict()
cj = set(cg.loc[cg.conj, 'contig'])
near, dist = np.zeros(len(h), bool), np.full(len(h), np.nan)
for i, (ct, a, b) in enumerate(zip(h.contig, h.start, h.end)):
    r = isd.get(ct)
    if r is None: continue
    d = np.maximum(0, np.maximum(r[:, 0] - b, a - r[:, 1]))
    dist[i] = d.min(); near[i] = d.min() <= 10000
h['IS_within_10kb'] = near; h['IS_dist'] = dist
h['conj_on_contig'] = h.contig.isin(cj)
# --- GC deviation (gene vs contig)
import os
seq, cur = {}, None
for line in (open(H + 'grex/hc_regions.fna') if os.path.exists(H + 'grex/hc_regions.fna') else []):
    if line[0] == '>': cur = line[1:].split()[0]; seq[cur] = []
    else: seq[cur].append(line.strip())
gc = lambda s: (s.upper().count('G') + s.upper().count('C')) / max(1, len(s))
reg = pd.read_csv(H + 'grex/regions.tsv', sep='\t', header=None, names=['region', 'tag'])
ggc = {t: gc(''.join(seq.get(r, ''))) for r, t in zip(reg.region, reg.tag) if not t.endswith('_flank')}
cgc = pd.read_csv(H + 'grex/hc_contig_gc.tsv', sep='\t')
cgc.columns = ['contig'] + list(cgc.columns[1:])
gccol = [c for c in cgc.columns if 'gc' in c.lower()][0]
h['gene_gc'] = h.gene_id.map(ggc) * 100
h['contig_gc'] = h.contig.map(cgc.set_index('contig')[gccol]).astype(float)
if h.contig_gc.max() <= 1: h['contig_gc'] *= 100
h['dGC'] = h.gene_gc - h.contig_gc
dep = pd.read_csv(H + 'grex/hc_contig_depths.tsv.gz', sep='\t', usecols=['contigName', 'contigLen', 'totalAvgDepth']).set_index('contigName')
h['contig_len'] = h.contig.map(dep.contigLen); h['contig_depth'] = h.contig.map(dep.totalAvgDepth).round(2)
cv = pd.read_csv(H + 'grex/hc_checkv.tsv', sep='\t').set_index('contig_id')
h['checkv_quality'] = h.contig.map(cv.checkv_quality).fillna('')
h['phage_contig'] &= ~(h.contig.map(cv.host_genes).fillna(0) > h.contig.map(cv.viral_genes).fillna(0))
h['mge_borne'] = h[['plasmid_strict', 'phage_contig', 'provirus', 'integron', 'genomic_island']].any(axis=1)
h['mobile_assoc'] = h.mge_borne | h.IS_within_10kb
h.to_csv(H + 'hc_context.tsv', sep='\t', index=False)
cols = ['plasmid_any', 'plasmid_strict', 'plasmid_conj', 'phage_contig', 'provirus', 'genomic_island', 'integron', 'IS_within_10kb', 'conj_on_contig', 'mge_borne', 'mobile_assoc']
print(h[cols].sum().to_string())
print(h.groupby('category')[cols].mean().mul(100).round(1).to_string())
print('dGC by mobile_assoc:', h.groupby('mobile_assoc').dGC.describe().round(2).to_string())
