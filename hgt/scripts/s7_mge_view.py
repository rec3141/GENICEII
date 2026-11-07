"""Catabolic MGEs: contigs carrying confident HC genes that are plasmid (strict / F- or T4SS-type), phage, provirus, integron or island.
Adds RNA time courses (summed HC CPM per pool/day) and DNA depth (CMO nanopore samples), plasmid:host copy ratio."""
import sys, re, json; sys.path.insert(0, '.'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import meta, nanopore_samples
from hcdef import *
H = '/data/dev/rnaseq-exploration/hgt/'; DN = '/data/dev/rnaseq-exploration/dna/'
h = pd.read_csv(H + 'hc_context.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
h = h.join(pd.read_csv(H + 'hc_host.tsv', sep='\t', keep_default_na=False).set_index('gene_id'))
for c in ['plasmid_any', 'plasmid_strict', 'plasmid_conj', 'phage_contig', 'provirus', 'genomic_island', 'integron', 'IS_within_10kb']:
    h[c] = h[c].astype(str).eq('True')
FALSE_PLASMID = {'contig_0030642': 'chromosomal: colinear with C. zancles 78-ME, flanked by fabA/fabB/flgE, copy ratio 1.1 (cyclo fork)'}
p = pd.read_csv(H + 'mge/plasmid_summary.tsv', sep='\t').set_index('seq_name')
E = pd.read_pickle(H + 'hc_cpm.pkl'); m = meta().loc[E.columns]
hc = h[h.tier == 'confident']
sel = hc[(hc.plasmid_strict | hc.phage_contig | hc.provirus | hc.integron | hc.genomic_island |
          (hc.plasmid_any & hc.contig.map(p.plasmid_score).ge(0.95))) & ~hc.contig.isin(FALSE_PLASMID)]
cg = pd.read_csv(H + 'grex/hc_contig_genes.tsv.gz', sep='\t', keep_default_na=False)
rx_is = re.compile(r'transposase|insertion element|\bIS\d|IS[0-9]+ family|integrase|resolvase|invertase', re.I)
dep = pd.read_csv(H + 'grex/hc_contig_depths.tsv.gz', sep='\t', index_col=0)
ns = nanopore_samples()
ns['col'] = ns.run_barcode + '.sorted.bam'
ns = ns[ns.col.isin(dep.columns) & ns.sample_name.str.match(r'CMO_(CTL|EXP)_')]
ns['pool'] = ns.sample_name.str.extract(r'CMO_(CTL|EXP)_')[0]
t = ns.sample_name.str.extract(r'_(PRE|DAY)(\d+)')
ns['day'] = np.where(t[0] == 'PRE', -t[1].astype(float), t[1].astype(float))
ns['season'] = np.where(ns.log.astype(str) == '2024', '2024-25', '2025-26')
bd = pd.read_csv(DN + 'bin_depths.tsv', sep='\t', index_col=0)
def tc(genes):
    """mean summed CPM per season/pool/day (core SW+SI samples)"""
    v = E.loc[genes].sum()
    d = m.assign(v=v)[m.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI'])]
    out = {}
    for (se, g), x in d.groupby(['season', 'group']):
        s = x.groupby('day').v.mean()
        out[f'{se} {g}'] = [[float(k), round(float(val), 3)] for k, val in s.items()]
    return out
rows, courses = [], {}
for ct, g in sel.groupby('contig'):
    allg = cg[cg.contig == ct]
    kinds = [k for k, f in [('plasmid (geNomad strict)', g.plasmid_strict.any()), ('plasmid (score>=0.95, no hallmark)', (~g.plasmid_strict & g.plasmid_any).any()),
                            ('phage', g.phage_contig.any()), ('provirus', g.provirus.any()), ('integron', g.integron.any()), ('genomic island', g.genomic_island.any())] if f]
    pr = p.loc[ct] if ct in p.index else None
    conj = str(pr.conjugation_genes) if pr is not None and pr.conjugation_genes == pr.conjugation_genes else ''
    mob = 'F-type T4SS' if 'F_tra' in conj else 'P/virB-type T4SS' if 'virB' in conj else 'T4SS/other' if conj else ''
    bin_ = g.bin.iloc[0]
    dna = {}
    if ct in dep.index:
        for se in ['2024-25', '2025-26']:
            for pool in ['CTL', 'EXP']:
                cols = ns[(ns.season == se) & (ns.pool == pool)].sort_values('day')
                dna[f'{se} {pool}'] = [[float(d), round(float(dep.loc[ct, c]), 3)] for d, c in zip(cols.day, cols.col)]
    copy = np.nan
    if bin_ and bin_ in bd.index and ct in dep.index:
        cols = [c for c in bd.columns if c in dep.columns]
        hb, pb = bd.loc[bin_, cols].astype(float), dep.loc[ct, cols].astype(float)
        ok = hb > 1
        if ok.sum() >= 3: copy = round(float((pb[ok] / hb[ok]).median()), 2)
    ex = E.loc[g.index]
    rows.append(dict(contig=ct, kind='; '.join(kinds), length=int(dep.contigLen.get(ct, 0)) if ct in dep.index else None,
                     plasmid_score=float(pr.plasmid_score) if pr is not None else None, hallmarks=int(pr.n_hallmarks) if pr is not None else None,
                     mobility=mob, conjugation_genes=conj, host_bin=bin_, host_genus=g.genus.iloc[0], host_family=g.family.iloc[0], host_tax_source=g.tax_source.iloc[0],
                     n_hc=len(g), hc_categories='; '.join(sorted(set(g.category))), hc_products='; '.join(g['product']),
                     n_IS_on_contig=int(allg['product'].str.contains(rx_is).sum()), n_genes_on_contig=len(allg),
                     dna_mean_depth=round(float(dep.totalAvgDepth.get(ct, np.nan)), 2), copies_per_host=copy,
                     rna_reads=int(0), max_cpm=round(float(ex.sum().max()), 2)))
    courses[ct] = dict(rna=tc(list(g.index)), dna=dna)
M = pd.DataFrame(rows)
oiled = {se: m.index[(m.season == se) & m.design.isin(['EXP_SW', 'EXP_SI']) & (m.phase == 'oiled')] for se in ['2024-25', '2025-26']}
ctl = {se: m.index[(m.season == se) & m.design.isin(['CTL_SW', 'CTL_SI']) & (m.phase == 'control')] for se in ['2024-25', '2025-26']}
for se in ['2024-25', '2025-26']:
    M[f'cpm_oiled_{se}'] = [round(float(E.loc[sel.index[sel.contig == c]].sum()[oiled[se]].mean()), 3) for c in M.contig]
    M[f'cpm_ctl_{se}'] = [round(float(E.loc[sel.index[sel.contig == c]].sum()[ctl[se]].mean()), 3) for c in M.contig]
    M[f'l2fc_{se}'] = np.log2((M[f'cpm_oiled_{se}'] + .05) / (M[f'cpm_ctl_{se}'] + .05)).round(2)
M = M.sort_values('cpm_oiled_2024-25', ascending=False)
M.to_csv(H + 'catabolic_mges.tsv', sep='\t', index=False)
json.dump(courses, open(H + 'catabolic_mge_courses.json', 'w'))
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 60)
print(len(M), M.kind.value_counts().to_string())
print(M[['contig', 'kind', 'length', 'mobility', 'host_genus', 'host_bin', 'hc_categories', 'n_IS_on_contig', 'copies_per_host', 'cpm_oiled_2024-25', 'cpm_ctl_2024-25', 'cpm_oiled_2025-26', 'cpm_ctl_2025-26']].head(30).to_string())
