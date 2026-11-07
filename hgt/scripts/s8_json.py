"""Assemble website-ready hgt.json."""
import sys, json; sys.path.insert(0, '.'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import meta
from hcdef import *
H = '/data/dev/rnaseq-exploration/hgt/'
h = pd.read_csv(H + 'hc_context.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
h = h.join(pd.read_csv(H + 'hc_host.tsv', sep='\t', keep_default_na=False).set_index('gene_id'))
h = h.join(pd.read_csv(H + 'hc_expr_summary.tsv', sep='\t', index_col=0))
B = ['plasmid_any', 'plasmid_strict', 'plasmid_conj', 'phage_contig', 'provirus', 'genomic_island', 'integron', 'IS_within_10kb', 'mge_borne', 'mobile_assoc']
for c in B: h[c] = h[c].astype(str).eq('True')
E = pd.read_pickle(H + 'hc_cpm.pkl'); m = meta().loc[E.columns]
conf = h[h.tier == 'confident']
J = {'generated': pd.Timestamp.now().isoformat(timespec='minutes'),
     'design': {'treatment': 'diesel both seasons (2024-25 added DAY05 16:30; 2025-26 ~136 L DAY01)', 'pools': 'one oiled (EXP) and one control (CTL) pool per season; samples are pseudoreplicates'},
     'definitions': {
         'confident': 'Bakta (extra) product match and/or CANT-HYD HMM (trusted cutoff); putative = coordinate-remapped eggNOG KO only',
         'mge_borne': 'geNomad plasmid with >=1 hallmark or conjugation genes (score>=0.9), geNomad virus contig (>=2 hallmarks, score>=0.95, CheckV host<viral), geNomad provirus, IslandPath island, or IntegronFinder integron (+-1 kb)',
         'IS_within_10kb': 'transposase/IS/integrase/resolvase (Bakta product, excl. RecA/Xer/Ruv) within 10 kb',
         'mobile_assoc': 'mge_borne or IS_within_10kb',
         'hgt_candidate': '>=95% aa id & >=80% cov (mmseqs), >=99% nt identity, different genus (GTDB/Kaiju harmonised), DNA depth-profile rho<0.2 across 276 nanopore samples (same-MAG contigs median 0.60, random -0.52); flank test: identical block in otherwise unaligned +-20 kb flanks'}}
cnt = {}
for tier, d in [('confident', conf), ('all', h)]:
    cnt[tier] = {'genes': int(len(d)), 'contigs': int(d.contig.nunique()), 'binned_frac': round(float((d.bin != '').mean()), 3),
                 **{c: int(d[c].sum()) for c in B}}
J['counts'] = cnt
J['by_category'] = (conf.groupby('category')[B].mean().mul(100).round(1).assign(n=conf.groupby('category').size(), pathway=lambda x: x.index.map(PATHWAY))
                    .reset_index().to_dict(orient='records'))
J['mobile_vs_chrom'] = pd.read_csv(H + 'mobile_vs_chrom_summary.tsv', sep='\t', index_col=[0, 1, 2]).reset_index().rename(
    columns={'level_0': 'tier', 'level_1': 'season', 'level_2': 'flag'}).replace({np.nan: None}).to_dict(orient='records')
J['mobile_share_by_category'] = pd.read_csv(H + 'mobile_share_by_category.tsv', sep='\t').replace({np.nan: None}).to_dict(orient='records')
# mobile HC gene table (expressed)
mg = pd.read_csv(H + 'mobile_hc_genes.tsv', sep='\t', index_col=0, keep_default_na=False)
keep = ['category', 'tier', 'product', 'contig', 'bin', 'genus', 'family', 'tax_source', 'plasmid_any', 'plasmid_strict', 'phage_contig', 'provirus', 'genomic_island',
        'integron', 'IS_within_10kb', 'IS_dist', 'contig_len', 'contig_depth', 'oiled_2024-25', 'ctl_2024-25', 'l2fc_2024-25', 'oiled_2025-26', 'ctl_2025-26', 'l2fc_2025-26', 'reads_total']
mg = mg[keep].head(150).reset_index().rename(columns={'index': 'gene_id'})
for c in ['oiled_2024-25', 'ctl_2024-25', 'l2fc_2024-25', 'oiled_2025-26', 'ctl_2025-26', 'l2fc_2025-26', 'contig_depth']:
    mg[c] = pd.to_numeric(mg[c], errors='coerce').round(3)
J['mobile_genes'] = mg.replace({np.nan: None, '': None}).to_dict(orient='records')
# catabolic MGEs + courses
cm = pd.read_csv(H + 'catabolic_mges.tsv', sep='\t', keep_default_na=False)
J['catabolic_mges'] = cm.head(80).replace({np.nan: None, '': None}).to_dict(orient='records')
courses = json.load(open(H + 'catabolic_mge_courses.json'))
J['catabolic_mge_courses'] = {c: courses[c] for c in cm.contig.head(15) if c in courses}
# top mobile gene time courses
def tc(genes):
    v = E.loc[genes].sum(); d = m.assign(v=v)[m.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI'])]
    return {f'{se} {g}': [[float(k), round(float(x), 3)] for k, x in s.groupby('day').v.mean().items()] for (se, g), s in d.groupby(['season', 'group'])}
top = mg[(mg.reads_total.astype(float) >= 50)].head(12)
J['mobile_gene_courses'] = {r.gene_id: {'label': f"{r.category} | {r.genus or '?'} | {r.contig}", 'course': tc([r.gene_id])} for r in top.itertuples()}
# HGT network
F = pd.read_csv(H + 'hgt_candidates_flanks.tsv', sep='\t', keep_default_na=False)
C = pd.read_csv(H + 'hgt_candidates.tsv', sep='\t', keep_default_na=False)
nodes, edges = {}, []
def node(genus, fam, cls, src):
    k = genus or f'?({fam})'
    nodes.setdefault(k, {'id': k, 'family': fam, 'class': cls, 'tax_source': src}); return k
for r in C.itertuples():
    f = F[((F.q == r.q) & (F.t == r.t))]
    interp = f.interpretation.iloc[0] if len(f) else F[(F.q_contig == r.q_contig) & (F.t_contig == r.t_contig)].interpretation.iloc[0] if len(F[(F.q_contig == r.q_contig) & (F.t_contig == r.t_contig)]) else ''
    a = node(r.q_genus, r.q_family, r.q_cls, r.q_tax_source); b = node(r.t_genus, r.t_family, r.t_cls, r.t_tax_source)
    edges.append({'source': a, 'target': b, 'category': r.q_category, 'level': r.level, 'aa_ident': float(r.fident), 'nt_ident': float(r.nt_ident),
                  'depth_rho': None if r.depth_rho == '' else round(float(r.depth_rho), 2), 'q_gene': r.q, 't_gene': r.t, 'q_contig': r.q_contig, 't_contig': r.t_contig,
                  'flank_call': interp,
                  'l2fc_2024_25_q': round(float(h.loc[r.q, 'l2fc_2024-25']), 2) if r.q in h.index else None,
                  'l2fc_2025_26_q': round(float(h.loc[r.q, 'l2fc_2025-26']), 2) if r.q in h.index else None})
J['hgt_network'] = {'nodes': list(nodes.values()), 'edges': edges}
J['hgt_flanks'] = F.drop(columns=['block_products']).replace({np.nan: None, '': None}).to_dict(orient='records')
sp = pd.read_csv(H + 'shared_pairs.tsv', sep='\t')
J['sharing_funnel'] = {'pairs_95aa_diff_contigs': int(len(sp)), 'nt_ge_99': int((sp.nt_ident >= 99).sum()),
                       'cross_genus_nt99': int(len(C) + 0), 'flank_confirmed_HGT': int(F.interpretation.str.startswith('HGT').sum())}
json.dump(J, open(H + 'hgt.json', 'w'), indent=None, separators=(',', ':'), default=lambda o: None if o != o else (o.item() if hasattr(o, 'item') else str(o)))
print('wrote', H + 'hgt.json', round(len(open(H + 'hgt.json').read()) / 1e3), 'kB')
# --- extras: induction by genus, gene maps of key catabolic plasmids
J['induction_by_genus'] = pd.read_csv(H + 'induction_by_genus.tsv', sep='\t', keep_default_na=False).to_dict(orient='records')
import re
cg = pd.read_csv(H + 'grex/hc_contig_genes.tsv.gz', sep='\t', keep_default_na=False)
rx_is = re.compile(r'transposase|insertion element|\bIS\d|IS[0-9]+ family|integrase|resolvase|invertase', re.I)
rx_conj = re.compile(r'conjugal|conjugative|type IV secretion|relaxase|\bTra[A-Z]\b|\bTrb[A-Z]\b|mobilization protein|\bMob[A-Z]\b|pilin|pilus', re.I)
hcat = h['category'].to_dict()
maps = {}
for ct in ['contig_0009679', 'contig_0014557', 'contig_0009104', 'contig_0003430', 'contig_0022367', 'contig_0025080']:
    x = cg[cg.contig == ct]
    maps[ct] = [{'gene_id': r.gene_id, 'start': int(r.start), 'end': int(r.end), 'strand': r.strand, 'product': r['product'],
                 'role': ('HC:' + hcat[r.gene_id]) if r.gene_id in hcat else 'IS/integrase' if rx_is.search(r['product']) else 'conjugation' if rx_conj.search(r['product']) else ''}
                for _, r in x.iterrows()]
J['plasmid_maps'] = maps
json.dump(J, open(H + 'hgt.json', 'w'), separators=(',', ':'), default=lambda o: None if o != o else (o.item() if hasattr(o, 'item') else str(o)))
print('rewrote with extras', round(len(open(H + 'hgt.json').read()) / 1e3), 'kB')
