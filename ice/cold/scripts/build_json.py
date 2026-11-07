"""Collect the cold-adaptation results into ../cold.json for the website."""
import sys, json; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
from scipy.stats import spearmanr, fisher_exact
C = C_DIR
r = lambda x, n=3: None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
m = design(); S = samples()
E = pd.read_csv(C + 'cold_cpm.tsv', sep='\t', index_col=0)
P = pd.read_csv(C + 'cold_paired.tsv', sep='\t')
W = pd.read_csv(C + 'cold_within_summary.tsv', sep='\t', index_col=0)
x = m[(m.medium == 'SI') & (m.season == '2024-25')].copy()
x['thick'] = pd.to_numeric(x.ice_thickness_cm.str.extract(r'(\d+)')[0], errors='coerce'); x = x.dropna(subset=['thick'])
cats = []
for c in E.columns:
    d = {'name': c, 'n_genes': int((pd.read_csv(C + 'cold_gene_members.tsv', sep='\t', usecols=['category']).category == c).sum()) if False else None}
    for kind in ['ice', 'top']:
        q = P[(P.category == c) & (P.contrast == kind)].iloc[0]
        d[kind] = {'median_log2': r(q.median_log2), 'frac_up': r(q.frac_up), 'n': int(q.n_blocks), 'q25': r(q.q25), 'q75': r(q.q75)}
    if c in W.index: d['within_mag'] = {'n_mags': int(W.loc[c, 'n_mags']), 'median_log2': r(W.loc[c, 'median']), 'frac_up': r(W.loc[c, 'frac_up'])}
    rho, p = spearmanr(x.thick, E.loc[x.index, c]); d['thickness_rho'] = r(rho); d['thickness_p'] = r(p, 4)
    cats.append(d)
mem = pd.read_csv(C + 'cold_gene_members.tsv', sep='\t', usecols=['category'])
for d in cats: d['n_genes'] = int((mem.category == d['name']).sum())
K = pd.read_csv(C + 'cold_ko_paired.tsv', sep='\t', keep_default_na=False)
K = K[K.mean_cpm >= 2]
ko = [{k: (r(v) if isinstance(v, float) else v) for k, v in row.items()} for row in K.to_dict('records')]
smp = [{'id': s, 'name': m.loc[s, 'sample_name'], 'season': m.loc[s, 'season'], 'group': m.loc[s, 'group'], 'medium': m.loc[s, 'medium'],
        'layer': m.loc[s, 'layer'], 'day': None if pd.isna(m.loc[s, 'day']) else int(m.loc[s, 'day'])} for s in S]
cpm = {c: [r(v, 2) for v in E.loc[S, c]] for c in E.columns}
# IBP
I = pd.read_csv(C + 'ibp_context.tsv', sep='\t', keep_default_na=False)
dom = pd.read_csv(C + 'ibp_domain.tsv', sep='\t', index_col=0)
tx = pd.read_csv(C + 'ibp_taxa.tsv', sep='\t', keep_default_na=False).head(20)
eu = pd.read_csv(C + 'ibp_euk_taxa.tsv', sep='\t', keep_default_na=False)
eut = eu.assign(group=eu.group.replace('', 'no MarFERReT hit')).groupby(['group', 'cls']).agg(n=('gene_id', 'size'), cpm_SI=('cpm_SI', 'sum'), cpm_SW=('cpm_SW', 'sum')).reset_index().sort_values('cpm_SI', ascending=False).head(10)
nearest = pd.read_csv(C + 'ibp_euk_nearest.tsv', sep='\t')
pairs = pd.read_csv(C + 'ibp_pairs.tsv', sep='\t', keep_default_na=False)
hgt = pairs[pairs.interpretation.str.startswith('HGT')]
pc = pd.read_csv(C + 'ibp_context.tsv', sep='\t')
mob = {'n_genes': len(I), 'mobile_assoc_pct': r(pc.mobile_assoc.mean() * 100, 1), 'mge_borne_pct': r(pc.mge_borne.mean() * 100, 1),
       'IS_within_10kb_pct': r(pc.IS_within_10kb.mean() * 100, 1), 'plasmid_strict_n': int(pc.plasmid_strict.sum()), 'phage_n': int(pc.phage_contig.sum()),
       'ice_transcripts_mobile_pct': r(pc.cpm_SI[pc.mobile_assoc].sum() / pc.cpm_SI.sum() * 100, 1),
       'ice_transcripts_mge_pct': r(pc.cpm_SI[pc.mge_borne].sum() / pc.cpm_SI.sum() * 100, 1)}
funnel = {'pairs_ge90aa_diff_contig': int(len(pairs) + 0), 'by_cross': pairs.cross.value_counts().to_dict(),
          'by_interpretation': pairs.interpretation.value_counts().to_dict()}
A = pd.read_csv(C + 'grex/ibp_ava.m8', sep='\t', header=None, names=['q', 't', 'fid', 'alnlen', 'qcov', 'tcov', 'evalue', 'bits'])
topg = I.sort_values('cpm_SI', ascending=False).head(20)
ibp = {'domain': [{'domain': k, **{c: r(v, 2) for c, v in row.items()}} for k, row in dom.iterrows()],
       'taxa': [{c: (r(v, 2) if isinstance(v, float) else v) for c, v in row.items()} for row in tx.to_dict('records')],
       'top_genes': [{'gene_id': q.gene_id, 'bin': q.bin, 'genus': q.genus, 'domain': q.domain, 'plen': int(q.plen) if q.plen != '' else None, 'product': q['product'],
                      'cpm_SI': r(float(q.cpm_SI), 2), 'cpm_SW': r(float(q.cpm_SW), 2), 'ice_log2': r(float(q.ice_log2), 2), 'mobile_assoc': q.mobile_assoc == 'True' or q.mobile_assoc is True} for _, q in topg.iterrows()],
       'euk_taxa': [{c: (r(v, 2) if isinstance(v, float) else v) for c, v in row.items()} for row in eut.to_dict('records')],
       'euk_nearest_bacterial_aa': {'n_euk': int(len(nearest)), 'with_hit': int(nearest.bac.notna().sum()), 'median_ident': r(nearest.fid.median()), 'max_ident': r(nearest.fid.max())},
       'mobility': mob, 'sharing_funnel': funnel,
       'transfers': [{c: (r(v, 2) if isinstance(v, float) else v) for c, v in row.items() if c not in ('q_contig', 't_contig')} for row in hgt.to_dict('records')]}
# MpIBP homolog
mp = pd.read_csv(C + 'grex/mpibp.m8', sep='\t', header=None, names=['q', 's', 'pid', 'len', 'qs', 'qe', 'ss', 'se', 'ev', 'bits', 'slen'])
mpibp = {'mag': 'dastool-lorbin_161', 'species': 'Marinomonas primoryensis', 'contig': 'contig_0024173',
         'hits': [{'query': q.q.split('|')[1], 'gene': q.s, 'pident': r(q.pid, 1), 'aln_len': int(q.len), 'q_range': f'{q.qs}-{q.qe}', 'gene_len_aa': int(q.slen)}
                  for q in mp.itertuples() if q.s.startswith('JPFHOD_0022676')],
         'note': 'MpIBP/MpAFP homolog fragmented into 8 short ORFs (JPFHOD_002267650-657, ~1.8 kaa total vs ~1.5 MDa native), next to a retention-module protein; RII repeats collapsed and frameshifted in the assembly. Transcripts ~0.3 CPM in both ice and water.'}
# MAGs
D = pd.read_csv(C + 'mag_cold_table.tsv', sep='\t', index_col=0)
D = D.dropna(subset=['ice_delta'])
D['hab'] = np.where((D.ice_delta >= 0.75) & (D.ice_q < 0.05), 'ice', np.where((D.ice_delta <= -0.75) & (D.ice_q < 0.05), 'water', 'neutral'))
t = pd.crosstab(D.hab == 'ice', D.ibp_n > 0); orr, pp = fisher_exact(t)
prev = {h: r((D[D.hab == h].ibp_n > 0).mean() * 100, 1) for h in ['ice', 'neutral', 'water']}
prev.update({'n': D.hab.value_counts().to_dict(), 'odds_ratio_ice_vs_rest': r(orr, 2), 'p': float(f'{pp:.3g}')})
mcols = ['g', 's', 'f', 'p', 'ice_delta', 'top_delta', 'hab', 'ibp_n', 'Cold-shock proteins', 'Compatible solute synthesis', 'Compatible solute transport',
         'Fatty-acid desaturases', 'PUFA synthase (pfa)', 'EPS export', 'Carotenoid synthesis', 'arg_lys', 'gln_thr', 'polar_unch', 'gly', 'pro', 'acidic', 'gc', 'genome_mb', 'SCG_completeness']
mags = [{'bin': b, **{c: (r(v, 4) if isinstance(v, (float, np.floating)) else (None if v != v else v)) for c, v in row.items()}} for b, row in D[mcols].iterrows()]
T = pd.read_csv(C + 'mag_test.tsv', sep='\t')
IP = pd.read_csv(C + 'index_paired.tsv', sep='\t')
out = {'categories': cats, 'ko': ko, 'samples': smp, 'cpm': cpm, 'ibp': ibp, 'mpibp': mpibp, 'mags': mags, 'mag_tests': [{c: (r(v, 4) if isinstance(v, float) else v) for c, v in row.items()} for row in T.to_dict('records')],
       'ibp_prevalence': prev, 'expressed_proteome': [{c: (r(v, 5) if isinstance(v, float) else v) for c, v in row.items()} for row in IP.to_dict('records')],
       'genus_pairs': pd.read_csv(C + 'genus_pairs.tsv', sep='\t')[['genus', 'ice_mag', 'water_mag', 'ice_delta_ice', 'ice_delta_water', 'd_arg_lys', 'd_gln_thr', 'd_gc', 'd_Ice-binding proteins']].round(4).to_dict('records')}
json.dump(out, open(C + 'cold.json', 'w'), separators=(',', ':'), default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print({k: (len(v) if hasattr(v, '__len__') else v) for k, v in out.items()}, os.path.getsize(C + 'cold.json') if (os := __import__('os')) else '')
