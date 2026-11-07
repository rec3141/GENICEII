"""site/islands.json: island table (MAG-binned islands), function enrichment, and gene maps for oil-gene, induced and highly mobile islands."""
import json, os, numpy as np, pandas as pd
I = '/data/dev/rnaseq-exploration/islands/out/'
GI = pd.read_csv(I + 'islands.tsv', sep='\t', keep_default_na=False, low_memory=False)
G = pd.read_csv(I + 'island_genes.tsv', sep='\t', keep_default_na=False, low_memory=False)
num = lambda d, c: pd.to_numeric(d[c], errors='coerce')
for c in GI.columns:
    if c.startswith(('cpm', 'log2', 'host', 'activity', 'dGC', 'plasmid')): GI[c] = num(GI, c)
for c in G.columns:
    if c.startswith('cpm'): G[c] = num(G, c)
GI['core_like'] = GI.core_like.astype(str).eq('True')
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
CLS = ['hydrocarbon degradation', 'integrase', 'transposase/IS', 'phage', 'conjugation/T4SS', 'restriction-modification', 'CRISPR/defence',
       'metal resistance', 'antibiotic resistance', 'EPS/cell surface', 'secretion/adhesion', 'transport', 'regulation', 'other enzyme/function', 'hypothetical']
S = ('2024-25', '2025-26')
def induced(x, s):
    return (x[f'cpm_oil_{s}'] >= 10) & (x[f'host_cpm_ctl_{s}'] >= 50) & (x[f'log2_share_oil_ctl_{s}'] >= 2) & ~x.core_like
sel = (GI.hc_genes > 0) | induced(GI, S[0]) | induced(GI, S[1]) | ((GI.mobility_score >= 3) & (GI.bin != 'unbinned') & ~GI.core_like) | \
      ((GI.bin != 'unbinned') & ~GI.core_like & ((GI.get('phage', 0) >= 5) | (GI.get('metal resistance', 0) >= 2) | (GI.get('CRISPR/defence', 0) >= 3)))
detail = set(GI.island_id[sel])
rows = []
for _, x in GI[(GI.bin != 'unbinned') | GI.island_id.isin(detail)].iterrows():
    rows.append({'id': x.island_id, 'contig': x.contig, 'start': int(x.start), 'end': int(x.end), 'len': int(x.len), 'bin': x.bin.replace('dastool-', ''),
                 'taxon': x.taxon, 'family': x.family, 'n': int(x.n_cds), 'dGC': r(x.dGC, 1), 'tRNA': str(x.tRNA_edge) == 'True',
                 'plasmid': r(x.plasmid_score), 'provirus': str(x.provirus) == 'True', 'integron': str(x.integron) == 'True', 'mob': int(x.mobility_score),
                 'core': bool(x.core_like), 'hc': int(x.hc_genes), 'hc_cat': x.hc_categories,
                 'cls': {k: int(float(x[k])) for k in CLS if k in GI.columns and float(x[k] or 0) > 0},
                 'oil': [r(x[f'cpm_oil_{s}'], 1) for s in S], 'ctl': [r(x[f'cpm_ctl_{s}'], 1) for s in S],
                 'lfc': [r(x[f'log2_oil_ctl_{s}'], 1) for s in S], 'share': [r(x[f'log2_share_oil_ctl_{s}'], 1) for s in S],
                 'act': [r(x[f'activity_oil_{s}'], 1) for s in S], 'detail': x.island_id in detail})
genes = {}
for gid, g in G[G.island_id.isin(detail)].groupby('island_id'):
    genes[gid] = [{'s': int(z['start']), 'e': int(z['stop']), 'st': z['strand'], 't': z['type'], 'g': z['gene'], 'p': z['product'], 'c': z['class'], 'hc': z['hc_category'] or z['canthyd'],
                   'ko': z['KO'], 'oil': [r(z[f'cpm_oil_{s}'], 1) for s in S], 'ctl': [r(z[f'cpm_ctl_{s}'], 1) for s in S]} for z in g.to_dict('records')]
en = pd.read_csv(I + 'function_enrichment.tsv', sep='\t', index_col=0)
J = {'stats': {'islands': len(GI), 'in_mags': int((GI.bin != 'unbinned').sum()), 'mags': int(GI.bin[GI.bin != 'unbinned'].nunique()),
               'core_like': int(GI.core_like.sum()), 'hc_islands': int((GI.hc_genes > 0).sum()), 'hc_genes': int(GI.hc_genes.sum()),
               'hc_total_confident': None, 'median_len': int(GI.len.median())},
     'enrichment': [{'cls': k, 'inside': int(x.inside), 'outside': int(x.outside), 'enrich': r(x.enrich, 2)} for k, x in en.iterrows()],
     'islands': rows, 'genes': genes}
hc = pd.read_csv('/data/dev/rnaseq-exploration/hgt/hc_genes.tsv', sep='\t', keep_default_na=False)
J['stats']['hc_total_confident'] = int((hc.tier == 'confident').sum())
json.dump(J, open('/data/dev/rnaseq-exploration/site/islands.json', 'w'), separators=(',', ':'), allow_nan=False)
print(len(rows), 'islands in json;', len(genes), 'with gene maps;', os.path.getsize('/data/dev/rnaseq-exploration/site/islands.json') // 1000, 'kB')
