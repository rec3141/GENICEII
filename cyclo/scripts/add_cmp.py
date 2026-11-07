"""Add reference comparison (ANI, AAI, isolation source, HC homology, cassette synteny vs C. zancles) to cyclo.json."""
import json, pandas as pd, numpy as np
C = '/data/dev/rnaseq-exploration/cyclo/'
j = json.load(open(C + 'cyclo.json'))
refs = {'GCA_047634535.1': ('Cycloclasticus sp047634535 (MAG)', 'Central Arctic Ocean sea-ice core, 87.3N (2020)'),
        'GCA_033791505.1': ('Cycloclasticus sp033791505 (MAG)', 'Kivalliq, Nunavut marine sediment + crude-oil microcosm (2020)'),
        'GCA_019090295.1': ('Cycloclasticus sp019090295 (MAG)', 'Douglas Channel BC, dilbit-coated beads in seawater, 360 m (2014)'),
        'GCA_913060865.1': ('Cycloclasticus sp001566825 (MAG)', 'Canadian marine water metagenome (2013)'),
        'GCF_000442595.1': ('C. zancles 78-ME (isolate)', 'complete genome + plasmid p7ME01'),
        'GCF_000384415.1': ('C. pugetii PS-1 (type strain)', 'Puget Sound sediment isolate'),
        'GCF_000299965.1': ('Cycloclasticus sp. P1 (isolate)', 'deep-sea sediment isolate'),
        'GCF_000444935.1': ('C. pugetii PY97M (isolate)', 'isolate'), 'GCF_002443175.1': ('Cycloclasticus sp. PY97N (isolate)', 'isolate')}
sk = pd.read_csv(C + 'skani_refs.tsv', sep='\t', skiprows=1, header=None) if False else None
import subprocess
rows = [l.rstrip('\n').split('\t') for l in open(C + 'skani_refs_local.tsv')][1:]
names = [r[0].split('/')[-1].replace('.fa', '').replace('.fna', '') for r in rows]
M = pd.DataFrame([[float(x) for x in r[1:]] for r in rows], index=names, columns=names)
aai = pd.read_csv(C + 'cmp/aai.tsv', sep='\t', header=None, names=['bin', 'ref', 'n', 'aai'])
hv = pd.read_csv(C + 'cmp/hc_vs_refs.tsv', sep='\t')
cz = pd.read_csv(C + 'cmp/cass_zancles.tsv', sep='\t', header=None, names=['gene_id', 'sseq', 'pid', 'len', 'qlen', 'sstart', 'send', 'evalue'])
cz = cz.sort_values(['gene_id', 'pid'], ascending=[True, False]).drop_duplicates('gene_id')
j['references'] = dict(
    refs=[dict(acc=a, name=v[0], source=v[1]) for a, v in refs.items()],
    ani=dict(names=M.index.tolist(), values=M.round(2).values.tolist(), note='skani triangle; 0 = below ~80% ANI detection'),
    aai=aai.to_dict('records'),
    hc_homology=hv.round(1).replace({np.nan: None}).to_dict('records'),
    cassette_vs_zancles=cz[['gene_id', 'pid', 'sstart', 'send']].to_dict('records'))
json.dump(j, open(C + 'cyclo.json', 'w'), separators=(',', ':'))
print('ok', M.shape)

# operon onset after diesel + alkane-oxidation gene expression
j = json.load(open(C + 'cyclo.json'))
j['operon_onset'] = pd.read_csv(C + 'operon_onset.tsv', sep='\t').replace({np.nan: None}).to_dict('records')
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import samples, meta
G = pd.read_csv(C + 'cyclo_genes_annot.tsv', sep='\t', keep_default_na=False); cpm = np.load(C + 'cyclo_cpm.npy')
S = samples(); m = meta().loc[S]
sel = G[G['product'].str.contains('Alkane 1-mono|Cytochrome P450|Methane monooxygenase', regex=True) & (G.cbin != 'unbinned_kaiju_cyclo')]
al = []
for i, r in sel.iterrows():
    d = dict(gene_id=r.gene_id, bin=r.cbin, contig=r.contig, product=r['product'], ko=r.KO)
    for se in ['2024-25', '2025-26']:
        o = ((m.season == se) & (m.phase == 'oiled') & m.medium.isin(['SW', 'SI'])).values; c = ((m.season == se) & (m.phase == 'control')).values
        d[f'cpm_oiled_{se}'] = round(float(cpm[i, o].mean()), 2); d[f'cpm_ctl_{se}'] = round(float(cpm[i, c].mean()), 2)
    al.append(d)
j['alkane_genes'] = al
json.dump(j, open(C + 'cyclo.json', 'w'), separators=(',', ':'))
print('onset+alkane added')
