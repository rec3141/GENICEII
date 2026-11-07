"""Aggregate Illumina vOTU read counts (map/counts/*.tsv; first 2 M pairs = 4 M reads per library) into viral share of reads,
per-class share, vOTU richness (>= 2 reads) and Shannon. Out: out/vdna_illumina.tsv, out/votu_il_counts.pkl"""
import glob, os, numpy as np, pandas as pd
V = '/data/dev/rnaseq-exploration/viruses/'
S = pd.read_csv(V + 'map/samples.tsv', sep='\t').set_index('sample')
vc = pd.read_csv(V + 'out/viral_confident.tsv', sep='\t'); info = vc[vc.region == vc.votu_rep].set_index('votu_rep')
C = {}
for f in glob.glob(V + 'map/counts/*.tsv'):
    s = os.path.basename(f)[:-4]
    x = pd.read_csv(f, sep='\t', header=None, names=['v', 'n'])
    if len(x): C[s] = x.set_index('v').n
C = pd.DataFrame(C).fillna(0)
S = S.loc[[s for s in S.index if s in C.columns]]
cls = info.cls.reindex(C.index).fillna('').replace('', 'unclassified')
out = S[['name', 'season', 'pool', 'medium', 'day']].copy()
out['viral_pct'] = C[S.index].sum() / 4e6 * 100
for k, v in C.groupby(cls.values).sum().iterrows(): out[k] = v[S.index] / 4e6 * 100
P = C[S.index] / C[S.index].sum()
out['richness'] = (C[S.index] >= 2).sum(); out['shannon'] = -(P * np.log(P.where(P > 0))).sum()
out.to_csv(V + 'out/vdna_illumina.tsv', sep='\t'); pd.to_pickle(C, V + 'out/votu_il_counts.pkl')
pd.set_option('display.width', 220)
print(out.groupby(['season', 'pool', 'medium', pd.cut(out.day, [-10, 0, 5, 12, 25, 70])])[['viral_pct', 'Caudoviricetes', 'Megaviricetes', 'richness', 'shannon']].median().round(2).dropna(how='all').to_string())
