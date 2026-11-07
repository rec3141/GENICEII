"""vOTU-level viral transcription: per-library viral CPM (by class), active vOTU richness (vOTUs with >= 1 CPM), Shannon of the viral
transcript pool, Bray-Curtis turnover, lytic/lysogeny marker classes; diesel vs control (SW, after dosing) and 2025-26 BACI per vOTU;
ice vs water and top vs bottom pairs (ice/out/pairs.tsv). Out: out/vrna_samples.tsv, out/votu_rna.tsv, out/votu_rna_cpm.pkl"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
V = '/data/dev/rnaseq-exploration/viruses/'
d = pd.read_pickle(V + 'out/vrna_gene_cpm.pkl'); Gv = d['genes']; S = d['samples']
C = pd.DataFrame(d['cpm'], columns=S)
vc = pd.read_csv(V + 'out/viral_confident.tsv', sep='\t')
rep = vc.set_index('region').votu_rep
U = C.groupby(Gv.contig.map(rep).values).sum()          # vOTU x sample CPM
info = vc[vc.region == vc.votu_rep].set_index('votu_rep')[['cls', 'fam', 'length', 'checkv_quality', 'completeness', 'bin', 'taxonomy']]
m = meta().loc[S]
P = U / U.sum()
H = -(P * np.log(P.where(P > 0))).sum()
rec = pd.DataFrame({'viral_cpm': U.sum(), 'richness': (U >= 1).sum(), 'shannon': H, 'season': m.season, 'group': m.group,
                    'medium': m.medium, 'day': m.day, 'phase': m.phase, 'name': m.sample_name, 'layer': m.ice_layer})
cls = info.cls.reindex(U.index).fillna('').replace('', 'unclassified')
for k, v in U.groupby(cls.values).sum().iterrows(): rec[k] = v
mk = pd.read_csv(V + 'out/vrna_classes.tsv', sep='\t', index_col=0)
for k in mk.index: rec['mk|' + k] = mk.loc[k, S].values
rec.to_csv(V + 'out/vrna_samples.tsv', sep='\t')
pd.to_pickle(U.astype(np.float32), V + 'out/votu_rna_cpm.pkl')
# per-vOTU diesel vs control and BACI
oil, ctl, pre = m.phase.eq('oiled'), m.phase.eq('control'), m.phase.eq('pre')
out = info.reindex(U.index).copy(); out['mean_cpm'] = U.mean(1); out['max_cpm'] = U.max(1)
for s in ('2024-25', '2025-26'):
    a = m.index[oil & (m.season == s) & (m.medium == 'SW')]; c = m.index[ctl & (m.season == s) & (m.medium == 'SW')]
    out[f'oil_{s}'] = U[a].median(1); out[f'ctl_{s}'] = U[c].median(1)
    out[f'lfc_{s}'] = np.log2((out[f'oil_{s}'] + 0.5) / (out[f'ctl_{s}'] + 0.5))
pe = m.index[pre & (m.season == '2025-26') & (m.medium == 'SW') & (m.group == 'EXP')]; pc = m.index[pre & (m.season == '2025-26') & (m.medium == 'SW') & (m.group == 'CTL')]
out['baci_2025-26'] = np.log2((out['oil_2025-26'] + 0.5) / (U[pe].median(1) + 0.5)) - np.log2((out['ctl_2025-26'] + 0.5) / (U[pc].median(1) + 0.5))
out.to_csv(V + 'out/votu_rna.tsv', sep='\t')
pd.set_option('display.width', 250)
x = rec[(rec.group != 'other') & rec.medium.isin(['SW', 'SI'])]
print(x.groupby(['season', 'group', 'medium', 'phase'])[['viral_cpm', 'richness', 'shannon']].median().round(2).to_string())
print('vOTUs active (>=1 CPM in >=1 library):', (U.max(1) >= 1).sum(), 'of', len(U))
c = ['cls', 'fam', 'length', 'checkv_quality', 'bin', 'oil_2025-26', 'ctl_2025-26', 'baci_2025-26', 'oil_2024-25', 'ctl_2024-25', 'lfc_2024-25']
print(out[(out['oil_2025-26'] >= 20)].sort_values('baci_2025-26', ascending=False)[c].head(12).round(1).to_string())
print(out[(out['oil_2024-25'] >= 20)].sort_values('lfc_2024-25', ascending=False)[c].head(12).round(1).to_string())
