"""site/viruses.json: viral DNA share (nanopore + Illumina), viral transcription by class and lytic/lysogeny markers, richness/Shannon,
phage-host time courses, diesel-bloom vOTUs with candidate hosts, prophage table, ice-vs-water contrasts, virus-concentrate method test."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
import contrasts as IC
V = '/data/dev/rnaseq-exploration/viruses/'
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
J = {}
vc = pd.read_csv(V + 'out/viral_confident.tsv', sep='\t')
U = pd.read_pickle(V + 'out/votu_rna_cpm.pkl'); R = pd.read_csv(V + 'out/vrna_samples.tsv', sep='\t', index_col=0)
info = pd.read_csv(V + 'out/votu_rna.tsv', sep='\t', index_col=0)
J['stats'] = {'genomad_all': 272410, 'confident': len(vc), 'votus': int(vc.votu_rep.nunique()), 'active_votus': int((U.max(1) >= 1).sum()),
              'proviruses': 441, 'classes': vc.drop_duplicates('votu_rep').cls.replace('', 'unclassified').value_counts().head(8).to_dict(),
              'quality': vc.drop_duplicates('votu_rep').checkv_quality.value_counts().to_dict()}
# RNA per library (pool samples only; concentrates separately)
CL = ['Caudoviricetes', 'Megaviricetes', 'Maveriviricetes', 'Malgrandaviricetes', 'Faserviricetes']
MK = [k for k in R.columns if k.startswith('mk|')]
pool = R[(R.group != 'other') & R.medium.isin(['SW', 'SI'])]
J['rna'] = [{'id': i, 'name': x['name'], 'season': x.season, 'pool': x.group, 'medium': x.medium, 'day': r(x.day, 0), 'phase': x.phase,
             'layer': x.layer if isinstance(x.layer, str) else '', 'viral': r(x.viral_cpm, 0), 'rich': int(x.richness), 'H': r(x.shannon, 2),
             'cls': {c: r(x[c], 0) for c in CL if c in R.columns}, 'mk': {k[3:]: r(x[k], 0) for k in MK}} for i, x in pool.iterrows()]
# nanopore DNA
nd = pd.read_csv(V + 'out/vdna_samples.tsv', sep='\t', index_col=0)
nd = nd[(nd.group != 'other') & nd.day.notna()]
J['dna_np'] = [{'name': x.sample_name, 'season': x.season, 'pool': x.group, 'medium': x.medium, 'day': int(x.day), 'viral': r(x.viral_pct, 2),
                'cls': {c: r(x[c], 2) for c in CL if c in nd.columns}} for _, x in nd.iterrows()]
# Illumina DNA (if mapped)
if os.path.exists(V + 'out/vdna_illumina.tsv'):
    il = pd.read_csv(V + 'out/vdna_illumina.tsv', sep='\t', index_col=0)
    J['dna_il'] = [{'name': x['name'], 'season': x.season, 'pool': x.pool, 'medium': x.medium, 'day': r(x.day, 0), 'viral': r(x.viral_pct, 2),
                    'rich': int(x.richness), 'H': r(x.shannon, 2), 'cls': {c: r(x.get(c), 2) for c in CL if c in il.columns}} for _, x in il.iterrows()]
# phage-host courses for the key bloom vOTUs
H = pd.read_csv(V + 'out/votu_hosts.tsv', sep='\t', index_col=0)
B = table('bin_tpm')
m = meta().loc[U.columns]
key = [('contig_0040828', 'dastool-lorbin_722_sub', 'Cycloclasticus phage (2025-26)'), ('contig_0069512', 'dastool-lorbin_722_sub', 'Cycloclasticus phage 2, Schitoviridae (2025-26)'),
       ('contig_0010334', 'dastool-comebin_1134', 'Diesel-pool phage of GCA-2708415 (2024-25)'), ('contig_0016704', 'dastool-comebin_1134', 'Zobellviridae phage (2024-25)'),
       ('contig_0049668', None, 'Virophage, control pool (2025-26)'), ('contig_0000564', None, 'Giant virus contig (2025-26)')]
x = m[(m.group != 'other') & (m.medium == 'SW')]
J['pairs'] = []
for v, h, lab in key:
    if v not in U.index: continue
    J['pairs'].append({'votu': v, 'label': lab, 'host': h.replace('dastool-', '') if h else None, 'len': r(info.loc[v, 'length'], 0) if v in info.index else None,
                       'pts': [{'season': x.season[i], 'pool': x.group[i], 'day': int(x.day[i]), 'phage': r(U.loc[v, i], 1), 'host': r(B.loc[h, i], 0) if h else None} for i in x.index]})
J['hosts'] = [{'votu': v, **{k: (r(y[k], 1) if isinstance(y[k], (float, np.floating)) else (y[k] if isinstance(y[k], str) else None)) for k in
               ['cls', 'length', 'checkv_quality', 'oil_2024-25', 'ctl_2024-25', 'oil_2025-26', 'ctl_2025-26', 'baci_2025-26', 'host1', 'rho1', 'host2', 'rho2', 'binned']}}
              for v, y in H.iterrows()]
# prophages
P = pd.read_csv(V + 'out/prophage.tsv', sep='\t', keep_default_na=False)
for c in P.columns:
    if c.startswith(('cpm', 'share', 'activity', 'baci', 'host_cpm')): P[c] = pd.to_numeric(P[c], errors='coerce')
P = P[P.bin != ''].sort_values('cpm_max', ascending=False)
J['prophage'] = [{'id': y.provirus, 'host': y.host, 'bin': y.bin.replace('dastool-', ''), 'kb': r(y.length / 1000, 1), 'n': int(y.n_genes),
                  'oil': [r(y['cpm_oil_2024-25'], 1), r(y['cpm_oil_2025-26'], 1)], 'ctl': [r(y['cpm_ctl_2024-25'], 1), r(y['cpm_ctl_2025-26'], 1)],
                  'baci': r(y['baci_share_2025-26'], 1), 'act': [r(y['activity_oil_2024-25'], 1), r(y['activity_oil_2025-26'], 1)]} for _, y in P.head(40).iterrows()]
J['prophage_n'] = {'all': 441, 'in_mags': int(len(P))}
# ice vs water / top vs bottom for viral metrics (pairs from ice/out/pairs.tsv)
L = pd.DataFrame({'all viral': np.log2(R.viral_cpm + 1), 'active vOTUs': np.log2(R.richness + 1), **{c: np.log2(R[c] + 1) for c in CL[:3]}}).T
L = L[[c for c in L.columns if c in U.columns]]
iw, _ = IC.contrast(L, IC.PI); tb, _ = IC.contrast(L, IC.PT)
J['ice'] = [{'k': k, 'iw': r(iw.loc[k, 'delta']), 'iw_up': int(iw.loc[k, 'n_up']), 'iw_n': int(iw.loc[k, 'n']), 'tb': r(tb.loc[k, 'delta']),
             'tb_up': int(tb.loc[k, 'n_up']), 'tb_n': int(tb.loc[k, 'n'])} for k in L.index]
# virus-concentrate method test
vi = R.index[R.name.str.startswith('VIRUS', na=False)]
ref = R.index[(R.season == '2025-26') & (R.group == 'CTL') & (R.medium == 'SW') & (R.day >= 49)]
enr = (np.log2(U[vi].median(1) + 1) - np.log2(U[ref].median(1) + 1))
top = info.assign(lfc=enr, conc=U[vi].median(1), sw=U[ref].median(1)).query('conc >= 10').sort_values('lfc', ascending=False).head(10)
from scipy.spatial.distance import pdist
P2 = U[vi] / U[vi].sum()
J['conc'] = {'libs': [{'name': R.name[i], 'viral': r(R.viral_cpm[i], 0), 'rich': int(R.richness[i]), 'H': r(R.shannon[i], 2), 'cls': {c: r(R.loc[i, c], 0) for c in CL if c in R.columns}} for i in vi],
             'ref': [{'name': R.name[i], 'viral': r(R.viral_cpm[i], 0), 'rich': int(R.richness[i]), 'H': r(R.shannon[i], 2), 'cls': {c: r(R.loc[i, c], 0) for c in CL if c in R.columns}} for i in ref],
             'bc_within': r(float(np.median(pdist(P2.T.values, 'braycurtis'))), 2),
             'bc_ref': r(float(np.median([pdist(np.vstack([P2[a].values, (U[b] / U[b].sum()).values]), 'braycurtis')[0] for a in vi for b in ref])), 2),
             'enriched': [{'votu': v, 'cls': y.cls if isinstance(y.cls, str) else '', 'fam': y.fam if isinstance(y.fam, str) else '', 'kb': r(y.length / 1000, 1),
                           'conc': r(y.conc, 0), 'sw': r(y.sw, 1), 'lfc': r(y.lfc, 1)} for v, y in top.iterrows()]}
json.dump(J, open('/data/dev/rnaseq-exploration/site/viruses.json', 'w'), separators=(',', ':'), allow_nan=False)
print({k: (len(v) if hasattr(v, '__len__') else v) for k, v in J.items()}, os.path.getsize('/data/dev/rnaseq-exploration/site/viruses.json') // 1000, 'kB')
