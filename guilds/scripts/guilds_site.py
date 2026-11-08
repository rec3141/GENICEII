"""site/guilds.json from guilds/out (guilds.py): genome fingerprints, both embeddings, factor scores and driver links."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
E0 = '/data/dev/rnaseq-exploration/'; O = E0 + 'guilds/out/'
def r(x, n=3):
    try: x = float(x)
    except Exception: return None
    return round(x, n) if np.isfinite(x) else None
F = pd.read_csv(O + 'fingerprints.tsv', sep='\t', index_col=0)
FS = pd.read_csv(O + 'factors_samples.tsv', sep='\t', index_col=0); FD = pd.read_csv(O + 'factor_drivers.tsv', sep='\t', index_col=0)
FN = pd.read_csv(O + 'function_pc1.tsv', sep='\t', index_col=0); MD = json.load(open(O + 'model.json'))
DR = ['winter', 'ice', 'diesel', 'time', 'light']; FAC = [c for c in FS.columns if c.startswith('F') and c[1:].isdigit()]
def fname(f):
    r = FD[f].sort_values(key=abs, ascending=False)
    a = [f"{k} ({'+' if v > 0 else '−'})" for k, v in r.items() if abs(v) >= 0.4][:2]
    return ' and '.join(a) if a else f"mixed, weak (strongest: {r.index[0]}, |r| {abs(r.iloc[0]):.2f})"
fac_name = {f: fname(f) for f in FAC}
F['guild_name'] = [f'{g}: {fac_name[g]}' for g in F.guild]
m = meta()
# phylum palette: the 7 commonest phyla in fixed order, the rest folded into 'other'
ph = F.phylum.value_counts(); top = list(ph.index[:7])
F['phylum_c'] = F.phylum.where(F.phylum.isin(top), 'other')
G = []
for b, x in F.iterrows():
    fn = FN.loc[b].dropna().sort_values(key=abs, ascending=False) if b in FN.index else pd.Series(dtype=float)
    G.append({'bin': b.replace('dastool-', ''), 'taxon': x.taxon, 'phylum': x.phylum, 'phylum_c': x.phylum_c, 'tpm': r(x.tpm, 1), 'n': int(x.n), 'single': bool(x.single_winter),
              'r2': r(x.r2_full), 'u': {d: r(x[f'u_{d}']) for d in DR}, 'q': {d: r(x[f'q_{d}'], 4) for d in DR}, 'seasonal': r(x.seasonal), 'other': r(x.other_shared),
              'dom': x.dominant, 'pc1_driver': x.pc1_driver if isinstance(x.pc1_driver, str) else None, 'fc': x.fcluster_name, 'guild': x.guild_name,
              'fx': r(x.fx), 'fy': r(x.fy), 'gx': r(x.gx), 'gy': r(x.gy), 'load': {f: r(x[f'load_{f}']) for f in FAC},
              'fn': [[k, r(v, 2)] for k, v in fn.head(6).items()]})
S = [s for s in FS.index]

J = {'model': {'K': MD['K'], 've': {k: r(v) for k, v in MD['ve'].items()}, 'stab': {k: r(v, 2) for k, v in MD['stab'].items()}, 'fmethod': MD['fmethod'], 'gmethod': MD['gmethod'],
               'n_genomes': MD['n_genomes'], 'n_samples': MD['n_samples'], 'nperm': MD['nperm']},
     'genomes': G, 'phyla': top + ['other'],
     'clusters': list(F.fcluster_name.value_counts().index), 'guilds': list(F.guild_name.value_counts().index), 'dominants': list(F.dominant.value_counts().index),
     'factors': [{'id': f, 'top': fac_name[f], 'r': {d: r(FD.loc[d, f], 2) for d in FD.index}} for f in FAC],
     'samples': [{'id': s, 'name': m.sample_name.get(s, s), 'season': FS.season[s], 'pool': FS.pool[s], 'medium': FS.medium[s], 'day': None if pd.isna(FS.day[s]) else int(FS.day[s])} for s in S],
     'scores': {f: [r(v, 3) for v in FS[f]] for f in FAC},
     'sig': {d: [int((F[f'q_{d}'] < 0.05).sum()), int(F[f'q_{d}'].notna().sum())] for d in DR},
     'med': {d: r(F[f'u_{d}'].median()) for d in DR} | {'seasonal': r(F.seasonal.median()), 'r2': r(F.r2_full.median())},
     'xtab': pd.crosstab(F.fcluster_name, F.guild_name).to_dict()}
# notable pairs: (a) same fingerprint cluster, different guild; (b) same guild, different phylum, near-identical loadings
L = F[[f'load_{f}' for f in FAC]].values; nrm = L / np.linalg.norm(L, axis=1, keepdims=True)
cos = nrm @ nrm.T; idx = list(F.index); ab = F.tpm.fillna(0).values
pa, pb = [], []
for i in range(len(idx)):
    for j in range(i + 1, len(idx)):
        a, b = F.iloc[i], F.iloc[j]
        if a.guild_name == b.guild_name and a.phylum.replace(' (Kaiju)', '') != b.phylum.replace(' (Kaiju)', '') and 'unclassified' not in (a.phylum, b.phylum) and cos[i, j] >= 0.97:
            pb.append((min(ab[i], ab[j]), idx[i], idx[j], cos[i, j]))
        if a.fcluster_name == b.fcluster_name and a.guild_name != b.guild_name:
            dfp = np.abs(F.loc[idx[i], ['u_' + d for d in DR]].fillna(0).values - F.loc[idx[j], ['u_' + d for d in DR]].fillna(0).values).sum()
            if dfp < 0.08: pa.append((min(ab[i], ab[j]), idx[i], idx[j], dfp))
def pair(t, kind):
    _, a, b, v = t; A, B = F.loc[a], F.loc[b]
    return {'a': a.replace('dastool-', ''), 'b': b.replace('dastool-', ''), 'ta': A.taxon, 'tb': B.taxon, 'pa': A.phylum, 'pb': B.phylum, 'ga': A.guild_name, 'gb': B.guild_name,
            'ca': A.fcluster_name, 'cb': B.fcluster_name, 'v': r(v, 3), 'kind': kind}
J['pairs_same_fp'] = [pair(t, 'fp') for t in sorted(pa, reverse=True)[:8]]
J['pairs_same_guild'] = [pair(t, 'guild') for t in sorted(pb, reverse=True)[:8]]
json.dump(J, open(E0 + 'site/guilds.json', 'w'), separators=(',', ':'), allow_nan=False, ensure_ascii=False)
print(len(G), 'genomes;', os.path.getsize(E0 + 'site/guilds.json') // 1000, 'kB;', len(J['pairs_same_fp']), 'fp pairs,', len(J['pairs_same_guild']), 'guild pairs')
print(json.dumps(J['factors'], ensure_ascii=False)); print(J['sig'], J['med'])
