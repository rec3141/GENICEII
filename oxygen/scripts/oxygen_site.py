"""site/oxygen.json from oxygen/out (oxygen.py) plus per-genome oxidase/anaerobic-marker shares for Cycloclasticus lorbin_184 and
the genomes behind the 2025-26 high-affinity-oxidase rise."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
O = '/data/dev/rnaseq-exploration/oxygen/out/'; E0 = '/data/dev/rnaseq-exploration/'
def r(x, n=2):
    try: x = float(x)
    except (TypeError, ValueError): return None
    return round(x, n) if np.isfinite(x) else None
I = pd.read_csv(O + 'sample_index.tsv', sep='\t', index_col=0)
E = pd.read_csv(O + 'group_cpm.tsv', sep='\t', index_col=0)
CT = pd.read_csv(O + 'contrasts.tsv', sep='\t'); CV = pd.read_csv(O + 'covariates.tsv', sep='\t'); GS = pd.read_csv(O + 'genome_shift.tsv', sep='\t')
IDX = ['comm_HA_LA', 'comm_cbb3_aa3', 'comm_bd_aa3', 'comm_hemN_hemF', 'comm_nrdDG_nrdAB', 'comm_anaer_cpm', 'within_HA_LA', 'within_HA_LA_noCyclo', 'regulation_HA_logodds', 'composition_HA_logodds']
S = list(I.index)
J = {'samples': [{'id': s, 'name': I.sample_name[s], 'season': I.season[s], 'pool': I.group[s], 'medium': I.medium[s], 'phase': I.phase[s],
                  'day': None if pd.isna(I.day[s]) else int(I.day[s]), **{k: r(I.loc[s, k]) for k in IDX + ['cyclo_tpm']}} for s in S],
     'contrasts': [{'contrast': x.contrast, 'index': x['index'], 'n': None if pd.isna(x.n) else int(x.n), 'mean': r(x['mean']), 'n_up': None if pd.isna(x.n_up) else int(x.n_up), 'p': r(x.p, 4)} for _, x in CT.iterrows()],
     'covariates': [{'subset': x.subset, 'covariate': x.covariate, 'index': x['index'], 'n': int(x.n), 'rho': r(x.rho), 'p': r(x.p, 4)} for _, x in CV[CV.p < 0.05].iterrows()],
     'genome_shift': [{'season': x.season, 'medium': x.medium, 'bin': x.bin.replace('dastool-', ''), 'taxon': x.taxon, 'shift': r(x['shift']), 'q': r(x.q, 3), 'n_oil': int(x.n_oil), 'n_ctl': int(x.n_ctl)}
                      for _, x in GS.sort_values('q').iterrows() if x.q < 0.25]}
m = meta().loc[S]
# gene groups: medians per season x pool phase (seawater)
gg = []
for g in E.index:
    row = {'group': g}
    for s in ('2024-25', '2025-26'):
        for k, ph, grp in (('pre_exp', 'pre', 'EXP'), ('pre_ctl', 'pre', 'CTL'), ('oil', 'oiled', 'EXP'), ('ctl', 'control', 'CTL')):
            ix = m.index[(m.season == s) & (m.medium == 'SW') & (m.phase == ph) & (m.group == grp)]
            row[f'{k}{s[2:4]}'] = r(E.loc[g, ix].median(), 1)
    gg.append(row)
J['groups'] = gg
# who carries the high-affinity oxidase transcripts (2025-26 seawater)
X, G = genes(); cds = G['type'].values == 'CDS'; tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel(); kos = G['KO'].fillna('').str.split(',')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def blab(b):
    if not b: return 'unbinned'
    t = [y.split('__', 1)[1] for y in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in y and y.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')
SETS = {'HA': {'K00404', 'K00405', 'K00406', 'K00407', 'K00425', 'K00426', 'K22501'}, 'LA': {'K02274', 'K02275', 'K02276', 'K02297', 'K02298', 'K02299', 'K02300'},
        'narGHI': {'K00370', 'K00371', 'K00374'}, 'hox': {'K00436', 'K18005', 'K18006', 'K18007'}, 'cbb3': {'K00404', 'K00405', 'K00406', 'K00407'}, 'aa3': {'K02274', 'K02275', 'K02276'}}
o = m.index[(m.season == '2025-26') & (m.medium == 'SW') & (m.phase == 'oiled')]; c = m.index[(m.season == '2025-26') & (m.medium == 'SW') & (m.phase == 'control')]
J['carriers'] = {}
for nm in ('HA', 'LA'):
    ix = np.where(kos.apply(lambda x: bool(SETS[nm].intersection(x))).values & cds)[0]
    Cm = pd.DataFrame(X[ix].toarray() / tot * 1e6, columns=m.index); b = G['bin'].iloc[ix].fillna('').values
    t = pd.DataFrame({'oil': Cm[o].groupby(b).sum().median(1), 'ctl': Cm[c].groupby(b).sum().median(1)}).sort_values('oil', ascending=False)
    J['carriers'][nm] = {'total_oil': r(Cm[o].sum().median(), 0), 'total_ctl': r(Cm[c].sum().median(), 0),
                         'rows': [{'bin': k.replace('dastool-', '') or 'unbinned', 'taxon': blab(k), 'oil': r(v.oil, 0), 'ctl': r(v.ctl, 0)} for k, v in t.head(8).iterrows()]}
# Cycloclasticus lorbin_184: per-genome shares (per mille of its own transcripts)
BIN = pd.read_pickle(E0 + 'nitrogen/out/bin_cds_cpm.pkl'); bn = 'dastool-lorbin_184'
cy = []
for nm in ('cbb3', 'aa3', 'narGHI', 'hox'):
    ix = np.where(kos.apply(lambda x: bool(SETS[nm].intersection(x))).values & cds & (G['bin'].values == bn))[0]
    v = pd.Series(np.asarray(X[ix].sum(0)).ravel() / tot * 1e6, index=m.index) / BIN.loc[bn, m.index] * 1e3
    row = {'gene': nm}
    for k, ph in (('pre', 'pre'), ('oil', 'oiled'), ('ctl', 'control')):
        kk = m.index[(m.season == '2025-26') & (m.medium == 'SW') & (m.phase == ph) & (BIN.loc[bn, m.index] >= 2000)]
        row[k] = r(v[kk].median(), 2); row['n_' + k] = len(kk)
    cy.append(row)
J['cyclo'] = cy
mos = I.index[I.sample_name == 'CMO_EXP_DAY04_MOS'][0]; same = I.index[(I.sample_name.isin(['CMO_EXP_DAY04_SW', 'CMO_EXP_DAY04_SI']))]
J['mos'] = {'rows': [{'group': g, 'mos': r(E.loc[g, mos], 1), **{I.sample_name[s].split('_')[-1]: r(E.loc[g, s], 1) for s in same}} for g in E.index]}
J['genome_count'] = int(len(pd.read_csv(O + 'genome_ha.tsv', sep='\t')))
json.dump(J, open(E0 + 'site/oxygen.json', 'w'), separators=(',', ':'), allow_nan=False)
print(os.path.getsize(E0 + 'site/oxygen.json') // 1000, 'kB'); print(json.dumps(J['carriers'])[:800]); print(J['cyclo'])
