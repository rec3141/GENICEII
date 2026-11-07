"""site/nitrogen.json from nitrogen/out (nitrogen.py): category table with hosts (MAG or Kaiju genus), time courses,
host-normalised N-acquisition response, nitrifier genomes (TPM time courses), denitrification step balance."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
O = '/data/dev/rnaseq-exploration/nitrogen/out/'; E0 = '/data/dev/rnaseq-exploration/'
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
summ = pd.read_csv(O + 'cat_summary.tsv', sep='\t', index_col=0)
E = pd.read_csv(O + 'cat_cpm.tsv', sep='\t', index_col=0)
G = pd.read_csv(O + 'genes.tsv', sep='\t', keep_default_na=False)
M = pd.read_pickle(O + 'gene_cpm.pkl'); mem = pd.read_pickle(O + 'members.pkl')
H = pd.read_csv(O + 'host_nstress.tsv', sep='\t'); H[['taxon', 'phylum']] = H[['taxon', 'phylum']].fillna('')
S = list(E.columns); m = meta().loc[S]
def kgen(l):
    l = [x.strip() for x in l.split(';')] if l else []
    if not l: return 'unbinned, unclassified'
    g = l[5] if len(l) > 5 and l[5] not in ('', 'NA') else (l[4] if len(l) > 4 else l[-1])
    return f'unbinned ({g}, Kaiju)'
bid = G.bin.str.replace('dastool-', '', regex=False)
G['host'] = np.where(G.bin == '', G.kaiju.map(kgen), np.where(G.taxon == bid, G.taxon, G.taxon + ' ' + bid))
oil = list(m.index[m.phase == 'oiled']); ctl = list(m.index[m.phase == 'control'])
J = {'groups': list(dict.fromkeys(summ.group)), 'categories': [], 'timecourse': {}}
for c in summ.index:
    x = summ.loc[c].copy()
    if x.n_genes == 0: continue
    idx = mem[c]; sub = G.iloc[idx]; mm = M.iloc[idx]
    hosts = {}
    for w, cols in (('diesel', oil), ('control', ctl), ('all', S)):
        v = mm[cols].sum(1).groupby(sub.host.values).sum(); t = v.sum()
        hosts[w] = [{'host': h, 'share': r(s / t, 3)} for h, s in v.sort_values(ascending=False).head(5).items()] if t > 0 else []
    z = lambda *v: all((not np.isfinite(t)) or t < 0.5 for t in v)   # both sides ~0 CPM: no ratio
    if z(x['oil_2024-25_SW'], x['ctl_2024-25_SW']): x['log2_2024-25_SW'] = np.nan
    if z(x['oil_2025-26_SW'], x['ctl_2025-26_SW']): x['log2_2025-26_SW'] = np.nan
    for se in ('2024-25', '2025-26'):
        if z(x[f'oil_{se}_SW'], x[f'ctl_{se}_SW'], x[f'pre_exp_{se}'], x[f'pre_ctl_{se}']): x[f'baci_{se}'] = np.nan
    J['categories'].append({'name': c, 'group': x.group, 'n_genes': int(x.n_genes), 'mean': r(x.mean_cpm, 1),
        'lfc24': r(x['log2_2024-25_SW']), 'lfc25': r(x['log2_2025-26_SW']), 'baci24': r(x['baci_2024-25']), 'baci25': r(x['baci_2025-26']),
        'oil24': r(x['oil_2024-25_SW'], 1), 'ctl24': r(x['ctl_2024-25_SW'], 1), 'oil25': r(x['oil_2025-26_SW'], 1), 'ctl25': r(x['ctl_2025-26_SW'], 1),
        'pre_exp24': r(x['pre_exp_2024-25'], 1), 'pre_ctl24': r(x['pre_ctl_2024-25'], 1), 'pre_exp25': r(x['pre_exp_2025-26'], 1), 'pre_ctl25': r(x['pre_ctl_2025-26'], 1),
        'ice': r(x.ice_log2), 'ice_up': int(x.ice_up), 'ice_n': int(x.ice_n), 'top': r(x.top_log2), 'top_up': int(x.top_up), 'top_n': int(x.top_n), 'hosts': hosts})
    J['timecourse'][c] = [r(v, 2) for v in E.loc[c, S]]
J['samples'] = [{'id': s, 'name': m.sample_name[s], 'season': m.season[s], 'pool': m.group[s], 'medium': m.medium[s],
                 'day': None if pd.isna(m.day[s]) else int(m.day[s])} for s in S]
# host-normalised N acquisition
J['nstress'] = [{'season': x.season, 'bin': x.bin.replace('dastool-', ''), 'taxon': x.taxon, 'phylum': x.phylum, 'host_oil': r(x.host_oil, 0), 'host_ctl': r(x.host_ctl, 0),
                 'share_oil': r(x.share_oil * 1e3, 2), 'share_ctl': r(x.share_ctl * 1e3, 2), 'lfc': r(x.log2_share)} for x in H.itertuples()
                if max(x.host_oil, x.host_ctl) >= 1000 and np.isfinite(x.log2_share)]
json.dump(J, open(E0 + 'site/nitrogen.json', 'w'), separators=(',', ':'), allow_nan=False)
print(len(J['categories']), 'categories;', len(J['nstress']), 'host rows;', os.path.getsize(E0 + 'site/nitrogen.json') // 1000, 'kB')

# genomes: nitrifiers and copper-monooxygenase (CuMMO) hydrocarbon degraders
GB = [('metabat_9833', 'Nitrosopumilus (ammonia-oxidising archaeon)'), ('comebin_369', 'Nitrosomonadaceae GCA-2721545 (ammonia-oxidising bacterium)'),
      ('maxbin_568', 'LS-NOB (nitrite oxidiser)'), ('comebin_2576', 'Porticoccaceae 50-400-T64 (CuMMO)'), ('metabat_9342', 'Cycloclasticus (CuMMO)')]
tpm = table('bin_tpm'); BIN = pd.read_pickle(O + 'bin_cds_cpm.pkl'); baci = pd.read_csv(E0 + 'out/mag_baci.tsv', sep='\t')
core = m[m.group.isin(['CTL', 'EXP']) & (m.medium == 'SW')]
J['genomes'] = []
for b, lab in GB:
    d = 'dastool-' + b; cu = G.gene_id[(G.bin == d) & G.KO.str.contains('K1094[456]')]
    row = {'bin': b, 'label': lab, 'tpm': [r(v, 1) for v in tpm.loc[d, S]], 'summary': []}
    for s in ('2024-25', '2025-26'):
        bb = baci[(baci.bin == d) & (baci.season == s)]
        o = {'season': s, 'baci': r(bb.interaction.mean()) if len(bb) else None, 'q': r(bb.q.min(), 4) if len(bb) else None}
        for k, sel in (('pre_exp', (core.phase == 'pre') & (core.group == 'EXP')), ('pre_ctl', (core.phase == 'pre') & (core.group == 'CTL')),
                       ('oil', core.phase == 'oiled'), ('ctl', core.phase == 'control')):
            ix = core.index[(core.season == s) & sel]
            o[k] = r(tpm.loc[d, ix].median(), 0)
            sh = (M.loc[cu, ix].sum() / BIN.loc[d, ix].replace(0, np.nan)).median() * 1e3 if len(cu) else np.nan
            o[k + '_mo'] = r(sh, 1)
        row['summary'].append(o)
    J['genomes'].append(row)
json.dump(J, open(E0 + 'site/nitrogen.json', 'w'), separators=(',', ':'), allow_nan=False)
print(json.dumps([g['summary'] for g in J['genomes']])[:1500])
