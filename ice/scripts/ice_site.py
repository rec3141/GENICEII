"""Assemble site/ice.json (ice vs water, top vs bottom) from ice/out, ice/phys and the two fork outputs (cold.json, popgen.json)."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
I = '/data/dev/rnaseq-exploration/ice/'; O = I + 'out/'
def r(x, n=3):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
J = {}
GT = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
def blab(b):
    if not isinstance(b, str) or not b or b == 'unbinned': return 'unbinned contigs'
    if plastid_label(b): return plastid_label(b)
    t = [x.split('__', 1)[1] for x in str(GT.get(b.replace('dastool-', ''), '')).split(';') if '__' in x and x.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')
# physics: per-core summary + temperature profiles
C = pd.read_csv(I + 'phys/cores.tsv', sep='\t'); P = pd.read_csv(I + 'phys/profiles.tsv', sep='\t')
J['cores'] = [{k: (r(v) if isinstance(v, (int, float, np.floating, np.integer)) and k not in ('season', 'pool') else v) for k, v in x.items()} for x in C.to_dict('records')]
prof = []
for (s, p, d), g in P[P['var'] == 'T'].dropna(subset=['lo']).groupby(['season', 'pool', 'day']):
    if s == '2025-26' and p == 'CTL' and d == 12: continue  # profile inverted (colder downward); likely logged bottom-up
    pts = sorted({((a + b) / 2, v) for a, b, v in zip(g.lo, g.hi, g.val)})
    if len(pts) >= 3: prof.append({'season': s, 'pool': p, 'day': int(d), 'pts': [[r(z, 1), r(v, 2)] for z, v in pts],
                                   'H': r(C[(C.season == s) & (C.pool == p) & (C.day == d)].thickness_cm.max(), 1)})
J['profiles'] = prof
# MAG contrasts (ice-water vs top-bottom), only MAGs with q<0.1 in either or top-200 by max TPM
iw = pd.read_csv(O + 'mag_icewater.tsv', sep='\t', index_col=0); tb = pd.read_csv(O + 'mag_topbot.tsv', sep='\t', index_col=0)
j = iw.join(tb[['delta', 'q', 'n_up', 'n_down']], rsuffix='_tb')
j['genus'] = j.taxon.fillna('').str.replace(r' \(Kaiju\)', '', regex=True).str.split(' ').str[0]
j.loc[j.taxon.fillna('').str.startswith('Plastid'), 'genus'] = 'Plastid'
keep = (j.q < 0.1) | (j.q_tb < 0.1) | (j.max_tpm.rank(ascending=False) <= 150)
J['mags'] = [{'bin': b.replace('dastool-', ''), 'taxon': x.taxon if isinstance(x.taxon, str) else '', 'family': x.family if isinstance(x.family, str) else '',
              'genus': x.genus, 'iw': r(x.delta, 2), 'iw24': r(x['delta_2024-25'], 2), 'iw25': r(x['delta_2025-26'], 2), 'iw_q': r(x.q, 4),
              'iw_up': int(x.n_up), 'iw_n': int(x.n), 'tb': r(x.delta_tb, 2), 'tb_q': r(x.q_tb, 4), 'tb_up': int(x.n_up_tb), 'tb_n': int(x.n_up_tb + x.n_down_tb),
              'tpm': r(x.max_tpm, 0)} for b, x in j[keep].iterrows()]
# per-pair deltas for a few named MAGs (for the strip plot)
J['mag_pairs'] = {k: pd.read_csv(O + f'mag_{k}_pairs.tsv', sep='\t', index_col=0).round(2).pipe(
    lambda d: {'cols': d.columns.tolist(), 'rows': {b.replace('dastool-', ''): d.loc[b].tolist() for b in d.index if b in j[keep].index}}) for k in ('icewater', 'topbot')}
# rRNA taxonomy contrasts
def rr(f, n=14):
    d = pd.read_csv(O + f, sep='\t', index_col=0)
    d = pd.concat([d.head(n), d.tail(n)]).drop_duplicates()
    return [{'lineage': i, 'short': ';'.join(i.split(';')[-3:]), 'delta': r(x.delta, 2), 'up': int(x.n_up), 'n': int(x.n), 'q': r(x.q, 4), 'share': r(x.mean_share, 5)} for i, x in d.iterrows()]
J['rrna'] = {f'{dom}_{k}': rr(f'rrna_{dom}_{k}.tsv') for dom in ('euk', 'prok') for k in ('icewater', 'topbot')}
# KO tops
def ko(f, n=25):
    d = pd.read_csv(O + f, sep='\t', index_col=0); d = d[(d.mean_cpm >= 20) & (d.q < 0.05)]
    d = pd.concat([d.head(n), d.tail(n)])
    return [{'ko': i, 'def': x['def'] if isinstance(x['def'], str) else '', 'delta': r(x.delta, 2), 'up': int(x.n_up), 'n': int(x.n), 'q': r(x.q, 4), 'cpm': r(x.mean_cpm, 1)} for i, x in d.iterrows()]
J['ko'] = {'icewater': ko('ko_icewater.tsv'), 'topbot': ko('ko_topbot.tsv')}
if os.path.exists(O + 'panel.tsv'):
    p = pd.read_csv(O + 'panel.tsv', sep='\t', index_col=0)
    pp = {k: pd.read_csv(O + f'panel_{k}_pairs.tsv', sep='\t', index_col=0) for k in ('icewater', 'topbot')}
    H = pd.read_csv(O + 'panel_hosts.tsv', sep='\t') if os.path.exists(O + 'panel_hosts.tsv') else None
    J['panel'] = [{'cat': c, 'iw': r(x.icewater_delta, 2), 'iw_up': int(x.icewater_n_up), 'iw_n': int(x.icewater_n), 'iw_q': r(x.icewater_q, 4),
                   'tb': r(x.topbot_delta, 2), 'tb_up': int(x.topbot_n_up), 'tb_n': int(x.topbot_n), 'tb_q': r(x.topbot_q, 4),
                   'cpm': {'SW': r(x.cpm_SW, 1), 'TOP': r(x.cpm_TOP, 1), 'BOT': r(x.cpm_BOT, 1)},
                   'iw_pairs': [r(v, 2) for v in pp['icewater'].loc[c]], 'tb_pairs': [r(v, 2) for v in pp['topbot'].loc[c]],
                   'hosts': [] if H is None else [{'taxon': blab(h.bin), 'bin': h.bin if isinstance(h.bin, str) else '', 'share': r(h.share, 3)} for _, h in H[H.category == c].iterrows()]} for c, x in p.iterrows()]
    J['panel_cols'] = {k: v.columns.tolist() for k, v in pp.items()}
# community divergence from paired seawater
dv = pd.read_csv(O + 'divergence.tsv', sep='\t')
J['divergence'] = [{'season': x.season, 'pool': x.pool, 'day': int(x.day), 'layer': x.layer, 'name': x['name'], 'd': r(x.d_water, 2), 'H': r(x.thickness, 1)}
                   for _, x in dv.dropna(subset=['d_water']).iterrows()]
for f in ('cold', 'popgen'):
    p = I + f'{f}/{f}.json'
    if os.path.exists(p): J[f] = json.load(open(p))
if 'popgen' in J:
    # paired same-day FST (ice vs water, top vs bottom) against consecutive-day water-vs-water FST, pairs with >= 1000 sites (2024-25),
    # and pool differentiation (control vs diesel pool, pooled) for every MAG with enough sites
    pf = pd.read_csv(I + 'popgen/paired_fst_vs_null.tsv', sep='\t')
    J['popgen']['paired'] = [{k: (r(v, 4) if isinstance(v, (float, np.floating)) else (int(v) if isinstance(v, (np.integer,)) else v)) for k, v in x.items()} for x in pf.to_dict('records')]
    ps = pd.read_csv(I + 'popgen/popgen_summary.tsv', sep='\t')
    pool = []
    for _, x in ps.iterrows():
        for med in ('water', 'ice'):
            f, n = x.get(f'fst|2024-25|{med}_ctl_vs_exp'), x.get(f'sites|2024-25|{med}_ctl_vs_exp')
            iw, ni = x.get('fst|2024-25|ice_vs_water'), x.get('sites|2024-25|ice_vs_water')
            if pd.notna(f) and n >= 200: pool.append({'bin': x.bin.replace('dastool-', ''), 'label': x.label, 'medium': med, 'fst': r(f, 3), 'sites': int(n),
                                                      'iw': r(iw, 3) if pd.notna(iw) and ni >= 200 else None})
    J['popgen']['pool'] = pool
    J['popgen'].pop('popgen', None)   # per-bin raw blocks are large; the site uses the summaries above
    J['popgen'].pop('samples', None)
json.dump(J, open('/data/dev/rnaseq-exploration/site/ice.json', 'w'), separators=(',', ':'), allow_nan=False)
print({k: len(v) if hasattr(v, '__len__') else v for k, v in J.items()}, os.path.getsize('/data/dev/rnaseq-exploration/site/ice.json'))
