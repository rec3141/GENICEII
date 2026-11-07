"""Paired ice-vs-water and top-vs-bottom contrasts for MAGs (bin TPM) and KOs (corrected KO CPM).
A 'pair' = same season, pool and day (2024-25 TOP/BOT ice + seawater; 2025-26 ice + seawater). Days are repeated measures
within one pool, so the sign test across pairs describes consistency, not pool-level replication.
Out: ice/out/{mag,ko}_{icewater,topbot}.tsv, ice/out/pairs.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import wilcoxon
O = '/data/dev/rnaseq-exploration/ice/out/'
import os; os.makedirs(O, exist_ok=True)
m = meta(); m = m[(m.group != 'other') & m.medium.isin(['SI', 'SW']) & m.day.notna()]
m['layer'] = m.ice_layer.fillna('').where(m.medium == 'SI', 'SW')

def pairs_iw():
    P = []
    for (s, g, d), x in m.groupby(['season', 'group', 'day']):
        si, sw = x.index[x.medium == 'SI'].tolist(), x.index[x.medium == 'SW'].tolist()
        if si and sw: P.append(dict(season=s, group=g, day=d, a=si, b=sw, kind='ice-water'))
    return P
def pairs_tb():
    P = []
    for (s, g, d), x in m.groupby(['season', 'group', 'day']):
        t, b = x.index[x.layer == 'TOP'].tolist(), x.index[x.layer == 'BOT'].tolist()
        if t and b: P.append(dict(season=s, group=g, day=d, a=t, b=b, kind='top-bottom'))
    return P

def contrast(L, P):
    """L: features x samples log-scale matrix. Returns per-feature mean paired delta, n pairs up/down, Wilcoxon p, per-season means."""
    D = pd.DataFrame({f"{p['season']}|{p['group']}|{int(p['day'])}": L[p['a']].mean(axis=1) - L[p['b']].mean(axis=1) for p in P})
    out = pd.DataFrame({'delta': D.mean(axis=1), 'n_up': (D > 0).sum(axis=1), 'n_down': (D < 0).sum(axis=1), 'n': D.shape[1]})
    for s in ('2024-25', '2025-26'):
        c = [k for k in D.columns if k.startswith(s)]
        if c: out[f'delta_{s}'] = D[c].mean(axis=1)
    out['p'] = [wilcoxon(r).pvalue if np.any(r != 0) else 1.0 for r in D.values]
    out['q'] = bh(out.p.values)
    return out, D

PI, PT = pairs_iw(), pairs_tb()
if __name__ == '__main__':
    pd.DataFrame([dict(p, a=','.join(p['a']), b=','.join(p['b'])) for p in PI + PT]).to_csv(O + 'pairs.tsv', sep='\t', index=False)
    print('ice-water pairs', len(PI), ' top-bottom pairs', len(PT))

    # MAGs: CLR of bin TPM (bins with max TPM >= 1000 somewhere in ice or water pairs)
    B = table('bin_tpm'); B = B[B.max(axis=1) >= 1000]
    Lb = clr(B, pseudo=1.0) / np.log(2)  # log2 units
    tax = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', index_col=0)['classification']
    def lab(b):
        if plastid_label(b): return plastid_label(b)
        t = str(tax.get(b.replace('dastool-', ''), '')).split(';')
        g = [x[3:] for x in t if len(x) > 3]
        return (g[-1] if g else '') or kaiju_bin_label(b)
    for name, P in (('icewater', PI), ('topbot', PT)):
        r, Dm = contrast(Lb, P)
        r['taxon'] = [lab(b) for b in r.index]
        r['family'] = [next((x[3:] for x in str(tax.get(b.replace('dastool-', ''), '')).split(';') if x.startswith('f__')), '') for b in r.index]
        r['max_tpm'] = B.max(axis=1)
        r.sort_values('delta').to_csv(O + f'mag_{name}.tsv', sep='\t'); Dm.to_csv(O + f'mag_{name}_pairs.tsv', sep='\t')
    # KOs: log2 CPM, KOs with mean CPM >= 5
    K = pd.read_csv(D + 'ko_cpm.tsv', sep='\t', index_col=0); K = K[K.mean(axis=1) >= 5]
    Lk = np.log2(K + 1)
    kd = pd.read_csv(D + 'ko_list.tsv', sep='\t', header=None, index_col=0)[1]
    for name, P in (('icewater', PI), ('topbot', PT)):
        r, Dm = contrast(Lk, P)
        r['def'] = kd.reindex(r.index).values
        r['mean_cpm'] = K.mean(axis=1)
        r.sort_values('delta').to_csv(O + f'ko_{name}.tsv', sep='\t'); Dm.to_csv(O + f'ko_{name}_pairs.tsv', sep='\t')
    print('done', B.shape, K.shape)
