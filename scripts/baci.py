"""BACI for 2025-26 seawater (5 pre-oil samples per pool): per-MAG interaction of pool × period on CLR, permutation test."""
import sys; sys.path.insert(0, '.')
from common import *
rng = np.random.default_rng(1)
m = meta(); bc = table('bin_counts'); S = bc.columns; m = m.loc[S]
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
C = pd.DataFrame(clr(bc.values.astype(float)), index=bc.index, columns=S)
out = []
for se, me in [('2025-26', 'SW'), ('2024-25', 'SW')]:
    sel = m[(m.season == se) & m.design.isin([f'EXP_{me}', f'CTL_{me}'])]
    pool = (sel.group == 'EXP').values; post = (sel.day >= OIL_DAY[f'CMO {se}']).values
    Y = C[sel.index].values
    keep = bc[sel.index].sum(1).values >= 500
    def inter(pool, post):
        f = lambda a: Y[:, a].mean(1)
        return (f(pool & post) - f(pool & ~post)) - (f(~pool & post) - f(~pool & ~post))
    obs = inter(pool, post)
    # permute pool labels within period (keeps the pre/post design)
    perm = np.zeros((2000, len(obs)))
    for i in range(2000):
        pp = pool.copy()
        for k in (False, True):
            ix = np.where(post == k)[0]; pp[ix] = rng.permutation(pool[ix])
        perm[i] = inter(pp, post)
    p = ((np.abs(perm) >= np.abs(obs)).sum(0) + 1) / 2001
    r = pd.DataFrame({'season': se, 'medium': me, 'bin': bc.index, 'interaction': obs, 'p': p})[keep]
    r['q'] = bh(r.p); out.append(r)
    print(se, me, 'n pre EXP/CTL', (pool & ~post).sum(), (~pool & ~post).sum(), 'post', (pool & post).sum(), (~pool & post).sum())
R = pd.concat(out)
R['taxon'] = [gt.get(b.replace('dastool-', ''), '').replace('d__Bacteria;', '') for b in R.bin]
R['genus'] = R.taxon.str.extract(r'g__([^;]*)')[0].fillna('') + '|' + R.taxon.str.extract(r'f__([^;]*)')[0].fillna('')
R.to_csv('../out/mag_baci.tsv', sep='\t', index=False)
for se, d in R.groupby('season'):
    print(f'\n{se}: q<0.1 up {((d.q<0.1)&(d.interaction>0)).sum()}, down {((d.q<0.1)&(d.interaction<0)).sum()} of {len(d)}')
    print(d.sort_values('interaction', ascending=False).head(15)[['bin', 'genus', 'interaction', 'p', 'q']].round(3).to_string(index=False))
    print(d.sort_values('interaction').head(8)[['bin', 'genus', 'interaction', 'p', 'q']].round(3).to_string(index=False))
