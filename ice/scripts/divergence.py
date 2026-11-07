"""How far ice communities sit from their paired seawater (and top from bottom), against ice age, thickness, temperature
and brine salinity. Aitchison distance on CLR of MAG TPM (MAGs >=1000 TPM somewhere). Out: ice/out/divergence.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
import contrasts as C
from scipy.stats import spearmanr
B = table('bin_tpm'); B = B[B.max(axis=1) >= 1000]; L = clr(B, 1.0)
cores = pd.read_csv('/data/dev/rnaseq-exploration/ice/phys/cores.tsv', sep='\t')
def core(season, pool, day):
    c = cores[(cores.season == season) & (cores.pool == pool)]
    if not len(c): return {}
    c = c.iloc[(c.day - day).abs().argsort()[:1]]
    return c.iloc[0].to_dict() if abs(c.day.iloc[0] - day) <= 2 else {}
rows = []
m = C.m
for (s, g, d), x in m.groupby(['season', 'group', 'day']):
    sw = x.index[x.medium == 'SW']; ph = core(s, g, d)
    for sid in x.index[x.medium == 'SI']:
        r = dict(season=s, pool=g, day=d, sample=sid, name=x.sample_name[sid], layer=x.layer[sid] or 'SI',
                 thickness=ph.get('thickness_cm'), T_top=ph.get('T_top'), T_bot=ph.get('T_bot'),
                 Sb_top=ph.get('Sbrine_top'), Sb_bot=ph.get('Sbrine_bot'), Vb_top=ph.get('Vb_top'), Vb_bot=ph.get('Vb_bot'))
        if len(sw): r['d_water'] = float(np.sqrt(((L[sid] - L[sw].mean(axis=1)) ** 2).sum()))
        rows.append(r)
R = pd.DataFrame(rows)
# replicate-water baseline: distance between seawater replicates of the same day (a, b, ...) for scale
base = [float(np.sqrt(((L[a] - L[b]) ** 2).sum())) for (_, _, _), x in m[m.medium == 'SW'].groupby(['season', 'group', 'day'])
        for a, b in zip(x.index[:-1], x.index[1:])]
R.to_csv(C.O + 'divergence.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print('SW replicate distance median %.1f (n=%d)' % (np.median(base), len(base)))
print(R.groupby(['season', 'layer']).d_water.describe().round(1))
for s in ('2024-25', '2025-26'):
    for lay in ('TOP', 'BOT', 'SI'):
        x = R[(R.season == s) & (R.layer == lay)].dropna(subset=['d_water'])
        if len(x) > 4:
            print(s, lay, 'rho(day)=%.2f p=%.3f' % spearmanr(x.day, x.d_water), ' rho(thick)=%.2f' % spearmanr(x.thickness, x.d_water, nan_policy='omit')[0])
