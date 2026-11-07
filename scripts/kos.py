"""KOs consistently up/down in the oiled pool (vs control, both seasons, SW and SI), on log2 TPM."""
import sys; sys.path.insert(0, '.')
from common import *
from scipy.stats import mannwhitneyu
m = meta(); k = pd.read_csv(D + 'ko_cpm.tsv', sep='\t', index_col=0); S = k.columns; m = m.loc[S]
names = pd.read_csv(D + 'ko_list.tsv', sep='\t', header=None, index_col=0)[1]
L = np.log2(k + 1)
res = {}
for se in ['2024-25', '2025-26']:
    for me in ['SW', 'SI']:
        a = m.index[(m.season == se) & (m.design == f'EXP_{me}') & (m.phase == 'oiled')]
        b = m.index[(m.season == se) & (m.design == f'CTL_{me}') & (m.phase == 'control')]
        lfc = L[a].mean(1) - L[b].mean(1)
        p = np.array([mannwhitneyu(L.loc[i, a], L.loc[i, b]).pvalue for i in L.index])
        res[(se, me, 'lfc')] = lfc; res[(se, me, 'q')] = pd.Series(bh(p), index=L.index)
R = pd.DataFrame(res)
expr = k.max(1) >= 20
lf = R.xs('lfc', axis=1, level=2); qq = R.xs('q', axis=1, level=2)
score = lf.min(1).where(lf.min(1) > 0, 0) + lf.max(1).where(lf.max(1) < 0, 0)
sig = (qq < 0.05).sum(1)
out = lf.round(2).copy(); out.columns = [f'{a} {b}' for a, b in out.columns]
out['n_sig'] = sig; out['min_tpm_max'] = k.max(1).round(0); out['name'] = names.reindex(out.index).str.slice(0, 70)
out = out[expr]
out.to_csv('../out/ko_oil_response.tsv', sep='\t')
pd.set_option('display.width', 250)
u = out[(lf[expr] > 1).all(1) & (sig[expr] >= 3)].copy(); u['m'] = lf[expr].loc[u.index].mean(1)
print(f'KOs up >2x in all 4 comparisons (q<0.05 in >=3): {len(u)}')
print(u.sort_values('m', ascending=False).head(45).drop(columns='m').to_string())
d = out[(lf[expr] < -1).all(1) & (sig[expr] >= 3)].copy(); d['m'] = lf[expr].loc[d.index].mean(1)
print(f'\nKOs down >2x in all 4: {len(d)}')
print(d.sort_values('m').head(25).drop(columns='m').to_string())
