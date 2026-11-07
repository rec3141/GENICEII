"""Post-processing of states.py output (no recount): what the clusters/axes align with (season, pool, medium, day),
PC and growth-index trends with day inside each pool after dosing. Out: states/out/{assoc.tsv, trends.tsv}"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score
O = '/data/dev/rnaseq-exploration/states/out/'
SC = pd.read_pickle(O + 'scores.pkl'); GR = pd.read_pickle(O + 'growth.pkl'); m = meta()
def r2(Y, lab):
    lab = np.asarray(lab); gm = Y.mean(0); sst = ((Y - gm) ** 2).sum()
    ssw = sum(((Y[lab == g] - Y[lab == g].mean(0)) ** 2).sum() for g in np.unique(lab)); return 1 - ssw / sst
A, T = [], []
for b, s in SC.items():
    mm = m.loc[s['samples']]; Y = s['pc'][:, :3]
    row = dict(bin=b, n=len(mm))
    for f in ('season', 'group', 'medium'):
        if mm[f].nunique() > 1: row[f'R2_{f}'] = r2(Y, mm[f].values); row[f'ari_{f}'] = adjusted_rand_score(mm[f].values, s['cluster'])
    A.append(row)
    for se in ('2024-25', '2025-26'):
        for P, ph in (('EXP', 'oiled'), ('CTL', 'control')):
            ix = ((mm.season == se) & (mm.phase.isin([ph, 'pre'])) & (mm.group == P)).values
            if ix.sum() < 6: continue
            d = mm.day.values[ix].astype(float); g = np.array([GR.get(b, {}).get(x, np.nan) for x in np.array(s['samples'])[ix]])
            T.append(dict(bin=b, season=se, pool=P, n=int(ix.sum()), r_pc1_day=spearmanr(d, Y[ix, 0])[0], r_pc2_day=spearmanr(d, Y[ix, 1])[0],
                          r_growth_day=spearmanr(d, g, nan_policy='omit')[0] if np.isfinite(g).sum() >= 6 else np.nan,
                          growth_first=np.nanmedian(g[d <= np.nanpercentile(d, 33)]) * 1e3 if np.isfinite(g).any() else np.nan,
                          growth_last=np.nanmedian(g[d >= np.nanpercentile(d, 67)]) * 1e3 if np.isfinite(g).any() else np.nan))
A = pd.DataFrame(A); T = pd.DataFrame(T)
A.to_csv(O + 'assoc.tsv', sep='\t', index=False); T.to_csv(O + 'trends.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(A.describe().round(2).to_string())
best = A[['R2_season', 'R2_group', 'R2_medium']].idxmax(1).value_counts(); print(best)
KEY = ['lorbin_184', 'metabat_9342', 'metabat_5432', 'comebin_3787_sub', 'comebin_2576', 'semibin_088', 'vamb_549', 'metabat_5376', 'vamb_4149', 'lorbin_095']
print(A[A.bin.str.replace('dastool-', '').isin(KEY)].round(2).to_string())
print(T[T.bin.str.replace('dastool-', '').isin(KEY)].round(2).to_string())
