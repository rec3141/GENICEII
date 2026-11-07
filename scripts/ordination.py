"""Aitchison PCA of MAG transcript counts + PERMANOVA (pseudo-F, 999 perms) for season, medium, pool."""
import sys; sys.path.insert(0, '.')
from common import *
rng = np.random.default_rng(0)
m = meta(); bc = table('bin_counts'); S = bc.columns; m = m.loc[S]
keep = bc.sum(1) >= 1000
C = clr(bc[keep].values.astype(float)).T          # samples × MAGs
def pca(Y):
    Y = Y - Y.mean(0); U, s, Vt = np.linalg.svd(Y, full_matrices=False)
    return U * s, s ** 2 / (s ** 2).sum()
pc, ev = pca(C)
pd.DataFrame(pc[:, :4], index=S, columns=['PC1', 'PC2', 'PC3', 'PC4']).join(m[['season', 'design', 'phase', 'day', 'sample_name']]).to_csv('../out/pca_mags.tsv', sep='\t')
print('PC var %:', (ev[:5] * 100).round(1))
def permanova(Y, labels, n=999):
    Yc = Y - Y.mean(0); sst = (Yc ** 2).sum()
    def F(l):
        ssw = sum(((Y[l == k] - Y[l == k].mean(0)) ** 2).sum() for k in np.unique(l))
        g = len(np.unique(l)); return ((sst - ssw) / (g - 1)) / (ssw / (len(l) - g)), 1 - ssw / sst
    f, r2 = F(labels); null = [F(rng.permutation(labels))[0] for _ in range(n)]
    return r2, (np.sum(np.array(null) >= f) + 1) / (n + 1)
core = m.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI']).values
Y = C[core]; mm = m[core]
for name, lab in [('season', mm.season), ('medium', mm.medium), ('pool', mm.group), ('phase(pre/ctl/oiled)', mm.phase),
                  ('season×pool×medium', mm.season + mm.group + mm.medium)]:
    r2, p = permanova(Y, lab.values); print(f'  {name:22s} R2={r2:.3f} p={p:.3f}')
for se in ['2024-25', '2025-26']:
    s = (mm.season == se).values
    for name, lab in [('medium', mm.medium), ('pool', mm.group), ('oiled vs not', mm.phase.eq('oiled'))]:
        r2, p = permanova(Y[s], lab.values[s]); print(f'  {se} {name:14s} R2={r2:.3f} p={p:.3f}')
