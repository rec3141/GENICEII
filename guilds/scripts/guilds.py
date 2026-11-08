"""Response guilds: embed all profiled genomes (states/) in one space by what drives their within-genome state.

Layer 1, driver fingerprints. For each genome, Y = its within-genome PC1-3 scores (states/out/scores.pkl; full-depth RNA,
depth-aware PCA). Multivariate (RDA-style) R² of Y on five driver blocks, all samples weighted equally:
  winter  : 2025-26 indicator (only for genomes profiled in both winters; otherwise NA)
  ice     : sea ice vs seawater, plus top (+1) / bottom (-1) ice layer
  diesel  : diesel pool (EXP) indicator + oiled (EXP after dosing) indicator  -> pool history + dosing (one pool per treatment)
  time    : experiment day within winter (z-scored per winter)
  light   : log modelled under-ice PAR (light/scripts/underice.py)
Commonality: unique R² per block = R²(all) - R²(all but block). 'seasonal' = the part shared among winter, time and light:
R²(all) - R²(all but winter,time,light) - their unique parts. 'other shared' = the rest. Season is treated as a meaningful
bundle of drivers (light, temperature, ice growth, nutrients), not as a nuisance. Unique fractions are tested by permuting
the block's rows within strata that hold the other design factors fixed (winter: within medium; ice: within winter;
diesel and light: within winter x medium; time: within winter x medium x pool; 199 permutations), BH across genomes per block.
Embedding: the six fractions (as shares of each genome's total variance) -> UMAP; clusters: sklearn HDBSCAN, falling back
to Gaussian mixtures if HDBSCAN leaves most genomes unclustered. Single-winter genomes: winter fraction imputed with the median
of the other genomes for the embedding only, and flagged.

Layer 2, co-response guilds. Shared-factor model over the common sample set: stack every genome's PC1-3 scores
(samples x 3*genomes; missing where a genome was not profiled), scale each genome block to unit total variance, and fit a
rank-k PCA with missing values by EM (iterative SVD imputation). k by variance explained and split-half stability (genomes
split in two halves, factor sample scores matched by absolute correlation). PC signs are arbitrary per genome, so a genome's
'loading' on a factor is the R² of its PC1-3 block on that factor's sample scores (observed samples only). Loadings -> UMAP,
clusters = guilds. Factors are interpreted by correlating their sample scores with the drivers.
Out: guilds/out/{fingerprints.tsv, factors_samples.tsv, factor_drivers.tsv, guild_loadings.tsv, model.json}"""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr
import umap
from sklearn.cluster import HDBSCAN
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
rng = np.random.default_rng(0)
E0 = '/data/dev/rnaseq-exploration/'; O = E0 + 'guilds/out/'; os.makedirs(O, exist_ok=True)
SC = pd.read_pickle(E0 + 'states/out/scores.pkl')
SC = {b: v for b, v in SC.items() if is_prokaryote(b)}
m = meta(); L = pd.read_csv(E0 + 'light/out/underice_samples.tsv', sep='\t', index_col=0)
allS = sorted(set(s for v in SC.values() for s in v['samples']))
X = pd.DataFrame(index=allS)
X['season'] = m.season.reindex(allS); X['pool'] = m.group.reindex(allS); X['medium'] = m.medium.reindex(allS)
X['layer'] = m.ice_layer.reindex(allS).fillna(''); X['phase'] = m.phase.reindex(allS); X['day'] = m.day.reindex(allS).astype(float)
X['light'] = np.log10(L.par_under.reindex(allS).astype(float).clip(lower=1e-3))
X['dayz'] = X.groupby('season').day.transform(lambda d: (d - d.mean()) / d.std())
lay = X.layer.str.upper(); X['topbot'] = np.select([lay.isin(['TOP', 'HIGH']), lay.isin(['BOT', 'LOW'])], [1.0, -1.0], 0.0)
BLOCKS = {'winter': lambda x: [(x.season == '2025-26').astype(float)],
          'ice': lambda x: [(x.medium == 'SI').astype(float), x.topbot],
          'diesel': lambda x: [(x.pool == 'EXP').astype(float), (x.phase == 'oiled').astype(float)],
          'time': lambda x: [x.dayz],
          'light': lambda x: [x.light]}
DR = list(BLOCKS)

def r2(Y, cols):
    if not cols: return 0.0
    A = np.column_stack([np.ones(len(Y))] + cols)
    B, *_ = np.linalg.lstsq(A, Y, rcond=None); R = Y - A @ B
    Yc = Y - Y.mean(0); return 1 - (R ** 2).sum() / (Yc ** 2).sum()
def design(x, blocks):
    cols = {}
    for b in blocks:
        cc = [np.asarray(c, float) for c in BLOCKS[b](x)]
        cc = [c for c in cc if np.nanstd(c) > 0]
        if cc: cols[b] = cc
    return cols

rows, NP = [], 199
CACHE = O + 'fingerprint_rows.pkl'
for b, v in ({} if os.path.exists(CACHE) else SC).items():
    x = X.loc[v['samples']]; Y = np.asarray(v['pc'], float)[:, :3]
    ok = x[['day', 'light']].notna().all(1).values; x = x[ok]; Y = Y[ok]
    blocks = [d for d in DR if not (d == 'winter' and x.season.nunique() < 2)]
    cols = design(x, blocks); blocks = [d for d in blocks if d in cols]
    flat = lambda bs: [c for d in bs for c in cols[d]]
    full = r2(Y, flat(blocks)); out = dict(bin=b, n=len(Y), seasons=x.season.nunique(), r2_full=full)
    # permutation strata per block: keep the other design factors fixed while shuffling the tested one
    ST = {'winter': x.medium.values, 'ice': x.season.values, 'diesel': (x.season + '|' + x.medium).values,
          'time': (x.season + '|' + x.medium + '|' + x.pool).values, 'light': (x.season + '|' + x.medium).values}
    for d in DR:
        if d not in blocks: out[f'u_{d}'] = np.nan; out[f'p_{d}'] = np.nan; continue
        rest = [k for k in blocks if k != d]; u = full - r2(Y, flat(rest)); out[f'u_{d}'] = u
        cnt = 0
        for _ in range(NP):
            perm = np.arange(len(Y))
            for s in np.unique(ST[d]):
                ix = np.where(ST[d] == s)[0]; perm[ix] = rng.permutation(ix)
            cp = flat(rest) + [c[perm] for c in cols[d]]
            if r2(Y, cp) - r2(Y, flat(rest)) >= u - 1e-12: cnt += 1
        out[f'p_{d}'] = (cnt + 1) / (NP + 1)
    seas = [k for k in ('winter', 'time', 'light') if k in blocks]
    out['seasonal'] = full - r2(Y, flat([k for k in blocks if k not in seas])) - sum(out[f'u_{k}'] for k in seas)
    out['other_shared'] = full - sum(np.nan_to_num(out[f'u_{d}']) for d in DR) - out['seasonal']
    rows.append(out)
if os.path.exists(CACHE): rows = pd.read_pickle(CACHE)
else: pd.to_pickle(rows, CACHE)
F = pd.DataFrame(rows).set_index('bin')
for d in DR: F[f'q_{d}'] = np.nan; ok = F[f'p_{d}'].notna(); F.loc[ok, f'q_{d}'] = bh(F.loc[ok, f'p_{d}'].values)
FR = ['u_winter', 'u_ice', 'u_diesel', 'u_time', 'u_light', 'seasonal']
# dominant driver: largest of the unique fractions and the shared seasonal part
lab = {'u_winter': 'winter', 'u_ice': 'ice vs water', 'u_diesel': 'diesel pool', 'u_time': 'time', 'u_light': 'light', 'seasonal': 'seasonal (shared)'}
F['dominant'] = F[FR].fillna(-1).idxmax(1).map(lab)
F['single_winter'] = F.seasons < 2
# taxonomy, abundance, function directions
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
lin = json.load(open(E0 + 'site/lineage.json'))
def phylum(b):
    t = str(gt.get(b.replace('dastool-', ''), '')); p = [y[3:] for y in t.split(';') if y.startswith('p__') and len(y) > 3]
    if p: return p[0]
    l = lin.get(b.replace('dastool-', ''), ''); parts = l.split(' › ')
    return (parts[1].split(' (')[0] + ' (Kaiju)') if len(parts) > 1 else 'unclassified'
F['phylum'] = [phylum(b) for b in F.index]
F['taxon'] = [special_label(b) or (lambda t: [y.split('__', 1)[1] for y in t.split(';') if '__' in y and y.split('__', 1)[1]])(str(gt.get(b.replace('dastool-', ''), '')))[-1:] and [y.split('__', 1)[1] for y in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in y and y.split('__', 1)[1]][-1] or kaiju_bin_label(b) or b.replace('dastool-', '') for b in F.index]
tpm = table('bin_tpm'); core = [s for s in allS if s in tpm.columns]
F['tpm'] = [float(tpm.loc[b, core].mean()) if b in tpm.index else np.nan for b in F.index]
fn = pd.read_csv(E0 + 'states/out/functions.tsv', sep='\t')
FN = fn.pivot_table(index='bin', columns='category', values='r_pc1')
# PC1 alignment with each driver (which driver the main axis follows)
al = {}
for b, v in SC.items():
    x = X.loc[v['samples']]; y = np.asarray(v['pc'], float)[:, :1]; ok = x[['day', 'light']].notna().all(1).values
    cols = design(x[ok], DR); al[b] = {d: r2(y[ok], cols[d]) if d in cols else np.nan for d in DR}
AL = pd.DataFrame(al).T; F['pc1_driver'] = AL.idxmax(1).reindex(F.index)
# embedding + clusters on fractions of total variance
Z = F[FR].copy(); Z['u_winter'] = Z.u_winter.fillna(Z.u_winter.median())
Z = Z.fillna(0)   # other blocks are NA only when the driver does not vary in that genome's samples: no information, no effect
Zs = (Z - Z.mean()) / Z.std()
emb = umap.UMAP(n_neighbors=15, min_dist=0.25, random_state=1).fit_transform(Zs.values)
F['fx'], F['fy'] = emb[:, 0], emb[:, 1]
def clusters(M, name):
    lab_h = HDBSCAN(min_cluster_size=6).fit_predict(M)
    if (lab_h < 0).mean() <= 0.35 and len(set(lab_h) - {-1}) >= 3: return lab_h, 'HDBSCAN'
    best = None
    for k in range(3, 9):
        g = GaussianMixture(k, covariance_type='full', random_state=0, n_init=5).fit(M); l = g.predict(M)
        if len(set(l)) < 2: continue
        s = silhouette_score(M, l)
        if best is None or s > best[0]: best = (s, l, k)
    print(name, 'HDBSCAN noise', round((lab_h < 0).mean(), 2), '-> GMM k', best[2], 'silhouette', round(best[0], 2)); return best[1], f'Gaussian mixture (k={best[2]}, silhouette {best[0]:.2f})'
cl, cmeth = clusters(Zs.values, 'fingerprint')
F['fcluster'] = cl
name = {}
for c in sorted(set(cl)):
    if c < 0: name[c] = 'unclustered'; continue
    mz = Zs[F.fcluster == c].mean(); top = mz.idxmax()      # driver most over-represented relative to all genomes
    lev = {'u_winter': 'winter-switchers', 'u_ice': 'ice-switchers', 'u_diesel': 'diesel-pool responders', 'u_time': 'time drifters',
           'u_light': 'light followers', 'seasonal': 'seasonal (shared drivers)'}[top]
    if F.r2_full[F.fcluster == c].median() < F.r2_full.median() - 0.1: lev += ', weakly structured'
    name[c] = lev
# make names unique
seen = {}
for c in sorted(name):
    n = name[c]; seen[n] = seen.get(n, 0) + 1
    if seen[n] > 1: name[c] = f'{n} {seen[n]}'
F['fcluster_name'] = F.fcluster.map(name)

# ---------- layer 2: shared factors ----------
gl = list(SC)
M = pd.DataFrame(index=allS, dtype=float)
for b in gl:
    pc = pd.DataFrame(np.asarray(SC[b]['pc'], float)[:, :3], index=SC[b]['samples'], columns=[f'{b}|{i}' for i in range(3)])
    pc = pc - pc.mean(); pc = pc / np.sqrt((pc ** 2).sum().sum() / len(pc))      # unit total variance per genome block
    M = M.join(pc, how='left')
obs = M.notna().values; A0 = M.values.copy()
def empca(A, mask, k, it=300):
    Z = np.where(mask, A, 0.0)
    for _ in range(it):
        U, s, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
        R = (U[:, :k] * s[:k]) @ Vt[:k] + Z.mean(0)
        Zn = np.where(mask, A, R)
        if np.abs(Zn - Z).max() < 1e-6: Z = Zn; break
        Z = Zn
    mu = Z.mean(0); U, s, Vt = np.linalg.svd(Z - mu, full_matrices=False)
    return U[:, :k] * s[:k], Vt[:k], mu
# choose k: observed-entry variance explained + split-half stability of sample scores
tot = np.nansum((A0 - np.nanmean(A0, 0)) ** 2)
ve, stab = {}, {}
srng = np.random.default_rng(42); SPLITS = [srng.permutation(len(gl)) for _ in range(20)]   # fixed splits: deterministic K
def halves(perm):
    h1 = set(np.array(gl)[perm[: len(gl) // 2]])
    return [i for i, c in enumerate(M.columns) if c.split('|')[0] in h1], [i for i, c in enumerate(M.columns) if c.split('|')[0] not in h1]
HV = [halves(p) for p in SPLITS]
for k in range(1, 9):
    S_, V_, mu_ = empca(A0, obs, k); R = S_ @ V_ + mu_
    ve[k] = 1 - np.nansum(np.where(obs, (A0 - R) ** 2, np.nan)) / tot
    st_ = []
    for c1, c2 in HV:
        Sa, _, _ = empca(A0[:, c1], obs[:, c1], k, it=150); Sb, _, _ = empca(A0[:, c2], obs[:, c2], k, it=150)
        C = np.abs(np.corrcoef(Sa.T, Sb.T)[:k, k:]); st_.append(float(np.mean(C.max(1))))
    stab[k] = float(np.mean(st_))
K = 1
while K + 1 in stab and stab[K + 1] >= 0.7: K += 1   # stop at the first K whose factors no longer replicate across halves (mean |r| < 0.7 over 20 splits)
K = max(K, 2)
print('variance explained', {k: round(v, 3) for k, v in ve.items()}); print('split-half stability', {k: round(v, 2) for k, v in stab.items()}, '-> K', K)
S_, V_, sv = empca(A0, obs, K)
FS = pd.DataFrame(S_ / S_.std(0), index=allS, columns=[f'F{i + 1}' for i in range(K)])
# genome loadings: R² of each genome's PC1-3 block on each factor's sample scores (observed samples)
GLd = {}
for b in gl:
    cols = [c for c in M.columns if c.startswith(b + '|')]; ok = M[cols].notna().all(1).values; Y = M.loc[ok, cols].values
    GLd[b] = {f: r2(Y, [FS.loc[ok, f].values]) for f in FS.columns}
GL = pd.DataFrame(GLd).T
# factor interpretation: correlation of sample scores with drivers
dv = pd.DataFrame({'winter (2025-26)': (X.season == '2025-26').astype(float), 'ice vs water': (X.medium == 'SI').astype(float), 'top vs bottom ice': X.topbot,
                   'diesel pool': (X.pool == 'EXP').astype(float), 'after dosing (diesel)': (X.phase == 'oiled').astype(float), 'day (within winter)': X.dayz, 'under-ice light': X.light})
FD = pd.DataFrame({f: {d: spearmanr(FS[f], dv[d], nan_policy='omit')[0] for d in dv} for f in FS.columns})
# guild embedding + clusters
emb2 = umap.UMAP(n_neighbors=15, min_dist=0.25, random_state=2).fit_transform(((GL - GL.mean()) / GL.std()).values)
F['gx'], F['gy'] = emb2[:, 0], emb2[:, 1]
# guilds: each genome joins the factor it loads on most, relative to all genomes (z-scored loadings). Soft clustering of the
# loadings (GMM/HDBSCAN) only split the two largest factors and hid the third, so the assignment rule is used instead.
GZ = (GL - GL.mean()) / GL.std()
def fname(f):
    r = FD[f].sort_values(key=abs, ascending=False); a = [f"{k} ({'+' if v > 0 else '−'})" for k, v in r.items() if abs(v) >= 0.4][:2]
    return ' and '.join(a) if a else r.index[0]
FNAME = {f: fname(f) for f in FS.columns}
F['guild'] = GZ.idxmax(1).reindex(F.index).values
gmeth = 'assignment to the factor with the highest standardised loading'
F['guild_name'] = [f'{f}: {FNAME[f]}' for f in F.guild]

for f in FS.columns: F[f'load_{f}'] = GL[f]
F.to_csv(O + 'fingerprints.tsv', sep='\t'); FS.join(X).to_csv(O + 'factors_samples.tsv', sep='\t'); FD.to_csv(O + 'factor_drivers.tsv', sep='\t')
GL.to_csv(O + 'guild_loadings.tsv', sep='\t'); FN.to_csv(O + 'function_pc1.tsv', sep='\t'); AL.to_csv(O + 'pc1_alignment.tsv', sep='\t')
json.dump({'fnames': FNAME, 'K': K, 've': ve, 'stab': stab, 'fmethod': cmeth, 'gmethod': gmeth, 'n_genomes': len(F), 'n_samples': len(allS), 'nperm': NP}, open(O + 'model.json', 'w'))
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(F.fcluster_name.value_counts()); print(F.guild_name.value_counts()); print(FD.round(2))
print(F[[*FR, 'r2_full']].describe().round(3))
for d in DR: print(d, 'q<0.05:', int((F[f'q_{d}'] < 0.05).sum()), 'of', int(F[f'q_{d}'].notna().sum()))
print(pd.crosstab(F.fcluster_name, F.guild_name))
