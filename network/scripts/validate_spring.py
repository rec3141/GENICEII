"""Validation of the Python SPRING port (spring.py). Out: out/spring_validation.tsv and printed summary.
  1. Bridge grid vs the exact truncated/truncated bridge function (random r, z1, z2): inversion error in r.
  2. Simulation from the model SPRING assumes (truncated Gaussian copula): latent X ~ N(0, Sigma) with a sparse precision
     (random graph); observed = exp(X) where X exceeds a per-variable cut giving 0-70% zeros, else 0; then made compositional
     by dividing by the row total. Recovery of Sigma (latent correlation) and of the true graph (precision, recall).
  3. Comparison with the R SPRING package (spring_ref.R) on the same simulated data: mclr, Kendall tau, latent correlation
     ('original' and 'approx'), the lambda path on the same correlation matrix, and the StARS-refitted graph.
Needs Rscript with latentcor (R_LIBS_USER=/data/scratch/grid-envs/Rlib), huge, pulsar, and a SPRING checkout (SPRING_SRC)."""
import os, sys, subprocess, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import spring
from netcommon import CACHE, OUT

SPRING_SRC = os.environ.get('SPRING_SRC', '/data/scratch/grid-cache/src_spring/SPRING')
R_LIBS = os.environ.get('R_LIBS_USER', '/data/scratch/grid-envs/Rlib')
rows = []
def rec(test, metric, value, note=''):
    rows.append(dict(test=test, metric=metric, value=value, note=note)); print(f'{test:28s} {metric:42s} {value}  {note}')

# ---------------------------------------------------------------- 1. grid vs exact
grid = spring.bridge_grid()
rng = np.random.default_rng(1)
z = rng.uniform(0, 0.9, (60, 2)); r = rng.uniform(-0.95, 0.95, 60)
K = np.array([spring.bridge_tt(a, b, c) for a, (b, c) in zip(r, z)])
rg = spring._invert_grid(K, z[:, 0], z[:, 1], grid)
rec('bridge grid', 'max |r_grid - r_true| (60 random points)', round(float(np.abs(rg - r).max()), 5))
rec('bridge grid', 'median |r_grid - r_true|', round(float(np.median(np.abs(rg - r))), 6))
rec('bridge (z = 0, 0)', 'max |F - (2/pi) asin r|', round(float(max(abs(spring.bridge_tt(x, 0.0, 0.0) - 2 / np.pi * np.arcsin(x)) for x in (-0.9, -0.3, 0.2, 0.7))), 6))

# ---------------------------------------------------------------- 2. simulation
def simulate(n, p, seed, edge_prob=0.06):
    g = np.random.default_rng(seed)
    A = np.triu(g.random((p, p)) < edge_prob, 1); A = A | A.T
    W = A * g.choice([-1, 1], (p, p)) * g.uniform(0.3, 0.6, (p, p)); W = np.triu(W, 1); W = W + W.T
    Om = W + np.eye(p) * (np.abs(np.linalg.eigvalsh(W)).max() + 0.2)
    Sig = np.linalg.inv(Om); d = np.sqrt(np.diag(Sig)); Sig = Sig / np.outer(d, d)
    Xl = g.multivariate_normal(np.zeros(p), Sig, n)
    zf = g.uniform(0, 0.7, p); cut = np.quantile(Xl, zf, axis=0).diagonal()
    Y = np.where(Xl > cut, np.exp(Xl + g.normal(0, 0.0, Xl.shape) + g.uniform(0, 3, p)), 0.0)
    Y = Y / Y.sum(1, keepdims=True)
    return Y, Sig, A

for n, p, seed in ((100, 40, 11), (200, 60, 12)):
    Y, Sig, A = simulate(n, p, seed)
    t0 = time.time(); fit = spring.spring(Y, rep_num=20, verbose=False); dt = time.time() - t0
    iu = np.triu_indices(p, 1)
    est, true = fit['adj'][iu], A[iu]
    tp = (est & true).sum()
    tag = f'simulation n={n} p={p}'
    rec(tag, 'true edges / estimated edges', f'{true.sum()} / {est.sum()}', f'StARS index {fit["opt_index"] + 1} of 20, {dt:.0f} s')
    rec(tag, 'precision', round(tp / max(est.sum(), 1), 3))
    rec(tag, 'recall', round(tp / max(true.sum(), 1), 3))
    rec(tag, 'mean |R_latent - Sigma| (off-diagonal)', round(float(np.abs(fit['R'] - Sig)[iu].mean()), 4))
    sgn = np.sign(-np.linalg.inv(Sig))[iu]                       # partial-correlation sign
    ok = est & true
    rec(tag, 'sign agreement with true partial corr.', round(float((np.sign(fit['beta'][iu][ok]) == sgn[ok]).mean()), 3))
    # naive Pearson on clr at the same edge count, for scale
    lc = np.log(Y + 1e-6); lc = lc - lc.mean(1, keepdims=True); C = np.corrcoef(lc.T)[iu]
    top = np.argsort(-np.abs(C))[:est.sum()]; tp2 = true[top].sum()
    rec(tag, 'precision, top-|Pearson clr| same count', round(tp2 / max(est.sum(), 1), 3))

# ---------------------------------------------------------------- 3. against the R package
Y, Sig, A = simulate(120, 30, 21)
p = Y.shape[1]
names = [f'v{i:02d}' for i in range(p)]
f = os.path.join(CACHE, 'spring_sim.tsv'); pref = os.path.join(CACHE, 'spring_ref')
pd.DataFrame(Y, columns=names, index=[f's{i}' for i in range(len(Y))]).to_csv(f, sep='\t')
env = dict(os.environ, R_LIBS_USER=R_LIBS)
res = subprocess.run(['Rscript', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'spring_ref.R'), SPRING_SRC, f, pref, '20'],
                     env=env, capture_output=True, text=True)
print(res.stdout[-400:], res.stderr[-800:])
rd = lambda s: pd.read_csv(pref + s, sep='\t', index_col=0).to_numpy()
iu = np.triu_indices(p, 1)
Q = spring.mclr(Y)
rec('vs R', 'max |mclr_py - mclr_R|', float(np.abs(Q - rd('_mclr.tsv')).max()))
R_py_o, K_py, _, Rpw_py_o = spring.latent_cor(Q, method='original', return_all=True)
rec('vs R', 'max |Kendall tau-a py - R| (off-diagonal)', float(np.abs(K_py - rd('_K.tsv'))[iu].max()))
rec('vs R', "max |R pointwise, 'original'| py - R", round(float(np.abs(Rpw_py_o - rd('_Rpw_original.tsv')).max()), 5), 'both root-find the exact bridge')
rec('vs R', "max |R after nearPD, 'original'| py - R", round(float(np.abs(R_py_o - rd('_R_original.tsv')).max()), 5))
R_py_a, _, _, Rpw_py_a = spring.latent_cor(Q, method='approx', return_all=True)
rec('vs R', "max |R pointwise| py grid - R 'approx'", round(float(np.abs(Rpw_py_a - rd('_Rpw_approx.tsv')).max()), 5), 'grid vs latentcor multilinear interpolant')
rec('vs R', "max |R pointwise| py grid - R 'original'", round(float(np.abs(Rpw_py_a - rd('_Rpw_original.tsv')).max()), 5))
S = rd('_R_original.tsv')
lam = spring.default_lambdas()
_, Apy = spring.mb_path(S, lam)
agree = []
for i in range(len(lam)):
    w = pd.read_csv(f'{pref}_path_{i + 1}.tsv', sep='\t')
    ar = np.zeros((p, p), bool); ar[w.iloc[:, 0] - 1, w.iloc[:, 1] - 1] = True; ar = ar | ar.T
    agree.append((int(Apy[i][iu].sum()), int(ar[iu].sum()), int((Apy[i][iu] != ar[iu]).sum())))
rec('vs R', 'lambda path on the same R: edges py / R / differing', '; '.join(f'{a}/{b}/{c}' for a, b, c in agree[::3]), 'every 3rd lambda')
rec('vs R', 'lambda path: total differing edge calls', sum(c for *_, c in agree), f'over {len(lam)} lambdas')
st = pd.read_csv(pref + '_stars.tsv', sep='\t')
fit = spring.spring(Y, rep_num=20, verbose=False, full_path=True)
w = pd.read_csv(pref + '_refit.tsv', sep='\t'); ar = np.zeros((p, p), bool); ar[w.iloc[:, 0] - 1, w.iloc[:, 1] - 1] = True; ar = ar | ar.T
ia = int(np.argmax(st['summary'].to_numpy() >= 0.1)) - 1
rec('vs R', 'StARS chosen lambda index py / R', f'{fit["opt_index"] + 1} / {max(ia, 0) + 1}', 'subsamples differ (numpy vs R RNG)')
a, b = fit['adj'][iu], ar[iu]
rec('vs R', 'refitted graph: edges py / R / shared', f'{a.sum()} / {b.sum()} / {(a & b).sum()}', f'Jaccard {(a & b).sum() / max((a | b).sum(), 1):.2f}')
rec('vs R', 'max |StARS instability py - R| over the path', round(float(np.abs(fit['summary'] - st['summary'].to_numpy()).max()), 4))
pd.DataFrame(rows).to_csv(os.path.join(OUT, 'spring_validation.tsv'), sep='\t', index=False)
