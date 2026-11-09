"""SPRING (Semi-Parametric Rank-based approach for INference in Graphical models) in Python.

Port of the R package SPRING (Yoon, Gaynanova & Mueller 2019, Front. Genet. 10:516, doi:10.3389/fgene.2019.00516;
https://github.com/GraceYoon/SPRING, R/SPRING.R and R/helpers.R) and the pieces it calls:
  - latentcor (Huang, Mueller & Gaynanova 2021, https://github.com/mingzehuang/latentcor; formerly mixedCCA,
    https://github.com/irinagain/mixedCCA): Kendall tau-a -> latent Gaussian-copula correlation through the
    truncated/truncated ("22", "tt") bridge function of Yoon, Carroll & Gaynanova (2020, Biometrika 107:609);
  - Matrix::nearPD (Higham 2002 alternating projections with Dykstra correction, corr = TRUE) plus the nu shrinkage;
  - huge::huge.mb (Meinshausen-Buhlmann neighbourhood selection on a correlation-matrix input, sym = "or");
  - pulsar::pulsar / stars.stability / refit (StARS: Liu, Roeder & Wasserman 2010).

Steps, as in SPRING.default(data, quantitative = FALSE):
  1. mclr: log of the non-zero entries centred by the mean log of the sample's non-zeros, then every non-zero shifted by
     |min| + 1 so non-zeros are strictly positive and zeros stay 0 (helpers.R mclr()).
  2. latent correlation per pair from Kendall tau-a and both variables' zero fractions (type "tru"): F_TT(r; z1, z2) = tau.
     The inverse is read off a precomputed grid of F_TT (r in [-0.999, 0.999] x z1 x z2; bilinear in z, linear inverse in r).
     latentcor's "approx" switches to root finding (r_sol) when |tau| > ratio x bound (bound = 1 - max(z1, z2)^2) because its
     stored interpolant only covers that range; the same switch is kept here (exact_tail=True, the default), and pairs with a
     zero fraction above the grid's 0.95 also go to root finding. method='original' uses root finding for every pair, as
     latentcor method = "original".
  3. nearest positive-definite correlation matrix, then R = (1 - nu) R + nu I (nu = 0.001).
  4. neighbourhood-selection lasso along a decreasing lambda path (default exp(seq(log 0.6, log 0.006, length 20))),
     coordinate descent on the correlation matrix, edges symmetrised with "or".
  5. StARS: rep_num subsamples of floor(0.8 n) samples; edge instability 4 sum theta(1 - theta) / (p (p - 1)) per lambda,
     monotonised with cummax; the chosen index is the last lambda before instability first reaches thresh (0.1); the graph is
     refitted on all samples at that lambda. Lambdas past the first crossing cannot change the choice, so they are skipped
     (full_path=True computes them, as pulsar does).

Differences from R: the subsamples come from numpy's RNG (R's set.seed stream cannot be reproduced), so the StARS curve
differs by subsampling noise; the latent correlation and the full-data path are deterministic and are compared directly in
validate_spring.py. The 4-variate normal CDF comes from scipy (Genz) instead of mnormt (also Genz); the interpolation grid
replaces latentcor's stored multilinear interpolant.
"""
import os
import numpy as np
from scipy.stats import norm, multivariate_normal
from scipy.optimize import brentq
import numba

GRID_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cache', 'bridge_tt_grid.npz')
S2 = 1 / np.sqrt(2)


# ---------------------------------------------------------------- 1. mclr
def mclr(X, tol=1e-16, atleast=1.0):
    """SPRING helpers.R mclr(): rows = samples. Zeros stay zero; non-zeros get clr over non-zeros, then a common positive shift."""
    X = np.asarray(X, float)
    nz = X >= tol
    L = np.where(nz, np.log(np.where(nz, X, 1.0)), 0.0)
    gm = L.sum(1, keepdims=True) / np.maximum(nz.sum(1, keepdims=True), 1)
    C = np.where(nz, L - gm, 0.0)
    shift = abs(C.min()) + atleast if C.min() < 0 else 0.0
    return np.where(nz, C + shift, 0.0)


# ---------------------------------------------------------------- 2. Kendall tau-a and the truncated/truncated bridge
def kendall_tau_a(X):
    """Kendall tau-a for all column pairs: sum over sample pairs of sign(dx) sign(dy) / (n (n - 1) / 2) (latentcor Kendalltau)."""
    X = np.asarray(X, float); n = X.shape[0]
    i, j = np.triu_indices(n, 1)
    Sg = np.sign(X[i] - X[j])                                  # (n0, p)
    return (Sg.T @ Sg) / len(i)


def bridge_tt(r, z1, z2):
    """latentcor internal.R comb "22" (mixedCCA bridgeF_tt): Kendall tau implied by latent correlation r for two truncated variables
    with zero fractions z1, z2."""
    d1, d2 = norm.ppf(z1), norm.ppf(z2)
    m1 = np.array([[1, 0, S2, -r * S2], [0, 1, -r * S2, S2], [S2, -r * S2, 1, -r], [-r * S2, S2, -r, 1]])
    m2 = np.array([[1, r, S2, r * S2], [r, 1, r * S2, S2], [S2, r * S2, 1, r], [r * S2, S2, r, 1]])
    up = np.array([-d1, -d2, 0, 0])
    kw = dict(abseps=1e-5, releps=1e-5)                        # scipy's Genz defaults; QMC error about 1e-5 in tau
    return -2 * multivariate_normal.cdf(up, np.zeros(4), m1, allow_singular=True, **kw) + \
        2 * multivariate_normal.cdf(up, np.zeros(4), m2, allow_singular=True, **kw)


def bound_tt(z1, z2):
    return 1 - np.maximum(z1, z2) ** 2


def _grid_axes():
    zg = np.round(np.linspace(0, 0.95, 39), 6)
    rg = np.unique(np.round(np.concatenate([np.linspace(-0.999, 0.999, 121), np.linspace(-0.999, -0.9, 12), np.linspace(0.9, 0.999, 12)]), 6))
    return zg, rg


def _grid_row(a, b, zg, rg):
    return a, b, np.array([bridge_tt(r, zg[a], zg[b]) for r in rg])


def bridge_grid(n_jobs=14, path=GRID_FILE):
    """F_TT on (z1, z2, r); computed once (about a minute on 14 cores) and cached. F is symmetric in (z1, z2) and increasing in r."""
    if os.path.exists(path):
        g = np.load(path); return g['zg'], g['rg'], g['F']
    from joblib import Parallel, delayed
    zg, rg = _grid_axes(); nz = len(zg)
    jobs = [(a, b) for a in range(nz) for b in range(a, nz)]
    out = Parallel(n_jobs=n_jobs)(delayed(_grid_row)(a, b, zg, rg) for a, b in jobs)
    F = np.zeros((nz, nz, len(rg)))
    for a, b, f in out:
        f = np.maximum.accumulate(f)                           # QMC noise must not break monotonicity in r
        F[a, b] = F[b, a] = f
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(path, zg=zg, rg=rg, F=F)
    return zg, rg, F


def _invert_grid(K, z1, z2, grid, chunk=20000):
    """Bilinear interpolation of the F(r) curve at (z1, z2), then linear inversion in r (vectorised over pairs, in chunks)."""
    if len(K) > chunk:
        return np.concatenate([_invert_grid(K[i:i + chunk], z1[i:i + chunk], z2[i:i + chunk], grid) for i in range(0, len(K), chunk)])
    zg, rg, F = grid
    z1 = np.clip(z1, zg[0], zg[-1]); z2 = np.clip(z2, zg[0], zg[-1])
    a = np.clip(np.searchsorted(zg, z1) - 1, 0, len(zg) - 2); b = np.clip(np.searchsorted(zg, z2) - 1, 0, len(zg) - 2)
    ta = ((z1 - zg[a]) / (zg[a + 1] - zg[a]))[:, None]; tb = ((z2 - zg[b]) / (zg[b + 1] - zg[b]))[:, None]
    Fc = (1 - ta) * (1 - tb) * F[a, b] + ta * (1 - tb) * F[a + 1, b] + (1 - ta) * tb * F[a, b + 1] + ta * tb * F[a + 1, b + 1]
    out = np.empty(len(K))
    for i in range(len(K)):                                   # curves differ per pair; np.interp per pair is fast enough
        out[i] = np.interp(K[i], Fc[i], rg)
    return out


def _invert_exact(k, z1, z2, tol=1e-6, r0=None, half=0.03):
    """latentcor r_sol: solve F_TT(r) = K on (-0.999, 0.999); a K beyond the attainable range returns the end point. With a starting
    value r0 (the grid estimate) the root is first bracketed in r0 +/- half, which needs a handful of bridge evaluations."""
    f = lambda r: bridge_tt(r, z1, z2) - k
    if max(z1, z2) >= 1: return 0.0                            # a variable that is all zeros carries no rank information
    if r0 is not None:
        a, b = max(-0.999, r0 - half), min(0.999, r0 + half)
        fa, fb = f(a), f(b)
        if fa <= 0 <= fb: return brentq(f, a, b, xtol=tol) if fa < 0 < fb else (a if fa == 0 else b)
    lo, hi = f(-0.999), f(0.999)
    if lo >= 0: return -0.999
    if hi <= 0: return 0.999
    return brentq(f, -0.999, 0.999, xtol=tol)


def latent_cor(X, method='approx', ratio=0.9, use_nearpd=True, nu=0.001, grid=None, return_all=False, exact_tail=True, n_jobs=14):
    """latentcor(X, types = "tru") for non-negative data whose zeros are the truncation point (mclr output).
    Returns R (after nearPD and nu shrinkage when use_nearpd) and optionally K (tau-a), zratios, R pointwise."""
    X = np.asarray(X, float); p = X.shape[1]
    X = X - X.min(0)                                           # encodeX: truncated columns shifted to a minimum of 0
    z = (X == 0).mean(0)
    K = kendall_tau_a(X)
    iu = np.triu_indices(p, 1); k = K[iu]; z1, z2 = z[iu[0]], z[iu[1]]
    r = np.empty(len(k))
    r0 = np.full(len(k), np.nan)
    if method == 'original':
        exact = np.ones(len(k), bool)
    else:
        grid = grid or bridge_grid()
        exact = np.maximum(z1, z2) > grid[0][-1]                   # sparser than the grid covers: root finding
        if exact_tail: exact |= np.abs(k) > ratio * bound_tt(z1, z2)  # latentcor 'approx': r_sol above ratio x bound
        inside = ~(np.maximum(z1, z2) > grid[0][-1])
        r0[inside] = _invert_grid(k[inside], z1[inside], z2[inside], grid)   # grid value, and the start for root finding
        r[~exact] = r0[~exact]
    ie = np.where(exact)[0]
    st = [None if np.isnan(r0[i]) else float(r0[i]) for i in ie]
    if len(ie) > 50 and n_jobs > 1:
        from joblib import Parallel, delayed
        r[ie] = Parallel(n_jobs=n_jobs, batch_size=16)(delayed(_invert_exact)(k[i], z1[i], z2[i], 1e-6, s0) for i, s0 in zip(ie, st))
    else:
        for i, s0 in zip(ie, st): r[i] = _invert_exact(k[i], z1[i], z2[i], 1e-6, s0)
    R = np.eye(p); R[iu] = r; R[(iu[1], iu[0])] = r
    Rp = R.copy()
    if use_nearpd:
        R = near_pd_corr(R)
        R = (1 - nu) * R + nu * np.eye(p)
    return (R, K, z, Rp) if return_all else R


# ---------------------------------------------------------------- 3. nearest positive-definite correlation matrix
def near_pd_corr(x, eig_tol=1e-6, conv_tol=1e-7, posd_tol=1e-8, maxit=100):
    """Matrix::nearPD(x, corr = TRUE) defaults: Higham (2002) alternating projections with Dykstra's correction, unit diagonal,
    convergence on the infinity norm, then the 'do2eigen' step that lifts eigenvalues below posd_tol x the largest."""
    X = np.array(x, float); n = len(X); DS = np.zeros_like(X)
    for _ in range(maxit):
        Y = X; R = Y - DS
        d, Q = np.linalg.eigh(R); d, Q = d[::-1], Q[:, ::-1]
        keep = d > eig_tol * d[0]
        X = (Q[:, keep] * d[keep]) @ Q[:, keep].T
        DS = X - R
        np.fill_diagonal(X, 1.0)
        conv = np.abs(Y - X).sum(1).max() / np.abs(Y).sum(1).max()
        if conv <= conv_tol: break
    d, Q = np.linalg.eigh(X); d, Q = d[::-1], Q[:, ::-1]
    eps = posd_tol * abs(d[0])
    if d[-1] < eps:
        d = np.maximum(d, eps); o = np.diag(X).copy()
        X = (Q * d) @ Q.T
        D = np.sqrt(np.maximum(eps, o) / np.diag(X))
        X = D[:, None] * X * D[None, :]
    np.fill_diagonal(X, 1.0)
    return (X + X.T) / 2


# ---------------------------------------------------------------- 4. Meinshausen-Buhlmann neighbourhood selection
@numba.njit(cache=True)
def _cd_node(S, j, b, g, lam, tol, maxit):
    """Lasso for node j at one lambda, warm-started from b (g = S b kept current). Active-set coordinate descent: sweep the
    non-zero coefficients until they settle, then one full sweep; stop when a full sweep changes nothing by more than tol."""
    p = S.shape[0]
    for it in range(maxit):
        dmax = 0.0
        for k in range(p):                                     # full sweep
            if k == j: continue
            rho = S[k, j] - (g[k] - S[k, k] * b[k])
            if rho > lam: nb = (rho - lam) / S[k, k]
            elif rho < -lam: nb = (rho + lam) / S[k, k]
            else: nb = 0.0
            dd = nb - b[k]
            if dd != 0.0:
                for m in range(p): g[m] += S[m, k] * dd
                b[k] = nb
                if abs(dd) > dmax: dmax = abs(dd)
        if dmax < tol: return
        act = np.where(b != 0)[0]
        for it2 in range(maxit):                               # active-set sweeps
            dmax2 = 0.0
            for k in act:
                rho = S[k, j] - (g[k] - S[k, k] * b[k])
                if rho > lam: nb = (rho - lam) / S[k, k]
                elif rho < -lam: nb = (rho + lam) / S[k, k]
                else: nb = 0.0
                dd = nb - b[k]
                if dd != 0.0:
                    for m in range(p): g[m] += S[m, k] * dd
                    b[k] = nb
                    if abs(dd) > dmax2: dmax2 = abs(dd)
            if dmax2 < tol: break


@numba.njit(parallel=True, cache=True)
def _mb_step(S, B, lam, tol, maxit):
    """All p neighbourhood regressions at one lambda, in place: B[k, j] = coefficient of node k in node j's regression, warm start
    from the incoming B (huge's SPMBgraph solves the same problem: 0.5 b' S_{-j,-j} b - S_{-j,j}' b + lam |b|_1)."""
    p = S.shape[0]
    for j in numba.prange(p):
        b = B[:, j].copy()
        g = S @ b
        _cd_node(S, j, b, g, lam, tol, maxit)
        B[:, j] = b


def mb_path(S, lambdas, tol=1e-7, maxit=100000):
    """huge.mb(S, lambda, sym = "or"): coefficient arrays and symmetrised adjacency per lambda (decreasing, warm starts)."""
    S = np.ascontiguousarray(S, dtype=np.float64); p = len(S)
    B = np.zeros((p, p)); Bs = np.zeros((len(lambdas), p, p))
    for i, lam in enumerate(lambdas):
        _mb_step(S, B, float(lam), tol, maxit); Bs[i] = B
    nzb = Bs != 0
    return Bs, nzb | nzb.transpose(0, 2, 1)


# ---------------------------------------------------------------- 5. StARS (pulsar) and the SPRING driver
def default_lambdas(nlambda=20, lambda_min_ratio=1e-2, lmax=0.6):
    return np.exp(np.linspace(np.log(lmax), np.log(lmax * lambda_min_ratio), nlambda))


def stars_index(summary, thresh):
    """pulsar .starsind: last index before the monotonised instability first reaches thresh (at least the first); if it never
    does, the last (densest) lambda, with a warning."""
    hit = np.where(summary >= thresh)[0]
    if len(hit) == 0:
        print('SPRING/StARS warning: optimal lambda may be outside the supplied values'); return len(summary) - 1
    return max(hit[0] - 1, 0)


def spring(data, quantitative=False, lambdas=None, nlambda=20, lambda_min_ratio=1e-2, thresh=0.1, subsample_ratio=0.8, rep_num=20,
           seed=10010, method='approx', nu=0.001, ratio=0.9, verbose=True, full_path=False, tol=1e-7, maxit=100000):
    """SPRING.default. data: n samples x p taxa (counts or compositions). Returns a dict with the refitted adjacency, the
    symmetrised coefficients at the chosen lambda (mean of b_kj and b_jk; their sign is the sign of the partial association),
    the full-data latent correlation R, the lambda path, StARS summary and edge selection frequencies at the chosen lambda."""
    X = np.asarray(data, float)
    if (X < 0).any(): raise ValueError('negative values: counts or compositions expected')
    Q = X if quantitative else mclr(X)
    n, p = Q.shape
    lam = default_lambdas(nlambda, lambda_min_ratio) if lambdas is None else np.asarray(lambdas, float)
    grid = bridge_grid() if method == 'approx' else None
    rng = np.random.default_rng(seed); m = int(np.floor(n * subsample_ratio))
    subs = [rng.choice(n, m, replace=False) for _ in range(rep_num)]
    Rs = [np.ascontiguousarray(latent_cor(Q[idx], method=method, nu=nu, ratio=ratio, grid=grid)) for idx in subs]
    Bsub = [np.zeros((p, p)) for _ in range(rep_num)]
    # lambda by lambda over all subsamples (warm starts per subsample). pulsar's choice depends only on the instability up to its
    # first crossing of thresh (cummax), so the path stops there unless full_path=True; the result is the same.
    merge = np.full((len(lam), p, p), np.nan); summary = np.full(len(lam), np.nan)
    for i, l in enumerate(lam):
        acc = np.zeros((p, p))
        for s in range(rep_num):
            _mb_step(Rs[s], Bsub[s], float(l), tol, maxit)
            nz = Bsub[s] != 0; acc += nz | nz.T
        merge[i] = acc / rep_num
        summary[i] = 4 * (merge[i] * (1 - merge[i])).sum() / (p * (p - 1))
        if verbose: print(f'  SPRING lambda {i + 1}/{len(lam)} ({l:.4f}): instability {summary[i]:.4f}', flush=True)
        if not full_path and np.nanmax(summary[:i + 1]) >= thresh: break
    done = ~np.isnan(summary)
    summary = np.where(done, np.maximum.accumulate(np.where(done, summary, 0)), np.nan)
    opt = stars_index(summary[done], thresh) if done.all() or np.nanmax(summary) >= thresh else int(done.sum()) - 1
    R = latent_cor(Q, method=method, nu=nu, ratio=ratio, grid=grid)
    B, A = mb_path(R, lam[:max(opt + 1, int(done.sum()))], tol=tol, maxit=maxit)
    Bs = (B[opt] + B[opt].T) / 2
    return {'adj': A[opt], 'beta': Bs, 'beta_raw': B[opt], 'R': R, 'lambdas': lam, 'summary': summary, 'opt_index': opt,
            'stability': merge[opt], 'path': A}
