"""Presence Pearson: for every node pair, Pearson r of log(share of the library) over only the libraries where BOTH nodes are
detected (> 0). Zeros are left out rather than replaced, so the method asks "when both are there, do they rise and fall
together?". Input: the same pseudo-counts as the other methods (netcommon.load_layer), as shares of the library total.
A pair is tested when the two nodes are co-present in at least N_MIN (netcommon.PP_N_MIN = 10) libraries. p-value: t test with n - 2 degrees of freedom;
BH-FDR over every tested pair of the layer. Edge call (netcommon): q < PP_Q (0.05) AND |r| >= PP_R (0.85) AND co-present in
>= PP_N (20) libraries; q alone is far too loose (see netcommon). It does not know the design (winter, pool, medium, day), so like SPRING
and SPIEC-EASI its links can reflect a shared response to those.
Out: out/edges_pearson_{layer}.tsv (the edge call: a, b, sign, weight = r, n = co-present libraries, q);
cache/pearson_q05_{layer}.tsv (every pair at q < 0.05, for out/pearson_thresholds.tsv, the edge counts at other thresholds)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
from scipy.stats import t as tdist
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import bh
N_MIN, Q = PP_N_MIN, PP_Q


def presence_pearson(S):
    """S: libraries x nodes shares. Returns r, n (nodes x nodes), pairwise over the libraries where both are > 0."""
    P = (S > 0).astype(float)
    L = np.where(S > 0, np.log(np.where(S > 0, S, 1.0)), 0.0)
    n = P.T @ P                                   # co-present libraries
    sx = L.T @ P                                  # sx[i, j]: sum of log i over libraries where j is present (i is present too, or L = 0)
    sxx = (L * L).T @ P
    sxy = L.T @ L
    with np.errstate(invalid='ignore', divide='ignore'):
        cov = n * sxy - sx * sx.T
        var = (n * sxx - sx ** 2) * (n * sxx - sx ** 2).T
        r = cov / np.sqrt(var)
    return r, n


TH = []
for layer in (sys.argv[1:] or LAYERS):
    X, lib, N = load_layer(layer)
    r, n = presence_pearson(rel_abund(X, lib, layer).to_numpy(float))
    iu = np.triu_indices(len(N), 1); ri, ni = r[iu], n[iu]
    test = (ni >= N_MIN) & np.isfinite(ri)
    rr = np.clip(ri[test], -0.999999, 0.999999); nn = ni[test]
    p = 2 * tdist.sf(np.abs(rr) * np.sqrt((nn - 2) / (1 - rr ** 2)), nn - 2)
    q = bh(p)
    ok = q < Q
    a, b = iu[0][test][ok], iu[1][test][ok]
    A = pd.DataFrame({'a': N.index[a], 'b': N.index[b], 'sign': np.sign(rr[ok]).astype(int), 'weight': rr[ok].round(4),
                      'n': nn[ok].astype(int), 'q': q[ok]})
    A[['a', 'b']] = [edge_key(x, y) for x, y in zip(A.a, A.b)]
    A.to_csv(os.path.join(CACHE, f'pearson_q05_{layer}.tsv'), sep='\t', index=False)
    E = A[(A.weight.abs() >= PP_R) & (A.n >= PP_N)]
    E.to_csv(os.path.join(OUT, f'edges_pearson_{layer}.tsv'), sep='\t', index=False)
    for nmin in (10, 20, 30):
        for r in (0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9):
            for qq in (0.05, 1e-3, 1e-5):
                k = int(((A.weight.abs() >= r) & (A.q < qq) & (A.n >= nmin)).sum())
                TH.append(dict(layer=layer, n_min=nmin, r_min=r, q_max=qq, edges=k, share_of_pairs=round(k / len(ri), 4)))
    print(layer, f'{len(N)} nodes, {len(lib)} libraries: {test.sum()} of {len(ri)} pairs co-present in >= {N_MIN}; {ok.sum()} at q < {Q}; '
          f'edges (|r| >= {PP_R}, n >= {PP_N}): {len(E)} ({(E.sign > 0).sum()} positive, {len(E) / len(ri):.1%} of pairs)')
pd.DataFrame(TH).to_csv(os.path.join(OUT, 'pearson_thresholds.tsv'), sep='\t', index=False)
