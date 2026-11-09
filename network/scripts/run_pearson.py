"""Presence Pearson: for every node pair, Pearson r of log(share of the library) over only the libraries where BOTH nodes are
detected (> 0). Zeros are left out rather than replaced, so the method asks "when both are there, do they rise and fall
together?". Input: the same pseudo-counts as the other methods (netcommon.load_layer), as shares of the library total.
A pair is tested when the two nodes are co-present in at least N_MIN (netcommon.PP_N_MIN = 10) libraries. p-value: t test with n - 2 degrees of freedom;
BH-FDR over every tested pair of the layer. Edges: q < Q (0.05). It does not know the design (winter, pool, medium, day), so like SPRING
and SPIEC-EASI its links can reflect a shared response to those.
Out: out/edges_pearson_{layer}.tsv (a, b, sign, weight = r, n = co-present libraries, q)"""
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
    E = pd.DataFrame({'a': N.index[a], 'b': N.index[b], 'sign': np.sign(rr[ok]).astype(int), 'weight': rr[ok].round(4),
                      'n': nn[ok].astype(int), 'q': q[ok]})
    E[['a', 'b']] = [edge_key(x, y) for x, y in zip(E.a, E.b)]
    E.to_csv(os.path.join(OUT, f'edges_pearson_{layer}.tsv'), sep='\t', index=False)
    print(layer, f'{len(N)} nodes, {len(lib)} libraries: {test.sum()} of {len(ri)} pairs co-present in >= {N_MIN}; {ok.sum()} edges at q < {Q} '
          f'({(E.sign > 0).sum()} positive); min |r| {np.abs(E.weight).min():.2f}; |r| >= 0.4: {(E.weight.abs() >= 0.4).sum()}, '
          f'strong (|r| >= 0.6, q < 0.001, n >= 20): {((E.weight.abs() >= 0.6) & (E.q < 1e-3) & (E.n >= 20)).sum()}')
