"""Design-residual Spearman baseline. Each node's centred log-ratio (log of reads or bases per million + 1, so zeros stay finite,
minus the library's mean over the layer's nodes; this removes library-wide shifts such as the share of reads on binned contigs) is
regressed by least squares on the design (netcommon.design: winter x pool x medium cells plus a day trend within each winter); the
residuals are ranked and correlated (Spearman) for every node pair. p-values from the t approximation with n - rank(design) - 2
degrees of freedom; BH-FDR over all pairs of the layer. Edges: q < Q. Out: out/edges_spearman_{layer}.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
from scipy.stats import t as tdist, rankdata
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import bh
Q = 0.05

for layer in LAYERS:
    X, lib, N = load_layer(layer)
    Y = np.log(X.div(lib['size'] / (1e6 if layer == 'rna' else 1e6 * DNA_UNIT), axis=0) + 1)   # per million (reads or bases)
    Y = Y.sub(Y.mean(axis=1), axis=0)                                                         # centred log-ratio per library
    Z = design(lib).to_numpy()
    beta, *_ = np.linalg.lstsq(Z, Y.to_numpy(), rcond=None)
    Rz = Y.to_numpy() - Z @ beta
    Rk = np.apply_along_axis(rankdata, 0, Rz)
    Rk = (Rk - Rk.mean(0)) / Rk.std(0)
    rho = Rk.T @ Rk / len(Rk)
    dof = len(Rk) - np.linalg.matrix_rank(Z) - 2
    iu = np.triu_indices(len(N), 1); r = np.clip(rho[iu], -0.999999, 0.999999)
    p = 2 * tdist.sf(np.abs(r) * np.sqrt(dof / (1 - r ** 2)), dof)
    q = bh(p)
    ok = q < Q
    E = pd.DataFrame({'a': N.index[iu[0][ok]], 'b': N.index[iu[1][ok]], 'sign': np.sign(r[ok]).astype(int), 'weight': r[ok].round(4), 'q': q[ok]})
    E[['a', 'b']] = [edge_key(a, b) for a, b in zip(E.a, E.b)]
    E.to_csv(os.path.join(OUT, f'edges_spearman_{layer}.tsv'), sep='\t', index=False)
    print(layer, f'{len(N)} nodes, {len(lib)} libraries, dof {dof}: {ok.sum()} edges at q < {Q} ({(E.sign > 0).sum()} positive); min |rho| {np.abs(E.weight).min():.2f}')
