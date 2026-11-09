"""Force-directed node positions per layer: one layout per method, built from that method's own edges as the page exports them
(netcommon.edge_table, so after the Spearman and Presence Pearson floors), and one for the consensus edges (consensus.py).
Only positive links pull, and of those each genome keeps its TOP_K strongest (a link stays when it is among the top of either
genome), so the dense edge sets give modules rather than one ball: a negative link says two genomes rarely peak together, which is no reason to draw them side by side.
Edge weight (attraction) = the method's score: FlashWeave |weight|, SPRING and SPIEC-EASI StARS stability, residual Spearman |rho|,
Presence Pearson |r|; consensus edges weigh the number of methods supporting them plus the mean percentile of those methods' scores. Each method's scores are divided by their
median, so every layout sees weights around 1.
Layout: igraph Fruchterman-Reingold (R, forcelayout.R), started from the t-SNE map, 2000 iterations, fixed seed, so the result is
the same on every run. Genomes without a positive link in a layout sit on rings around it (outside the drawn extent), ordered by
domain and then by abundance, instead of piling up in the middle.
Out: out/forcelayout_{layer}.tsv (<layout>_x, <layout>_y per node)"""
import os, sys, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *

DOM_ORDER = {'eukaryote': 0, 'virus': 1, 'organelle': 2, 'prokaryote': 3}
RING_GAP = 0.034            # spacing of ringed nodes, in units of the laid-out radius
TOP_K = 10                  # each genome keeps its TOP_K strongest positive links for the layout


def backbone(i, j, w, k=TOP_K):
    """Rows kept: links among the k highest-weight links of at least one of their two genomes. Dense edge sets (the residual
    Spearman, the Presence Pearson, the RNA consensus: tens of thousands of positive links) otherwise pull every genome to every
    other and the layout is one round ball; sparse ones (FlashWeave) are hardly changed."""
    D = pd.DataFrame({'n': np.r_[i, j], 'w': np.r_[w, w], 'row': np.r_[np.arange(len(i)), np.arange(len(i))]})
    top = D.sort_values(['n', 'w'], ascending=[True, False]).groupby('n').head(k)
    keep = np.zeros(len(i), bool); keep[top.row.to_numpy()] = True
    return keep


def score(T, m):
    if m in ('spring', 'spieceasi'): return T[f'{m}_stab'].to_numpy(float)
    return np.abs(T[f'{m}_w'].to_numpy(float))


def ring(ids, N, ab, R0):
    """Positions on concentric rings outside radius R0, ordered by domain then abundance (largest first)."""
    order = sorted(ids, key=lambda i: (DOM_ORDER.get(N.domain[i], 9), -ab[i]))
    out, k, rad = {}, 0, R0 * 1.10
    while k < len(order):
        cap = max(12, int(2 * np.pi * rad / (RING_GAP * R0)))
        chunk = order[k:k + cap]
        for t, i in enumerate(chunk):
            a = np.pi / 2 - 2 * np.pi * t / len(chunk)
            out[i] = (rad * np.cos(a), rad * np.sin(a))
        k += cap; rad += R0 * 0.06
    return out


for layer in LAYERS:
    X, lib, N = load_layer(layer)
    ab = rel_abund(X, lib, layer).median()
    ids = list(N.index); ix = {k: i + 1 for i, k in enumerate(ids)}
    T, _ = edge_table(layer)
    C = pd.read_csv(os.path.join(OUT, 'edges_consensus.tsv'), sep='\t').query('layer == @layer')
    rows = []
    for m in METHODS:
        E = T[T[f'{m}_sign'] > 0]
        w = np.maximum(score(E, m), 0.02); w = w / np.median(w)        # SPRING keeps a few refit edges with stability 0
        rows.append(pd.DataFrame({'layout': m, 'i': E.a.map(ix), 'j': E.b.map(ix), 'w': w.round(4)}))
    # consensus: weight = number of supporting methods plus the mean percentile of their scores (0-1), which breaks the many ties
    Cp = C[C.sign > 0][['a', 'b', 'support']].merge(T, on=['a', 'b'], how='left')
    pct = np.column_stack([np.where(Cp[f'{m}_sign'] > 0, pd.Series(score(T, m)).where(T[f'{m}_sign'] > 0).rank(pct=True)
                                    .set_axis(pd.MultiIndex.from_frame(T[['a', 'b']])).reindex(pd.MultiIndex.from_frame(Cp[['a', 'b']])).to_numpy(), np.nan)
                           for m in METHODS])
    w = Cp.support.to_numpy() + np.nanmean(pct, axis=1) * 0.99
    rows.append(pd.DataFrame({'layout': 'consensus', 'i': Cp.a.map(ix), 'j': Cp.b.map(ix), 'w': (w / np.median(w)).round(4)}))
    Ein = pd.concat(rows)
    Ein = pd.concat([g[backbone(g.i.to_numpy(), g.j.to_numpy(), g.w.to_numpy())] for _, g in Ein.groupby('layout', sort=False)])
    tsne = pd.read_csv(os.path.join(OUT, f'layout_{layer}.tsv'), sep='\t', index_col=0).loc[ids]
    fe, fi, fo = (os.path.join(CACHE, f'forcelayout_{layer}_{s}.tsv') for s in ('edges', 'init', 'out'))
    Ein.to_csv(fe, sep='\t', index=False); tsne[['x', 'y']].to_csv(fi, sep='\t', index=False)
    print(layer, flush=True)
    subprocess.run(['Rscript', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'forcelayout.R'), fe, fi, fo], check=True)
    R = pd.read_csv(fo, sep='\t')
    out = pd.DataFrame(index=ids)
    for k in list(METHODS) + ['consensus']:
        r = R[R.layout == k]
        xy = r[['x', 'y']].to_numpy(float); xy = xy - np.median(xy, axis=0)
        R0 = np.quantile(np.hypot(xy[:, 0], xy[:, 1]), 0.99)       # a stray component far out should not set the ring radius
        pos = {ids[n - 1]: tuple(p) for n, p in zip(r.node, xy)}
        iso = [i for i in ids if i not in pos]
        pos |= ring(iso, N, ab, R0)
        out[f'{k}_x'] = [round(pos[i][0] / R0, 4) for i in ids]; out[f'{k}_y'] = [round(pos[i][1] / R0, 4) for i in ids]
        print(f'  {k}: {len(r)} genomes linked, {len(iso)} on the rings')
    out.index.name = 'node'
    out.to_csv(os.path.join(OUT, f'forcelayout_{layer}.tsv'), sep='\t')
