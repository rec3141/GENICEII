"""Consensus of the four methods per layer, method comparison and the known-pair check.
Methods: FlashWeave (HE-S with meta-variables; meta-variable edges are kept apart as 'drivers'), SPRING (Python port), SPIEC-EASI
(mb), residual Spearman (q < 0.05). An edge is a consensus edge when at least MIN_SUPPORT methods report it with the same sign,
no method reports it with the opposite sign (those are counted as sign conflicts and left out), and at least one of the supporting
methods accounts for the design (DESIGN_AWARE: FlashWeave with the meta-variables as nodes, or the residual Spearman). SPRING and
SPIEC-EASI do not condition on winter, pool or medium, so two of them agreeing can still be a shared response to those.
Known pairs (validation only; the methods are not tuned to them): the giant viruses contig_0000564 and contig_0000623 and the
virophage contig_0049668 with EE07 (choanoflagellate); contig_0000362 with EE04; contig_0000610 with EE22; the plastid
metabat_7332_sub, whose contigs all sit in plastid group P01, with EE05.
Out: out/edges_all.tsv (every edge any method reports except Spearman-only edges), out/edges_consensus.tsv, out/drivers.tsv,
out/method_summary.tsv, out/method_overlap.tsv, out/known_pairs.tsv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *

MIN_SUPPORT = 2
DESIGN_AWARE = ('flashweave', 'spearman')
KNOWN = [('contig_0000564', 'EE07', 'giant virus and choanoflagellate'), ('contig_0000623', 'EE07', 'giant virus and choanoflagellate'),
         ('contig_0049668', 'EE07', 'virophage and choanoflagellate'), ('contig_0000362', 'EE04', 'giant virus and cercozoan'),
         ('contig_0000610', 'EE22', 'giant virus and Micromonas'), ('P01', 'EE05', 'plastid group (holds metabat_7332_sub) and Minidiscus')]


def read_flashweave(layer, meta_cols):
    f = os.path.join(OUT, f'flashweave_{layer}.edgelist')
    E = pd.read_csv(f, sep='\t', comment='#', header=None, names=['a', 'b', 'weight'])
    meta = E.a.isin(meta_cols) | E.b.isin(meta_cols)
    Dr = E[meta].copy()
    Dr['driver'] = np.where(Dr.a.isin(meta_cols), Dr.a, Dr.b); Dr['node'] = np.where(Dr.a.isin(meta_cols), Dr.b, Dr.a)
    Dr = Dr[~(Dr.a.isin(meta_cols) & Dr.b.isin(meta_cols))]
    E = E[~meta].copy(); E['sign'] = np.sign(E.weight).astype(int)
    E[['a', 'b']] = [edge_key(a, b) for a, b in zip(E.a, E.b)]
    return E, Dr[['node', 'driver', 'weight']]


ALL, CONS, DRV, SUMM, OVL, KP = [], [], [], [], [], []
for layer in LAYERS:
    N = layer_nodes(layer)
    meta_cols = list(pd.read_csv(os.path.join(CACHE, f'{layer}_meta.tsv'), sep='\t', index_col=0).columns)
    M = {}
    M['flashweave'], dr = read_flashweave(layer, meta_cols)
    DRV.append(dr.assign(layer=layer))
    for m in ('spring', 'spieceasi', 'spearman'):
        M[m] = pd.read_csv(os.path.join(OUT, f'edges_{m}_{layer}.tsv'), sep='\t')
    for m, E in M.items():
        bad = set(E.a) | set(E.b)
        assert bad <= set(N.index), (m, layer, sorted(bad - set(N.index))[:5])
        E['key'] = list(zip(E.a, E.b))
    keys = sorted(set().union(*[set(E.key) for E in M.values()]))
    T = pd.DataFrame(index=pd.MultiIndex.from_tuples(keys, names=['a', 'b']))
    for m, E in M.items():
        s = E.set_index(['a', 'b'])
        T[f'{m}_sign'] = s['sign'].reindex(T.index).fillna(0).astype(int)
        T[f'{m}_w'] = s['weight'].reindex(T.index)
    S = T[[f'{m}_sign' for m in METHODS]]
    pos, neg = (S > 0).sum(1), (S < 0).sum(1)
    T['sign'] = np.where(pos >= neg, 1, -1)
    T['support'] = np.maximum(pos, neg)
    T['conflict'] = (pos > 0) & (neg > 0)
    T['methods'] = [','.join(m for m in METHODS if row[f'{m}_sign'] == sg) for sg, (_, row) in zip(T.sign, S.iterrows())]
    T['design_aware'] = np.column_stack([S[f'{m}_sign'].to_numpy() == T.sign.to_numpy() for m in DESIGN_AWARE]).any(1)
    T['consensus'] = (T.support >= MIN_SUPPORT) & ~T.conflict & T.design_aware
    T['layer'] = layer
    T = T.reset_index()
    ALL.append(T[(T.support >= 2) | (T.spearman_sign == 0)])          # Spearman-only edges are too many to carry
    CONS.append(T[T.consensus])
    # per-method counts and pairwise overlap (unsigned Jaccard; sign agreement on shared edges)
    for m, E in M.items():
        SUMM.append(dict(layer=layer, method=m, nodes=len(N), edges=len(E), positive=int((E.sign > 0).sum()),
                         in_consensus=int(T.consensus[T[f'{m}_sign'] != 0].sum())))
    SUMM.append(dict(layer=layer, method='consensus', nodes=len(N), edges=int(T.consensus.sum()), positive=int((T.consensus & (T.sign > 0)).sum()),
                     in_consensus=int(T.consensus.sum()), sign_conflicts=int(T.conflict.sum())))
    for i, m1 in enumerate(METHODS):
        for m2 in METHODS[i + 1:]:
            a, b = set(M[m1].key), set(M[m2].key); sh = a & b
            s1 = M[m1].set_index('key').sign; s2 = M[m2].set_index('key').sign
            agree = np.mean([s1[k] == s2[k] for k in sh]) if sh else np.nan
            OVL.append(dict(layer=layer, m1=m1, m2=m2, shared=len(sh), jaccard=round(len(sh) / max(len(a | b), 1), 4), sign_agreement=round(float(agree), 3) if sh else None))
    # known pairs
    Ti = T.set_index(['a', 'b'])
    for v, h, what in KNOWN:
        k = edge_key(v, h); row = dict(layer=layer, a=v, b=h, pair=what, in_layer=bool(v in N.index and h in N.index))
        for m in METHODS:
            sg = int(Ti.loc[k, f'{m}_sign']) if k in Ti.index else 0
            row[m] = sg
        row['support'] = int(Ti.loc[k, 'support']) if k in Ti.index else 0
        row['consensus'] = bool(Ti.loc[k, 'consensus']) if k in Ti.index else False
        # rank of the host among the virus's (or plastid's) consensus partners, by support then |weight| of SPRING
        if k in Ti.index or True:
            Cn = T[T.consensus & ((T.a == v) | (T.b == v))].copy()
            Cn['other'] = np.where(Cn.a == v, Cn.b, Cn.a)
            Cn = Cn.sort_values(['support', 'sign'], ascending=False)
            row['consensus_partners'] = len(Cn)
            row['eukaryote_partners'] = ', '.join(o for o in Cn.other if str(o).startswith('EE'))
        KP.append(row)
    print(layer, {m: len(E) for m, E in M.items()}, 'consensus', int(T.consensus.sum()), 'conflicts', int(T.conflict.sum()))

pd.concat(ALL).to_csv(os.path.join(OUT, 'edges_all.tsv'), sep='\t', index=False)
pd.concat(CONS).to_csv(os.path.join(OUT, 'edges_consensus.tsv'), sep='\t', index=False)
pd.concat(DRV).to_csv(os.path.join(OUT, 'drivers.tsv'), sep='\t', index=False)
pd.DataFrame(SUMM).to_csv(os.path.join(OUT, 'method_summary.tsv'), sep='\t', index=False)
pd.DataFrame(OVL).to_csv(os.path.join(OUT, 'method_overlap.tsv'), sep='\t', index=False)
K = pd.DataFrame(KP); K.to_csv(os.path.join(OUT, 'known_pairs.tsv'), sep='\t', index=False)
pd.set_option('display.width', 250)
print(pd.DataFrame(SUMM).to_string(index=False)); print(pd.DataFrame(OVL).to_string(index=False)); print(K.drop(columns=['eukaryote_partners']).to_string(index=False))
