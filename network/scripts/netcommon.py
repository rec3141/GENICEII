"""Shared definitions for the network scripts: layer tables after the prevalence rule, meta-variables, design matrix, paths.

Prevalence rule, per layer (nanopore DNA, Illumina DNA, RNA): a node enters a layer's network if it is detected (> 0) in at least PREV of that layer's
libraries and averages at least MIN_PPM per million (RNA and Illumina DNA: contig-assigned read pairs; nanopore DNA: mapped bases). EE MAGs, giant viruses,
virophages and organelle groups are always kept where data exist (detected in at least PREV_FORCED of the libraries).
Nodes without contigs (plastid bins wholly inside an organelle group, e.g. metabat_7332_sub inside P01) are dropped."""
import os
import numpy as np, pandas as pd

NET = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
CACHE, OUT = os.path.join(NET, 'cache'), os.path.join(NET, 'out')
REPO = os.path.abspath(os.path.join(NET, '..'))
PREV, MIN_PPM, PREV_FORCED = 0.25, 1.0, 0.10
DNA_UNIT = 1000.0                      # DNA pseudo-counts: kb of aligned nanopore read bases (coverage depth x length / 1000)
# Layers. Each is one abundance table in cache/ (nodes x libraries, written by nodes.py) plus a library table with 'size' (the
# library total in the same unit) and the design columns season, group, medium, day, layer, par_under. 'unit' turns the stored
# values into the pseudo-counts the methods get (value / unit). To add a layer (e.g. Illumina DNA), write its two tables and add
# an entry here; every script, the R and Julia runs (they read cache/*_X.tsv) and the page follow. 'page' sets the page order.
LAYER_CFG = {
    'dna': dict(counts='dna_cov.tsv', lib='dna_lib.tsv', unit=DNA_UNIT, label='Nanopore DNA', long='Nanopore DNA (presence)', page=2,
                what='nanopore metagenomes'),
    # Illumina metagenomes (nodes.py, 2026-11-09): read pairs on the node's contigs, CMO pool libraries; pseudo-counts = read pairs
    'mgdna': dict(counts='mgdna_counts.tsv', lib='mgdna_lib.tsv', unit=1.0, label='Illumina DNA', long='Illumina DNA (presence)', page=3,
                  what='Illumina metagenomes'),
    'rna': dict(counts='rna_counts.tsv', lib='rna_lib.tsv', unit=1.0, label='RNA', long='RNA (activity)', page=1,
                what='Illumina RNA libraries'),
}
LAYERS = tuple(LAYER_CFG)
METHODS = ('flashweave', 'spring', 'spieceasi', 'spearman', 'pearson')
METHOD_LABEL = {'flashweave': 'FlashWeave', 'spring': 'SPRING', 'spieceasi': 'SPIEC-EASI', 'spearman': 'Residual Spearman',
                'pearson': 'Presence Pearson'}
METHOD_ABBR = {'flashweave': 'F', 'spring': 'S', 'spieceasi': 'E', 'spearman': 'R', 'pearson': 'P'}
DESIGN_AWARE = ('flashweave', 'spearman')     # methods that account for winter, pool, medium and time


def nodes_all():
    return pd.read_csv(os.path.join(OUT, 'nodes_all.tsv'), sep='\t', index_col=0, keep_default_na=False, na_values=[''])


def forced(N):
    return (N.domain == 'eukaryote') | N.kind.isin(['giant virus', 'virophage', 'organelle group'])


def raw_table(layer):
    """Stored abundance (nodes x libraries) and library table of a layer."""
    c = LAYER_CFG[layer]
    return (pd.read_csv(os.path.join(CACHE, c['counts']), sep='\t', index_col=0),
            pd.read_csv(os.path.join(CACHE, c['lib']), sep='\t', index_col=0))


def layer_nodes(layer, N=None):
    N = nodes_all() if N is None else N
    N = N[N.contigs > 0].copy()
    if f'{layer}_prev' not in N:                      # a layer added after nodes.py: prevalence from its table
        R, lib = raw_table(layer); R = R.reindex(N.index).fillna(0)[lib.index]
        N[f'{layer}_prev'] = (R > 0).mean(1).round(3); N[f'{layer}_share_ppm'] = (R.div(lib['size'], axis=1).mean(1) * 1e6).round(2)
    keep = (forced(N) & (N[f'{layer}_prev'] >= PREV_FORCED)) | ((N[f'{layer}_prev'] >= PREV) & (N[f'{layer}_share_ppm'] >= MIN_PPM))
    return N[keep]


def load_layer(layer):
    """X: libraries x nodes (counts; DNA in kb of aligned bases, rounded), lib: library table, N: node table of the layer."""
    N = layer_nodes(layer)
    R, lib = raw_table(layer)
    X = (R / LAYER_CFG[layer]['unit']).loc[N.index, lib.index].T.round()
    return X, lib, N


def rel_abund(X, lib, layer):
    """Share of the library (libraries x nodes) from load_layer's pseudo-counts."""
    return X.mul(LAYER_CFG[layer]['unit']).div(lib['size'], axis=0)


def meta_vars(lib):
    """Numeric meta-variables for FlashWeave: winter (2025-26 = 1), diesel pool, control pool, seawater, sea ice, top ice, bottom
    ice (0/1 each; 'other' material is 0 on both pool and medium), experiment day and modelled under-ice PAR (light/underice.py).
    Day and PAR are missing for some non-pool libraries; they get the median of their winter. Water temperature is logged for
    only 40 of 128 RNA libraries (none of the ice samples), so it is not used."""
    L = lib.layer.fillna('').astype(str).str.upper()
    M = pd.DataFrame({'winter_2025_26': (lib.season == '2025-26').astype(int), 'diesel_pool': (lib.group == 'EXP').astype(int),
                      'control_pool': (lib.group == 'CTL').astype(int), 'seawater': (lib.medium == 'SW').astype(int),
                      'sea_ice': (lib.medium == 'SI').astype(int), 'top_ice': L.isin(['TOP', 'HIGH']).astype(int),
                      'bottom_ice': L.isin(['BOT', 'LOW']).astype(int)}, index=lib.index)
    for c, src in (('day', 'day'), ('under_ice_par', 'par_under')):
        v = pd.to_numeric(lib[src], errors='coerce')
        M[c] = v.fillna(v.groupby(lib.season).transform('median')).round(4)
    return M.loc[:, M.std() > 0]


def design(lib):
    """Design for the residual-Spearman baseline: one level per winter x pool x medium cell (full interaction) plus a linear day
    trend within each winter; missing days are set to the winter median."""
    cell = (lib.season + '|' + lib.group + '|' + lib.medium).astype(str)
    Z = pd.get_dummies(cell, dtype=float)
    d = pd.to_numeric(lib.day, errors='coerce'); d = d.fillna(d.groupby(lib.season).transform('median'))
    for s in sorted(lib.season.unique()):
        Z[f'day|{s}'] = np.where(lib.season == s, d - d[lib.season == s].mean(), 0.0)
    return Z


def edge_key(a, b):
    return (a, b) if a < b else (b, a)


def meta_columns(layer):
    return list(pd.read_csv(os.path.join(CACHE, f'{layer}_meta.tsv'), sep='\t', index_col=0).columns)


def read_flashweave(layer):
    """FlashWeave edges of a layer: (genome-genome edges with sign, meta-variable 'driver' edges)."""
    meta = meta_columns(layer)
    E = pd.read_csv(os.path.join(OUT, f'flashweave_{layer}.edgelist'), sep='\t', comment='#', header=None, names=['a', 'b', 'weight'])
    ism = E.a.isin(meta) | E.b.isin(meta)
    Dr = E[ism & ~(E.a.isin(meta) & E.b.isin(meta))].copy()
    Dr['driver'] = np.where(Dr.a.isin(meta), Dr.a, Dr.b); Dr['node'] = np.where(Dr.a.isin(meta), Dr.b, Dr.a)
    E = E[~ism].copy(); E['sign'] = np.sign(E.weight).astype(int)
    E[['a', 'b']] = [edge_key(a, b) for a, b in zip(E.a, E.b)]
    return E, Dr[['node', 'driver', 'weight']]


def read_methods(layer):
    """Every method's edges of a layer, keyed the same way: {method: DataFrame(a, b, sign, weight[, stability, q, n])}, drivers.
    One key order for every method: R (SPIEC-EASI) sorts ids case-insensitively, Python by code point, so 'contig_..' vs 'EE04'
    came out in opposite orders and the same pair was two different edges (fixed 2026-11-08)."""
    M = {}
    M['flashweave'], dr = read_flashweave(layer)
    for m in METHODS[1:]:
        M[m] = pd.read_csv(os.path.join(OUT, f'edges_{m}_{layer}.tsv'), sep='\t')
        M[m][['a', 'b']] = [edge_key(a, b) for a, b in zip(M[m].a, M[m].b)]
    return M, dr


# ---------------------------------------------------------------- the edge table the page gets (network_site.py, forcelayout.py)
# The residual Spearman is dense, so edges that ONLY it reports are kept only at |rho| >= SP_FLOOR (RNA: 69,711 Spearman edges
# at q < 0.05). The Presence Pearson needs no floor: its edge call (PP_R, PP_N) is strict.
# 'Strong' scores, per method (also used by the permissivity steps; the page recomputes them from the exported values):
#   SPRING, SPIEC-EASI: StARS stability >= STRONG_STAB; residual Spearman: |rho| >= STRONG_RHO and q < STRONG_Q;
#   Presence Pearson: |r| >= STRONG_PP_R, q < STRONG_Q and co-present in >= STRONG_PP_N libraries (stricter than its edge call on
#   all three); FlashWeave: every edge.
SP_FLOOR = 0.4
# Presence Pearson edge call (run_pearson.py): pairs co-present in >= PP_N_MIN libraries are tested (BH over all of them); an edge
# needs q < PP_Q AND |r| >= PP_R AND co-presence in >= PP_N. q alone is far too loose: with up to 128 libraries, |r| ~ 0.25 passes
# q < 0.05 and 61% of RNA pairs (40% of DNA) would be edges. The rule is the least strict one (raising |r| before lowering q) that
# keeps the edge count at or below SPRING's in both layers (DNA 5,694, RNA 10,051): DNA 4,479 edges (3.7% of pairs), RNA 1,994 (1.2%).
# |r| >= 0.8 would pass RNA (4,019) but not DNA (7,470 at n >= 20; 6,692 even at q < 1e-5). See NOTES.md for the full table.
PP_N_MIN, PP_Q, PP_R, PP_N = 10, 0.05, 0.85, 20
STRONG_STAB, STRONG_RHO, STRONG_Q = 0.9, 0.6, 1e-3
STRONG_PP_R, STRONG_PP_N = 0.9, 30
BIT = {m: 1 << i for i, m in enumerate(METHODS)}
# exported values, as the page gets them: score x100 (|beta| x1000), q as -log10 q x10
q10 = lambda q: np.where(np.isfinite(q), np.minimum(np.round(-np.log10(np.maximum(q, 1e-300)) * 10), 3000), 0).astype(int)
i100 = lambda v: np.round(np.abs(np.nan_to_num(np.asarray(v, float))) * 100).astype(int)


def strong_bits(T):
    """Bit mask of the methods with a strong score on each edge, on the exported (rounded) values, exactly as the page does it."""
    qt = -np.log10(STRONG_Q) * 10
    st = BIT['flashweave'] * (T.flashweave_sign.to_numpy() != 0)
    st = st + BIT['spring'] * (i100(T.spring_stab) >= STRONG_STAB * 100) + BIT['spieceasi'] * (i100(T.spieceasi_stab) >= STRONG_STAB * 100)
    st = st + BIT['spearman'] * ((i100(T.spearman_w) >= STRONG_RHO * 100) & (q10(T.spearman_q.to_numpy(float)) > qt))
    st = st + BIT['pearson'] * ((i100(T.pearson_w) >= STRONG_PP_R * 100) & (q10(T.pearson_q.to_numpy(float)) > qt) &
                                (np.nan_to_num(T.pearson_n.to_numpy(float)) >= STRONG_PP_N))
    return st.astype(int)


def edge_table(layer):
    """All methods' edges of a layer joined on the pair: <m>_sign, <m>_w (signed score), <m>_stab (SPRING, SPIEC-EASI), spearman_q,
    pearson_q, pearson_n; sign = majority sign (ties positive); mask / anti = bit masks of the methods with / against that sign.
    Single-method Spearman and Presence Pearson edges below their floors are dropped. Returns (table, {method: dropped})."""
    M, _ = read_methods(layer)
    T = None
    for m, E in M.items():
        cols = {'sign': f'{m}_sign', 'weight': f'{m}_w'}
        if 'stability' in E: cols['stability'] = f'{m}_stab'
        for c in ('q', 'n'):
            if c in E: cols[c] = f'{m}_{c}'
        e = E.set_index(['a', 'b'])[list(cols)].rename(columns=cols)
        T = e if T is None else T.join(e, how='outer')
    for m in METHODS: T[f'{m}_sign'] = T[f'{m}_sign'].fillna(0).astype(int)
    S = T[[f'{m}_sign' for m in METHODS]]
    pos, neg = (S > 0).sum(axis=1), (S < 0).sum(axis=1)
    T['sign'] = np.where(pos >= neg, 1, -1)
    T['mask'] = sum(np.where(S[f'{m}_sign'].to_numpy() == T.sign.to_numpy(), BIT[m], 0) for m in METHODS)
    T['anti'] = sum(np.where(S[f'{m}_sign'].to_numpy() == -T.sign.to_numpy(), BIT[m], 0) for m in METHODS)
    only = lambda m: T['mask'].eq(BIT[m]) & T.anti.eq(0)
    drop_sp = only('spearman') & (T.spearman_w.abs() < SP_FLOOR)
    drop_pp = only('pearson') & False                 # its edge call is already strict (PP_R, PP_N); no extra floor
    T = T[~(drop_sp | drop_pp)]
    return T.reset_index(), {'spearman': int(drop_sp.sum()), 'pearson': int(drop_pp.sum())}
