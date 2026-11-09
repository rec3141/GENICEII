"""Shared definitions for the network scripts: layer tables after the prevalence rule, meta-variables, design matrix, paths.

Prevalence rule, per layer (DNA, RNA): a node enters a layer's network if it is detected (> 0) in at least PREV of that layer's
libraries and averages at least MIN_PPM per million (RNA: contig-assigned read pairs; DNA: mapped bases). EE MAGs, giant viruses,
virophages and organelle groups are always kept where data exist (detected in at least PREV_FORCED of the libraries).
Nodes without contigs (plastid bins wholly inside an organelle group, e.g. metabat_7332_sub inside P01) are dropped."""
import os
import numpy as np, pandas as pd

NET = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
CACHE, OUT = os.path.join(NET, 'cache'), os.path.join(NET, 'out')
REPO = os.path.abspath(os.path.join(NET, '..'))
LAYERS = ('dna', 'rna')
PREV, MIN_PPM, PREV_FORCED = 0.25, 1.0, 0.10
DNA_UNIT = 1000.0                      # DNA pseudo-counts: kb of aligned nanopore read bases (coverage depth x length / 1000)
METHODS = ('flashweave', 'spring', 'spieceasi', 'spearman')
METHOD_LABEL = {'flashweave': 'FlashWeave', 'spring': 'SPRING', 'spieceasi': 'SPIEC-EASI', 'spearman': 'Residual Spearman'}


def nodes_all():
    return pd.read_csv(os.path.join(OUT, 'nodes_all.tsv'), sep='\t', index_col=0, keep_default_na=False, na_values=[''])


def forced(N):
    return (N.domain == 'eukaryote') | N.kind.isin(['giant virus', 'virophage', 'organelle group'])


def layer_nodes(layer, N=None):
    N = nodes_all() if N is None else N
    N = N[N.contigs > 0]
    keep = (forced(N) & (N[f'{layer}_prev'] >= PREV_FORCED)) | ((N[f'{layer}_prev'] >= PREV) & (N[f'{layer}_share_ppm'] >= MIN_PPM))
    return N[keep]


def load_layer(layer):
    """X: libraries x nodes (counts; DNA in kb of aligned bases, rounded), lib: library table, N: node table of the layer."""
    N = layer_nodes(layer)
    if layer == 'rna':
        X = pd.read_csv(os.path.join(CACHE, 'rna_counts.tsv'), sep='\t', index_col=0)
        lib = pd.read_csv(os.path.join(CACHE, 'rna_lib.tsv'), sep='\t', index_col=0)
    else:
        X = pd.read_csv(os.path.join(CACHE, 'dna_cov.tsv'), sep='\t', index_col=0) / DNA_UNIT
        lib = pd.read_csv(os.path.join(CACHE, 'dna_lib.tsv'), sep='\t', index_col=0)
    X = X.loc[N.index, lib.index].T.round()
    return X, lib, N


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
