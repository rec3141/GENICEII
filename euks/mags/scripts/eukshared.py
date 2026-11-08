"""Shared rules for the eukaryotic/viral group pipeline.
ids: one stable ID allocator for every step (final.py splits, register_blobs.py blobs). out/id_registry.tsv maps a stable key
     ('split:E03:2', 'blob:57,-47', ...) to an ID; a key keeps its ID across reruns; new keys take max(existing E or V IDs) + 1.
split_parts: a group splits into its two components of a 2-component Gaussian mixture on contig t-SNE coordinates (contigs
     replicated per 5 kb, random_state 0) when Ashman's D > 3 (clearly bimodal), the component centres are > 5 t-SNE units apart
     (~2% of the map width), and the minor part is >= 2 Mbp and >= 10% of the group. Parts are returned largest first."""
import os, re
import numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
REG = O + 'id_registry.tsv'
D_MIN, SEP_MIN, MINOR_BP, MINOR_FRAC = 3.0, 5.0, 2e6, 0.10
def load_ids():
    return pd.read_csv(REG, sep='\t').set_index('key').id.to_dict() if os.path.exists(REG) else {}
def save_ids(R):
    pd.DataFrame({'key': list(R), 'id': list(R.values())}).to_csv(REG, sep='\t', index=False)
def allocate(R, key, prefix, taken=()):
    if key in R: return R[key]
    nums = [int(v[1:]) for v in list(R.values()) + list(taken) if isinstance(v, str) and re.fullmatch(prefix + r'\d+', v)]
    R[key] = f'{prefix}{(max(nums) if nums else 0) + 1:02d}'; return R[key]
def split_parts(g):
    """g: contig rows (length, tsne_x, tsne_y) of one group. Returns [index-list, ...] (one element = no split) and stats."""
    L = g.length.values; xy = g[['tsne_x', 'tsne_y']].values
    if L.sum() < 2 * MINOR_BP or len(g) < 20: return [list(g.index)], {}
    rep = np.repeat(xy, np.maximum(1, np.round(L / 5000).astype(int)), axis=0)
    if len(rep) > 40000: rep = rep[np.random.default_rng(0).choice(len(rep), 40000, replace=False)]
    gm = GaussianMixture(2, random_state=0, n_init=3).fit(rep)
    sep = float(np.linalg.norm(gm.means_[0] - gm.means_[1])); sd = np.sqrt([np.trace(c) / 2 for c in gm.covariances_])
    D = float(np.sqrt(2) * sep / np.sqrt((sd ** 2).sum())); lab = gm.predict(xy)
    bp = np.array([L[lab == k].sum() for k in (0, 1)])
    st = dict(ashman_D=round(D, 2), sep=round(sep, 1), minor_bp=int(bp.min()), minor_frac=round(bp.min() / bp.sum(), 2))
    if D > D_MIN and sep > SEP_MIN and bp.min() >= MINOR_BP and bp.min() / bp.sum() >= MINOR_FRAC:
        order = np.argsort(-bp); return [list(g.index[lab == k]) for k in order], st
    return [list(g.index)], st
