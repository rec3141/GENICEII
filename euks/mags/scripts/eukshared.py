"""Shared rules for the eukaryotic/viral group pipeline.
ids: one stable ID allocator for every step (final.py splits, register_blobs.py blobs). out/id_registry.tsv maps a stable key
     ('split:E03:2', 'blob:57,-47', ...) to an ID; a key keeps its ID across reruns; new keys take max(existing E or V IDs) + 1.
split_parts: a group splits into its two components of a 2-component Gaussian mixture on contig t-SNE coordinates (contigs
     replicated per 5 kb, random_state 0) when Ashman's D > 3 (clearly bimodal), the component centres are > 5 t-SNE units apart
     (~2% of the map width), and the minor part is >= 2 Mbp and >= 10% of the group, AND (RNA guard, user rule 2026-11-08) the
     parts' RNA profiles differ: Spearman rho < RHO_SPLIT (0.9) between the parts' summed per-library RNA (reads per million
     contig-assigned reads over the analysed RNA libraries; the metric of expr.py / hull_expr.py). Two lobes of one t-SNE group that
     rise and fall together across the libraries are one genome drawn twice by t-SNE, so they are not cut (E03/E31 before this guard).
     Parts are returned largest first.
rna_profile / dna_profile: the shared per-library profiles of a contig set (RNA as above; DNA = share of each CMO nanopore
     library's depth x length, data from dna/depths.txt.gz), used by the split guard and by merge_mags.py."""
import os, re
import numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
REG = O + 'id_registry.tsv'
D_MIN, SEP_MIN, MINOR_BP, MINOR_FRAC = 3.0, 5.0, 2e6, 0.10
RHO_SPLIT = 0.9
_RNA = {}
def _rna():
    if not _RNA:
        import sys, scipy.io; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
        from common import samples, all_samples, D
        S = samples(); ALL = all_samples(); kc = [ALL.index(x) for x in S]
        C = scipy.io.mmread(D + 'marine_coassembly.contig_counts.mtx.gz').tocsr()[:, kc]
        allids = np.sort(pd.read_pickle(O + 'contigs_full.pkl').id.values)       # matrix rows = contig ids ascending
        assert len(allids) == C.shape[0]
        _RNA.update(C=C, tot=np.asarray(C.sum(0)).ravel(), pos=pd.Series(np.arange(len(allids)), index=allids), S=S)
    return _RNA
_MG = {}
def _mgdna():
    """Illumina DNA (metagenomes; common.illumina_dna_contig_matrix: lanes summed per library, QC-flagged and < 1 M pair libraries
    left out), in the same layout as _rna(): C contigs x libraries (read pairs), tot = contig-assigned pairs per library, pos, S."""
    if not _MG:
        import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
        from common import illumina_dna_contig_matrix
        C, S = illumina_dna_contig_matrix()
        allids = np.sort(pd.read_pickle(O + 'contigs_full.pkl').id.values)
        assert len(allids) == C.shape[0]
        _MG.update(C=C, tot=np.asarray(C.sum(0)).ravel(), pos=pd.Series(np.arange(len(allids)), index=allids), S=S)
    return _MG
def mgdna_profile(ids):
    """Share of each Illumina metagenome's contig-assigned read pairs on a contig set (= depth x length share with fixed-length reads)."""
    r = _mgdna(); ids = [c for c in ids if c in r['pos'].index]
    return pd.Series(np.asarray(r['C'][r['pos'][ids].values].sum(0)).ravel() / r['tot'], index=r['S'])
def rna_profile(ids):
    r = _rna(); return pd.Series(np.asarray(r['C'][r['pos'][list(ids)].values].sum(0)).ravel() / r['tot'] * 1e6, index=r['S'])
_DNA = {}
def dna_matrix():
    """contig x CMO nanopore library depth*length (DataFrame) and per-library totals; read once."""
    if not _DNA:
        import sys, gzip; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
        from common import nanopore_samples
        f = O + '../dna/depths.txt.gz'; hdr = gzip.open(f, 'rt').readline().rstrip('\n').split('\t'); ns = set(nanopore_samples().run_barcode)
        use = [c for c in hdr[3:] if not c.endswith('-var') and c.replace('.sorted.bam', '') in ns]
        X = pd.read_csv(f, sep='\t', usecols=[hdr[0], hdr[1]] + use, index_col=0)
        V = X[use].mul(X[hdr[1]], axis=0); V.columns = [c.replace('.sorted.bam', '') for c in use]
        _DNA.update(V=V, tot=V.sum(0))
    return _DNA['V'], _DNA['tot']
def dna_profile(ids):
    V, tot = dna_matrix(); return V.reindex(list(ids)).fillna(0).sum(0) / tot
def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).correlation)
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
        order = np.argsort(-bp); parts = [list(g.index[lab == k]) for k in order]
        rho = spearman(rna_profile(parts[0]), rna_profile(parts[1])); st['rna_rho'] = round(rho, 3)
        if rho < RHO_SPLIT: return parts, st
        st['kept_whole'] = f'parts co-vary in RNA (rho {rho:.3f} >= {RHO_SPLIT})'
    return [list(g.index)], st

def resolve_source(gid, valid, R=None):
    """A hull's source group that the RNA split guard retired (registry 'retired:<gid> ...' and 'split:<parent>:k' -> gid) resolves
    to the parent group; any other id is returned unchanged."""
    R = load_ids() if R is None else R
    if gid in valid or not any(k.startswith(f'retired:{gid} ') for k in R): return gid
    par = [k.split(':')[1] for k, v in R.items() if k.startswith('split:') and v == gid]
    return par[0] if par else gid
MERGE_RX = re.compile(r'^retired:(EE\d+) \(merged into (EE\d+), \d{4}-\d\d-\d\d\)$')
def merged_into(R=None):
    """{retired EE: kept EE} from merge_mags.py registry keys 'retired:EEnn (merged into EEmm, <date>)', resolved transitively."""
    R = load_ids() if R is None else R
    m = {a: b for a, b in (MERGE_RX.match(k).groups() for k in R if MERGE_RX.match(k))}
    def root(x):
        while x in m: x = m[x]
        return x
    return {a: root(a) for a in m}
