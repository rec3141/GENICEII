"""Outlier-resistant RNA profile of a contig set (eukaryotic and organelle MAGs; user-approved 2026-11-08).
A MAG's summed RNA (reads on all its contigs per million) can be dominated by a few contigs: an 18S/28S rRNA contig (EE04, EE01),
conserved tubulin/ubiquitin genes that take cross-mapped reads from many eukaryotes (EE08), or a few contigs a hull took in from a
neighbour (EE31). The robust profile follows the bulk of the genome instead (the heatmap.py method):
  1. per contig and library: reads per kb per million contig-assigned reads (RPKM; data/marine_coassembly.contig_counts.mtx.gz,
     the analysed RNA libraries, common.samples());
  2. keep contigs detected (> 0 reads) in >= MIN_DET (25%) of the libraries;
  3. scale each kept contig's profile to mean 1;
  4. MAG profile = per-library STAT of the scaled profiles x the MAG's median contig level (median over kept contigs of their
     mean RPKM). Units stay reads per kb per million.
     STAT = 10% trimmed mean (TRIM: the highest and lowest 10% of the kept contigs' values in each library are dropped, the rest
     averaged). The approved plan said 'median'; but eukaryotic contigs are sparse, and a contig detected in 25% of libraries is
     absent from most of them, so the per-library median was 0 in 45% of MAG x library cells (more than half of the MAGs had a median
     profile of 0), which leaves no time course and makes rank correlations mostly ties at 0. The trimmed mean is 0 in 9% of cells.
     Outlier resistance comes from step 3 (every contig weighs the same, whatever its read count) plus the trimming of per-library
     spikes. Set TRIM = None for the plain median. The merge decisions of 2026-11-08 hold under either (EE03+EE31 rho 0.982 median,
     0.996 trimmed; EE06+EE13 0.978 / 0.986).
MAGs with < MIN_N (5) kept contigs fall back to the length-weighted mean RPKM (total reads / total kb), and are flagged.
A profile built from < THIN_N (20) kept contigs is flagged 'thin': with so few contigs the few dominant ones can still be half of
them (EE08: 9 kept contigs, 5 of them the conserved tubulin/ubiquitin/HSP90 contigs), so it is not clearly the genome's own expression.
data='mgdna' applies the same steps to the Illumina metagenomes (eukshared._mgdna; reads per kb per million contig-assigned read pairs
of each metagenome), for the Illumina DNA abundance of the EE MAGs and the Illumina DNA network layer.
API: robust_profiles(groups, data='rna'|'mgdna') with groups = {name: contig ids} or a Series contig -> group; returns (profiles DataFrame
groups x libraries, info DataFrame: contigs, kb, n_used, method, thin, top5_share). robust_profile(ids) -> one Series.
summed_profiles(groups): the old metric (reads per million contig-assigned reads, summed over the contigs), for reference."""
import numpy as np, pandas as pd
from scipy.stats import trim_mean
import eukshared as E
MIN_DET, MIN_N, TOP_N, THIN_N = 0.25, 5, 5, 20
TRIM = 0.10
_L = {}
def _len():
    if not _L: _L['s'] = pd.read_pickle(E.O + 'contigs_full.pkl').set_index('id').length
    return _L['s']
def _groups(groups):
    if isinstance(groups, pd.Series): return {k: list(v.index) for k, v in groups.groupby(groups)}
    return {k: list(v) for k, v in groups.items()}
def robust_profiles(groups, data='rna'):
    r = E._rna() if data == 'rna' else E._mgdna(); C, tot, pos, S = r['C'], r['tot'], r['pos'], r['S']; L = _len()
    P, info = {}, []
    for g, ids in _groups(groups).items():
        ids = [c for c in dict.fromkeys(ids) if c in pos.index]
        reads = C[pos[ids].values].toarray().astype(float)                        # contigs x libraries
        kb = L.reindex(ids).values / 1e3
        rpkm = reads / tot * 1e6 / kb[:, None]
        ok = (reads > 0).mean(axis=1) >= MIN_DET; n = int(ok.sum())
        rs = reads.sum(axis=1); top = float(np.sort(rs)[::-1][:TOP_N].sum() / max(rs.sum(), 1))
        if n >= MIN_N:
            lvl = rpkm[ok].mean(axis=1); Z = rpkm[ok] / lvl[:, None]
            prof = (np.median(Z, axis=0) if TRIM is None else trim_mean(Z, TRIM, axis=0)) * np.median(lvl); how = 'robust'
        else:
            prof = reads.sum(axis=0) / tot * 1e6 / kb.sum(); how = 'length-weighted mean (fallback)'
        P[g] = prof; info.append(dict(group=g, contigs=len(ids), kb=round(kb.sum(), 1), n_used=n, method=how, thin=n < THIN_N, top5_share=round(top, 3)))
    return pd.DataFrame(P, index=S).T, pd.DataFrame(info).set_index('group')
def robust_profile(ids, data='rna'):
    P, _ = robust_profiles({'x': list(ids)}, data); return P.loc['x']
def summed_profiles(groups):
    return pd.DataFrame({g: E.rna_profile([c for c in ids if c in E._rna()['pos'].index]) for g, ids in _groups(groups).items()}).T
