"""Per-library RNA expression of the hull MAGs (EE##/VV##; run after hull_mags.py).
Same metric as hull_mags.py (expr_ppm_median) and expr.py (E## groups): read pairs on the MAG's contigs
(data/marine_coassembly.contig_counts.mtx.gz) per million contig-assigned reads in each of the analysed RNA libraries
(common.samples(), QC-excluded libraries removed). The P## organelle groups (out/contig_organelle_group.tsv) are included too.
Out: out/hull_mags_expr_samples.tsv (rows EE##/VV##/P##, columns = libraries)"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
import pandas as pd, numpy as np, scipy.io, scipy.sparse as sp
from common import samples, all_samples, D
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
S = samples(); ALL = all_samples(); kc = [ALL.index(x) for x in S]
C = scipy.io.mmread(D + 'marine_coassembly.contig_counts.mtx.gz').tocsr()[:, kc]
tot = np.asarray(C.sum(0)).ravel()
allids = np.sort(pd.read_pickle(O + 'contigs_full.pkl').id.values); pos = pd.Series(np.arange(len(allids)), index=allids)
assert len(allids) == C.shape[0]
M = pd.concat([pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')[['contig', 'mag']],
               pd.read_csv(O + 'contig_organelle_group.tsv', sep='\t').rename(columns={'oid': 'mag'})[['contig', 'mag']]])
codes, names = pd.factorize(M.mag)
A = sp.csr_matrix((np.ones(len(M)), (codes, pos[M.contig].values)), shape=(len(names), C.shape[0]))
E = pd.DataFrame((A @ C).toarray() / tot * 1e6, index=names, columns=S).sort_index()
E.to_csv(O + 'hull_mags_expr_samples.tsv', sep='\t', float_format='%.4g')
H = pd.read_csv(O + 'hull_mags.tsv', sep='\t').set_index('mag')
chk = pd.DataFrame({'median': E.median(axis=1).round(1), 'hull_mags.tsv': H.expr_ppm_median.reindex(E.index)})
print(chk.to_string()); print('max |diff| vs hull_mags.tsv:', float((chk['median'] - chk['hull_mags.tsv']).abs().max()))
