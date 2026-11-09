"""Illumina DNA (metagenome) abundance of the EE## MAGs and P## organelle groups, and its differential abundance. Run after
hull_mags.py / merge_mags.py (it reads out/contig_hull_mag.tsv and out/contig_organelle_group.tsv); rerun when the MAGs or the
metagenome run change.
Data: common.illumina_dna_contig_matrix() (danaSeq illumina_rna in DNA mode against the marine coassembly; lanes summed per library;
QC-flagged libraries and libraries with < 1 M read pairs on contigs left out: 162 of 169).
Metric: length-weighted mean = read pairs on the MAG's contigs per kb of MAG per million contig-assigned read pairs of the library
(reads per kb per million, the unit of the RNA profile). With fixed-length reads this is the MAG's mean depth as a share of the
library, i.e. the nanopore 'depth x length share' divided by MAG length, so the two DNA platforms use one metric.
Why not the outlier-resistant RNA profile (robustexpr.py)? DNA is not spiky the way RNA is: a contig's DNA level follows the genome's
copy number, not how busy its genes are. The 5 contigs with most reads hold a median of 9% of a MAG's Illumina DNA reads (EE04: 6%,
against 94% of its RNA on one rRNA contig). The robust profile computed on the same DNA agrees with the length-weighted mean (median
Spearman rho 0.94 over the 162 libraries), and the length-weighted mean agrees slightly better with nanopore DNA on the shared
samples. MAGs whose DNA does rest on few contigs (top-5 share >= TOP5_FLAG) are flagged.
Agreement with nanopore DNA: Spearman rho per MAG over the samples with both an Illumina and a nanopore metagenome, matched by winter
and corrected sample name (nanopore resequenced runs averaged, duplicate Illumina libraries of one sample averaged).
Differential abundance: dacore.run_mg (the da.py contrasts; ice vs water from ice and water libraries of the same winter, pool and day),
data label 'Illumina DNA'.
Out: out/hull_mags_mgdna_samples.tsv (MAG x library, reads per kb per million), out/hull_mags_mgdna_method.tsv (contigs, kb,
top-5 read share, flag, rho with the robust profile, rho with nanopore DNA and the number of shared samples, medians),
out/hull_mags_da_mgdna.tsv (long table, as hull_mags_da.tsv)."""
import sys
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from common import illumina_dna_samples, nanopore_samples
import eukshared as E, robustexpr as RX, dacore
O = E.O
TOP5_FLAG = 0.3

RET = {v for k, v in E.load_ids().items() if k.startswith('retired:')}
H = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')[['contig', 'mag']]
H = H[~H.mag.isin(RET) & H.mag.str.startswith('EE')]
G = pd.concat([H, pd.read_csv(O + 'contig_organelle_group.tsv', sep='\t').rename(columns={'oid': 'mag'})[['contig', 'mag']]]).set_index('contig').mag

r = E._mgdna(); C, tot, pos, S = r['C'], r['tot'], r['pos'], r['S']; L = RX._len()
P, rows = {}, []
for m, ids in G.groupby(G):
    ids = [c for c in ids.index if c in pos.index]
    reads = C[pos[ids].values]; rs = np.asarray(reads.sum(1)).ravel(); kb = L.reindex(ids).sum() / 1e3
    P[m] = np.asarray(reads.sum(0)).ravel() / tot * 1e6 / kb
    top = float(np.sort(rs)[::-1][:5].sum() / max(rs.sum(), 1))
    rows.append(dict(mag=m, contigs=len(ids), kb=round(kb, 1), top5_share=round(top, 3), few_contigs=top >= TOP5_FLAG))
P = pd.DataFrame(P, index=S).T.sort_index(); I = pd.DataFrame(rows).set_index('mag').sort_index()
RB, RI = RX.robust_profiles(G, 'mgdna')
I['robust_contigs_used'] = RI.n_used.reindex(I.index)
I['rho_lw_vs_robust'] = [round(spearmanr(P.loc[m], RB.loc[m]).correlation, 3) for m in I.index]

# agreement with nanopore DNA (eukshared.dna_profile: depth x length share of the CMO nanopore libraries)
lib = illumina_dna_samples().loc[S]
V, vt = E.dna_matrix()
ns = nanopore_samples().set_index('run_barcode').reindex(V.columns)
nkey = (ns.log.astype(str).map({'2024': '2024-25', '2025': '2025-26'}) + '|' + ns.sample_name).values
ND = pd.DataFrame({m: V.reindex(list(ids.index)).fillna(0).sum(axis=0) / vt for m, ids in G.groupby(G)}).T
ND = ND.T.groupby(nkey).mean().T
ikey = (lib.season + '|' + lib.sample_name).values
IL = P.T.groupby(ikey).mean().T
sh = sorted(set(IL.columns) & set(ND.columns))
I['n_shared_nanopore'] = len(sh)
I['rho_vs_nanopore'] = [round(spearmanr(IL.loc[m, sh], ND.loc[m, sh]).correlation, 3) for m in I.index]
pool = lib.index[lib.pool]
I['median_rpkm_pool'] = P[pool].median(axis=1).round(3)
I['det_pool'] = (P[pool] > 0).mean(axis=1).round(3)

P.to_csv(O + 'hull_mags_mgdna_samples.tsv', sep='\t', float_format='%.4g')
I.to_csv(O + 'hull_mags_mgdna_method.tsv', sep='\t')
pd.DataFrame({'key': sh}).to_csv(O + 'mgdna_nanopore_shared_samples.tsv', sep='\t', index=False)
DA = dacore.run_mg(P.loc[P.index.str.startswith('EE')])
DA.to_csv(O + 'hull_mags_da_mgdna.tsv', sep='\t', index=False)

pd.set_option('display.width', 220)
print(len(S), 'Illumina metagenomes,', len(pool), 'pool libraries;', len(sh), 'samples shared with nanopore')
print(I.to_string())
print('EE median rho vs nanopore:', I.loc[I.index.str.startswith('EE'), 'rho_vs_nanopore'].median(), '; lw vs robust:', I.rho_lw_vs_robust.median())
print(DA[DA.q < 0.05].sort_values(['contrast', 'q']).to_string(index=False))
print(DA.groupby('contrast')[['n_a', 'n_b']].first())
