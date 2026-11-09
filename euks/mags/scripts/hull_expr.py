"""Per-library RNA expression of the hull MAGs (EE##) and organelle groups (P##); run after hull_mags.py / merge_mags.py.
Main table (since 2026-11-08): the outlier-resistant profile of robustexpr.py: reads per kb per million contig-assigned reads,
per-library median over the mean-scaled profiles of the MAG's contigs detected in >= 25% of the analysed RNA libraries
(common.samples(), QC-excluded libraries removed), times the MAG's median contig level. MAGs with < 5 such contigs use the
length-weighted mean (flagged in the method table). A few rRNA / conserved-gene / neighbour contigs therefore cannot set the profile.
Reference table: the old summed metric (read pairs on all the MAG's contigs per million contig-assigned reads; its median equals
hull_mags.tsv:expr_ppm_median, which the status and caution rules still use).
Out: out/hull_mags_expr_samples.tsv (robust; rows EE##/P##, columns = libraries), out/hull_mags_expr_samples_summed.tsv (old),
out/hull_mags_expr_method.tsv (contigs, kb, contigs used, method, share of reads on the 5 most-read contigs, rho robust vs summed)."""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
import pandas as pd, numpy as np
import eukshared as E, robustexpr as RX
O = E.O
M = pd.concat([pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')[['contig', 'mag']],
               pd.read_csv(O + 'contig_organelle_group.tsv', sep='\t').rename(columns={'oid': 'mag'})[['contig', 'mag']]])
g = M.set_index('contig').mag
P, I = RX.robust_profiles(g); P = P.sort_index(); I = I.sort_index()
Sm = RX.summed_profiles(g).reindex(P.index)
P.to_csv(O + 'hull_mags_expr_samples.tsv', sep='\t', float_format='%.4g')
Sm.to_csv(O + 'hull_mags_expr_samples_summed.tsv', sep='\t', float_format='%.4g')
I['rho_robust_vs_summed'] = [round(E.spearman(P.loc[m], Sm.loc[m]), 3) for m in I.index]
I['median_robust_rpkm'] = P.median(axis=1).round(2); I['median_summed_ppm'] = Sm.median(axis=1).round(1)
I.to_csv(O + 'hull_mags_expr_method.tsv', sep='\t')
H = pd.read_csv(O + 'hull_mags.tsv', sep='\t').set_index('mag')
pd.set_option('display.width', 200); print(I.to_string())
print('max |summed median - hull_mags.tsv expr_ppm_median|:', float((I.median_summed_ppm - H.expr_ppm_median.reindex(I.index)).abs().max()))
