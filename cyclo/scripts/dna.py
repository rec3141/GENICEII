"""Cycloclasticus DNA abundance in nanopore metagenomes (CMO samples) from per-contig depths; plasmid copy number."""
import sys, re, json
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
C = '/data/dev/rnaseq-exploration/cyclo/'
d = pd.read_csv(C + 'cyclo_depths.txt', sep='\t')
cols = [c for c in d.columns if c.endswith('.sorted.bam')]
cg = pd.read_csv(C + 'cyclo_contigs.tsv', sep='\t', header=None, names=['contig', 'bin']).set_index('contig')['bin']
d['bin'] = d.contigName.map(cg)
ns = nanopore_samples()
key = {c: c.replace('.sorted.bam', '') for c in cols}
ns = ns[ns.run_barcode.isin(key.values())]
print('CMO nanopore samples in depths:', len(ns))
def wmean(x):
    w = x.contigLen.values[:, None]; return pd.Series((x[cols].values * w).sum(0) / w.sum(), index=cols)
B = d.groupby('bin').apply(wmean).rename(columns=key)
B.T.to_csv(C + 'cyclo_dna_depth_all.tsv', sep='\t')
cm = B[ns.run_barcode].T; cm.index = ns.sample_name.values
cm.insert(0, 'run_barcode', ns.run_barcode.values)
# total depth per sample to normalise (sum over all contigs unknown here) -> report raw depth
pd.set_option('display.width', 250)
print(cm.round(2).sort_values('dastool-lorbin_184', ascending=False).head(40).to_string())
# plasmid copy number: plasmid contig depth / chromosome median contig depth of the same bin
pl = {'contig_0030642': 'dastool-lorbin_184', 'contig_0003999': 'dastool-lorbin_184', 'contig_0011601': 'dastool-lorbin_184', 'contig_0047870': 'dastool-lorbin_184',
      'contig_0030704': 'dastool-comebin_3787_sub', 'contig_0380606': 'dastool-metabat_5432'}
rows = []
for c, b in pl.items():
    chrom = d[(d.bin == b) & (~d.contigName.isin(pl)) & (d.contigLen > 20000)]
    med = chrom[cols].median()
    pc = d.set_index('contigName').loc[c, cols]
    for s in ns.itertuples():
        col = s.run_barcode + '.sorted.bam'
        if med[col] >= 1: rows.append((c, b, s.sample_name, round(pc[col], 2), round(med[col], 2), round(pc[col] / med[col], 2)))
P = pd.DataFrame(rows, columns=['contig', 'bin', 'sample', 'contig_depth', 'chrom_median_depth', 'copy_ratio'])
P.to_csv(C + 'plasmid_copy.tsv', sep='\t', index=False)
print(P.groupby(['contig', 'bin']).copy_ratio.describe().round(2).to_string())
print(P[P.contig == 'contig_0030642'].to_string())
