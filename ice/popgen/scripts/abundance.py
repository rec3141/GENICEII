"""Per-bin read counts / coverage / % of reads from the competitive BAMs -> abund_reads.tsv, abund_cov.tsv, abund_pct.tsv (bins x samples)."""
import pandas as pd, subprocess, io, os, sys
D = '/data/dev/rnaseq-exploration/ice/popgen'
S = pd.read_csv(f'{D}/samples_map.tsv', sep='\t', index_col=0)
cb = pd.read_csv(f'{D}/ref/contig_bin_len.tsv', sep='\t').set_index('contig')
blen = cb.groupby('bin').len.sum()
R = {}
for s in S.index:
    b = f'{D}/bam/{s}.bam'
    if not os.path.exists(b + '.bai'): continue
    x = pd.read_csv(io.StringIO(subprocess.run(['samtools', 'idxstats', b], capture_output=True, text=True).stdout), sep='\t', header=None, names=['contig', 'len', 'm', 'u'])
    x = x[x.contig != '*'].set_index('contig')
    R[s] = x.m.groupby(cb.bin.reindex(x.index).values).sum()
R = pd.DataFrame(R).fillna(0)
R.to_csv(f'{D}/abund_reads.tsv', sep='\t')
(R.mul(150).div(blen.reindex(R.index), axis=0)).round(4).to_csv(f'{D}/abund_cov.tsv', sep='\t')
(R / 4e6 * 100).round(5).to_csv(f'{D}/abund_pct.tsv', sep='\t')   # % of all 4M input reads (2M pairs)
print(R.shape, (R.sum() / 4e6).describe().round(3).to_dict())
