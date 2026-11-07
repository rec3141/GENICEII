"""Cross-platform identity check for early 2024-25 CMO samples: per-library signature of pool-discriminating genomes."""
import sys; sys.path.insert(0, 'scripts')
import pandas as pd, numpy as np
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
SIG = ['dastool-lorbin_032', 'dastool-comebin_2596', 'dastool-semibin_088']
# nanopore per run
D = pd.read_csv('dna/bin_depths.tsv', sep='\t', index_col=0)
L = D.pop('length'); D.columns = [c.replace('.sorted.bam', '') for c in D.columns]
rel = D.mul(L, axis=0); rel = 100 * rel / rel.sum()
ns = pd.read_csv('dna/nanopore_samples.tsv', sep='\t')
print(ns.columns.tolist())
import re
from common import meta
def lab(n):
    m_ = re.search(r'CMO_(CTL|EXP)_(PRE\d+|DAY\d+)_(SW|SI\w*)', n or '')
    return m_ and m_.group(0)
EARLY = lambda s: bool(re.search(r'CMO_(CTL|EXP)_(PRE0\d|DAY0\d|DAY1[0-2])_', s or ''))
rows = []
nsx = ns.set_index('run_barcode')
for rb in rel.columns:
    if rb in nsx.index and nsx.loc[rb, 'log'] == 2024 and EARLY(nsx.loc[rb, 'sample_name']):
        rows.append(('nanopore', nsx.loc[rb, 'sample_name'], rb, *[rel.loc[b, rb] for b in SIG]))
# illumina (k-mer %; different scale)
I = pd.read_csv('out/mg_kmer.tsv', sep='\t', index_col=0)
for s, x in I[I.season == '2024-25'].iterrows():
    if EARLY(x['name']): rows.append(('illumina', x['name'], s, *[x[b.replace('dastool-', '')] for b in SIG]))
# RNA bin TPM
T = pd.read_csv('data/marine_coassembly.bin_tpm.tsv', sep='\t', index_col=0); M = meta()
for lib in T.columns:
    if lib in M.index and M.loc[lib, 'season'] == '2024-25' and EARLY(M.loc[lib, 'sample_name']):
        rows.append(('rna', M.loc[lib, 'sample_name'], lib, *[T.loc[b, lib] / 1e4 for b in SIG]))
R = pd.DataFrame(rows, columns=['platform', 'sample', 'lib', 'sulfito', 'casp', 'sulfurim'])
R['day'] = R['sample'].str.extract(r'_(PRE|DAY)(\d+)').apply(lambda r: -int(r[1]) if r[0] == 'PRE' else int(r[1]), axis=1)
R = R.sort_values(['day', 'sample', 'platform'])
print(R.round(3).to_string(index=False))
R.to_csv('swapcheck/signature.tsv', sep='\t', index=False)
