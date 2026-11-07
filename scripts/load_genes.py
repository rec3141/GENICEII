"""Gene × sample counts (CSR) + gene map -> cached npz/parquet for fast reuse."""
import gzip, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp
D = '/data/dev/rnaseq-exploration/data/'
X = scipy.io.mmread(gzip.open(D + 'marine_coassembly.gene_counts.mtx.gz')).tocsr()
sp.save_npz(D + 'gene_counts.npz', X)
g = pd.read_csv(D + 'genemap.tsv.gz', sep='\t', dtype={'KO': str, 'gene': str, 'product': str}, keep_default_na=False)
c2b = pd.read_csv(D + 'contig2bin.tsv', sep='\t', header=None, names=['contig', 'bin']).drop_duplicates('contig')
g = g.merge(c2b, on='contig', how='left')
g['bin'] = g['bin'].fillna('')
g.to_pickle(D + 'genemap.pkl')
print(X.shape, X.nnz, len(g), (g['bin'] != '').mean())
