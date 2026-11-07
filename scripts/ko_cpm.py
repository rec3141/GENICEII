"""KO × sample CPM of CDS read pairs from gene counts with the coordinate-remapped KOs (replaces the pipeline ko_tpm)."""
import sys; sys.path.insert(0, '.')
from common import *
import scipy.sparse as sp
X, g = genes(); S = samples()
cds = (g['type'] == 'CDS').values
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
k = g['KO'].fillna('')
idx = np.where(cds & k.ne('').values)[0]
rows, cols = [], []
kid = {}
for i, ks in zip(idx, k.values[idx]):
    for x in ks.split(','):
        rows.append(kid.setdefault(x, len(kid))); cols.append(i)
A = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(kid), X.shape[0]))
K = (A @ X).toarray() / tot * 1e6
pd.DataFrame(K, index=list(kid), columns=S).round(3).to_csv(D + 'ko_cpm.tsv', sep='\t', index_label='ko')
print(len(kid))
