"""Match each RNA library to nanopore DNA runs by whole-community MAG profile (Spearman on log abundance)."""
import sys; sys.path.insert(0, 'scripts')
import pandas as pd, numpy as np
from common import meta
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 400)
D = pd.read_csv('dna/bin_depths.tsv', sep='\t', index_col=0); L = D.pop('length')
D.columns = [c.replace('.sorted.bam', '') for c in D.columns]
rel = D.mul(L, axis=0); rel = rel / rel.sum()
ns = pd.read_csv('dna/nanopore_samples.tsv', sep='\t').set_index('run_barcode')
runs = [r for r in rel.columns if r in ns.index and str(ns.loc[r, 'sample_name']).startswith('CMO_') and rel[r].sum() > 0]
T = pd.read_csv('data/marine_coassembly.bin_tpm.tsv', sep='\t', index_col=0); M = meta()
libs = [l for l in T.columns if l in M.index and M.loc[l, 'group'] in ('CTL', 'EXP')]
bins = rel.index.intersection(T.index)
A = np.log10(rel.loc[bins, runs] + 1e-5).rank()
B = np.log10(T.loc[bins, libs] / 1e6 + 1e-5).rank()
C = pd.DataFrame(np.corrcoef(B.T.values, A.T.values)[:len(libs), len(libs):], index=libs, columns=runs)
out = []
for l in libs:
    nm = M.loc[l, 'sample_name']; c = C.loc[l].sort_values(ascending=False)
    same = [r for r in runs if ns.loc[r, 'sample_name'] == nm]
    out.append(dict(lib=l, name=nm, season=M.loc[l, 'season'], best=ns.loc[c.index[0], 'sample_name'], best_r=c.iloc[0],
                    same_r=C.loc[l, same].max() if same else np.nan, same_rank=(1 + int((c > C.loc[l, same].max()).sum())) if same else np.nan))
O = pd.DataFrame(out).sort_values(['season', 'name'])
O.to_csv('swapcheck/rna_vs_nanopore.tsv', sep='\t', index=False)
print(O.round(3).to_string(index=False))
