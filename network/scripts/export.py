"""Write the per-layer input tables for the R (SPIEC-EASI) and Julia (FlashWeave) runs: cache/{layer}_X.tsv (libraries x nodes,
counts) and cache/{layer}_meta.tsv (meta-variables, netcommon.meta_vars)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *

for layer in LAYERS:
    X, lib, N = load_layer(layer)
    X.index.name = 'library'
    X.astype(int).to_csv(os.path.join(CACHE, f'{layer}_X.tsv'), sep='\t')
    M = meta_vars(lib); M.index.name = 'library'
    M.to_csv(os.path.join(CACHE, f'{layer}_meta.tsv'), sep='\t')
    print(layer, X.shape, 'meta:', list(M.columns), '; nodes by domain:', N.domain.value_counts().to_dict())
