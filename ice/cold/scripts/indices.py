"""Proteome cold-adaptation indices.
Per MAG (pooled composition of all its proteins), per sample (transcript-weighted composition of the expressed CDS proteome),
and a within/between decomposition (fixed MAG proteomes vs expression shifts).
-> ../mag_indices.tsv, ../sample_indices.tsv, ../index_paired.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
from cats import AA
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
ix = {a: i for i, a in enumerate(AA)}

def indices(F):
    """F: (n, 20) amino-acid fractions -> DataFrame of indices"""
    F = np.asarray(F, float); c = lambda s: F[:, [ix[a] for a in s]].sum(1)
    R, K = c('R'), c('K')
    return pd.DataFrame({'arg_lys': R / (R + K), 'pro': c('P'), 'acidic': c('DE'), 'charged': c('DEKR'), 'gly': c('G'),
                         'ser': c('S'), 'gln_thr': c('QT'), 'aliphatic': 100 * (c('A') + 2.9 * c('V') + 3.9 * c('IL')),
                         'gravy': F @ np.array([KD[a] for a in AA]), 'polar_unch': c('STNQ')})

if __name__ == '__main__':
    cnt = np.load(C_DIR + 'grex/counts.npy').astype(np.float32)
    ids = open(C_DIR + 'grex/ids.txt').read().split()
    X, g = genes(); m = design(); S = samples()
    pos = pd.Series(np.arange(len(ids)), index=ids).reindex(g.gene_id.values)
    ok = pos.notna().values & (g['type'].values == 'CDS')
    print('CDS with protein:', ok.sum(), 'of', (g['type'].values == 'CDS').sum())
    A = np.zeros((len(g), 20), np.float32); A[ok] = cnt[pos.values[ok].astype(int)]
    L = A.sum(1)
    # --- per MAG
    q = pd.read_csv('/data/dev/rnaseq-exploration/data/summary.tsv', sep='\t', index_col=0)
    bins = g['bin'].fillna('').values
    sel = ok & (bins != '')
    df = pd.DataFrame(A[sel], columns=list(AA)); df['bin'] = bins[sel]
    comp = df.groupby('bin').sum(); comp = comp.div(comp.sum(1), axis=0)
    MI = indices(comp.values); MI.index = comp.index
    gc = pd.read_csv(C_DIR + 'grex/bin_gc.tsv', sep='\t', header=None, names=['bin', 'gc', 'size_bp']).set_index('bin')
    MI = MI.join(gc).join(q[['SCG_completeness', 'SCG_redundancy']])
    MI.round(5).to_csv(C_DIR + 'mag_indices.tsv', sep='\t')
    # --- per sample: transcript-weighted (reads / length = molecules) mean of per-protein fractions
    Fp = np.zeros_like(A); Fp[ok] = A[ok] / np.maximum(L[ok, None], 1)
    w = X.multiply(1.0 / np.maximum(L * 3, 1)[:, None]).tocsc()      # molecules ~ reads per nt
    w = w[np.where(ok)[0]]
    Wt = np.asarray(w.sum(0)).ravel()
    comp_s = (w.T @ Fp[ok]) / Wt[:, None]
    SI = indices(comp_s); SI.index = S
    # within-MAG fixed-composition counterfactual: each MAG's transcript weight x its pooled proteome composition
    binned = ok & (bins != '')
    codes, uniq = pd.factorize(pd.Series(bins[binned]))
    import scipy.sparse as sp
    B = sp.csr_matrix((np.ones(binned.sum()), (codes, np.where(binned)[0])), shape=(len(uniq), len(g)))
    Wm = (B @ X.multiply(1.0 / np.maximum(L * 3, 1)[:, None]).tocsr()).toarray()   # MAG x sample molecules
    Fb = comp.reindex(uniq).values
    comp_fixed = (Wm.T @ Fb) / Wm.sum(0)[:, None]
    # actual composition restricted to binned genes
    wb = X.multiply(1.0 / np.maximum(L * 3, 1)[:, None]).tocsr()[np.where(binned)[0]]
    comp_bin = (wb.T @ Fp[binned]) / np.asarray(wb.sum(0)).ravel()[:, None]
    SF = indices(comp_fixed); SF.index = S
    SB = indices(comp_bin); SB.index = S
    out = pd.concat({'all': SI, 'binned': SB, 'fixed': SF}, axis=1)
    out.columns = [f'{a}_{b}' for a, b in out.columns]
    out.round(6).to_csv(C_DIR + 'sample_indices.tsv', sep='\t')
    rows = []
    for kind in ['ice', 'top']:
        for col in out.columns:
            d = [out.loc[a, col].mean() - out.loc[b, col].mean() for _, a, b in blocks(m, kind)]
            rows.append({'contrast': kind, 'index': col, 'median_diff': np.median(d), 'frac_up': np.mean(np.array(d) > 0), 'n': len(d),
                         'sd_samples': out[col].std()})
    R = pd.DataFrame(rows); R.round(6).to_csv(C_DIR + 'index_paired.tsv', sep='\t', index=False)
    pd.set_option('display.width', 200); pd.set_option('display.max_rows', 100)
    print(R.round(4).to_string(index=False))
