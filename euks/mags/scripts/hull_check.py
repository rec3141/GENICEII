"""Zoomed re-check panel for one hull MAG (user request 2026-11-08: "could EE30 be expanded ~2-fold?").
Usage: hull_check.py <hull id in hulls.tsv> x0 x1 y0 y1 [--ref out/hulls/check/<hull>_ref.tsv]
Panel 1: all contigs >= 2 kb in the window sized by cumulative RNA (out/cum_ppm.pkl, as plot_hulls.py), coloured by Tiara;
         every curated polygon in the window outlined (this hull thick), 1-unit grid.
Panel 2: same points coloured by the mean of the RNA and DNA Pearson r (log1p profiles, out/kernel_*.pkl as kernel.py) with the
         REFERENCE profile = bp-weighted mean log profile of the reference contigs (default: the hull's members now in
         out/contig_hull_mag.tsv; the first run writes them to --ref so later runs after a polygon edit keep the old core as reference).
         Contigs without signal (< 10 RNA reads summed / no kernel profile) are grey.
Also prints, for the current polygon: Mbp, contigs, Tiara eukarya share, and the bp-weighted r of the contigs that are new
relative to the reference set.
Out: out/hulls/check/<hull>.png, out/hulls/check/<hull>_ref.tsv"""
import sys, os, argparse
import numpy as np, pandas as pd
from matplotlib.path import Path
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
O = '/data/dev/rnaseq-exploration/euks/mags/out/'; HF = '/data/dev/rnaseq-exploration/euks/mags/hulls.tsv'
ap = argparse.ArgumentParser(); ap.add_argument('hull'); ap.add_argument('win', type=float, nargs=4); ap.add_argument('--ref')
a = ap.parse_args(); x0, x1, y0, y1 = a.win
os.makedirs(O + 'hulls/check', exist_ok=True); REF = a.ref or O + f'hulls/check/{a.hull}_ref.tsv'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id'); T = T[T.tsne_x.notna() & (T.length >= 1000)]
T['cum_ppm'] = pd.read_pickle(O + 'cum_ppm.pkl').cum_ppm.reindex(T.index).fillna(0)
H = pd.read_csv(HF, sep='\t', comment='#'); CM = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')
if not os.path.exists(REF):
    CM[CM.hull == a.hull][['contig']].to_csv(REF, sep='\t', index=False); print('reference written:', REF)
ref = pd.read_csv(REF, sep='\t').contig
RN = np.log1p(pd.read_pickle(O + 'kernel_rna_cpm.pkl')); DN = np.log1p(pd.read_pickle(O + 'kernel_dna_share.pkl') * 1e6)
W = T[T.tsne_x.between(x0, x1) & T.tsne_y.between(y0, y1)]
def prof(M, ids):
    ids = [i for i in ids if i in M.index and M.loc[i].sum() > 0]; w = T.length.reindex(ids).values
    return (M.loc[ids].values * w[:, None]).sum(0) / w.sum()
def corr(M, p, ids):
    ids = [i for i in ids if i in M.index]; X = M.loc[ids].values; X = X - X.mean(1, keepdims=True); q = p - p.mean()
    den = np.sqrt((X ** 2).sum(1) * (q ** 2).sum()); return pd.Series(np.where(den > 0, X @ q / np.where(den > 0, den, 1), np.nan), index=ids)
pr, pd_ = prof(RN, ref), prof(DN, ref)
rr, rd = corr(RN, pr, W.index), corr(DN, pd_, W.index)
reads = pd.read_pickle(O + 'kernel_rna_cpm.pkl').reindex(W.index).sum(axis=1)
rr[reads.reindex(rr.index) < 1] = np.nan       # too few reads to judge (cpm summed < 1)
r = pd.concat([rr, rd], axis=1).mean(axis=1).reindex(W.index)
h = H[H.id == a.hull].sort_values('vertex_order')[['x', 'y']].values
inpoly = Path(h).contains_points(T[['tsne_x', 'tsne_y']].values); P = T[inpoly & ~T.tiara.isin(['bacteria', 'archaea']).values]
new = P.index[~P.index.isin(set(ref))]; old = P.index[P.index.isin(set(ref))]
def wr(ids):
    s = r.reindex(ids).dropna(); w = T.length[s.index]; return (s * w).sum() / w.sum() if len(s) else np.nan
tl = P.length[P.tiara.notna()].sum()
print(f'{a.hull} polygon now: {P.length.sum() / 1e6:.2f} Mbp, {len(P)} contigs, Tiara eukarya {P.length[P.tiara == "eukarya"].sum() / tl:.3f}; '
      f'reference {T.length.reindex(ref).sum() / 1e6:.2f} Mbp / {len(ref)} contigs; '
      f'bp-weighted mean r with reference profile: kept reference {wr(old):.2f}, new {wr(new):.2f} ({T.length[new].sum() / 1e6:.2f} Mbp, {len(new)} contigs)')
SR = W[W.length >= 2000].sort_values('cum_ppm'); q = T.cum_ppm.quantile(.999); size = lambda v: 1 + 120 * np.sqrt(v / q)
TC = {'eukarya': '#2f9e44', 'organelle': '#f08c00', 'bacteria': '#1971c2', 'archaea': '#9c36b5'}
fig, axs = plt.subplots(1, 2, figsize=(20, 10.5))
for k, ax in enumerate(axs):
    if k == 0:
        c = [TC.get(t, '#adb5bd') for t in SR.tiara]; ax.scatter(SR.tsne_x, SR.tsne_y, s=size(SR.cum_ppm), c=c, lw=0, alpha=.55)
        for t, cc in TC.items(): ax.scatter([], [], c=cc, label=t)
        ax.scatter([], [], c='#adb5bd', label='unknown/none'); ax.legend(loc='lower left', fontsize=8)
        ax.set_title(f'{a.hull} window: contigs >= 2 kb sized by cumulative RNA, coloured by Tiara', fontsize=11)
    else:
        rv = r.reindex(SR.index); g = rv.isna()
        ax.scatter(SR.tsne_x[g], SR.tsne_y[g], s=size(SR.cum_ppm[g]), c='#dee2e6', lw=0, alpha=.5)
        sc = ax.scatter(SR.tsne_x[~g], SR.tsne_y[~g], s=size(SR.cum_ppm[~g]), c=rv[~g], cmap='RdYlBu_r', vmin=-.2, vmax=1, lw=0, alpha=.75)
        plt.colorbar(sc, ax=ax, shrink=.6, label='mean RNA/DNA r with reference profile')
        ax.set_title(f'coloured by profile r with the {a.hull} reference core ({len(ref)} contigs)', fontsize=11)
    for hid, g in H.groupby('id'):
        v = g.sort_values('vertex_order')[['x', 'y']].values
        if not ((v[:, 0].max() > x0) & (v[:, 0].min() < x1) & (v[:, 1].max() > y0) & (v[:, 1].min() < y1)): continue
        ax.plot(np.r_[v[:, 0], v[0, 0]], np.r_[v[:, 1], v[0, 1]], c='k', lw=2.2 if hid == a.hull else .9)
        ax.annotate(hid, (v[:, 0].mean(), v[:, 1].max()), fontsize=9, weight='bold', ha='center', va='bottom')
    st = 1 if (x1 - x0) <= 30 else 2
    for ax_ in (ax.xaxis, ax.yaxis): ax_.set_major_locator(MultipleLocator(st))
    ax.grid(True, color='#495057', lw=.35, alpha=.5); ax.tick_params(labelsize=7); plt.setp(ax.get_xticklabels(), rotation=90)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect('equal')
plt.tight_layout(); plt.savefig(O + f'hulls/check/{a.hull}.png', dpi=90); plt.close(); print('wrote', O + f'hulls/check/{a.hull}.png')
# length-matched comparison (short contigs have noisier profiles, so a raw mean r understates halo membership)
lb = pd.cut(T.length, [1000, 5000, 10000, 30000, 1e9], labels=['1-5 kb', '5-10 kb', '10-30 kb', '>30 kb'])
cmp_ = pd.DataFrame({'reference (in polygon)': r.reindex(old).groupby(lb.reindex(old), observed=True).median(),
                     'new': r.reindex(new).groupby(lb.reindex(new), observed=True).median(),
                     'outside polygon, window': r.reindex(W.index.difference(P.index)).groupby(lb.reindex(W.index.difference(P.index)), observed=True).median()})
print('median r by contig length:'); print(cmp_.round(2).to_string())
