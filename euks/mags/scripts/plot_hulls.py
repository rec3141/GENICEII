"""Figures for the hull MAGs (run after hull_mags.py).
out/euk_hull_mags_tsne.png: 3 panels as euk_groups_tsne.png (full map, eukaryotic-region zoom, zoom sized by cumulative RNA) with
EE/VV member contigs coloured and the curated polygons drawn as outlines labelled EE##/VV## (retired ids in out/id_registry.tsv left out).
out/hulls/verify_expr_sheet_##.png: per MAG, hull outline over all contigs sized by cumulative RNA (expression check)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id'); T = T[T.tsne_x.notna()]
T['cum_ppm'] = pd.read_pickle(O + 'cum_ppm.pkl').cum_ppm.reindex(T.index).fillna(0)
H = pd.read_csv('/data/dev/rnaseq-exploration/euks/mags/hulls.tsv', sep='\t', comment='#')
S = pd.read_csv(O + 'hull_mags.tsv', sep='\t'); CM = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts'); from eukshared import load_ids
RET = {v for k, v in load_ids().items() if k.startswith('retired:')}         # e.g. VV01, a tailed-phage cluster, not a genome
S = S[~S.mag.isin(RET)]; CM = CM[CM.mag.isin(S.mag)]
PAL = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#F0E442', '#332288', '#117733', '#882255', '#44AA99', '#999933', '#AA4499', '#88CCEE', '#DDCC77', '#CC6677']
col = {m: ('#000000' if m.startswith('VV') else PAL[k % len(PAL)]) for k, m in enumerate(S.sort_values('bp', ascending=False).mag)}
poly = {x.mag: [H[H.id == h].sort_values('vertex_order')[['x', 'y']].values for h in x.hull.split('+')] for x in S.itertuples()}     # merged MAGs: several hulls
Z = (-50, 60, -118, -25)
fig = plt.figure(figsize=(16, 17.5)); gsp = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1.05])
axes = [fig.add_subplot(gsp[0, 0]), fig.add_subplot(gsp[0, 1]), fig.add_subplot(gsp[1, :])]
bg = T[T.length >= 2000]; bgs = bg.sample(min(len(bg), 400000), random_state=0)
q = T.cum_ppm.quantile(.999); size = lambda v: 0.3 + 60 * np.sqrt(v / q)
for i, ax in enumerate(axes):
    zoom = i > 0; b = bgs if not zoom else bg[bg.tsne_x.between(Z[0], Z[1]) & bg.tsne_y.between(Z[2], Z[3])]
    if i < 2: ax.scatter(b.tsne_x, b.tsne_y, s=.3 if not zoom else .6, c='#e3e6e8', lw=0, rasterized=True)
    else:
        bb = b.sort_values('cum_ppm'); ax.scatter(bb.tsne_x, bb.tsne_y, s=size(bb.cum_ppm), c='#b8bec4', lw=0, alpha=.45, rasterized=True)
    for m, g in CM.groupby('mag'):
        s = T.loc[g.contig]
        if zoom: s = s[s.tsne_x.between(Z[0], Z[1]) & s.tsne_y.between(Z[2], Z[3])]
        if len(s) == 0: continue
        ax.scatter(s.tsne_x, s.tsne_y, s=(size(s.cum_ppm) if i == 2 else (1.2 if not zoom else 2.5)), color=col[m], lw=0, alpha=.7 if i == 2 else .8, rasterized=True)
        for v in poly[m]: ax.plot(np.r_[v[:, 0], v[0, 0]], np.r_[v[:, 1], v[0, 1]], c='k' if not m.startswith('VV') else '#c92a2a', lw=.9 if zoom else .5)
        for k, v in enumerate(poly[m]):      # a merged MAG labels each of its hulls; the total goes on the first
            lab = m if (i < 2 or k > 0) else f"{m} {s.cum_ppm.sum() / 1e3:.0f}k"
            ax.annotate(lab, (v[:, 0].max(), v[:, 1].max()), fontsize=7.5 if i == 0 else 8.5, weight='bold', ha='left', va='bottom',
                        bbox=dict(boxstyle='round,pad=0.12', fc='white' if not m.startswith('VV') else '#fff3bf', ec='none', alpha=.8))
    st = 20 if not zoom else 10
    ax.xaxis.set_major_locator(MultipleLocator(st)); ax.yaxis.set_major_locator(MultipleLocator(st)); ax.grid(True, color='#b9c0c6', lw=.5, alpha=.7)
    ax.tick_params(labelsize=8); ax.set_xlabel('t-SNE 1', fontsize=9); ax.set_ylabel('t-SNE 2', fontsize=9)
    if zoom: ax.set_xlim(Z[0], Z[1]); ax.set_ylim(Z[2], Z[3])
    else: ax.add_patch(plt.Rectangle((Z[0], Z[2]), Z[1] - Z[0], Z[3] - Z[2], fill=False, ec='#444', lw=.8, ls='--'))
axes[0].set_title('Contig t-SNE, marine-coassembly-flye', fontsize=11); axes[1].set_title('Eukaryotic region (dashed box)', fontsize=11)
axes[2].set_title('Eukaryotic region, contigs ≥ 2 kb sized by cumulative RNA (per-library reads per million summed over 128 libraries); label = MAG total', fontsize=10)
fig.suptitle('Hull MAGs: EE## eukaryotic MAGs (coloured, black outline). Polygons were drawn visually\n(euks/mags/hulls.tsv); a contig belongs to the MAG whose polygon contains it (eukaryotic hulls exclude Tiara bacteria/archaea).', fontsize=10.5)
plt.tight_layout(rect=[0, 0, 1, 0.965]); plt.savefig(O + 'euk_hull_mags_tsne.png', dpi=100); plt.close()
# expression check sheets
mags = list(S.sort_values('mag').mag)
for k in range(0, len(mags), 9):
    fig, axs = plt.subplots(3, 3, figsize=(15, 15))
    for ax, m in zip(axs.ravel(), mags[k:k + 9]):
        v = np.vstack(poly[m]); x0, x1 = v[:, 0].min() - 4, v[:, 0].max() + 4; y0, y1 = v[:, 1].min() - 4, v[:, 1].max() + 4
        b = T[T.tsne_x.between(x0, x1) & T.tsne_y.between(y0, y1) & (T.length >= 1000)].sort_values('cum_ppm')
        ins = b.index.isin(CM.contig[CM.mag == m])
        ax.scatter(b.tsne_x[~ins], b.tsne_y[~ins], s=size(b.cum_ppm[~ins]), c='#adb5bd', lw=0, alpha=.5)
        ax.scatter(b.tsne_x[ins], b.tsne_y[ins], s=size(b.cum_ppm[ins]), c=col[m] if not m.startswith('VV') else '#1c7ed6', lw=0, alpha=.6)
        for vv in poly[m]: ax.plot(np.r_[vv[:, 0], vv[0, 0]], np.r_[vv[:, 1], vv[0, 1]], c='k', lw=1)
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect('equal'); ax.grid(alpha=.3); ax.tick_params(labelsize=7)
        x = S.set_index('mag').loc[m]; ax.set_title(f'{m} ({x.hull}) {x.bp / 1e6:.1f} Mbp, median {x.expr_ppm_median:.0f}/M', fontsize=9)
    for ax in axs.ravel()[len(mags[k:k + 9]):]: ax.axis('off')
    fig.suptitle('Expression check: size = cumulative RNA; coloured = MAG members; grey = other contigs; black = hull', fontsize=11)
    plt.tight_layout(); plt.savefig(O + f'hulls/verify_expr_sheet_{k // 9 + 1:02d}.png', dpi=75); plt.close()
print('ok', len(mags))
