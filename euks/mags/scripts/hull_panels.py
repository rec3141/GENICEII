"""Zoomed t-SNE panels for visual hull drawing (one per eukaryotic group >= 2 Mbp and viral group >= 1 Mbp).
All contigs >= 1 kb light grey; the group's contigs coloured (size ~ sqrt(length)); integer grid every 1 unit (2 or 5 for wide
windows) with labelled ticks; window = the group's 0.5-99.5% coordinate range plus a 3-unit margin.
Out: out/hulls/panels/<id>.png and 2x2 contact sheets out/hulls/panels/sheet_##.png"""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
O = '/data/dev/rnaseq-exploration/euks/mags/out/'; P = O + 'hulls/panels/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id'); T = T[T.tsne_x.notna()]
E = pd.read_csv(O + 'euk_groups.tsv', sep='\t'); CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t')
V = pd.read_csv(O + 'virus_groups.tsv', sep='\t'); CV = pd.read_csv(O + 'contig_virus_group.tsv', sep='\t')
groups = [(g, CE.contig[CE.euk_group == g].tolist(), '#d9480f') for g, b in zip(E.group, E.bp) if b >= 2e6] + \
         [(v, CV.contig[CV.vid == v].tolist(), '#1c7ed6') for v, b in zip(V.vid, V.bp) if b >= 1e6]
skipped = [g for g, b in zip(E.group, E.bp) if b < 2e6] + [v for v, b in zip(V.vid, V.bp) if b < 1e6]
print('panels:', len(groups), '; skipped (< 2 Mbp euk / < 1 Mbp virus):', skipped)
bg = T[T.length >= 1000]
def draw(ax, gid, ids, col):
    s = T.loc[[i for i in ids if i in T.index]]
    x0, x1 = np.quantile(s.tsne_x, [.005, .995]); y0, y1 = np.quantile(s.tsne_y, [.005, .995])
    x0, x1, y0, y1 = np.floor(x0) - 3, np.ceil(x1) + 3, np.floor(y0) - 3, np.ceil(y1) + 3
    w = max(x1 - x0, y1 - y0); cxm, cym = (x0 + x1) / 2, (y0 + y1) / 2; x0, x1, y0, y1 = cxm - w / 2, cxm + w / 2, cym - w / 2, cym + w / 2
    b = bg[bg.tsne_x.between(x0, x1) & bg.tsne_y.between(y0, y1)]
    ax.scatter(b.tsne_x, b.tsne_y, s=.6, c='#ced4da', lw=0, rasterized=True)
    ax.scatter(s.tsne_x, s.tsne_y, s=np.clip(np.sqrt(s.length) / 25, 1, 12), c=col, lw=0, alpha=.55, rasterized=True)
    st = 1 if w <= 16 else (2 if w <= 40 else 5)
    for a in (ax.xaxis, ax.yaxis): a.set_major_locator(MultipleLocator(st))
    ax.grid(True, color='#495057', lw=.35, alpha=.5); ax.tick_params(labelsize=7)
    plt.setp(ax.get_xticklabels(), rotation=90)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect('equal')
    ax.set_title(f'{gid}: {s.length.sum() / 1e6:.1f} Mbp, {len(s)} contigs (grid {st})', fontsize=10)
for gid, ids, col in groups:
    fig, ax = plt.subplots(figsize=(7, 7)); draw(ax, gid, ids, col); plt.tight_layout(); plt.savefig(P + f'{gid}.png', dpi=110); plt.close()
for k in range(0, len(groups), 4):
    fig, axs = plt.subplots(2, 2, figsize=(13, 13))
    for ax, (gid, ids, col) in zip(axs.ravel(), groups[k:k + 4]): draw(ax, gid, ids, col)
    for ax in axs.ravel()[len(groups[k:k + 4]):]: ax.axis('off')
    plt.tight_layout(); plt.savefig(P + f'sheet_{k // 4 + 1:02d}.png', dpi=85); plt.close()
print([g for g, _, _ in groups])
