"""Figure of eukaryotic (E##) and viral (V##) groups from the final tables (run after register_blobs.py).
Reads out/euk_groups.tsv, out/contig_euk_group.tsv, out/virus_groups.tsv, out/contig_virus_group.tsv, out/cum_ppm.pkl.
Panels: full map, zoom on the eukaryotic region, zoom sized by cumulative RNA. Viruses = black outline + label.
Out: out/euk_groups_tsne.png"""
import pandas as pd, numpy as np, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import samples
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id'); T['cum_ppm'] = pd.read_pickle(O + 'cum_ppm.pkl').cum_ppm.reindex(T.index).fillna(0)
S = pd.read_csv(O + 'euk_groups.tsv', sep='\t').fillna({'taxonomy_hint': ''}); CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t')
VG = pd.read_csv(O + 'virus_groups.tsv', sep='\t'); CV = pd.read_csv(O + 'contig_virus_group.tsv', sep='\t')
if 'status' not in S.columns: S['status'] = np.where(S.network_ok, 'network', 'doubtful')
def fmt_ppm(v):
    return f'{v / 1e3:.0f}k' if v >= 1e3 else f'{v:.0f}'
# figure
PAL = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#F0E442', '#332288', '#117733', '#882255', '#44AA99', '#999933', '#AA4499', '#88CCEE', '#DDCC77', '#CC6677']
ing = set(CE.contig)
ue = T[(T.tiara == 'eukarya') & (T.length >= 5000) & ~T.index.isin(ing)]
bg = T[T.length >= 2000]; bg = bg.sample(min(len(bg), 400000), random_state=0)
fig = plt.figure(figsize=(16, 17.5))
gsp = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1, 1.05])
axes = [fig.add_subplot(gsp[0, 0]), fig.add_subplot(gsp[0, 1])]; ax3 = fig.add_subplot(gsp[1, :])
Z = (-50, 60, -118, -25)
for ax, zoom in zip(axes, (False, True)):
    b = bg if not zoom else bg[bg.tsne_x.between(Z[0], Z[1]) & bg.tsne_y.between(Z[2], Z[3])]
    ax.scatter(b.tsne_x, b.tsne_y, s=.25 if not zoom else .5, c='#e3e6e8', lw=0, rasterized=True)
    u = ue if not zoom else ue[ue.tsne_x.between(Z[0], Z[1]) & ue.tsne_y.between(Z[2], Z[3])]
    ax.scatter(u.tsne_x, u.tsne_y, s=.35 if not zoom else .9, c='#7d8388', lw=0, alpha=.55, rasterized=True, label='uncaptured Tiara-eukarya contigs ≥ 5 kb')
    k = 0
    for _, x in S.sort_values('bp', ascending=False).iterrows():
        s = T.loc[CE.contig[CE.euk_group == x.group]]
        if zoom: s = s[s.tsne_x.between(Z[0], Z[1]) & s.tsne_y.between(Z[2], Z[3])]
        if len(s) == 0: continue
        lab = f"{x.group}{' ' + x.taxonomy_hint if x.taxonomy_hint else ''}"
        if x.status == 'network':
            col = PAL[k % len(PAL)]; k += 1
            ax.scatter(s.tsne_x, s.tsne_y, s=1.5 if not zoom else 3, color=col, lw=0, rasterized=True)
        elif x.status == 'doubtful':
            ax.scatter(s.tsne_x, s.tsne_y, s=3 if not zoom else 6, facecolors='none', edgecolors='#9aa0a5', linewidths=.25, alpha=.6, rasterized=True)
            lab += ' (doubtful)'
        else:   # pure, compact, but low RNA or < 4 Mbp
            col = PAL[k % len(PAL)]; k += 1
            ax.scatter(s.tsne_x, s.tsne_y, s=1.2 if not zoom else 2.5, color=col, lw=0, alpha=.35, rasterized=True)
            lab += f' ({x.status})'
        cx, cy = np.median(s.tsne_x), np.median(s.tsne_y)
        if not zoom and x.status != 'network': continue
        if zoom and x.status != 'network' and x.bp < 3e6: continue      # avoid label clutter from small non-network fragments
        ax.annotate(lab, (cx, cy), fontsize=7.5 if not zoom else 8.5, weight='bold', ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=.75))
    if zoom: ax.set_xlim(Z[0], Z[1]); ax.set_ylim(Z[2], Z[3])
    else: ax.add_patch(plt.Rectangle((Z[0], Z[2]), Z[1] - Z[0], Z[3] - Z[2], fill=False, ec='#444', lw=.8, ls='--'))
    from matplotlib.ticker import MultipleLocator
    st = 10 if zoom else 20
    ax.xaxis.set_major_locator(MultipleLocator(st)); ax.yaxis.set_major_locator(MultipleLocator(st))
    ax.grid(True, which='major', color='#b9c0c6', lw=.5, alpha=.7); ax.set_axisbelow(False)
    ax.tick_params(labelsize=8); ax.set_xlabel('t-SNE 1', fontsize=9); ax.set_ylabel('t-SNE 2', fontsize=9)
    ax.set_title('Contig t-SNE, marine-coassembly-flye' if not zoom else 'Eukaryotic region (dashed box)', fontsize=11)
axes[0].legend(loc='lower left', fontsize=8, markerscale=8, frameon=False)
# panel 3: eukaryotic region, contigs sized by cumulative expression
from matplotlib.ticker import MultipleLocator
W = T[T.tsne_x.between(Z[0], Z[1]) & T.tsne_y.between(Z[2], Z[3]) & (T.length >= 2000)]
size = lambda v: 0.3 + 60 * np.sqrt(v / W.cum_ppm.quantile(0.999))
grp_of = CE.set_index('contig').euk_group
gcol, k = {}, 0
for _, x in S.sort_values('bp', ascending=False).iterrows():
    if x.status != 'doubtful': gcol[x.group] = PAL[k % len(PAL)]; k += 1
w0 = W[~W.index.isin(grp_of.index)].sort_values('cum_ppm')
ax3.scatter(w0.tsne_x, w0.tsne_y, s=size(w0.cum_ppm), c='#aab0b5', lw=0, alpha=.45, rasterized=True)
w1 = W[W.index.isin(grp_of.index)].copy(); w1['g'] = grp_of.reindex(w1.index).values; w1 = w1.sort_values('cum_ppm')
ax3.scatter(w1.tsne_x, w1.tsne_y, s=size(w1.cum_ppm), c=[gcol.get(g, '#6b7176') for g in w1.g], lw=0, alpha=.75, rasterized=True)
for _, x in S.iterrows():
    if x.status == 'doubtful' or x.bp < 3e6: continue
    s_ = T.loc[CE.contig[CE.euk_group == x.group]]; s_ = s_[s_.tsne_x.between(Z[0], Z[1]) & s_.tsne_y.between(Z[2], Z[3])]
    if len(s_): ax3.annotate(f"{x.group} {fmt_ppm(s_.cum_ppm.sum())}", (np.median(s_.tsne_x), np.median(s_.tsne_y)), fontsize=8.5, weight='bold', ha='center', va='center', bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=.75))
for v in (100, 1000, 10000, 100000):
    ax3.scatter([], [], s=size(v), c='#6b7176', label=f'{v:,} cumulative per-million')
ax3.legend(loc='lower right', fontsize=8, frameon=True, title='contig expression, summed over libraries', title_fontsize=8, labelspacing=1.2, borderpad=1)
ax3.set_xlim(Z[0], Z[1]); ax3.set_ylim(Z[2], Z[3])
ax3.xaxis.set_major_locator(MultipleLocator(10)); ax3.yaxis.set_major_locator(MultipleLocator(10)); ax3.grid(True, color='#b9c0c6', lw=.5, alpha=.7)
ax3.tick_params(labelsize=8); ax3.set_xlabel('t-SNE 1', fontsize=9); ax3.set_ylabel('t-SNE 2', fontsize=9)
ax3.set_title(f'Eukaryotic region, contigs ≥ 2 kb sized by cumulative RNA expression ({len(samples())} libraries; per-library reads per million, summed). Colour = group; grey = other contigs. Label = group total.', fontsize=10)
fig.suptitle('Eukaryotic genome groups. Solid = network-ready; faded = pure but little CMO RNA (assembly is ~50% non-CMO samples) or < 4 Mbp;\ngrey outline = doubtful (impure or diffuse); black outline = virus (V##); grey dots = uncaptured Tiara-eukarya contigs ≥ 5 kb', fontsize=10.5)
for ax_, zoom in ((axes[0], False), (axes[1], True), (ax3, True)):
    for _, v in VG.iterrows():
        s_ = T.loc[CV.contig[CV.vid == v.vid]]
        if zoom: s_ = s_[s_.tsne_x.between(Z[0], Z[1]) & s_.tsne_y.between(Z[2], Z[3])]
        if len(s_) == 0: continue
        ax_.scatter(s_.tsne_x, s_.tsne_y, s=4 if zoom else 2, facecolors='none', edgecolors='k', linewidths=.3, alpha=.6, rasterized=True)
        ax_.annotate(f'{v.vid} (virus)', (np.median(s_.tsne_x), np.median(s_.tsne_y)), fontsize=8, weight='bold', ha='center', va='center', color='k',
                     bbox=dict(boxstyle='round,pad=0.15', fc='#fff3bf', ec='none', alpha=.85))
plt.tight_layout(rect=[0, 0, 1, 0.965]); plt.savefig(O + 'euk_groups_tsne.png', dpi=100)
