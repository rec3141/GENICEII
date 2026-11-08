"""Wiring tables + figure for the eukaryotic genome groups (run after final.py and expr.py).
out/contig_euk_group.tsv : contig, euk_group, source_bins (listed bins containing it), tiara, gc, length
out/euk_groups.tsv       : group summary with compactness, purity, taxonomy hint, expression, network_ok
   network_ok = Tiara-eukarya purity >= 0.9 AND >= 4 Mbp AND median contig-read share >= 10 per million (core libraries)
   AND t-SNE RMS radius <= 12 (single compact genome; clean eukaryotic genomes here are 3-6, diffuse multi-genome groups ~20).
   taxonomy_hint = lowest SILVA 18S rank shared by >= 2/3 of the group's 18S hits (>= 2 hits), else '' (no confident hint).
out/virus_groups.tsv, out/contig_virus_group.tsv: groups that are >= 50% geNomad-viral get V## IDs (out/virus_ids.tsv) and are
   kept out of the eukaryotic tables.
out/euk_groups_tsne.png  : all contigs light grey; groups coloured (colour-blind-safe Okabe-Ito + Tol muted); doubtful groups
   (network_ok False) as hollow grey outlines; unbinned Tiara-eukarya contigs >= 5 kb as small dark-grey dots; zoom panel."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
G = pd.read_csv(O + 'euk_groups_final.tsv', sep='\t').fillna('')
CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t')
X = pd.read_csv(O + 'group_expr_summary.tsv', sep='\t', index_col=0)
B = pd.read_csv(O + 'bin_stats.tsv', sep='\t'); COL = {'semibin': 'semibin_bin', 'comebin': 'comebin_bin', 'metabat': 'metabat_bin', 'vamb': 'vamb_bin', 'binette': 'binette_bin', 'dastool': 'bin'}
listed = B[B.user_flag == 'euk']
# source bins per contig
src = {}
for _, b in listed.iterrows():
    for c in T.index[T[COL[b.binner]] == b.bin]: src.setdefault(c, []).append(b.bin)
CE['source_bins'] = CE.contig.map(lambda c: ','.join(src.get(c, [])))
CE['tiara'] = T.tiara.reindex(CE.contig).fillna('').values; CE['gc'] = T.gc.reindex(CE.contig).values; CE['length'] = T.length.reindex(CE.contig).values
# taxonomy hint from 18S hits
RANKS = ['rrna_phylum', 'rrna_class', 'rrna_order', 'rrna_family', 'rrna_genus']
def hint(ids):
    h = T.loc[ids]; h = h[h.rrna_domain == 'Eukaryota'][RANKS].fillna('')
    if len(h) < 2: return ''
    best = ''
    for r in RANKS:
        v = h[r][h[r] != ''].value_counts()
        if v.empty or v.iloc[0] / len(h) < 2 / 3: break
        if not v.index[0].startswith(('uncultured', 'Incertae', 'Novel')): best = v.index[0]
    return best
rows = []
for _, x in G.iterrows():
    ids = CE.contig[CE.euk_group == x.group]
    e = X.loc[x.group] if x.group in X.index else None
    expr = float(e.contig_ppm_all_med) if e is not None else 0.0
    ok = (x.tiara_euk >= 0.9) and (x.bp >= 4e6) and (expr >= 10) and (x.rms <= 12)
    rows.append(dict(group=x.group, note=x.note, parts=x.parts, contigs=x.contigs, bp=x.bp, gc=x.gc, n50=x.n50, depth_med=x.depth_med,
                     tiara_euk=x.tiara_euk, tiara_bac=x.tiara_bac, who_euk=x.who_euk, tsne_rms=x.rms, cx=x.cx, cy=x.cy,
                     taxonomy_hint=hint(ids), rrna18S=x.rrna18S,
                     expr_ppm_median=round(expr, 1), ctl_2025_SW=e.get('ctl_2025-26_SW') if e is not None else np.nan, oil_2025_SW=e.get('oil_2025-26_SW') if e is not None else np.nan,
                     baci_2025=round(e.get('baci_2025-26'), 2) if e is not None else np.nan, ice_vs_water_log2=round(e.get('ice_vs_water_log2'), 2) if e is not None else np.nan,
                     status='network' if ok else ('doubtful' if (x.tiara_euk < 0.9 or x.rms > 12) else ('little CMO RNA' if expr < 10 else '<4 Mbp')),
                     network_ok=ok))
S = pd.DataFrame(rows)
# viral groups: >= 50% of bp geNomad 'Viruses' -> V## (stable via out/virus_ids.tsv), kept aside from the eukaryotic tables
S['virus_bp'] = [T.length.reindex(CE.contig[CE.euk_group == g]).where(T.genomad_domain.reindex(CE.contig[CE.euk_group == g]).eq('Viruses').values, 0).sum() /
                 T.length.reindex(CE.contig[CE.euk_group == g]).sum() for g in S.group]
import os
VID = pd.read_csv(O + 'virus_ids.tsv', sep='\t').set_index('from_group').vid.to_dict() if os.path.exists(O + 'virus_ids.tsv') else {}
for g in S.group[S.virus_bp >= 0.5]:
    if g not in VID: VID[g] = f'V{len(VID) + 1:02d}'
pd.DataFrame({'from_group': list(VID), 'vid': list(VID.values())}).to_csv(O + 'virus_ids.tsv', sep='\t', index=False)
V = S[S.group.isin(VID)].copy(); V.insert(0, 'vid', V.group.map(VID)); V['note'] = 'was ' + V.group; V['status'] = 'virus'; V['network_ok'] = False
V.drop(columns=['group']).to_csv(O + 'virus_groups.tsv', sep='\t', index=False)
CV = CE[CE.euk_group.isin(VID)].copy(); CV.insert(1, 'vid', CV.euk_group.map(VID)); CV.drop(columns=['euk_group']).to_csv(O + 'contig_virus_group.tsv', sep='\t', index=False)
S = S[~S.group.isin(VID)]; CE = CE[~CE.euk_group.isin(VID)]
CE[['contig', 'euk_group', 'source_bins', 'tiara', 'gc', 'length']].to_csv(O + 'contig_euk_group.tsv', sep='\t', index=False)
S.to_csv(O + 'euk_groups.tsv', sep='\t', index=False)
# cumulative expression per contig: sum over analysed RNA libraries of contig read share (per million)
import sys, scipy.io; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import samples, all_samples, D as DD
Cm = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsc()
ALL = all_samples(); keep = [ALL.index(x) for x in samples()]; Cm = Cm[:, keep]
ppm = Cm.multiply(1e6 / np.asarray(Cm.sum(0)).ravel()).tocsr()
allids = np.sort(T.index.values); CUM = pd.Series(np.asarray(ppm.sum(1)).ravel(), index=allids)
T['cum_ppm'] = CUM.reindex(T.index).fillna(0).values
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
fig.suptitle('Eukaryotic genome groups. Solid = network-ready; faded = pure but little CMO RNA (assembly is ~50% non-CMO samples) or < 4 Mbp;\ngrey outline = doubtful (impure or diffuse); grey dots = uncaptured Tiara-eukarya contigs ≥ 5 kb', fontsize=10.5)
plt.tight_layout(); plt.savefig(O + 'euk_groups_tsne.png', dpi=100)
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 40)
print(S[['group', 'note', 'parts', 'bp', 'gc', 'tiara_euk', 'who_euk', 'tsne_rms', 'taxonomy_hint', 'expr_ppm_median', 'baci_2025', 'ice_vs_water_log2', 'status', 'network_ok']].to_string(index=False))
