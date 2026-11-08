"""t-SNE compactness, fragmentation and footprint overlap of the listed eukaryotic bins and consolidated groups;
compactness-refined groups. Reference = all DAS Tool prokaryotic MAGs.
Rules (reproducible): spatial clusters = DBSCAN(eps=1.5, min_samples=20, sample_weight=contig kb) on contig t-SNE coords
(eukaryotic genomes are expected to span ~10% of the map width, ~24 units). Every DBSCAN(eps=1.5, min_samples=20 kb-weighted) cluster holding >= 5% and >= 2 Mbp of a group's bp and >= 50% Tiara-eukarya bp is its own genome
(a density gap separates genomes); smaller clusters join the nearest such cluster within 24 units AND <= 3 GC points; a cluster farther away holding >= 15% of the group's bp becomes its own group
(split); DBSCAN noise and small far clusters are trimmed. Candidate additions = unlisted Tiara-eukarya contigs >= 5 kb within
3 units of a kept contig, with GC within +-3 of the group and depth within 0.5-2x the group median. Footprint = 2x2-unit grid cells holding the top 90% of a bin's bp.
Out: out/compact_bins.tsv, out/compact_groups.tsv, out/footprint_overlap.tsv, out/groups_refined.tsv, out/contig_group_refined.tsv, out/*.png"""
import pandas as pd, numpy as np, itertools
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
# scale: eukaryotic genomes are expected to span ~10% of the map width (user), i.e. ~24 t-SNE units
WIDTH = None; EPS = 1.5; MINS = 20; SPLIT_GAP = 24.0; ADD_R = 3.0; GC_TOL = 3.0
T = pd.read_pickle(O + 'contigs_full.pkl'); T = T[T.tsne_x.notna()].set_index('id')
B = pd.read_csv(O + 'bin_stats.tsv', sep='\t'); CG = pd.read_csv(O + 'contig_group.tsv', sep='\t')
COL = {'semibin': 'semibin_bin', 'comebin': 'comebin_bin', 'metabat': 'metabat_bin', 'vamb': 'vamb_bin', 'binette': 'binette_bin', 'dastool': 'bin'}
mem = {b: T.index[T[COL[k]] == b] for b, k in zip(B.bin, B.binner)}
def geom(ids):
    g = T.loc[ids]; w = g.length.values.astype(float); x = g[['tsne_x', 'tsne_y']].values
    c = np.average(x, axis=0, weights=w); d = np.sqrt(((x - c) ** 2).sum(1)); o = np.argsort(d); cw = np.cumsum(w[o]) / w.sum()
    r80 = d[o][np.searchsorted(cw, 0.8)]
    lab = DBSCAN(eps=EPS, min_samples=MINS).fit_predict(x, sample_weight=w / 1000)
    bp = pd.Series(w).groupby(lab).sum(); main = bp.drop(-1, errors='ignore').idxmax() if (bp.index != -1).any() else -1
    return dict(rms=np.sqrt(np.average(d ** 2, weights=w)), d_med=np.median(d), d_p90=np.quantile(d, 0.9), r80=r80, cx=c[0], cy=c[1],
                n_clusters=int((bp.index != -1).sum()), main_share=bp.get(main, 0) / w.sum(), outlier_share=bp.get(-1, 0) / w.sum(),
                second_share=(bp.drop([-1, main], errors='ignore').max() / w.sum()) if len(bp.drop([-1, main], errors='ignore')) else 0.0), lab, main
# reference: DAS Tool prokaryotic MAGs
dom = pd.read_csv('/data/dev/rnaseq-exploration/data/bin_domain.tsv', sep='\t').set_index('bin')
ref = []
for b, g in T[T.bin.notna()].groupby('bin'):
    if dom.call.get(b) != 'prokaryote' or len(g) < 10: continue
    w = g.length.values; x = g[['tsne_x', 'tsne_y']].values; c = np.average(x, axis=0, weights=w)
    ref.append(np.sqrt(np.average(((x - c) ** 2).sum(1), weights=w)))
ref = np.array(ref); WIDTH = T.tsne_x.max() - T.tsne_x.min(); print('map width %.0f; 10%% = %.0f units' % (WIDTH, WIDTH / 10)); print('reference prokaryotic MAG RMS radius: median %.1f, IQR %.1f-%.1f, 90th %.1f (n=%d)' % (np.median(ref), *np.quantile(ref, [.25, .75, .9]), len(ref)))
rows = []
for _, x in B.iterrows():
    gm, lab, main = geom(mem[x.bin]); g = T.loc[mem[x.bin]]
    out = g[lab != main]
    other = []
    for b2, ids2 in mem.items():
        if b2 == x.bin: continue
        sh = out.length[out.index.isin(ids2)].sum() / max(out.length.sum(), 1)
        if sh > 0.2: other.append(f'{b2} {sh:.0%}')
    rows.append(dict(bin=x.bin, flag=x.user_flag, bp=x.bp, **{k: round(v, 2) if isinstance(v, float) else v for k, v in gm.items()},
                     rms_pctile_vs_prok=round((ref < gm['rms']).mean(), 2),
                     off_main_tiara_euk=round(out.length[out.tiara == 'eukarya'].sum() / max(out.length.sum(), 1), 2),
                     off_main_gc=round(np.average(out.gc, weights=out.length), 1) if len(out) else np.nan, off_main_depth=round(out.depth.median(), 2) if len(out) else np.nan,
                     off_main_claimed_by='; '.join(other)[:150]))
CB = pd.DataFrame(rows); CB.to_csv(O + 'compact_bins.tsv', sep='\t', index=False)
# footprints
def foot(ids):
    g = T.loc[ids]; cell = (np.floor(g.tsne_x / 2).astype(int).astype(str) + '_' + np.floor(g.tsne_y / 2).astype(int).astype(str))
    s = g.length.groupby(cell.values).sum().sort_values(ascending=False); keep = s.index[(s.cumsum() / s.sum()) <= 0.9].tolist() + [s.index[0]]
    return set(keep), cell
fp = {b: foot(ids) for b, ids in mem.items()}
ov = []
euk = list(B.bin[B.user_flag == 'euk'])
for a, b in itertools.permutations(euk, 2):
    A = T.loc[mem[a]]; inB = fp[a][1].isin(fp[b][0]).values
    f = A.length[inB].sum() / A.length.sum()
    sh = A.length[A.index.isin(mem[b])].sum() / A.length.sum()
    if f > 0.1 or sh > 0.1: ov.append(dict(a=a, b=b, footprint_frac_a_in_b=round(f, 2), shared_contig_frac_a=round(sh, 2)))
OV = pd.DataFrame(ov); OV.to_csv(O + 'footprint_overlap.tsv', sep='\t', index=False)
# refine groups by compactness
allk = set().union(*[set(v) for v in mem.values()])
tree_all = None
grows, cgr = [], []
for gid, g0 in CG.groupby('group'):
    ids = g0.contig.values; gm, lab, main = geom(ids); g = T.loc[ids]; w = g.length
    bp = w.groupby(lab).sum(); xy = g[['tsne_x', 'tsne_y']]
    cen = {c: np.average(xy[lab == c].values, axis=0, weights=w[lab == c]) for c in bp.index if c != -1}
    gcc = {c: np.average(g.gc[lab == c], weights=w[lab == c]) for c in cen}
    ok = lambda a, b: np.linalg.norm(cen[a] - cen[b]) <= SPLIT_GAP and abs(gcc[a] - gcc[b]) <= GC_TOL
    # merge clusters near a seed; seeds = clusters far (> SPLIT_GAP) from every bigger seed and holding >= 15% bp
    # every dense cluster holding >= 5% and >= 2 Mbp of the group's bp is its own genome (density gap = separate genome);
    # smaller clusters attach to the nearest such cluster within SPLIT_GAP units and GC_TOL points, else are trimmed
    teuk = {c: (w[(lab == c) & (g.tiara == 'eukarya').values].sum() / max(w[(lab == c) & g.tiara.notna().values].sum(), 1)) for c in cen}
    # a genome seed must be mostly eukaryotic by Tiara (drops prokaryotic parts of mixed bins such as semibin_1384)
    seeds = [c for c in sorted(cen, key=lambda c: -bp[c]) if bp[c] / w.sum() >= 0.05 and bp[c] >= 2e6 and teuk[c] >= 0.5] or \
            [c for c in [max(cen, key=lambda c: bp[c] * (teuk[c] >= 0.5))] if teuk[c] >= 0.5]
    assign = {}
    for c in cen:
        if c in seeds: assign[c] = c; continue
        d_ = {s_: np.linalg.norm(cen[c] - cen[s_]) for s_ in seeds if ok(c, s_)}
        if d_: assign[c] = min(d_, key=d_.get)
    keepc = seeds
    for j, c in enumerate(sorted(keepc, key=lambda c: -bp[c])):
        members_ = [k for k, v in assign.items() if v == c]
        sub = g[np.isin(lab, members_)]; nid = gid if len(keepc) == 1 else f'{gid}{"abcdefgh"[j]}'
        gcg = np.average(sub.gc, weights=sub.length); dg = sub.depth.median()
        cand = T[(T.tiara == 'eukarya') & (T.length >= 5000) & ~T.index.isin(allk) & T.gc.between(gcg - 3, gcg + 3) & T.depth.between(0.5 * dg, 2 * dg)]
        tr = cKDTree(sub[['tsne_x', 'tsne_y']].values); dd, _ = tr.query(cand[['tsne_x', 'tsne_y']].values, distance_upper_bound=ADD_R)
        add = cand[np.isfinite(dd)]
        tl = sub.length[sub.tiara.notna()].sum()
        grows.append(dict(group=nid, from_group=gid, bins=g0.contig.map(lambda _: '').iloc[0] or '', contigs=len(sub), bp=int(sub.length.sum()),
                          trimmed_bp=int(w.sum() - sum(bp[k] for k in assign)), purity=round(sub.length[sub.tiara == 'eukarya'].sum() / tl, 3) if tl else np.nan,
                          rms=round(geom(sub.index)[0]['rms'], 1), r80=round(geom(sub.index)[0]['r80'], 1), gc=round(np.average(sub.gc, weights=sub.length), 1), depth_med=round(sub.depth.median(), 2),
                          cand_add_contigs=len(add), cand_add_bp=int(add.length.sum()),
                          cand_add_who_euk=round(add.length[add.whokaryote == 'eukaryote'].sum() / max(add.length[add.whokaryote.notna()].sum(), 1), 2),
                          cx=round(sub.tsne_x.mean(), 1), cy=round(sub.tsne_y.mean(), 1)))
        for cid in sub.index: cgr.append(dict(contig=cid, group=nid, source='binned'))
        for cid in add.index: cgr.append(dict(contig=cid, group=nid, source='added_tiara_neighbour'))
    # panel
    fig, ax = plt.subplots(figsize=(6.5, 6))
    bg = T[(T.length >= 3000) & T.tsne_x.between(gm['cx'] - 40, gm['cx'] + 40) & T.tsne_y.between(gm['cy'] - 40, gm['cy'] + 40)]
    ax.scatter(bg.tsne_x, bg.tsne_y, s=0.5, c='#d0d4d8', lw=0, rasterized=True)
    cols = ['#d9480f', '#1971c2', '#2b8a3e', '#7048e8', '#e8590c', '#0c8599', '#c2255c']
    bins_here = [b for b in euk if len(set(mem[b]) & set(ids))]
    for k, b in enumerate(bins_here):
        s = T.loc[mem[b]]; ax.scatter(s.tsne_x, s.tsne_y, s=np.clip(s.length / 3000, 1, 20), c=cols[k % 7], lw=0, alpha=.6, label=f'{b} ({s.length.sum() / 1e6:.1f} Mbp)')
        s = s[s.tsne_x.between(gm['cx'] - 40, gm['cx'] + 40) & s.tsne_y.between(gm['cy'] - 40, gm['cy'] + 40)]; p = s[['tsne_x', 'tsne_y']].values
        if len(p) >= 3:
            h = ConvexHull(p); ax.plot(np.r_[p[h.vertices, 0], p[h.vertices[0], 0]], np.r_[p[h.vertices, 1], p[h.vertices[0], 1]], c=cols[k % 7], lw=.8)
    ax.legend(fontsize=7, loc='lower left'); ax.set_title(f'{gid}: rms {gm["rms"]:.1f}, main cluster {gm["main_share"]:.0%}, outliers {gm["outlier_share"]:.0%}', fontsize=9)
    ax.set_xlim(gm['cx'] - 40, gm['cx'] + 40); ax.set_ylim(gm['cy'] - 40, gm['cy'] + 40); ax.set_xticks([]); ax.set_yticks([])
    plt.tight_layout(); plt.savefig(O + f'tsne_{gid}.png', dpi=100); plt.close()
GR = pd.DataFrame(grows); GR.to_csv(O + 'groups_refined.tsv', sep='\t', index=False); pd.DataFrame(cgr).to_csv(O + 'contig_group_refined.tsv', sep='\t', index=False)
# overview
fig, ax = plt.subplots(figsize=(10, 9))
bg = T[T.length >= 5000].sample(200000, random_state=0); ax.scatter(bg.tsne_x, bg.tsne_y, s=.3, c='#d9dde0', lw=0, rasterized=True)
e = T[(T.tiara == 'eukarya') & (T.length >= 5000)]; ax.scatter(e.tsne_x, e.tsne_y, s=.3, c='#9bd3a5', lw=0, rasterized=True)
R = pd.read_csv(O + 'contig_group_refined.tsv', sep='\t'); cm = plt.get_cmap('tab20')
for k, (gid, g) in enumerate(R[R.source == 'binned'].groupby('group')):
    s = T.loc[g.contig]; ax.scatter(s.tsne_x, s.tsne_y, s=1.5, color=cm(k % 20), lw=0); ax.text(s.tsne_x.median(), s.tsne_y.median(), gid, fontsize=7, weight='bold')
ax.set_title('Eukaryotic genome groups (coloured) on the contig t-SNE; light green = all Tiara-eukarya contigs >= 5 kb'); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig(O + 'tsne_overview.png', dpi=110)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
print(CB[['bin', 'flag', 'bp', 'rms', 'd_med', 'd_p90', 'r80', 'rms_pctile_vs_prok', 'n_clusters', 'main_share', 'second_share', 'outlier_share', 'off_main_tiara_euk', 'off_main_gc', 'off_main_claimed_by']].to_string(index=False))
print(GR.drop(columns=['bins']).to_string(index=False))
print(OV[OV.footprint_frac_a_in_b >= 0.3].sort_values('footprint_frac_a_in_b', ascending=False).head(40).to_string(index=False))
