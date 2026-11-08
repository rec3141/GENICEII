"""Final eukaryotic genome groups and uncaptured eukaryotic blobs.
1) Merge refined groups when >= 50% of the smaller one's bp lies in the larger one's footprint (2x2 cells holding its top 90% bp)
   and GC differs by <= 2 points.
2) Uncaptured blobs: Tiara-eukarya contigs >= 3 kb outside all groups, DBSCAN(eps=1.5, min_samples=20 kb-weighted); blobs >= 2 Mbp.
3) Taxonomy hints: SILVA 18S/16S hits on the contigs, plastid bins overlapping (data/plastid_bins.tsv), sendsketch hits.
Split rule (eukshared.split_parts): a final group splits into its two 2-component GMM parts when clearly bimodal on the t-SNE
(Ashman D > 3, centres > 5 units apart, minor part >= 2 Mbp and >= 10%); the larger part keeps the ID, others get IDs from the
shared registry (out/id_registry.tsv).
Out: out/euk_groups_final.tsv, out/contig_euk_group.tsv (contig -> group, source), out/uncaptured_blobs.tsv, out/tsne_final.png"""
import pandas as pd, numpy as np, itertools
from sklearn.cluster import DBSCAN
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
R = pd.read_csv(O + 'contig_group_refined.tsv', sep='\t'); R = R[R.source == 'binned']
def fp(ids):
    g = T.loc[ids]; cell = np.floor(g.tsne_x / 2).astype(int).astype(str) + '_' + np.floor(g.tsne_y / 2).astype(int).astype(str)
    s = g.length.groupby(cell.values).sum().sort_values(ascending=False)
    return set(s.index[(s.cumsum() / s.sum()) <= 0.9]) | {s.index[0]}, cell
grp = {g: list(d.contig) for g, d in R.groupby('group')}
gc = {g: np.average(T.loc[i].gc, weights=T.loc[i].length) for g, i in grp.items()}
bp = {g: T.loc[i].length.sum() for g, i in grp.items()}
par = {g: g for g in grp}
def f(x):
    while par[x] != x: x = par[x]
    return x
F = {g: fp(i) for g, i in grp.items()}
for a, b in itertools.combinations(sorted(grp, key=lambda g: -bp[g]), 2):
    small, big = (a, b) if bp[a] < bp[b] else (b, a)
    inb = T.loc[grp[small]].length[F[small][1].isin(F[big][0]).values].sum() / bp[small]
    if inb >= 0.5 and abs(gc[a] - gc[b]) <= 2: par[f(small)] = f(big)
final = {}
for g in grp: final.setdefault(f(g), []).append(g)
rows, cmap = [], []
pl = pd.read_csv('/data/dev/rnaseq-exploration/data/plastid_bins.tsv', sep='\t')
# stable IDs: a new group inherits an earlier E## (out/v1/contig_euk_group.tsv) when it holds >= 50% of that group's bp;
# if it inherits several, it keeps the largest and notes the merge; groups inheriting nothing get the next free number
import os
V1 = pd.read_csv(O + 'v1/contig_euk_group.tsv', sep='\t') if os.path.exists(O + 'v1/contig_euk_group.tsv') else pd.DataFrame(columns=['contig', 'euk_group'])
v1bp = V1.assign(L=T.length.reindex(V1.contig).values).groupby('euk_group').L.sum()
order = sorted(final.items(), key=lambda kv: -sum(bp[g] for g in kv[1]))
names, notes, used = {}, {}, set()
for root, gs in order:
    ids = set(sum((grp[g] for g in gs), [])); hit = V1[V1.contig.isin(ids)].assign(L=lambda d: T.length.reindex(d.contig).values).groupby('euk_group').L.sum()
    newL = T.length[list(ids)].sum()
    inh = [e for e in hit.index if hit[e] / v1bp[e] >= 0.5 and e not in used]
    absorbed = [e for e in inh if v1bp[e] < newL / 3]      # an old group that is < 1/3 of the new one does not lend its ID
    inh = [e for e in inh if e not in absorbed]
    if inh:
        keep = max(inh, key=lambda e: v1bp[e]); names[root] = keep; used.update(inh)
        newbp = newL; notes[root] = ('merged ' + '+'.join(sorted(inh)) + '; ' if len(inh) > 1 else '') + (f'grew {v1bp[keep] / 1e6:.1f} -> {newbp / 1e6:.1f} Mbp' if newbp > 1.1 * v1bp[keep] else ('shrank' if newbp < 0.9 * v1bp[keep] else 'unchanged')) + (('; absorbed ' + '+'.join(absorbed)) if absorbed else '')
    elif absorbed: notes[root] = 'new; absorbed ' + '+'.join(absorbed); used.update(absorbed)
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
from eukshared import load_ids, save_ids, allocate
IDS = load_ids()
for root, gs in order:     # brand-new groups take IDs from the shared registry (key = the refined group they came from)
    if root not in names: names[root] = allocate(IDS, f'group:{root}', 'E', taken=list(v1bp.index)); notes.setdefault(root, 'new')
claimed = set()      # a contig claimed by two groups goes to the larger one (order = largest first)
for k, (root, gs) in enumerate(order, 1):
    ids = sorted(set(sum((grp[g] for g in gs), [])) - claimed); claimed.update(ids)
    from eukshared import split_parts, load_ids, save_ids, allocate
    parts, sst = split_parts(T.loc[ids])
    for j, ids in enumerate(parts):
        g = T.loc[ids]; L = g.length.sum(); tl = g.length[g.tiara.notna()].sum()
        name = names[root] if j == 0 else allocate(IDS, f'split:{names[root]}:{j + 1}', 'E', taken=list(names.values()) + list(v1bp.index))
        note_ = notes[root] if j == 0 else f'split from {names[root]} (Ashman D {sst.get("ashman_D")}, {sst.get("sep")} units apart)'
        if j == 0 and len(parts) > 1: note_ = f'{notes[root]}; split into {len(parts)} (Ashman D {sst.get("ashman_D")})'
        r18 = g[g.rrna_domain == 'Eukaryota'][['rrna_phylum', 'rrna_class', 'rrna_order', 'rrna_family', 'rrna_genus']].fillna('')
        r18 = r18.apply(lambda z: ' > '.join(x for x in z if x), axis=1).value_counts()
        r16 = g[g.rrna_domain == 'Bacteria'][['rrna_phylum', 'rrna_genus']].fillna('').apply(lambda z: ' > '.join(x for x in z if x), axis=1).value_counts()
        ss = g[g.sendsketch_domain.notna()][['sendsketch_phylum', 'sendsketch_genus']].fillna('').apply(lambda z: ' '.join(z), axis=1).value_counts()
        db = g.bin.dropna().value_counts()
        pb = [b for b in db.index if b in set(pl.bin)]
        rows.append(dict(group=name, note=note_, parts='+'.join(gs), contigs=len(g), bp=int(L), gc=round(np.average(g.gc, weights=g.length), 1), depth_med=round(g.depth.median(), 2),
                         n50=int(np.sort(g.length.values)[::-1][np.searchsorted(np.cumsum(np.sort(g.length.values)[::-1]), L / 2)]),
                         tiara_euk=round(g.length[g.tiara == 'eukarya'].sum() / tl, 3) if tl else np.nan, tiara_bac=round(g.length[g.tiara == 'bacteria'].sum() / tl, 3) if tl else np.nan,
                         tiara_org=round(g.length[g.tiara == 'organelle'].sum() / tl, 3) if tl else np.nan,
                         who_euk=round(g.length[g.whokaryote == 'eukaryote'].sum() / max(g.length[g.whokaryote.notna()].sum(), 1), 3),
                         cx=round(np.average(g.tsne_x, weights=g.length), 1), cy=round(np.average(g.tsne_y, weights=g.length), 1),
                         rms=round(np.sqrt(np.average((g.tsne_x - np.average(g.tsne_x, weights=g.length)) ** 2 + (g.tsne_y - np.average(g.tsne_y, weights=g.length)) ** 2, weights=g.length)), 1),
                         dastool_bins='; '.join(f'{b.replace("dastool-", "")} {v}' for b, v in db.head(3).items()),
                         rrna18S='; '.join(f'{t} ({n})' for t, n in r18.head(4).items()), rrna16S='; '.join(f'{t} ({n})' for t, n in r16.head(3).items()),
                         sendsketch='; '.join(f'{t} ({n})' for t, n in ss.head(3).items()), plastid_bins=','.join(pb)))
        cmap += [dict(contig=c, euk_group=name, source='binned') for c in ids]
save_ids(IDS)
G = pd.DataFrame(rows).sort_values('group'); G.to_csv(O + 'euk_groups_final.tsv', sep='\t', index=False); pd.DataFrame(cmap).to_csv(O + 'contig_euk_group.tsv', sep='\t', index=False)
# uncaptured eukaryotic blobs
ing = set(pd.DataFrame(cmap).contig)
U = T[(T.tiara == 'eukarya') & (T.length >= 3000) & ~T.index.isin(ing)]
lab = DBSCAN(eps=1.5, min_samples=20).fit_predict(U[['tsne_x', 'tsne_y']].values, sample_weight=U.length.values / 1000)
U = U.assign(blob=lab); ub = []
for b, g in U[U.blob >= 0].groupby('blob'):
    L = g.length.sum()
    if L < 2e6: continue
    sb = g.semibin_bin.fillna('none').groupby(g.semibin_bin.fillna('none')).size()
    sbbp = g.length.groupby(g.semibin_bin.fillna('none')).sum().sort_values(ascending=False)
    allc = T[(T.tsne_x - np.average(g.tsne_x, weights=g.length)) ** 2 + (T.tsne_y - np.average(g.tsne_y, weights=g.length)) ** 2 <= 25]
    r18 = g[g.rrna_domain == 'Eukaryota'][['rrna_phylum', 'rrna_class', 'rrna_genus']].fillna('').apply(lambda z: ' > '.join(x for x in z if x), axis=1).value_counts()
    ub.append(dict(blob=f'U{b}', contigs=len(g), bp=int(L), gc=round(np.average(g.gc, weights=g.length), 1), depth_med=round(g.depth.median(), 2),
                   who_euk=round(g.length[g.whokaryote == 'eukaryote'].sum() / max(g.length[g.whokaryote.notna()].sum(), 1), 2),
                   cx=round(np.average(g.tsne_x, weights=g.length), 1), cy=round(np.average(g.tsne_y, weights=g.length), 1),
                   semibin_bins='; '.join(f'{k} {v / 1e6:.1f}Mb' for k, v in sbbp.head(3).items()),
                   in_dastool_bp=int(g.length[g.bin.notna()].sum()), rrna18S='; '.join(f'{t} ({n})' for t, n in r18.head(3).items())))
UB = pd.DataFrame(ub).sort_values('bp', ascending=False); UB.to_csv(O + 'uncaptured_blobs.tsv', sep='\t', index=False)
# figure
fig, ax = plt.subplots(figsize=(11, 10))
bg = T[T.length >= 3000]; bg = bg.sample(min(len(bg), 300000), random_state=0); ax.scatter(bg.tsne_x, bg.tsne_y, s=.3, c='#d9dde0', lw=0, rasterized=True)
cm = plt.get_cmap('tab20')
for k, (gname, d) in enumerate(pd.DataFrame(cmap).groupby('euk_group')):
    s = T.loc[d.contig]; ax.scatter(s.tsne_x, s.tsne_y, s=1, color=cm(k % 20), lw=0)
    x = G.set_index('group').loc[gname]; ax.text(x.cx, x.cy, gname, fontsize=8, weight='bold', ha='center')
for _, x in UB.iterrows():
    s = U[U.blob == int(x.blob[1:])]; ax.scatter(s.tsne_x, s.tsne_y, s=1, c='k', lw=0, alpha=.5); ax.text(x.cx, x.cy, x.blob, fontsize=7, color='k', ha='center')
ax.set_title('Final eukaryotic genome groups (coloured, E##) and uncaptured Tiara-eukarya blobs >= 2 Mbp (black, U##)'); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig(O + 'tsne_final.png', dpi=105)
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 70)
print(G.drop(columns=['rrna16S', 'sendsketch']).to_string(index=False))
print(G[['group', 'rrna16S', 'sendsketch']].to_string(index=False))
print('uncaptured blobs >=2 Mbp:', len(UB), 'total %.0f Mbp' % (UB.bp.sum() / 1e6)); print(UB.to_string(index=False))

# split the dense multi-genome blob(s) (> 150 Mbp) at a finer density scale
fine = []
for _, x in UB[UB.bp > 150e6].iterrows():
    g = U[U.blob == int(x.blob[1:])]
    lab2 = DBSCAN(eps=0.8, min_samples=60).fit_predict(g[['tsne_x', 'tsne_y']].values, sample_weight=g.length.values / 1000)
    for b2, h in g.assign(b2=lab2)[lab2 >= 0].groupby('b2'):
        L = h.length.sum()
        if L < 5e6: continue
        sbbp = h.length.groupby(h.semibin_bin.fillna('none')).sum().sort_values(ascending=False)
        r18 = h[h.rrna_domain == 'Eukaryota'][['rrna_phylum', 'rrna_class', 'rrna_genus']].fillna('').apply(lambda z: ' > '.join(v for v in z if v), axis=1).value_counts()
        fine.append(dict(blob=f'{x.blob}.{b2}', contigs=len(h), bp=int(L), gc=round(np.average(h.gc, weights=h.length), 1), depth_med=round(h.depth.median(), 2),
                         who_euk=round(h.length[h.whokaryote == 'eukaryote'].sum() / max(h.length[h.whokaryote.notna()].sum(), 1), 2),
                         cx=round(np.average(h.tsne_x, weights=h.length), 1), cy=round(np.average(h.tsne_y, weights=h.length), 1),
                         rms=round(np.sqrt(np.average((h.tsne_x - np.average(h.tsne_x, weights=h.length)) ** 2 + (h.tsne_y - np.average(h.tsne_y, weights=h.length)) ** 2, weights=h.length)), 1),
                         semibin_bins='; '.join(f'{k} {v / 1e6:.1f}Mb' for k, v in sbbp.head(3).items()),
                         rrna18S='; '.join(f'{t} ({n})' for t, n in r18.head(3).items())))
FB = pd.DataFrame(fine).sort_values('bp', ascending=False); FB.to_csv(O + 'uncaptured_blobs_fine.tsv', sep='\t', index=False)
print('fine split of large blobs: %d sub-blobs >= 5 Mbp, %.0f Mbp' % (len(FB), FB.bp.sum() / 1e6)); print(FB.to_string(index=False))
