"""Survey of user-listed eukaryotic bins across binners on the marine-coassembly-flye contig table, plus a
reproducible consolidation into eukaryotic genome groups.
Rule: candidate bins = listed bins; contig-level union; groups = connected components of the bin-overlap graph
with containment(bp) >= 0.5 (shared bp / bp of the smaller bin), so sub-bins nest into the larger bin they came from. Purity = Tiara eukarya share of length among Tiara-called bp.
Out: out/bin_stats.tsv, out/overlap_jaccard.tsv, out/groups.tsv, out/contig_group.tsv, out/missed_euk.tsv"""
import pandas as pd, numpy as np, itertools
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl')
# ---- candidate bins (config) ----------------------------------------------------------------------------
# user-listed eukaryotic bins (2026-11-08) plus bins found by the uncaptured-blob search (semibin_1338/1315/1317/1384);
# 'semibin_bin_not' = bins the user ruled out (reported, never grouped). Large mixed bins (e.g. semibin_1384) are trimmed
# to their dense Tiara-eukaryotic part by the density-gap rule in compact.py.
LIST = {'semibin_bin': [1280, 1279, 1301, 1338, 1315, 1317, 1384], 'semibin_bin_not': [1397, 1323, 630, 1325, 1398, 1389, 1390],
        'comebin_bin': [1991, 751, 4693, 3836, 3802, 1009, 3707], 'metabat_bin': [4578, 4788, 6251, 10781, 11622],
        'vamb_bin': [2078, 3134, 1186, 3866], 'binette_bin': [830, 629, 978, 758, 825, 996, 1014, 952, 591, 842, 818],
        'bin': [3937, 4527]}
PFX = {'semibin_bin': 'semibin_', 'semibin_bin_not': 'semibin_', 'comebin_bin': 'comebin_', 'metabat_bin': 'metabat_', 'vamb_bin': 'vamb_', 'binette_bin': 'binette_', 'bin': 'dastool-comebin_'}
def n50(l):
    l = np.sort(l)[::-1]; c = np.cumsum(l); return int(l[np.searchsorted(c, c[-1] / 2)])
rows, members, missing = [], {}, []
for col, nums in LIST.items():
    c = 'semibin_bin' if col == 'semibin_bin_not' else col
    for n in nums:
        name = f"{PFX[col]}{n:03d}" if col in ('vamb_bin', 'binette_bin', 'semibin_bin', 'semibin_bin_not', 'comebin_bin', 'metabat_bin') and n < 100 else f"{PFX[col]}{n}"
        g = T[T[c] == name]
        if g.empty: missing.append(name); continue
        L = g.length.sum(); tc = g.tiara.notna()
        frac = lambda s, v: g.length[s == v].sum() / L
        r = dict(bin=name, binner=c.replace('_bin', '') if c != 'bin' else 'dastool', user_flag='NOT' if col == 'semibin_bin_not' else 'euk',
                 contigs=len(g), bp=int(L), n50=n50(g.length.values), gc=round(np.average(g.gc, weights=g.length), 1), depth_med=round(g.depth.median(), 2),
                 tiara_euk=frac(g.tiara, 'eukarya'), tiara_org=frac(g.tiara, 'organelle'), tiara_bac=frac(g.tiara, 'bacteria'), tiara_arc=frac(g.tiara, 'archaea'),
                 tiara_unk=1 - frac(g.tiara, 'eukarya') - frac(g.tiara, 'organelle') - frac(g.tiara, 'bacteria') - frac(g.tiara, 'archaea'),
                 who_euk=frac(g.whokaryote, 'eukaryote'), who_prok=frac(g.whokaryote, 'prokaryote'),
                 rrna18s='; '.join(sorted(set(' > '.join(x for x in z if isinstance(x, str) and x) for z in g[g.rrna_domain == 'Eukaryota'][['rrna_phylum', 'rrna_class', 'rrna_order', 'rrna_genus']].values)))[:300],
                 rrna16s_bac=int((g.rrna_domain == 'Bacteria').sum()), plastid_contigs=int((g.tiara == 'organelle').sum()),
                 tsne_x=round(np.average(g.tsne_x, weights=g.length), 1), tsne_y=round(np.average(g.tsne_y, weights=g.length), 1),
                 tsne_spread=round(np.sqrt(np.average((g.tsne_x - np.average(g.tsne_x, weights=g.length)) ** 2 + (g.tsne_y - np.average(g.tsne_y, weights=g.length)) ** 2, weights=g.length)), 1))
        rows.append(r); members[name] = set(g.id)
B = pd.DataFrame(rows); B.to_csv(O + 'bin_stats.tsv', sep='\t', index=False)
LEN = T.set_index('id').length
# overlap among user-euk bins (exclude the NOT list from grouping, but report their overlap)
euk = [b for b in B.bin if B.set_index('bin').loc[b, 'user_flag'] == 'euk']
allb = list(B.bin)
J = pd.DataFrame(0.0, index=allb, columns=allb)
for a, b in itertools.combinations(allb, 2):
    i = members[a] & members[b]
    if i:
        u = members[a] | members[b]; J.loc[a, b] = J.loc[b, a] = LEN[list(i)].sum() / LEN[list(u)].sum()
J.to_csv(O + 'overlap_jaccard.tsv', sep='\t')
C = pd.DataFrame(0.0, index=allb, columns=allb)
for a, b in itertools.combinations(allb, 2):
    i = members[a] & members[b]
    if i: C.loc[a, b] = C.loc[b, a] = LEN[list(i)].sum() / min(LEN[list(members[a])].sum(), LEN[list(members[b])].sum())
C.to_csv(O + 'overlap_containment.tsv', sep='\t')
# connected components among euk bins at Jaccard >= 0.2
par = {b: b for b in euk}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a, b in itertools.combinations(euk, 2):
    if C.loc[a, b] >= 0.5: par[f(a)] = f(b)
comp = {}
for b in euk: comp.setdefault(f(b), []).append(b)
grows, cg = [], []
for k, bs in enumerate(sorted(comp.values(), key=lambda bs: -sum(B.set_index('bin').loc[b, 'bp'] for b in bs)), 1):
    cnt = {}
    for b in bs:
        for c in members[b]: cnt[c] = cnt.get(c, 0) + 1
    ids = list(cnt); g = T.set_index('id').loc[ids]
    core = [c for c, n in cnt.items() if n >= 2]
    gid = f'EUK{k:02d}'
    for c, n in cnt.items(): cg.append(dict(contig=c, group=gid, n_bins=n, core=n >= 2))
    L = g.length.sum(); tl = g.length[g.tiara.notna()].sum()
    grows.append(dict(group=gid, bins=','.join(bs), n_bins=len(bs), contigs=len(ids), union_bp=int(L), core_bp=int(LEN[core].sum()),
                      gc=round(np.average(g.gc, weights=g.length), 1), depth_med=round(g.depth.median(), 2),
                      tiara_euk_purity=round(g.length[g.tiara == 'eukarya'].sum() / tl, 3) if tl else np.nan,
                      tiara_org=round(g.length[g.tiara == 'organelle'].sum() / tl, 3) if tl else np.nan,
                      who_euk=round(g.length[g.whokaryote == 'eukaryote'].sum() / max(g.length[g.whokaryote.notna()].sum(), 1), 3),
                      rrna18s='; '.join(sorted(set(' > '.join(x for x in z if isinstance(x, str) and x) for z in g[g.rrna_domain == 'Eukaryota'][['rrna_phylum', 'rrna_class', 'rrna_order', 'rrna_family', 'rrna_genus']].values)))[:400],
                      tsne_x=round(np.average(g.tsne_x, weights=g.length), 1), tsne_y=round(np.average(g.tsne_y, weights=g.length), 1)))
G = pd.DataFrame(grows); G.to_csv(O + 'groups.tsv', sep='\t', index=False)
CG = pd.DataFrame(cg); CG.to_csv(O + 'contig_group.tsv', sep='\t', index=False)
# missed eukaryotic content: Tiara eukarya contigs >= 5 kb not in any listed euk bin
inany = set(CG.contig)
E = T[(T.tiara == 'eukarya') & (T.length >= 5000)]
miss = E[~E.id.isin(inany)]
miss[['id', 'length', 'gc', 'depth', 'bin', 'semibin_bin', 'comebin_bin', 'metabat_bin', 'vamb_bin', 'binette_bin', 'whokaryote', 'tsne_x', 'tsne_y']].to_csv(O + 'missed_euk.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80)
print('missing bins:', missing)
print(B.drop(columns=['rrna18s']).round(2).to_string(index=False))
print(G.drop(columns=['rrna18s']).to_string(index=False))
for _, x in G.iterrows(): print(x.group, '18S:', x.rrna18s[:300])
print('Tiara-eukarya contigs >=5kb: total %.1f Mbp; in groups %.1f Mbp; missed %.1f Mbp (%d contigs)' % (E.length.sum() / 1e6, E[E.id.isin(inany)].length.sum() / 1e6, miss.length.sum() / 1e6, len(miss)))
print('missed: binned by dastool %.1f Mbp; by any binner %.1f Mbp' % (miss[miss.bin.notna()].length.sum() / 1e6, miss[miss[['semibin_bin', 'comebin_bin', 'metabat_bin', 'vamb_bin', 'binette_bin']].notna().any(axis=1)].length.sum() / 1e6))
print('missed top semibin bins:', (miss.groupby('semibin_bin').length.sum() / 1e6).sort_values(ascending=False).head(8).round(2).to_dict())
