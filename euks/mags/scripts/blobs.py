"""Characterise highly expressed blobs on the contig t-SNE at user-given seed points (config SEEDS).
Blob rule: DBSCAN(min_samples=20, kb-weighted) on contigs >= 2 kb in a 40x40 window around the seed; the blob is the cluster holding
the most bp within 3 units of the seed; eps starts at 1.5 and shrinks x0.75 until the blob is <= 80 Mbp (floor 0.4), so dense
multi-genome regions are cut at their density gaps. Reported per blob: size, GC, depth, cumulative expression (per-library
reads per million summed over analysed libraries), Tiara/Whokaryote fractions, DAS Tool MAG membership (prokaryote calls via
data/bin_domain.tsv), other binners' bins, eukaryotic-group overlap, rRNA hits, geNomad calls.
Call: 'prokaryote MAG' if >= 50% bp in prokaryotic DAS Tool bins; 'eukaryote, in E##' / 'eukaryote, uncaptured' if Tiara-eukarya
>= 50%; 'virus' if geNomad virus >= 50%; else 'unbinned/mixed'.
Out: out/blobs.tsv, out/blob_contigs.tsv"""
import pandas as pd, numpy as np, sys
from sklearn.cluster import DBSCAN
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
# ---- config: seed points (t-SNE 1, t-SNE 2) given by the user 2026-11-08 ----
MAX_BP = 80e6
SEEDS = [(55, -95), (45, -100), (57, -47), (27, -30), (19, -30), (9, -48), (14, -58), (-3, -53), (-2, -38), (-12, -60), (-42, -32), (-48, -28)]
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
import scipy.io
from common import samples, all_samples, D as DD
Cm = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsc()
ALL = all_samples(); keep = [ALL.index(x) for x in samples()]; Cm = Cm[:, keep]
ppm = Cm.multiply(1e6 / np.asarray(Cm.sum(0)).ravel()).tocsr()
T['cum_ppm'] = pd.Series(np.asarray(ppm.sum(1)).ravel(), index=np.sort(T.index.values)).reindex(T.index).fillna(0).values
dom = pd.read_csv(DD + 'bin_domain.tsv', sep='\t').set_index('bin').call
CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t').set_index('contig').euk_group
rows, cont = [], []
for (sx, sy) in SEEDS:
    W = T[T.tsne_x.between(sx - 20, sx + 20) & T.tsne_y.between(sy - 20, sy + 20) & (T.length >= 2000)]
    eps = 1.5
    while True:      # tighten the density scale until the blob at the seed is genome-sized (<= MAX_BP) or eps hits the floor
        lab = DBSCAN(eps=eps, min_samples=20).fit_predict(W[['tsne_x', 'tsne_y']].values, sample_weight=W.length.values / 1000)
        Wl = W.assign(lab=lab); near = Wl[((Wl.tsne_x - sx) ** 2 + (Wl.tsne_y - sy) ** 2 <= 9) & (Wl.lab >= 0)]
        if near.empty: break
        c = near.groupby('lab').length.sum().idxmax(); g = Wl[Wl.lab == c]; L = g.length.sum()
        if L <= MAX_BP or eps <= 0.4: break
        eps *= 0.75
    if near.empty: rows.append(dict(seed=f'{sx},{sy}', call='no dense blob at seed')); continue
    frac = lambda m: g.length[m].sum() / L
    tl = g.length[g.tiara.notna()].sum(); te = g.length[g.tiara == 'eukarya'].sum() / tl if tl else np.nan
    wl = g.length[g.whokaryote.notna()].sum(); we = g.length[g.whokaryote == 'eukaryote'].sum() / wl if wl else np.nan
    db = g.bin.fillna(''); prok = frac(db.map(lambda b: dom.get(b, '') == 'prokaryote').values)
    dbs = g.length[db != ''].groupby(db[db != '']).sum().sort_values(ascending=False)
    eg = g.index.map(lambda i: CE.get(i, '')); egs = g.length[eg != ''].groupby(eg[eg != '']).sum().sort_values(ascending=False)
    vir = frac((g.genomad_domain == 'Viruses').values)
    sb = g.length.groupby(g.semibin_bin.fillna('-')).sum().sort_values(ascending=False)
    r = g[g.rrna_domain.notna()][['rrna_domain', 'rrna_phylum', 'rrna_class', 'rrna_genus']].fillna('').apply(lambda z: ' > '.join(v for v in z if v), axis=1).value_counts()
    kj = g.length.groupby(g.kaiju_phylum.fillna('NA')).sum().sort_values(ascending=False)
    if prok >= 0.5: call = 'prokaryote MAG: ' + ', '.join(f'{b.replace("dastool-", "")} {v / L:.0%}' for b, v in dbs.head(2).items())
    elif egs.sum() / L >= 0.5: call = 'eukaryote, in ' + ', '.join(f'{e} {v / L:.0%}' for e, v in egs.head(2).items())
    elif te >= 0.5: call = 'eukaryote, uncaptured'
    elif vir >= 0.5: call = 'virus'
    else: call = 'unbinned/mixed prokaryotic' if (te < 0.5) else 'mixed'
    rows.append(dict(seed=f'{sx},{sy}', call=call, eps=round(eps, 2), contigs=len(g), bp=int(L), gc=round(np.average(g.gc, weights=g.length), 1), depth_med=round(g.depth.median(), 2),
                     cum_ppm=round(g.cum_ppm.sum()), tiara_euk=round(te, 2), who_euk=round(we, 2), prok_mag_bp=round(prok, 2), euk_group_bp=round(egs.sum() / L, 2), virus_bp=round(vir, 2),
                     cx=round(np.average(g.tsne_x, weights=g.length), 1), cy=round(np.average(g.tsne_y, weights=g.length), 1),
                     dastool=', '.join(f'{b.replace("dastool-", "")} {v / L:.0%}' for b, v in dbs.head(3).items()),
                     semibin=', '.join(f'{b} {v / 1e6:.1f}Mb' for b, v in sb.head(3).items()),
                     kaiju_phylum=', '.join(f'{k} {v / L:.0%}' for k, v in kj.head(3).items()),
                     rrna='; '.join(f'{t} ({n})' for t, n in r.head(3).items())))
    cont += [dict(contig=i, seed=f'{sx},{sy}') for i in g.index]
B = pd.DataFrame(rows); B.to_csv(O + 'blobs.tsv', sep='\t', index=False); pd.DataFrame(cont).to_csv(O + 'blob_contigs.tsv', sep='\t', index=False)
pd.set_option('display.width', 300); pd.set_option('display.max_colwidth', 70)
print(B[['seed', 'call', 'eps', 'bp', 'gc', 'depth_med', 'cum_ppm', 'tiara_euk', 'who_euk', 'prok_mag_bp', 'euk_group_bp', 'virus_bp']].to_string(index=False))
print(B[['seed', 'dastool', 'semibin', 'kaiju_phylum', 'rrna']].to_string(index=False))
