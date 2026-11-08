"""Blob capture with an expression/coverage kernel (user suggestion 2026-11-08).
For each seed (config SEEDS, t-SNE units): candidates = contigs >= 2 kb within R=10 of the seed.
Profiles: RNA = per-library contig reads per million (128 analysed libraries); DNA = per-sample contig depth x length as a share
of the sample's total (CMO nanopore libraries from mapping/depths.txt via nanopore_samples()). Both log1p-transformed.
Start: core = contigs within 1.5 of the seed. Blob profile = bp-weighted mean of members' log profiles.
Membership score = exp(-d^2 / 2 sigma^2) * r, with r = mean of the RNA and DNA Pearson correlations with the blob profile (each used
only if the contig has signal: >= 10 RNA reads summed, or >= 1x mean DNA depth). Member if score >= 0.5 (sigma = 3), iterated
until membership changes by < 2% bp (max 6 rounds). Contigs with no signal in either profile are never added (cannot be judged).
Out: out/kernel_blobs.tsv, out/kernel_blob_contigs.tsv, out/kernel_blobs.png"""
import pandas as pd, numpy as np, sys, scipy.io, gzip
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import samples, all_samples, nanopore_samples, D as DD
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
# ---- config ----
SEEDS = [(55, -95), (45, -100), (57, -47), (27, -30), (19, -30), (9, -48), (14, -58), (-3, -53), (-2, -38), (-12, -60), (-42, -32), (-48, -28)]
R, SIGMA, CORE_R, THR, MIN_RNA, MIN_DNA = 10.0, 3.0, 1.5, 0.5, 10, 1.0
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
win = np.zeros(len(T), bool)
for sx, sy in SEEDS: win |= (np.hypot(T.tsne_x - sx, T.tsne_y - sy) <= R).values
W = T[win & (T.length >= 2000).values].copy()
print('candidate contigs in seed windows:', len(W), flush=True)
# RNA
Cm = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsr()
ALL = all_samples(); S = samples(); kc = [ALL.index(x) for x in S]
allids = np.sort(T.index.values); pos = pd.Series(np.arange(len(allids)), index=allids)
tot = np.asarray(Cm[:, kc].sum(0)).ravel()
rr = Cm[pos[W.index].values][:, kc].toarray()
RNA = pd.DataFrame(rr / tot * 1e6, index=W.index, columns=S); RNAreads = pd.Series(rr.sum(1), index=W.index)
# DNA (CMO nanopore libraries)
ns = nanopore_samples(); di = pd.read_csv('/data/dev/rnaseq-exploration/out/dna_samples.tsv', sep='\t')
hdr = gzip.open(O + '../dna/depths.txt.gz', 'rt').readline().rstrip('\n').split('\t')
dcols = [c for c in hdr[3:] if not c.endswith('-var')]
rb = [c.replace('.sorted.bam', '') for c in dcols]
use = [c for c, b in zip(dcols, rb) if b in set(ns.run_barcode)]
print('DNA libraries (CMO nanopore):', len(use), flush=True)
need = set(W.index); rows = []; totals = np.zeros(len(use)); idx = [hdr.index(c) for c in use]
with gzip.open(O + '../dna/depths.txt.gz', 'rt') as fh:
    next(fh)
    for line in fh:
        f = line.rstrip('\n').split('\t'); L = float(f[1]); v = np.array([float(f[i]) for i in idx])
        totals += v * L
        if f[0] in need: rows.append((f[0], v * L, v))
DNA = pd.DataFrame(np.array([r[1] for r in rows]) / totals, index=[r[0] for r in rows], columns=[c.replace('.sorted.bam', '') for c in use])
DNAdepth = pd.Series(np.array([r[2] for r in rows]).mean(1), index=DNA.index)
DNA.to_pickle(O + 'kernel_dna_share.pkl'); RNA.to_pickle(O + 'kernel_rna_cpm.pkl')
lR = np.log1p(RNA * 10); lD = np.log1p(DNA.reindex(W.index).fillna(0) * 1e6)
okR = RNAreads >= MIN_RNA; okD = DNAdepth.reindex(W.index).fillna(0) >= MIN_DNA
def corr_rows(M, p):
    Mc = M.sub(M.mean(1), axis=0); pc = p - p.mean()
    return (Mc @ pc) / (np.sqrt((Mc ** 2).sum(1)) * np.sqrt((pc ** 2).sum()) + 1e-12)
out, mem = [], []
for sx, sy in SEEDS:
    d = np.hypot(W.tsne_x - sx, W.tsne_y - sy); cand = W.index[d <= R]
    members = set(W.index[(d <= CORE_R) & (okR | okD)])
    if not members: out.append(dict(seed=f'{sx},{sy}', note='no expressed/covered contig in core')); continue
    for it in range(6):
        m = list(members); w = W.length[m]
        pR = (lR.loc[m].mul(w, axis=0).sum() / w.sum()); pD = (lD.loc[m].mul(w, axis=0).sum() / w.sum())
        cR = corr_rows(lR.loc[cand], pR).where(okR[cand]); cD = corr_rows(lD.loc[cand], pD).where(okD[cand])
        r = pd.concat([cR, cD], axis=1).mean(1)
        score = np.exp(-d[cand] ** 2 / (2 * SIGMA ** 2)) * r
        new = set(score.index[score >= THR])
        ch = W.length[list(new ^ members)].sum() / max(W.length[list(new | members)].sum(), 1)
        members = new
        if ch < 0.02 or not members: break
    g = W.loc[list(members)]; L = g.length.sum()
    tl = g.length[g.tiara.notna()].sum()
    out.append(dict(seed=f'{sx},{sy}', iterations=it + 1, contigs=len(g), bp=int(L), gc=round(np.average(g.gc, weights=g.length), 1) if L else np.nan,
                    tiara_euk=round(g.length[g.tiara == 'eukarya'].sum() / tl, 2) if tl else np.nan, tiara_bac=round(g.length[g.tiara == 'bacteria'].sum() / tl, 2) if tl else np.nan,
                    virus_bp=round(g.length[g.genomad_domain == 'Viruses'].sum() / L, 2) if L else np.nan,
                    rms=round(np.sqrt(np.average((g.tsne_x - g.tsne_x.mean()) ** 2 + (g.tsne_y - g.tsne_y.mean()) ** 2, weights=g.length)), 1) if L else np.nan,
                    cum_ppm=round(RNA.loc[list(members)].sum().sum()) if L else 0))
    mem += [dict(contig=c, seed=f'{sx},{sy}') for c in members]
K = pd.DataFrame(out); K.to_csv(O + 'kernel_blobs.tsv', sep='\t', index=False); KM = pd.DataFrame(mem); KM.to_csv(O + 'kernel_blob_contigs.tsv', sep='\t', index=False)
cum = RNA.sum(1); q = cum.quantile(.999)
fig, axs = plt.subplots(3, 4, figsize=(16, 12))
for ax, (sx, sy) in zip(axs.ravel(), SEEDS):
    V = W[np.hypot(W.tsne_x - sx, W.tsne_y - sy) <= 10]
    col = np.where(V.tiara == 'eukarya', '#2b8a3e', np.where(V.tiara == 'bacteria', '#1971c2', '#868e96'))
    ax.scatter(V.tsne_x, V.tsne_y, s=0.5 + 60 * np.sqrt(cum.reindex(V.index).fillna(0) / q), c=col, lw=0, alpha=.3)
    m = KM.contig[KM.seed == f'{sx},{sy}'] if len(KM) else []
    g = W.loc[m] if len(m) else W.iloc[:0]
    ax.scatter(g.tsne_x, g.tsne_y, s=3, facecolors='none', edgecolors='crimson', linewidths=.4)
    k = K[K.seed == f'{sx},{sy}'].iloc[0]
    ax.set_title(f'({sx},{sy}) {k.get("bp", 0) / 1e6:.1f} Mbp, {k.get("contigs", 0)} ctg', fontsize=9)
    ax.plot(sx, sy, 'k+', ms=12); ax.set_xlim(sx - 10, sx + 10); ax.set_ylim(sy - 10, sy + 10); ax.tick_params(labelsize=7)
fig.suptitle('Kernel capture (red rings = members; size = cumulative RNA; green/blue/grey = Tiara eukarya/bacteria/other)'); plt.tight_layout(); plt.savefig(O + 'kernel_blobs.png', dpi=80)
pd.set_option('display.width', 250); print(K.to_string(index=False))
