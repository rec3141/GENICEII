"""Zoom panel of a t-SNE window coloured by per-contig MarFERReT protein lineage, for drawing hulls around lineage clusters
(2026-11-08, user: copepod EE12/EE14/EE17 belong to one larger cluster at x -20..0, y -65..-85).
Usage: lineage_panel.py NAME x0 x1 y0 y1 RANK TAXON [TAXON2 ...]
Per contig (>= 1 kb): MarFERReT best hits of its MetaEuk proteins (out/identify/grex/marferret_proteins.tsv), lineage via
identify.named(); the contig is 'TAXON' when >= 50% of its classified proteins (aln-length weighted) are TAXON at RANK,
'other euk' when it has classified proteins otherwise, 'no protein' when none. Points sized by cumulative RNA.
Panel 2: Tiara. All curated polygons (hulls.tsv) are outlined; 1-unit grid.
Out: out/hulls/check/<NAME>_lineage.png and out/hulls/check/<NAME>_lineage_contigs.tsv"""
import sys, gzip
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
import identify
O = '/data/dev/rnaseq-exploration/euks/mags/'
name, x0, x1, y0, y1, rank = sys.argv[1], *map(float, sys.argv[2:6]), sys.argv[6]; taxa = set(sys.argv[7:])

T = pd.read_pickle(O + 'out/contigs_full.pkl').set_index('id')
W = T[T.tsne_x.between(x0, x1) & T.tsne_y.between(y0, y1) & (T.length >= 1000)].copy()
P = pd.read_csv(O + 'out/identify/grex/marferret_proteins.tsv', sep='\t', usecols=['protein_id', 'best_hit', 'aln_length'])
P['contig'] = P.protein_id.str.split('|').str[1]; P = P[P.contig.isin(W.index)]
want = set(P.best_hit); tx = {}
with gzip.open(identify.MARF_TAX, 'rt') as fh:
    next(fh)
    for l in fh:
        f = l.split('\t')
        if f[1] in want: tx[f[1]] = int(f[2])
P['taxid'] = P.best_hit.map(tx); P = P[P.taxid.notna()]
lin = {t: identify.named(int(t)) for t in P.taxid.unique()}
P['hit'] = [lin[t].get(rank) in taxa for t in P.taxid]
s = P.groupby('contig').apply(lambda g: (g.aln_length * g.hit).sum() / g.aln_length.sum())
W['frac'] = s.reindex(W.index)
W['cls'] = np.where(W.frac >= 0.5, '/'.join(sorted(taxa)), np.where(W.frac.notna(), 'other euk', 'no protein'))
cum = pd.read_pickle(O + 'out/cum_ppm.pkl').cum_ppm.reindex(W.index).fillna(0)
W['size'] = 2 + 60 * (cum / cum.max()) ** 0.5
H = pd.read_csv(O + 'hulls.tsv', sep='\t', comment='#')

fig, ax = plt.subplots(1, 2, figsize=(22, 11), sharex=True, sharey=True)
col = {'/'.join(sorted(taxa)): '#d62728', 'other euk': '#1f77b4', 'no protein': '#cccccc'}
for k in ['no protein', 'other euk', '/'.join(sorted(taxa))]:
    d = W[W.cls == k]; ax[0].scatter(d.tsne_x, d.tsne_y, s=d['size'], c=col[k], alpha=.6, lw=0, label=f'{k} ({len(d)})')
tc = {'eukarya': '#2ca02c', 'bacteria': '#9467bd', 'archaea': '#8c564b', 'organelle': '#ff7f0e', 'unknown': '#bbbbbb'}
for k, c in tc.items():
    d = W[W.tiara.astype(str).str.contains(k)]; ax[1].scatter(d.tsne_x, d.tsne_y, s=d['size'], c=c, alpha=.5, lw=0, label=k)
for a in ax:
    for hid, g in H.groupby('id'):
        v = g.sort_values('vertex_order')[['x', 'y']].values
        if ((v[:, 0] > x0 - 2) & (v[:, 0] < x1 + 2) & (v[:, 1] > y0 - 2) & (v[:, 1] < y1 + 2)).any():
            v = np.vstack([v, v[:1]]); a.plot(v[:, 0], v[:, 1], 'k-', lw=1); a.text(v[:, 0].mean(), v[:, 1].max() + .3, hid, fontsize=9)
    a.set_xticks(np.arange(np.floor(x0), x1 + 1, 1), minor=True); a.set_yticks(np.arange(np.floor(y0), y1 + 1, 1), minor=True)
    a.grid(which='both', lw=.3, alpha=.5); a.set_xlim(x0, x1); a.set_ylim(y0, y1); a.legend(loc='upper right', fontsize=9)
ax[0].set_title(f'{rank} in {sorted(taxa)} (>= 50% of classified protein aln length), size = cumulative RNA')
ax[1].set_title('Tiara')
fig.savefig(O + f'out/hulls/check/{name}_lineage.png', dpi=90, bbox_inches='tight')
W[['tsne_x', 'tsne_y', 'length', 'tiara', 'frac', 'cls']].to_csv(O + f'out/hulls/check/{name}_lineage_contigs.tsv', sep='\t')
print(W.groupby('cls').length.agg(['size', 'sum']))
