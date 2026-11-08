"""Heatmap of the hull-curated eukaryotic MAGs (EE##) against DNA and RNA relative abundance (user request 2026-11-08).
Per contig (>= 1 kb, as in contig_hull_mag.tsv):
  RNA = reads per kb per million CDS-library reads (contig counts, 128 analysed libraries);
  DNA = depth x 1e9 / sample's total mapped bp (depth*length summed over all contigs), CMO nanopore libraries.
RNA: a per-sample median over contigs is 0 for most MAGs (eukaryotic contigs are sparse), so the median is taken over contig
profiles: contigs detected in >= MIN_DET of libraries; each contig's profile is scaled to mean 1; MAG profile = per-sample median of
those, times the MAG's median contig level. A few very highly expressed contigs (e.g. EE01's two) therefore cannot dominate.
MAGs with < MIN_N such contigs are not shown (the 'little CMO RNA' groups).
DNA: nanopore depth per contig per library is too sparse for any contig median; length-weighted mean depth per Gbp mapped.
Transform: 4th root. Rows: Ward clustering on both blocks, each block z-scored per row after the transform.
Columns: Ward clustering within each block. Colour: viridis, a separate scale per block.
Out: out/hull_mags_relabund_{rna,dna}.tsv, out/hull_mags_heatmap.png"""
import sys, gzip
import numpy as np, pandas as pd, scipy.io
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import samples, all_samples, meta, nanopore_samples, D as DD
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
MIN_DET, MIN_N = 0.25, 5

H = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')
H = H[H.mag.str.startswith('EE')]
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
L = T.length.reindex(H.contig).values / 1e3

# RNA
Cm = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsr()
ALL = all_samples(); S = samples(); kc = [ALL.index(x) for x in S]
pos = pd.Series(np.arange(len(T)), index=np.sort(T.index.values))
tot = np.asarray(Cm[:, kc].sum(0)).ravel()
rr = Cm[pos[H.contig].values][:, kc].toarray()
RNA_c = pd.DataFrame(rr / tot * 1e6 / L[:, None], index=H.contig.values, columns=S)

# DNA
ns = nanopore_samples()
hdr = gzip.open(O + '../dna/depths.txt.gz', 'rt').readline().rstrip('\n').split('\t')
dcols = [c for c in hdr[3:] if not c.endswith('-var') and c.replace('.sorted.bam', '') in set(ns.run_barcode)]
idx = [hdr.index(c) for c in dcols]; need = set(H.contig); rows = {}; totals = np.zeros(len(idx))
with gzip.open(O + '../dna/depths.txt.gz', 'rt') as fh:
    next(fh)
    for line in fh:
        f = line.rstrip('\n').split('\t'); v = np.array([float(f[i]) for i in idx])
        totals += v * float(f[1])
        if f[0] in need: rows[f[0]] = v
DNA_c = pd.DataFrame.from_dict(rows, orient='index', columns=[c.replace('.sorted.bam', '') for c in dcols]) / totals * 1e9
DNA_c = DNA_c.reindex(H.contig.values).fillna(0)

g = H.mag.values; kb = pd.Series(L, index=H.contig.values)
def robust(C, name):
    # contigs detected (>0) in >= MIN_DET of the block's libraries carry a usable per-sample profile
    det = (C > 0).mean(axis=1); ok = det >= MIN_DET
    pooled = C.mean(axis=1)
    level = pooled[ok].groupby(g[ok.values]).median()                       # typical contig level of the MAG
    rel = C[ok].div(pooled[ok], axis=0)                                      # each contig's profile, mean 1
    prof = rel.groupby(g[ok.values]).median()                                # per-sample median over contigs
    n = ok.groupby(g).sum(); print(name, 'contigs used per MAG (detected in >=', MIN_DET, 'of libraries):'); print(n.to_string())
    out = prof.mul(level, axis=0)
    return out[n.reindex(out.index) >= MIN_N]
RNA = robust(RNA_c, 'RNA')
# DNA: nanopore depth per contig per library is too sparse for a contig median (few contigs are covered in >= 25% of libraries);
# the length-weighted mean is used (no expression-style outlier contigs in DNA: e.g. EE01's top contigs carry ~8% of DNA signal)
DNA = DNA_c.mul(kb, axis=0).groupby(g).sum().div(kb.groupby(g).sum(), axis=0)
DNA = DNA.reindex(RNA.index)
print('MAGs shown:', len(RNA), '; dropped (< MIN_N RNA-detected contigs):', sorted(set(g) - set(RNA.index)))
RNA.to_csv(O + 'hull_mags_relabund_rna.tsv', sep='\t'); DNA.to_csv(O + 'hull_mags_relabund_dna.tsv', sep='\t')

# column annotations
m = meta().reindex(S)
rna_ann = pd.DataFrame({'season': m.season.values, 'group': m.group.values, 'medium': m.medium.values}, index=S)
nsx = ns.set_index('run_barcode').reindex(DNA.columns)
season = nsx.log.astype(str).map({'2024': '2024-25', '2025': '2025-26'}).fillna('other')
nm = nsx.sample_name.fillna('')
grp = np.where(nm.str.contains('_EXP_'), 'EXP', np.where(nm.str.contains('_CTL_'), 'CTL', 'other'))
med = np.where(nm.str.contains('_SW'), 'SW', np.where(nm.str.contains('_SI'), 'SI', 'other'))
dna_ann = pd.DataFrame({'season': season.values, 'group': grp, 'medium': med}, index=DNA.columns)

R4, D4 = RNA ** 0.25, DNA ** 0.25
z = lambda X: X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1).replace(0, 1), axis=0)
ro = leaves_list(linkage(pd.concat([z(R4), z(D4)], axis=1).values, 'ward'))
co_r = leaves_list(linkage(R4.T.values, 'ward')); co_d = leaves_list(linkage(D4.T.values, 'ward'))
rows_ord = R4.index[ro]

PAL = {'season': {'2024-25': '#9ecae1', '2025-26': '#08519c', 'other': '#d9d9d9'},
       'group': {'EXP': '#d95f02', 'CTL': '#1b9e77', 'other': '#d9d9d9'},
       'medium': {'SW': '#6baed6', 'SI': '#9e9ac8', 'other': '#f0f0f0'}}
fig = plt.figure(figsize=(22, 11))
gs = fig.add_gridspec(4, 5, width_ratios=[R4.shape[1], D4.shape[1], 3, 12, 3], height_ratios=[1, 1, 1, 26], wspace=0.06, hspace=0.05)
def block(col, X, ann, co, title, cb_col):
    for i, k in enumerate(['season', 'group', 'medium']):
        a = fig.add_subplot(gs[i, col]); c = [PAL[k].get(v, '#f0f0f0') for v in ann[k].values[co]]
        a.imshow([[plt.matplotlib.colors.to_rgb(x) for x in c]], aspect='auto'); a.set_xticks([]); a.set_yticks([])
        if col == 0: a.set_yticks([0]); a.set_yticklabels([k], fontsize=8)
        if i == 0: a.set_title(title, fontsize=11)
    ax = fig.add_subplot(gs[3, col]); M = X.loc[rows_ord].values[:, co]
    im = ax.imshow(M, aspect='auto', cmap='viridis', interpolation='nearest')
    ax.set_xticks([]); ax.set_xlabel(f'{X.shape[1]} libraries (Ward order)', fontsize=9)
    if col == 0: ax.set_yticks(range(len(rows_ord))); ax.set_yticklabels(rows_ord, fontsize=8)
    else: ax.set_yticks([])
    cax = fig.add_subplot(gs[3, cb_col]); cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=7)
    cb.set_label(title.split(':')[0] + ' (4th root)', fontsize=8)
block(0, R4, rna_ann, co_r, 'RNA: reads per kb per million', 2)
block(1, D4, dna_ann, co_d, 'DNA: depth per Gbp mapped', 4)
hs = [Patch(color=c, label=f'{k}: {v}') for k in PAL for v, c in PAL[k].items() if v != 'other'] + [Patch(color='#d9d9d9', label='other / unassigned')]
fig.legend(handles=hs, loc='lower center', ncol=7, fontsize=8, frameon=False)
fig.suptitle('Eukaryotic hull MAGs (EE##): RNA = median over contig profiles, DNA = length-weighted mean; 4th root, Ward clustering', fontsize=13)
fig.savefig(O + 'hull_mags_heatmap.png', dpi=130, bbox_inches='tight')
print('rows', len(rows_ord), 'RNA libs', R4.shape[1], 'DNA libs', D4.shape[1])
print(pd.DataFrame({'rna_median': RNA.median(axis=1).round(2), 'dna_median': DNA.median(axis=1).round(2)}).loc[rows_ord].to_string())
