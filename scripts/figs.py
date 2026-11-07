"""Figures for the CMO RNA-seq exploration report."""
import sys; sys.path.insert(0, '.')
from common import *
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': '#52514e',
                     'axes.labelcolor': '#0b0b0b', 'xtick.color': '#52514e', 'ytick.color': '#52514e', 'axes.grid': True,
                     'grid.color': '#e6e5e0', 'grid.linewidth': 0.6, 'figure.dpi': 150, 'savefig.bbox': 'tight', 'font.family': 'DejaVu Sans'})
CTL, EXP, GRY = '#2a78d6', '#eb6834', '#8a8984'
F = '/data/dev/rnaseq-exploration/fig/'
m = meta(); s = table('sample_summary'); S = s.index; m = m.loc[S]
core = m.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI'])

# 1. read fate per sample
pairs = s.reads_raw / 2
fr = pd.DataFrame({
    'rRNA (filtered)': (s.reads_raw - s.reads_norrna) / s.reads_raw,
    'rRNA (escaped filter)': s.pairs_rrna_features / pairs,
    'CDS': s.pairs_cds / pairs,
    'mapped, no feature/ambiguous': (s.reads_mapped / 2 - s.pairs_cds - s.pairs_rrna_features) / pairs,
    'unmapped / short local hit': (s.reads_nohuman - s.reads_mapped) / s.reads_raw}).clip(lower=0) * 100
order = m.assign(o=m.group.map({'CTL': 0, 'EXP': 1, 'other': 2})).sort_values(['season', 'o', 'medium', 'day']).index
cols = ['#9db8d9', '#d55181', '#1baf7a', '#c3c2b7', '#e6e5e0']
fig, axs = plt.subplots(2, 1, figsize=(11, 5.2), sharey=True)
for ax, se in zip(axs, ['2024-25', '2025-26']):
    o = [i for i in order if m.loc[i, 'season'] == se]; x = np.arange(len(o)); bot = np.zeros(len(o))
    for c, col in zip(fr.columns, cols):
        ax.bar(x, fr.loc[o, c], bottom=bot, color=col, width=0.85, label=c, linewidth=0); bot += fr.loc[o, c].values
    lab = [f"{m.loc[i,'group']} {m.loc[i,'medium']} {'' if np.isnan(m.loc[i,'day']) else int(m.loc[i,'day'])}" for i in o]
    ax.set_xticks(x, lab, rotation=90, fontsize=5.5); ax.set_xlim(-0.6, len(o) - 0.4); ax.grid(False)
    ax.set_ylabel('% of read pairs'); ax.set_title(f'CMO {se}', loc='left', fontsize=10)
axs[0].legend(ncol=5, frameon=False, fontsize=7.5, loc='lower left', bbox_to_anchor=(0, 1.08))
fig.tight_layout(); fig.savefig(F + 'f1_read_fate.png'); plt.close()

# 2. PCA
p = pd.read_csv('../out/pca_mags.tsv', sep='\t', index_col=0)
fig, ax = plt.subplots(figsize=(6.4, 4.6))
for (g, med), d in p.join(m[['group', 'medium']]).groupby(['group', 'medium']):
    col = {'CTL': CTL, 'EXP': EXP}.get(g, GRY); mk = {'SW': 'o', 'SI': '^'}.get(med, 's')
    pre = d.phase.eq('pre')
    ax.scatter(d.PC1[~pre], d.PC2[~pre], s=26, c=col, marker=mk, edgecolors='white', linewidths=0.6, label=f'{g} {med}')
    ax.scatter(d.PC1[pre], d.PC2[pre], s=26, facecolors='none', edgecolors=col, marker=mk, linewidths=1.1)
ax.text(0.02, 0.97, 'hollow = before oil added', transform=ax.transAxes, fontsize=7.5, color='#52514e', va='top')
ax.text(-14, 30, 'CMO 2024-25', fontsize=9.5, color='#0b0b0b'); ax.text(22, 10, 'CMO 2025-26', fontsize=9.5, color='#0b0b0b')
ax.set_xlabel('PC1 (46.7%)'); ax.set_ylabel('PC2 (13.1%)')
ax.set_title('MAG transcript composition (Aitchison PCA, 911 MAGs)', loc='left', fontsize=10)
ax.legend(frameon=False, fontsize=7.5, loc='lower right', ncol=2)
fig.savefig(F + 'f2_pca.png'); plt.close()

# 3. key MAG time courses
tpm = table('bin_tpm')
groups = {'Cycloclasticus (3 MAGs)': ['dastool-lorbin_184', 'dastool-metabat_9342', 'dastool-metabat_5432'],
          'Porticoccaceae: 50-400-T64 + Porticoccus (4)': ['dastool-comebin_2576', 'dastool-metabat_1191', 'dastool-comebin_3425', 'dastool-lorbin_288'],
          'Sulfitobacter lorbin_032': ['dastool-lorbin_032'],
          'Casp-alpha2 (Rhodospirillales) comebin_2596': ['dastool-comebin_2596'],
          'Oleispira lorbin_369': ['dastool-lorbin_369'],
          'Sulfurimonas (2 MAGs)': ['dastool-lorbin_167', 'dastool-semibin_088']}
def tc(series, fname, title, ylab):
    n = len(series)
    fig, axs = plt.subplots(n, 2, figsize=(9, 1.75 * n), sharex='col')
    for i, (name, v) in enumerate(series.items()):
        for j, se in enumerate(['2024-25', '2025-26']):
            ax = axs[i, j]
            for g, col in (('CTL', CTL), ('EXP', EXP)):
                for med, ls, mk in (('SW', '-', 'o'), ('SI', '--', '^')):
                    d = m[core & (m.season == se) & (m.group == g) & (m.medium == med)]
                    y = v.reindex(d.index); t = d.assign(y=y).groupby('day').y.mean()
                    ax.plot(t.index, t.values + 1, ls=ls, color=col, lw=1.6, marker=mk, ms=3.2, label=f'{g} {med}')
            ax.set_yscale('log'); ax.axvline(OIL_DAY[f'CMO {se}'], color='#0b0b0b', lw=0.8, ls=':')
            if i == 0: ax.set_title(f'CMO {se}', fontsize=10, loc='left')
            if j == 0: ax.set_ylabel(ylab[i] if isinstance(ylab, list) else ylab, fontsize=7.5)
            ax.text(0.01, 0.97, name, transform=ax.transAxes, fontsize=7.5, va='top', color='#0b0b0b')
    for ax in axs[-1]: ax.set_xlabel('experiment day (dotted = oil added)')
    axs[0, 1].legend(ncol=4, frameon=False, fontsize=7, loc='lower right', bbox_to_anchor=(1, 1.18))
    fig.suptitle(title, x=0.01, ha='left', fontsize=10.5, y=1.0); fig.tight_layout(); fig.savefig(F + fname); plt.close()
tc({k: tpm.loc[v].sum() for k, v in groups.items()}, 'f3_mag_timecourse.png', 'Oil-responsive MAGs, transcripts per million (+1)', 'TPM + 1')

# 4. hydrocarbon genes
hc = pd.read_csv('../out/hc_cpm.tsv', sep='\t', index_col=0)
ko = pd.read_csv(D + 'ko_cpm.tsv', sep='\t', index_col=0)  # coordinate-remapped KOs; the pipeline ko_tpm joins KOs to the wrong genes
genes_ = {'alkB alkane monooxygenase (K00496)': hc['alkB'], 'todC1/C2 benzene/toluene dioxygenase': hc['todC (benzene/toluene DO)'],
          'nahC/bphC extradiol ring-cleavage': hc['nahC/bphC (extradiol DO)'], 'catE/xylE catechol 2,3-dioxygenase': hc['catE/xylE (catechol 2,3-DO)'],
          'maleylacetate reductase (K00217)': ko.loc['K00217'], 'tmoA-F toluene monooxygenase': hc['tmoA-F (toluene MO)']}
tc(genes_, 'f4_hc_genes.png', 'Hydrocarbon-degradation transcripts per million CDS read pairs (+1)', 'CPM + 1')
print('ok')
