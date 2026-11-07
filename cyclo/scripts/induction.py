"""Genome-normalised induction of Cycloclasticus HC genes: within-MAG TPM (samples with >=3000 MAG CDS pairs),
EXP early-oiled (first 2 weeks) / late vs EXP pre-oil and vs CTL, per season."""
import sys, json
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
C = '/data/dev/rnaseq-exploration/cyclo/'
G = pd.read_csv(C + 'cyclo_genes_annot.tsv', sep='\t', keep_default_na=False)
pr = np.load(C + 'cyclo_pairs.npy'); wt = np.load(C + 'cyclo_wtpm.npy')
S = samples(); m = meta().loc[S]
CORE = ['dastool-lorbin_184', 'dastool-metabat_9342', 'dastool-metabat_5432', 'dastool-comebin_3787_sub']
cds = (G['type'] == 'CDS').values
rows = []
for b in CORE:
    ib = np.where((G.cbin == b).values & cds)[0]
    tot = pr[ib].sum(0)
    ok = tot >= 3000
    W = np.where(ok[None, :], wt[ib], np.nan)
    for se in ['2024-25', '2025-26']:
        sel = lambda f: ok & f
        e = (m.season == se).values & (m.group == 'EXP').values
        pre = sel(e & (m.phase == 'pre').values); early = sel(e & (m.phase == 'oiled').values & (m.day <= OIL_DAY[f'CMO {se}'] + 14).values)
        late = sel(e & (m.phase == 'oiled').values & (m.day > OIL_DAY[f'CMO {se}'] + 14).values)
        ctl = sel((m.season == se).values & (m.group == 'CTL').values & (m.phase == 'control').values)
        for j, gi in enumerate(ib):
            if G.hc.iloc[gi] in ('', 'chemotaxis'): continue
            f = lambda s: np.nanmean(W[j, s]) if s.sum() else np.nan
            rows.append(dict(bin=b, season=se, gene_id=G.gene_id.iloc[gi], contig=G.contig.iloc[gi], product=G['product'].iloc[gi], cat=G.hc.iloc[gi],
                             n_pre=int(pre.sum()), n_early=int(early.sum()), n_late=int(late.sum()), n_ctl=int(ctl.sum()),
                             w_pre=f(pre), w_early=f(early), w_late=f(late), w_ctl=f(ctl)))
R = pd.DataFrame(rows)
for a in ['early', 'late', 'ctl']:
    R[f'lfc_{a}_vs_pre'] = np.log2((R[f'w_{a}'] + 20) / (R.w_pre + 20))
R['lfc_early_vs_ctl'] = np.log2((R.w_early + 20) / (R.w_ctl + 20))
R.to_csv(C + 'induction.tsv', sep='\t', index=False)
pd.set_option('display.width', 250, 'display.max_colwidth', 50)
print(R.groupby(['bin', 'season'])[['n_pre', 'n_early', 'n_late', 'n_ctl']].first().to_string())
for (b, se), d in R.groupby(['bin', 'season']):
    if d.n_pre.iloc[0] == 0 or d.n_early.iloc[0] == 0: continue
    print(f'\n## {b} {se}: median lfc early/pre {d.lfc_early_vs_pre.median():.2f}  late/pre {d.lfc_late_vs_pre.median():.2f}')
    print(d.sort_values('lfc_early_vs_pre', ascending=False).head(8)[['product', 'cat', 'w_pre', 'w_early', 'w_late', 'w_ctl', 'lfc_early_vs_pre']].round(1).to_string(index=False))
