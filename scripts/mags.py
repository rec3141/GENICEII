"""MAG-level response: CLR of bin read pairs, EXP oiled vs CTL within season×medium (Mann-Whitney + BH)."""
import sys; sys.path.insert(0, '.')
from common import *
from scipy.stats import mannwhitneyu

m = meta()
bc = table('bin_counts'); tpm = table('bin_tpm')
S = bc.columns; m = m.loc[S]
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
q = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin')
tax = pd.Series({b: gt.get(b.replace('dastool-', ''), '') for b in bc.index})
gen = tax.str.extract(r'g__([^;]*)')[0].fillna('') ; fam = tax.str.extract(r'f__([^;]*)')[0].fillna('')
lab = (gen.where(gen != '', 'f:' + fam)).str.slice(0, 28)
C = pd.DataFrame(clr(bc.values.astype(float)), index=bc.index, columns=S)
res = []
for se in ['2024-25', '2025-26']:
    for me in ['SW', 'SI']:
        a = m.index[(m.season == se) & (m.design == f'EXP_{me}') & (m.phase == 'oiled')]
        b = m.index[(m.season == se) & (m.design == f'CTL_{me}') & (m.phase == 'control')]
        for bn in bc.index:
            if bc.loc[bn, a.union(b)].sum() < 500: continue
            p = mannwhitneyu(C.loc[bn, a], C.loc[bn, b]).pvalue
            res.append((se, me, bn, lab[bn], C.loc[bn, a].mean() - C.loc[bn, b].mean(), tpm.loc[bn, a].median(), tpm.loc[bn, b].median(), p))
R = pd.DataFrame(res, columns=['season', 'medium', 'bin', 'taxon', 'dCLR', 'tpm_oiled', 'tpm_ctl', 'p'])
R['q'] = R.groupby(['season', 'medium'])['p'].transform(lambda p: bh(p))
R = R.join(q[['SCG_completeness', 'SCG_redundancy', 'size']], on='bin')
R.to_csv('../out/mag_oil_response.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
for (se, me), d in R.groupby(['season', 'medium']):
    d = d[d.q < 0.05]
    print(f'\n{se} {me}: {len(d)} MAGs q<0.05 ({(d.dCLR>0).sum()} up, {(d.dCLR<0).sum()} down)')
    print(d.sort_values('dCLR', ascending=False).head(12)[['bin', 'taxon', 'dCLR', 'tpm_oiled', 'tpm_ctl', 'q', 'SCG_completeness']].round(3).to_string(index=False))
    print(d.sort_values('dCLR').head(6)[['bin', 'taxon', 'dCLR', 'tpm_oiled', 'tpm_ctl', 'q']].round(3).to_string(index=False))
# consistent responders across all four comparisons
w = R.pivot_table(index=['bin', 'taxon'], columns=['season', 'medium'], values='dCLR')
qq = R.pivot_table(index=['bin', 'taxon'], columns=['season', 'medium'], values='q')
up = w[(w > 1).all(1) & (qq < 0.05).sum(1).ge(3)]
print('\nUp in all 4 (dCLR>1, q<0.05 in >=3):'); print(up.round(2).to_string())
dn = w[(w < -1).all(1) & (qq < 0.05).sum(1).ge(3)]
print('\nDown in all 4:'); print(dn.round(2).to_string())
