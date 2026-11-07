"""Hydrocarbon-degradation gene expression (CPM of CDS read pairs) by pool, medium and phase."""
import sys; sys.path.insert(0, '.')
from common import *
from scipy.stats import mannwhitneyu
HC = {'alkB': ['K00496'], 'alkT/rubB': ['K05297'], 'ladA (product)': [], 'CYP153 (product)': [],
      'nahAc/Ad (naphthalene DO)': ['K14579', 'K14580'], 'nidAB (PAH DO)': ['K11943', 'K11944'],
      'bphA (biphenyl DO)': ['K08689', 'K15750'], 'todC (benzene/toluene DO)': ['K03268', 'K16268'],
      'tmoA-F (toluene MO)': ['K15760', 'K15761', 'K15762', 'K15763', 'K15764', 'K15765'],
      'dmpKLMNOP (phenol MO)': ['K16242', 'K16243', 'K16244', 'K16245', 'K16246', 'K16249'],
      'catA (catechol 1,2-DO)': ['K03381'], 'catE/xylE (catechol 2,3-DO)': ['K00446', 'K07104'],
      'pcaGH (protocatechuate 3,4-DO)': ['K00448', 'K00449'], 'gentisate 1,2-DO': ['K00450'],
      'nahC/bphC (extradiol DO)': ['K14583', 'K00462'], 'pmoA/mmoX (methane MO)': ['K10944', 'K16157'],
      'prmA (propane MO)': ['K18223']}
X, g = genes(); S = samples(); m = meta().loc[S]
cds = g['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
kos = g['KO'].str.split(',')
prod = g['product'].str.lower()
rows, members = {}, {}
for k, v in HC.items():
    if k.startswith('ladA'): sel = prod.str.contains('long-chain alkane monooxygenase|ladA', regex=True).values
    elif k.startswith('CYP153'): sel = prod.str.contains('alkane') & prod.str.contains('p450|cytochrome')
    else: sel = kos.apply(lambda x: any(i in v for i in x)).values
    sel = np.asarray(sel) & cds
    idx = np.where(sel)[0]; members[k] = idx
    rows[k] = np.asarray(X[idx].sum(0)).ravel() / tot * 1e6
E = pd.DataFrame(rows, index=S)
E.to_csv('../out/hc_cpm.tsv', sep='\t')
pd.to_pickle(members, '../out/hc_members.pkl')
x = E.join(m[['season', 'design', 'phase']])
core = x[x.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI'])]
print('genes per category:', {k: len(v) for k, v in members.items()})
print(core.groupby(['season', 'design', 'phase']).median().T.round(2).to_string())
print('\nEXP oiled vs CTL (same season+medium), Mann-Whitney:')
for (se, me), d in core.groupby(['season', core.design.str[-2:]]):
    a = d[d.phase == 'oiled']; b = d[d.phase == 'control']
    for k in E.columns:
        if a[k].sum() + b[k].sum() == 0: continue
        p = mannwhitneyu(a[k], b[k]).pvalue
        fc = (a[k].median() + .1) / (b[k].median() + .1)
        if p < 0.05: print(f'  {se} {me} {k:34s} n={len(a)}/{len(b)} median {a[k].median():8.2f} vs {b[k].median():8.2f}  fc={fc:6.1f}  p={p:.2g}')
