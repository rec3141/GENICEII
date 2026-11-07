"""Per-HC-gene expression: CPM (of CDS pairs) per sample; oiled vs control summaries per season."""
import sys; sys.path.insert(0, '.'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
h = pd.read_csv('../hc_genes.tsv', sep='\t', keep_default_na=False)
X, g = genes(); S = samples(); m = meta().loc[S]
cds = np.where(g['type'].values == 'CDS')[0]
tot = np.asarray(X[cds].sum(0)).ravel()
E = X[h.row.values].toarray() / tot * 1e6
E = pd.DataFrame(E, index=h.gene_id.values, columns=S)
E.to_pickle('../hc_cpm.pkl')
out = pd.DataFrame(index=h.gene_id.values)
for se in ['2024-25', '2025-26']:
    o = m.index[(m.season == se) & m.design.isin(['EXP_SW', 'EXP_SI']) & (m.phase == 'oiled')]
    c = m.index[(m.season == se) & m.design.isin(['CTL_SW', 'CTL_SI']) & (m.phase == 'control')]
    p = m.index[(m.season == se) & m.design.isin(['EXP_SW', 'EXP_SI', 'CTL_SW', 'CTL_SI']) & (m.phase == 'pre')]
    out[f'oiled_{se}'] = E[o].mean(1); out[f'ctl_{se}'] = E[c].mean(1); out[f'pre_{se}'] = E[p].mean(1)
    out[f'l2fc_{se}'] = np.log2((E[o].mean(1) + 0.05) / (E[c].mean(1) + 0.05))
out['reads_total'] = X[h.row.values].sum(1).A1
out.to_csv('../hc_expr_summary.tsv', sep='\t')
print(out.describe().round(2).to_string())
