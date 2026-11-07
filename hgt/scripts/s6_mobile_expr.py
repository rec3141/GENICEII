"""Is the diesel response carried by mobile HC genes? Compare MGE-associated vs chromosomal HC gene expression."""
import sys; sys.path.insert(0, '.'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import meta
from hcdef import *
from scipy.stats import mannwhitneyu
H = '/data/dev/rnaseq-exploration/hgt/'
h = pd.read_csv(H + 'hc_context.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
h = h.join(pd.read_csv(H + 'hc_host.tsv', sep='\t', keep_default_na=False).set_index('gene_id'))
h = h.join(pd.read_csv(H + 'hc_expr_summary.tsv', sep='\t', index_col=0))
E = pd.read_pickle(H + 'hc_cpm.pkl')
for c in ['mge_borne', 'mobile_assoc', 'IS_within_10kb', 'plasmid_strict', 'phage_contig', 'genomic_island', 'integron', 'provirus', 'plasmid_any']:
    h[c] = h[c].astype(str).eq('True')
res = {}
for tier in ['confident', 'all']:
    d = h if tier == 'all' else h[h.tier == 'confident']
    for se in ['2024-25', '2025-26']:
        o, c = d[f'oiled_{se}'], d[f'ctl_{se}']
        delta = (o - c).clip(lower=0)
        for flag in ['mobile_assoc', 'mge_borne', 'plasmid_strict', 'plasmid_any', 'IS_within_10kb']:
            m = d[flag]
            res[(tier, se, flag)] = dict(genes=int(m.sum()), frac_genes=round(m.mean() * 100, 1),
                                         frac_oiled_cpm=round(o[m].sum() / o.sum() * 100, 1), frac_ctl_cpm=round(c[m].sum() / c.sum() * 100, 1),
                                         frac_of_induction=round(delta[m].sum() / delta.sum() * 100, 1))
        ex = d[(d.reads_total >= 20)]
        a, b = ex.loc[ex.mobile_assoc, f'l2fc_{se}'], ex.loc[~ex.mobile_assoc, f'l2fc_{se}']
        res[(tier, se, 'l2fc_mobile_vs_chrom')] = dict(n_mobile=len(a), n_chrom=len(b), median_mobile=round(a.median(), 2), median_chrom=round(b.median(), 2),
                                                       p=float(f'{mannwhitneyu(a, b).pvalue:.3g}') if len(a) > 3 else None)
R = pd.DataFrame(res).T
print(R.to_string())
R.to_csv(H + 'mobile_vs_chrom_summary.tsv', sep='\t')
# per category: share of oiled-pool induction carried by mobile genes
rows = []
for (cat), d in h[h.tier == 'confident'].groupby('category'):
    for se in ['2024-25', '2025-26']:
        delta = (d[f'oiled_{se}'] - d[f'ctl_{se}']).clip(lower=0)
        rows.append(dict(category=cat, season=se, induced_cpm=round(delta.sum(), 2),
                         mobile_share=round(delta[d.mobile_assoc].sum() / delta.sum() * 100, 1) if delta.sum() else np.nan,
                         l2fc_total=round(np.log2((d[f'oiled_{se}'].sum() + .05) / (d[f'ctl_{se}'].sum() + .05)), 2)))
P = pd.DataFrame(rows); P.to_csv(H + 'mobile_share_by_category.tsv', sep='\t', index=False)
print(P.pivot(index='category', columns='season', values=['l2fc_total', 'mobile_share']).to_string())
# top induced mobile genes
mob = h[h.mobile_assoc & (h.reads_total >= 20)].copy()
mob['best_l2fc'] = mob[['l2fc_2024-25', 'l2fc_2025-26']].max(axis=1)
cols = ['category', 'tier', 'product', 'contig', 'bin', 'family', 'genus', 'tax_source', 'plasmid_any', 'plasmid_strict', 'phage_contig', 'provirus', 'genomic_island',
        'integron', 'IS_within_10kb', 'IS_dist', 'dGC', 'contig_len', 'contig_depth', 'oiled_2024-25', 'ctl_2024-25', 'l2fc_2024-25', 'oiled_2025-26', 'ctl_2025-26', 'l2fc_2025-26', 'reads_total']
top = mob.sort_values('best_l2fc', ascending=False)[cols]
top.to_csv(H + 'mobile_hc_genes.tsv', sep='\t')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 40)
print(top.head(30)[['category', 'tier', 'product', 'genus', 'family', 'plasmid_any', 'IS_within_10kb', 'genomic_island', 'phage_contig', 'l2fc_2024-25', 'l2fc_2025-26', 'oiled_2025-26', 'reads_total']].round(2).to_string())
