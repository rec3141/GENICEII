import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
SS=pd.read_csv(R+'data/marine_coassembly.sample_summary.tsv',sep='\t',index_col=0)
SS['rrna_of_genes']=SS.pairs_rrna_features/(SS.pairs_cds+SS.pairs_rrna_features)*100
SS['rrna_of_assigned']=SS.pairs_rrna_features/SS.pairs_assigned*100
x=S.join(SS[['rrna_of_genes','rrna_of_assigned','pct_rrna_total']])
print(x.groupby(['season','pool','medium']).rrna_of_genes.agg(['min','median','max','count']).round(1))
c=x[(x.season=='2025-26')&(x.pool=='CTL')]; print(c[['name','medium','day','rrna_of_genes']].round(1).sort_values('day').to_string())
u=pd.Series(D['unbinned_share'],index=S.index); print('unbinned share mean/median/range', u.mean().round(1), u.median(), u.min(), u.max())
BC=pd.read_csv(R+'data/marine_coassembly.bin_counts.tsv',sep='\t',index_col=0); ub=BC.loc['unbinned']/BC.sum()*100
print('unbinned of CDS pairs', ub.mean().round(1), ub.median().round(1), 'pooled', (BC.loc['unbinned'].sum()/BC.values.sum()*100).round(1))
print('n samples', len(S), 'n MAGs', BC.shape[0]-1)
