import pandas as pd
d='/data/dev/rnaseq-exploration/hgt/'
E=pd.read_csv(d+'hc_expr_summary.tsv',sep='\t',index_col=0); Hh=pd.read_csv(d+'hc_host.tsv',sep='\t',index_col=0); G=pd.read_csv(d+'hc_genes.tsv',sep='\t',index_col=1)
x=E.join(Hh.genus).join(G[['tier','category']])
x=x[(x.tier=='confident')&x.genus.fillna('').str.contains('Pelagibacter')]
print(x[['oiled_2024-25','ctl_2024-25','pre_2024-25','oiled_2025-26','ctl_2025-26','pre_2025-26']].sum().round(1))
