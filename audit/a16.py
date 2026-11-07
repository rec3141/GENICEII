import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
P=pd.DataFrame(D['families']['pct'],index=S.index)
x=P.join(S[['season','pool','medium','day','phase']])
x=x[(x.medium=='SW')&(x.pool!='other')]
t=x.groupby(['season','day','pool'])[['Cycloclasticaceae','Rhodobacteraceae','Sulfurimonadaceae','Porticoccaceae','Flavobacteriaceae']].mean().unstack('pool').round(1)
print(t.to_string())
