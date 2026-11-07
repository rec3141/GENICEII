import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
T=pd.DataFrame(C['mag_timecourse']).set_index('sample')
b=['lorbin_184','metabat_9342','metabat_5432','comebin_3787_sub','lorbin_722_sub']; T['cyclo']=T[b].sum(axis=1)
P=pd.read_pickle('genus_pct.pkl')
for se in ['2024-25','2025-26']:
  x=T[(T.season==se)&(T.pool!='other')].sort_values(['pool','medium','day'])
  x=x.join(P.T[['Sulfitobacter','Sulfurimonas','Colwellia']])
  print(x[['name','pool','medium','day','cyclo','metabat_9342','lorbin_184','Sulfitobacter','Sulfurimonas']].round(2).to_string())
