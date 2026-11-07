import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',250)
T=pd.DataFrame(C['mag_timecourse']).set_index('sample')
b=['lorbin_184','metabat_9342','metabat_5432','comebin_3787_sub','lorbin_722_sub']
T['cyclo']=T[b].sum(1)
x=T[(T.pool=='EXP')&(T.day>=50)].sort_values(['season','day'])
print(x[['name','season','medium','day']+b+['cyclo']].round(2).to_string())
print(x.groupby(['season','medium']).cyclo.agg(['min','max','mean']).round(1))
y=T[(T.pool=='CTL')&(T.day>=50)]; print('CTL late max', y.cyclo.max())
# cross-check from bin_counts
BC=pd.read_csv(R+'data/marine_coassembly.bin_counts.tsv',sep='\t',index_col=0)
print(BC.index[:5].tolist(), BC.shape)
