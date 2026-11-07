import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
BC=pd.read_csv(R+'data/marine_coassembly.bin_counts.tsv',sep='\t',index_col=0)
SS=pd.read_csv(R+'data/marine_coassembly.sample_summary.tsv',sep='\t',index_col=0)
T=pd.DataFrame(C['mag_timecourse']).set_index('sample')
s=T.index[T.name=='CMO_EXP_DAY62'][0]
b=['dastool-lorbin_184','dastool-metabat_9342','dastool-metabat_5432','dastool-comebin_3787_sub','dastool-lorbin_722_sub']
n=BC.loc[b,s].sum(); print('cyc counts',n, 'tc pct',T.loc[s,['lorbin_184','metabat_9342','metabat_5432','comebin_3787_sub','lorbin_722_sub']].sum())
print('implied denom', n/T.loc[s,'metabat_9342':'lorbin_722_sub'].sum()*100 if False else n/ T.loc[s,['lorbin_184','metabat_9342','metabat_5432','comebin_3787_sub','lorbin_722_sub']].sum()*100)
print(SS.loc[s].to_dict()); print('sum binned', BC[s].sum())
