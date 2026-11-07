import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
B=pd.read_csv(R+'out/mag_baci.tsv',sep='\t')
print(B[B.bin.isin(['dastool-lorbin_032','dastool-comebin_2596','dastool-metabat_1122','dastool-vamb_549'])].to_string())
BC=pd.read_csv(R+'data/marine_coassembly.bin_counts.tsv',sep='\t',index_col=0); P=BC/BC.sum()*100
for b in ['dastool-lorbin_032','dastool-comebin_2596']:
  x=S[(S.season=='2024-25')&(S.medium=='SW')&(S.pool!='other')].copy(); x['v']=P.loc[b,x.index]
  print(b, x.sort_values(['pool','day'])[['name','day','v']].round(2).to_string())
