import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
K=pd.read_csv(R+'out/ko_oil_response.tsv',sep='\t',index_col=0)
c=['2024-25 SW','2024-25 SI','2025-26 SW','2025-26 SI']
d=K[(K[c]<-1).all(axis=1)&(K.n_sig>=3)]; u=K[(K[c]>1).all(axis=1)&(K.n_sig>=3)]
print('down',len(d),'up',len(u)); print(d.to_string())
for k in ['K00169','K05927','K15230','K17229']: print(K.loc[k].to_dict() if k in K.index else k)
