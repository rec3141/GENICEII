import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300); pd.set_option('display.max_colwidth',70)
K=pd.read_csv(R+'out/ko_oil_response.tsv',sep='\t',index_col=0)
c=['2024-25 SW','2024-25 SI','2025-26 SW','2025-26 SI']
print(K.describe().round(2).to_string()); print(K.n_sig.value_counts().sort_index())
neg=K[(K[c]<0).all(axis=1)]; pos=K[(K[c]>0).all(axis=1)]
print('all-neg',len(neg),'all-pos',len(pos), 'total',len(K))
print('all neg & n_sig==4', ((K[c]<0).all(axis=1)&(K.n_sig==4)).sum(), 'all pos & n_sig==4', ((K[c]>0).all(axis=1)&(K.n_sig==4)).sum())
pat='notch|wif|protocadherin|hsp90|htpG|psbZ|heat shock'
print(K[K.name.str.contains(pat,case=False,na=False)].to_string())
print(neg[neg.n_sig==4].sort_values('2025-26 SW').head(40).to_string())
for k in ['K15730','K21966','K16506','K04079','K09903','K02725']:
    if k in K.index: print(k, K.loc[k].to_dict())
