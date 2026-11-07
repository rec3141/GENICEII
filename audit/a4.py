import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
K=pd.read_csv(R+'data/ko_cpm.tsv',sep='\t',index_col=0)
sets={'tmoA-F (toluene MO)':['K15760','K15761','K15762','K15763','K15764','K15765'],'nahAc/Ad (naphthalene DO)':['K14579','K14580'],'nahC/bphC (extradiol DO)':['K14583','K00462'],'alkB':['K00496']}
for k,v in sets.items():
    a=K.reindex(v).fillna(0).sum().reindex(HC.index); print(k, 'max diff vs ko_cpm sum', float((a-HC[k]).abs().max()), 'corr', np.corrcoef(a,HC[k])[0,1])
