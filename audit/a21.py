import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit"); sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
from lib import *
from common import genes, samples
X,g=genes(); Sm=samples()
G=pd.read_csv(R+'data/gtdbtk_taxonomy.tsv',sep='\t',index_col=0)['classification']
kj=pd.read_csv(R+'data/bin_kaiju_genus.tsv',sep='\t',index_col=0)['kaiju_genus']
def gen(b):
    if not isinstance(b,str) or not b or b=='unbinned': return 'unbinned'
    k=b.replace('dastool-','')
    if k in G.index and G[k].split(';')[5][3:]: return G[k].split(';')[5][3:]
    v=kj.get(b,''); return (str(v).split(' (')[0]+'(kaiju)') if isinstance(v,str) and v else 'NA'
kos={'porA':'K00169','hydA':'K05927','aclA':'K15230','fccB':'K17229','pmoA':'K10944','cydA':'K00425'}
KO=g['KO'].fillna('')
ix_s={s:i for i,s in enumerate(Sm)}
for se in ['2024-25','2025-26']:
  oil=[ix_s[s] for s in grp(se,'EXP',['SW','SI'],'oiled')]
  for n,k in kos.items():
    sel=np.where(KO.str.contains(k).values & (g['type'].values=='CDS'))[0]
    c=np.asarray(X[sel][:,oil].sum(1)).ravel()
    df=pd.DataFrame({'bin':g['bin'].values[sel],'c':c}); df['genus']=df.bin.map(gen)
    t=df.groupby('genus').c.sum().sort_values(ascending=False); t=t/t.sum()*100
    print(se,n,len(sel),'top:',t.head(4).round(1).to_dict())
KJ=pd.read_csv(R+'data/kaiju_contigs.tsv',sep='\t',header=None)
print(KJ.head(2))
