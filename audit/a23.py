import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit"); sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
from lib import *
from common import genes, samples
X,g=genes(); Sm=samples(); ixs={s:i for i,s in enumerate(Sm)}
print(g['type'].value_counts())
r=np.where(g['type'].values!='CDS')[0]
ctl=[ixs[s] for s in S.index[(S.season=='2025-26')&(S.pool=='CTL')]]
c=np.asarray(X[r][:,ctl].sum(1)).ravel()
df=pd.DataFrame({'gene':g['gene_id'].values[r],'contig':g['contig'].values[r],'type':g['type'].values[r],'prod':g['product'].values[r],'c':c}).sort_values('c',ascending=False)
df['pct']=df.c/df.c.sum()*100; print(df.head(6).to_string())
per=[]
top=df.index[0]
for s in ctl:
    tot=X[r][:,s].sum(); per.append(X[top, s]/tot*100 if False else float(X[g.index.get_loc(top) if False else top, s])/float(tot)*100)
print('per-sample share of top gene among rRNA', np.round(per,1))
