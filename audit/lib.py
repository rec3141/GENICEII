import sys, json; sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
import numpy as np, pandas as pd
from common import meta
R='/data/dev/rnaseq-exploration/'
D=json.load(open(R+'site/data.json')); C=json.load(open(R+'site/cyclo.json')); H=json.load(open(R+'site/hgt.json'))
S=pd.DataFrame(D['samples']).set_index('id')
FAM={k:pd.Series(v['cpm'],index=S.index) for k,v in D['pathways']['families'].items()}
FAM=pd.DataFrame(FAM)
HC=pd.read_csv(R+'out/hc_cpm.tsv',sep='\t',index_col=0)
m=meta()
def grp(se,pool=None,med=None,phase=None,src=S):
    x=src[src.season==se]
    if pool: x=x[x.pool==pool]
    if med: x=x[x.medium.isin(med if isinstance(med,list) else [med])]
    if phase: x=x[x.phase.isin(phase if isinstance(phase,list) else [phase])]
    return x.index
