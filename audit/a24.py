import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
T=pd.read_csv(R+'data/marine_coassembly.taxonomy_rrna.tsv',sep='\t')
print(T['rank'].value_counts())
ctl=list(S.index[(S.season=='2025-26')&(S.pool=='CTL')])
e=T[T.lineage.str.contains('Eugleno|Euglen',na=False)]
print(e[['rank','lineage']+ctl[:3]].to_string()[:2000])
