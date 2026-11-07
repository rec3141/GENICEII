import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
print(json.dumps(H['counts'],indent=0)[:1500]); print(H['sharing_funnel'])
mv=pd.DataFrame(H['mobile_vs_chrom']); print(mv[mv.tier=='confident'].drop(columns=['tier']).to_string())
print(pd.DataFrame(H['induction_by_genus']).to_string())
