import sys; sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
from common import meta
import pandas as pd; pd.set_option('display.width',250); pd.set_option('display.max_colwidth',90)
m=meta()
x=m[(m.day<=8)&(m.group!='other')].sort_values(['season','day'])

for i,r in x.iterrows(): print(r.season, r.sample_name, r.day, r.ice_thickness_cm, "|", str(r.log_comments)[:200])
