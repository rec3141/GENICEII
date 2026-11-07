import sys; sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
import pandas as pd
from common import meta
m=meta()
for s in m.index[m.sample_name.str.contains('DAY62|DAY63|MOS|MYSTERY|DAY05_SW|DAY60|DAY59|DAY55_SW')]:
    r=m.loc[s]; print(s, r.sample_name, '|', r.env_medium,'|',r.medium_code,'|', r.note,'|', r.log_comments,'|', r.site_description, '|', r.collection_date)
