import sys; sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
import pandas as pd
from common import meta
m=meta()
pd.set_option('display.width',250); pd.set_option('display.max_rows',300)
print(m.columns.tolist())
print(m[['sample_name','season','group','medium','day','phase']].sort_values(['season','group','medium','day']).to_string())
