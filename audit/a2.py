import sys; sys.path.insert(0,'/data/dev/rnaseq-exploration/scripts')
import pandas as pd
from common import meta
m=meta()
pd.set_option('display.max_colwidth',200)
print(m.loc[['V350357237_L02_F05','V350357237_L02_F06','V350357237_L02_F04']].T.to_string())
