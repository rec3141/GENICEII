import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
print(S.columns.tolist()); print(S.groupby(['season','pool','medium','phase']).size())
# compare FAM vs HC
com=[c for c in HC.columns if c in FAM.columns]
print('max abs diff fam vs hc_cpm', (FAM[com]-HC.loc[FAM.index,com]).abs().max().round(3).to_dict())
