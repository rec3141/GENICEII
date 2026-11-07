import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
for ps in [0.1,0.5]:
  print('pseudo',ps)
  for se in ['2024-25','2025-26']:
    for med in ['SW','SI']:
      a=HC.loc[grp(se,'EXP',med,'oiled')].median(); b=HC.loc[grp(se,'CTL',med,'control')].median()
      print(se,med, ((a+ps)/(b+ps)).round(1)[['tmoA-F (toluene MO)','todC (benzene/toluene DO)','nahC/bphC (extradiol DO)','nahAc/Ad (naphthalene DO)','catE/xylE (catechol 2,3-DO)','alkB','CYP153 (product)','prmA (propane MO)','pmoA/mmoX (methane MO)']].to_dict())
