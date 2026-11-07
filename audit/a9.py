import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
ar=['tmoA-F (toluene MO)','nahC/bphC (extradiol DO)','nahAc/Ad (naphthalene DO)','todC (benzene/toluene DO)','bphA (biphenyl DO)','catE/xylE (catechol 2,3-DO)','dmpKLMNOP (phenol MO)','alkB','CYP153 (product)']
mos=[i for i in S.index if S.loc[i,'name'] in ('CMO_EXP_DAY62','CMO_EXP_DAY63')]
print(HC.loc[mos,ar].round(0).to_string())
sw=grp('2024-25','EXP','SW','oiled'); late=[i for i in sw if S.loc[i,'day']>=55]
print('late SW', list(S.loc[late,'name'])); print(HC.loc[late,ar].round(0).to_string())
for lab,ix in [('all oiled SW mean',sw),('late SW mean',late)]:
    ref=HC.loc[ix,ar].mean(); r=HC.loc[mos,ar].div(ref+0.5)
    print(lab); print(ref.round(1).to_dict()); print(r.round(1).to_string())
ref=HC.loc[sw,ar].median(); print('median oiled SW', ref.round(1).to_dict()); print(HC.loc[mos,ar].div(ref+0.5).round(0).to_string())
