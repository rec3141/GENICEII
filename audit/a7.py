import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',250)
ar=['nahAc/Ad (naphthalene DO)','nahC/bphC (extradiol DO)','tmoA-F (toluene MO)','todC (benzene/toluene DO)','bphA (biphenyl DO)','catE/xylE (catechol 2,3-DO)','dmpKLMNOP (phenol MO)','alkB']
for se in ['2024-25','2025-26']:
    for p in ['EXP','CTL']:
        ix=grp(se,p,'SW','pre'); print(se,p,list(S.loc[ix,'name']))
        print(HC.loc[ix,ar].round(2).to_string())
out={}
for se in ['2024-25','2025-26']:
    ixe=grp(se,'EXP','SW','pre'); ixc=grp(se,'CTL','SW','pre'); ixb=ixe.append(ixc)
    ixe1=[i for i in ixe if S.loc[i,'day']<0]
    out[se+' EXPmean']=HC.loc[ixe,ar].mean(); out[se+' EXPmed']=HC.loc[ixe,ar].median(); out[se+' EXP(day<0)']=HC.loc[ixe1,ar].mean()
    out[se+' CTLmean']=HC.loc[ixc,ar].mean(); out[se+' bothmean']=HC.loc[ixb,ar].mean(); out[se+' bothmed']=HC.loc[ixb,ar].median()
O=pd.DataFrame(out).round(2); print(O.to_string())
print((O['2025-26 bothmean']/O['2024-25 bothmean']).round(1))
print((O['2025-26 EXPmean']/O['2024-25 EXPmean']).round(1))
