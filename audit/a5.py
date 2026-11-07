import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',250)
rows=[]
for k in HC.columns:
    r={'fam':k}
    for se in ['2024-25','2025-26']:
        for medn,med in [('SW','SW'),('SI','SI'),('both',['SW','SI'])]:
            a=HC.loc[grp(se,'EXP',med,'oiled'),k].mean(); b=HC.loc[grp(se,'CTL',med,'control'),k].mean()
            r[f'{se[2:]}{medn}']=round((a+.5)/(b+.5),1)
        # medians
        a=HC.loc[grp(se,'EXP','SW','oiled'),k].median(); b=HC.loc[grp(se,'CTL','SW','control'),k].median()
        r[f'{se[2:]}SWmed']=round((a+.5)/(b+.5),1)
        # BACI SW (pre available)
        ep=HC.loc[grp(se,'EXP','SW','pre'),k].mean(); cp=HC.loc[grp(se,'CTL','SW','pre'),k].mean()
        a=HC.loc[grp(se,'EXP','SW','oiled'),k].mean(); b=HC.loc[grp(se,'CTL','SW','control'),k].mean()
        r[f'{se[2:]}SWvsPre']=round((a+.5)/(ep+.5),1)
        r[f'{se[2:]}BACI']=round(((a+.5)/(ep+.5))/((b+.5)/(cp+.5)),1)
        r[f'{se[2:]}preE']=round(ep,2); r[f'{se[2:]}preC']=round(cp,2)
    rows.append(r)
T=pd.DataFrame(rows).set_index('fam'); print(T.to_string())
