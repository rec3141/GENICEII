import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
BC=pd.read_csv(R+'data/marine_coassembly.bin_counts.tsv',sep='\t',index_col=0)
SS=pd.read_csv(R+'data/marine_coassembly.sample_summary.tsv',sep='\t',index_col=0)
M=D['mags']; 
gen={}
for b,v in M.items(): gen['dastool-'+b]=v['tax'].get('g') or v['label']
# fall back via gtdbtk file for all bins
G=pd.read_csv(R+'data/gtdbtk_taxonomy.tsv',sep='\t',index_col=0)['classification']
def g_of(b):
    k=b.replace('dastool-','')
    if k in G.index:
        gg=G[k].split(';')[5][3:] if ';' in G[k] else ''
        if gg: return gg
    lab=M.get(k,{}).get('label','')
    return lab.replace(' (Kaiju)','') if lab else 'NA'
genus=pd.Series({b:g_of(b) for b in BC.index})
# denominator: total mapped pairs (reads_mapped/2)
den=BC.sum()
P=BC.groupby(genus).sum().div(den,axis=1)*100
# check vs mag_timecourse cyclo
T=pd.DataFrame(C['mag_timecourse']).set_index('sample')
print('cyclo check', (P.loc['Cycloclasticus',T.index]-T[['lorbin_184','metabat_9342','metabat_5432','comebin_3787_sub','lorbin_722_sub']].sum(axis=1)).abs().max())
for se in ['2024-25','2025-26']:
    sw=grp(se,'EXP','SW','oiled'); si=grp(se,'EXP','SI','oiled')
    r=pd.DataFrame({'SW':P[sw].mean(axis=1),'SI':P[si].mean(axis=1),'SWmed':P[sw].median(axis=1),'SImed':P[si].median(axis=1)})
    r['SI/SW']=r.SI/r.SW; r['SImed/SWmed']=r.SImed/(r.SWmed+1e-9)
    print(se); print(r.loc[['Cycloclasticus','Oleispira','Sulfurimonas','Sulfitobacter']].round(2).to_string())
    print(r.sort_values('SI',ascending=False).head(8).round(2).to_string())
    # time matched days
    dsw=set(S.loc[sw,'day']); dsi=set(S.loc[si,'day']); com=sorted(dsw&dsi)
    sw2=[i for i in sw if S.loc[i,'day'] in com]; si2=[i for i in si if S.loc[i,'day'] in com]
    print(' matched days',com); print(pd.DataFrame({'SW':P[sw2].mean(axis=1),'SI':P[si2].mean(axis=1)}).loc[['Cycloclasticus','Oleispira','Sulfurimonas']].round(2).to_string())
P.to_pickle('genus_pct.pkl')
