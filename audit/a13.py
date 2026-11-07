import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
A=pd.read_csv(R+'out/dna_rel_abundance.tsv',sep='\t',index_col=0)
print('colsum range', A.sum().min(), A.sum().max())
A=A/A.sum()*100
G=pd.read_csv(R+'data/gtdbtk_taxonomy.tsv',sep='\t',index_col=0)['classification']
def g_of(b):
    k=b.replace('dastool-','')
    if k in G.index and G[k].split(';')[5][3:]: return G[k].split(';')[5][3:]
    lab=D['mags'].get(k,{}).get('label','');  return lab.replace(' (Kaiju)','') if lab else 'NA'
gen=pd.Series({b:g_of(b) for b in A.index})
GA=A.groupby(gen).sum()
cols=[c for c in A.columns]
ds=pd.read_csv(R+'out/dna_samples.tsv',sep='\t').set_index('key')
pre=[c for c in cols if ('PRE' in c) or c.startswith('2024-25|CMO_') and ('DAY05' in c)]
sel=['Cycloclasticus','Sulfurimonas','Sulfitobacter','Casp-alpha2']
gen2=gen.copy(); fam=pd.Series({b:(G[b.replace('dastool-','')].split(';')[4][3:] if b.replace('dastool-','') in G.index else '') for b in A.index})
GF=A.groupby(fam).sum()
show=[c for c in cols if 'PRE' in c or 'DAY05' in c or 'SITE' in c or 'DAY01' in c or 'DAY02' in c or 'DAY03' in c or 'DAY04' in c]
X=pd.concat([GA.loc[['Cycloclasticus','Sulfurimonas','Sulfitobacter']],GF.loc[['Casp-alpha2']]]).T
X['lorbin_184']=A.loc['dastool-lorbin_184']; X['metabat_9342']=A.loc['dastool-metabat_9342']; X['lorbin032']=A.loc['dastool-lorbin_032']; X['comebin_2596']=A.loc['dastool-comebin_2596']
print(X.loc[sorted(show)].round(3).to_string())
cx=A.loc[['dastool-lorbin_184','dastool-comebin_3787_sub','dastool-lorbin_722_sub']].sum()
pc=[c for c in cols if c.startswith('2025-26') and 'PRE' in c]
print('complex %binned 2025-26 pre', cx[pc].round(2).to_dict())
print('lorbin184 %binned', A.loc['dastool-lorbin_184',pc].round(2).to_dict())
