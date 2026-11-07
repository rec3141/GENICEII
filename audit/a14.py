import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
X=pd.DataFrame(C['dna'])
x=X[(X.season=='2025-26')&X['sample'].str.contains('PRE')]
p=x.pivot_table(index=['sample','run_barcode'],columns='bin',values='dna_pct')
p['cyclo']=p.sum(axis=1); print(p.round(3).to_string())
print(p.groupby(p.index.get_level_values(0).str.contains('EXP')).agg(['min','max','mean']).round(2).T.to_string())
# RNA:DNA 9342 2024-25
y=X[(X.bin=='dastool-metabat_9342')&(X.season=='2024-25')&X.rna_dna.notna()]
print(y[['sample','run_barcode','dna_pct','rna_pct','rna_dna']].sort_values('sample').to_string())
