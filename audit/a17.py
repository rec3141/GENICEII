import sys; sys.path.insert(0,"/data/dev/rnaseq-exploration/audit")
from lib import *
pd.set_option('display.width',300)
O=pd.read_csv(R+'out/mag_oil_response.tsv',sep='\t')
B=pd.read_csv(R+'out/mag_baci.tsv',sep='\t')
for se,me in [('2024-25','SW'),('2024-25','SI'),('2025-26','SW'),('2025-26','SI')]:
    o=O[(O.season==se)&(O.medium==me)&(O.q<0.05)].sort_values('dCLR',ascending=False)
    print(se,me,'n sig',len(o)); print(o.head(8)[['bin','taxon','dCLR','tpm_oiled','tpm_ctl','q']].round(3).to_string(index=False))
print(O[O.taxon.str.contains('Sulfitobacter',na=False)&(O.season=='2024-25')].sort_values('tpm_oiled',ascending=False).head(6)[['medium','bin','taxon','dCLR','tpm_oiled','tpm_ctl','q']].to_string())
print(B.sort_values('interaction',ascending=False).head(10)[['season','medium','bin','genus','interaction','q']].to_string())
print(B[B.genus.str.contains('Sulfurimonas|Cycloclasticus',na=False)][['season','medium','bin','genus','interaction','q']].to_string())
