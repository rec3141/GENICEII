"""Relative DNA abundance of Cycloclasticus MAGs (share of mapped bases, nanopore) and RNA:DNA ratios for matched samples."""
import sys
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
C = '/data/dev/rnaseq-exploration/cyclo/'
bd = pd.read_csv('/data/dev/rnaseq-exploration/dna/bin_depths.tsv', sep='\t', index_col=0)
cols = [c for c in bd.columns if c.endswith('.sorted.bam')]
bases = bd[cols].mul(bd['length'], axis=0)
rel = bases / bases.sum() * 100
rel.columns = [c.replace('.sorted.bam', '') for c in cols]
ns = nanopore_samples()
ns = ns[ns.run_barcode.isin(rel.columns)]
CY = ['dastool-lorbin_184', 'dastool-metabat_9342', 'dastool-metabat_5432', 'dastool-comebin_3787_sub', 'dastool-lorbin_722_sub']
dna = rel.loc[CY, ns.run_barcode].T; dna.index = ns.sample_name.values; dna['season'] = np.where(ns.log.values == 2024, '2024-25', '2025-26')
dna['run_barcode'] = ns.run_barcode.values
print('unbinned share median', rel.loc['unbinned', ns.run_barcode].median().round(1) if 'unbinned' in rel.index else 'n/a')
# RNA share of mapped pairs
bc = table('bin_counts'); rshare = bc.loc[CY] / bc.sum() * 100
m = meta().loc[bc.columns]
out = []
for s, r in dna.iterrows():
    hit = m.index[(m.sample_name == s) & (m.season == r.season)]
    for b in CY:
        out.append(dict(sample=s, season=r.season, run_barcode=r.run_barcode, bin=b, dna_pct=r[b],
                        rna_pct=float(rshare.loc[b, hit].mean()) if len(hit) else np.nan, rna_sample=','.join(hit)))
O = pd.DataFrame(out)
O['rna_dna'] = O.rna_pct / O.dna_pct.where(O.dna_pct > 0.01)
O.to_csv(C + 'dna_rna_cyclo.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
w = O.pivot_table(index=['season', 'sample'], columns='bin', values='dna_pct').round(3)
print(w.to_string())
print(O.dropna(subset=['rna_pct'])[O.bin.isin(CY[:3])].round(3).to_string(index=False))
