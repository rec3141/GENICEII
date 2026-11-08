"""Coding fraction and BUSCO completeness for the EE MAGs, merged into out/ee_taxonomy.tsv (user request 2026-11-08).
coding_pct  = bp covered by MetaEuk exons (union of intervals, both strands) / MAG bp. Source: danaSeq mag_marine
              eukaryotic/metaeuk/metaeuk_headers.tsv (fetched by fetch_grex.sh). MetaEuk is homology-based (reference
              proteins), so the value is a lower bound and is low for lineages without close references.
genes_per_mb = MetaEuk predictions per Mbp; contigs_with_genes_pct = share of MAG bp on contigs with >= 1 prediction.
completeness = BUSCO 6.1.0, eukaryota_odb12, miniprot (out/busco/busco_eukaryota_odb12.tsv): 'C% [D%]'.
busco_dup_nt_id = median blastn identity between duplicated BUSCO copies (out/busco/dup_identity.tsv, if present):
                  >= 98% haplotypes, ~90-97% strains, < 90% different species."""
import re
import numpy as np, pandas as pd
O = '/data/dev/rnaseq-exploration/euks/mags/out/'

H = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')
H = H[H.mag.str.startswith('EE')]
mag_of = dict(zip(H.contig, H.mag))
size = H.groupby('mag').length.sum()

iv, ngenes = {}, {}
EX = re.compile(r'(\d+)\[\d+\]:(\d+)\[\d+\]:\d+\[\d+\]')
for line in open(O + 'identify/grex/metaeuk_headers.tsv'):
    h = line.rstrip('\n').split('\t')[-1].split('|')
    c = h[1]
    if c not in mag_of: continue
    ngenes[c] = ngenes.get(c, 0) + 1
    for a, b in EX.findall('|'.join(h[8:])):
        a, b = int(a), int(b)
        iv.setdefault(c, []).append((min(a, b), max(a, b)))

def union_len(x):
    x.sort(); tot = 0; s, e = x[0]
    for a, b in x[1:]:
        if a > e + 1: tot += e - s + 1; s, e = a, b
        else: e = max(e, b)
    return tot + e - s + 1

cod = pd.Series({c: union_len(v) for c, v in iv.items()})
C = H.set_index('contig')
C['coding'] = cod.reindex(C.index).fillna(0)
C['genes'] = pd.Series(ngenes).reindex(C.index).fillna(0)
g = C.groupby('mag')
R = pd.DataFrame({
    'coding_pct': (g.coding.sum() / size * 100).round(1),
    'genes_per_mb': (g.genes.sum() / size * 1e6).round(0),
    'contigs_with_genes_pct': (C.length.where(C.genes > 0, 0).groupby(C.mag).sum() / size * 100).round(1)})

B = pd.read_csv(O + 'busco/busco_eukaryota_odb12.tsv', sep='\t', index_col=0)
R['completeness'] = [f'{B.C[m]:.1f}% [D {B.D[m]:.1f}%]' if m in B.index else '' for m in R.index]
try:
    DI = pd.read_csv(O + 'busco/dup_identity.tsv', sep='\t', index_col=0)
    R['busco_dup_nt_id'] = DI.median_nt_id.reindex(R.index)
except FileNotFoundError:
    pass

T = pd.read_csv(O + 'ee_taxonomy.tsv', sep='\t')
T = T.drop(columns=[c for c in R.columns if c in T.columns]).merge(R, left_on='mag', right_index=True, how='left')
T.to_csv(O + 'ee_taxonomy.tsv', sep='\t', index=False)
print(T[['mag', 'Mbp', 'call', 'confidence'] + list(R.columns)].to_string(index=False))
