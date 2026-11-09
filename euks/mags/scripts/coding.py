"""Coding fraction and BUSCO completeness for the EE MAGs, merged into out/ee_taxonomy.tsv (user request 2026-11-08).
coding_pct  = bp covered by MetaEuk exons (union of intervals, both strands) / MAG bp. Source: danaSeq mag_marine
              eukaryotic/metaeuk/metaeuk_headers.tsv (fetched by fetch_grex.sh). MetaEuk is homology-based (reference
              proteins), so the value is a lower bound and is low for lineages without close references.
genes_per_mb = MetaEuk predictions per Mbp; contigs_with_genes_pct = share of MAG bp on contigs with >= 1 prediction.
completeness = BUSCO 6.1.0, eukaryota_odb12, miniprot (out/busco/busco_eukaryota_odb12.tsv): 'C% [D%]'.
busco_dup_nt_id = median blastn identity between duplicated BUSCO copies (out/busco/dup_identity.tsv, if present):
                  >= 98% haplotypes, ~90-97% strains, < 90% different species.
genome_status (user-approved relabel rules, 2026-11-08):
  'doubtful, not a genome'        BUSCO C < 1% AND coding_pct < 2% AND the MAG is not a metazoan call (protein kingdom vote
                                  Metazoa >= 0.5, or 'Metazoa'/a metazoan taxon in the call: copepods have few BUSCO/MetaEuk hits
                                  because their genomes are large and their reference proteins scarce, not because they are not genomes).
  'multi-genome (>= 2 species)'   BUSCO D >= 10% AND median duplicate-copy identity < 90%;
  'multi-genome (>= 2 strains)'   BUSCO D >= 10% AND median duplicate-copy identity 90-98% (>= 98% = haplotypes of one genome).
  ''                              otherwise.
Retired ids (out/id_registry.tsv 'retired:...') are dropped from ee_taxonomy.tsv."""
import re
import numpy as np, pandas as pd
O = '/data/dev/rnaseq-exploration/euks/mags/out/'

H = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t')
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts'); from eukshared import load_ids
RET = {v for k, v in load_ids().items() if k.startswith('retired:')}
H = H[H.mag.str.startswith('EE') & ~H.mag.isin(RET)]
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
R['busco_C'] = B.C.reindex(R.index); R['busco_D'] = B.D.reindex(R.index)
try:
    DI = pd.read_csv(O + 'busco/dup_identity.tsv', sep='\t', index_col=0)
    R['busco_dup_nt_id'] = DI.median_nt_id.reindex(R.index)
except FileNotFoundError:
    pass

T = pd.read_csv(O + 'ee_taxonomy.tsv', sep='\t')
T = T[~T.mag.isin(RET)]
T = T.drop(columns=[c for c in list(R.columns) + ['genome_status'] if c in T.columns]).merge(R, left_on='mag', right_index=True, how='left')
METAZOA = r'Metazoa|Arthropoda|Copepoda|Calanidae|Calanus|Cnidaria|Mollusca|Annelida|Chordata|Echinodermata|Nematoda|Rotifera'
def kingdom_metazoa(s):
    m = re.match(r'Metazoa ([\d.]+)', str(s)); return bool(m and float(m.group(1)) >= 0.5)
def status(x):
    meta = kingdom_metazoa(x.protein_kingdom) or bool(re.search(METAZOA, str(x.call)))
    if x.busco_C < 1 and x.coding_pct < 2 and not meta: return 'doubtful, not a genome'
    if x.busco_D >= 10 and x.busco_dup_nt_id == x.busco_dup_nt_id and x.busco_dup_nt_id < 98:
        return 'multi-genome (≥ 2 species)' if x.busco_dup_nt_id < 90 else 'multi-genome (≥ 2 strains)'
    return ''
T['genome_status'] = T.apply(status, axis=1)
T.to_csv(O + 'ee_taxonomy.tsv', sep='\t', index=False)
print(T[['mag', 'Mbp', 'call', 'confidence'] + list(R.columns) + ['genome_status']].to_string(index=False))
