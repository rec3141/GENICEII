"""Host taxonomy per HC gene: GTDB of DAS Tool bin, else Kaiju contig lineage -> ../hc_host.tsv."""
import sys; sys.path.insert(0, '.')
from hcdef import *
D = '/data/dev/rnaseq-exploration/data/'
h = pd.read_csv('../hc_genes.tsv', sep='\t', keep_default_na=False)
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
kj = pd.read_csv(D + 'kaiju_contigs.tsv', sep='\t', usecols=['contig_id', 'lineage'])
kj = kj[kj.contig_id.isin(set(h.contig))].set_index('contig_id')['lineage']
R = ['domain', 'phylum', 'class', 'order', 'family', 'genus']
rows = []
for b, c in zip(h.bin, h.contig):
    t = gt.get(b.replace('dastool-', ''), '') if b else ''
    if t:
        d = dict(x.split('__', 1) for x in t.split(';') if '__' in x)
        r = [d.get(k[0], '') for k in R]; src = 'GTDB(MAG)'
    else:
        l = str(kj.get(c, ''))
        parts = [x.strip() for x in l.split(';')] if l and l not in ('Unclassified', 'nan') and not l[:1].isupper() or ';' in l else []
        parts = [x if x != 'NA' else '' for x in parts][:6]
        r = parts + [''] * (6 - len(parts)); src = 'Kaiju(contig)' if parts else ''
    rows.append(r + [src])
T = pd.DataFrame(rows, columns=R + ['tax_source'])
T.insert(0, 'gene_id', h.gene_id.values)
T.to_csv('../hc_host.tsv', sep='\t', index=False)
print(T.tax_source.value_counts().to_string()); print(T.family.replace('', '?').value_counts().head(15).to_string())
