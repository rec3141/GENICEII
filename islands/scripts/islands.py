"""Genomic islands (IslandPath-DIMOB, mag_marine) annotated with mobility markers, gene functions, hydrocarbon genes and expression.
Inputs: islands/grex/{genomic_islands.tsv, island_contig_features.tsv.gz, island_gc.tsv} (Grex islands.sh), hgt/mge (geNomad, integrons),
data/genemap.pkl (corrected KO, bins), canthyd/canthyd_hits.tsv, hgt/hc_genes.tsv (hydrocarbon gene set), gene counts + corrected metadata.
Out: islands/out/islands.tsv (one row per island), island_genes.tsv (genes inside islands), function_enrichment.tsv, island_expr.tsv."""
import sys, re, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
I = '/data/dev/rnaseq-exploration/islands/'; H = '/data/dev/rnaseq-exploration/hgt/'; O = I + 'out/'
os.makedirs(O, exist_ok=True)

GI = pd.read_csv(I + 'grex/genomic_islands.tsv', sep='\t')
GI['len'] = GI.end - GI.start + 1
gc = pd.read_csv(I + 'grex/island_gc.tsv', sep='\t').set_index('island_id')
GI = GI.join(gc[['contig_len', 'island_gc', 'contig_gc_outside']], on='island_id')
GI['dGC'] = (GI.island_gc - GI.contig_gc_outside) * 100
F = pd.read_csv(I + 'grex/island_contig_features.tsv.gz', sep='\t', header=None, comment='#',
                names=['contig', 'type', 'start', 'stop', 'strand', 'gene_id', 'gene', 'product', 'xref'], keep_default_na=False)

# genes inside islands (feature midpoint inside the island)
F['mid'] = (F.start + F.stop) / 2
parts = []
for c, g in GI.groupby('contig'):
    f = F[F.contig == c]
    for _, r in g.iterrows():
        x = f[(f.mid >= r.start) & (f.mid <= r.end)].copy(); x['island_id'] = r.island_id; parts.append(x)
G = pd.concat(parts)
# tRNA at or near a boundary (classic island insertion site): any tRNA/tmRNA within 2 kb of either end
trna = F[F.type.isin(['tRNA', 'tmRNA'])]
def near_trna(r):
    t = trna[trna.contig == r.contig]
    return bool(((t.stop >= r.start - 2000) & (t.start <= r.start + 2000)).any() or ((t.stop >= r.end - 2000) & (t.start <= r.end + 2000)).any())
GI['tRNA_edge'] = GI.apply(near_trna, axis=1)

# function classes from Bakta products (coarse, mobility first)
CLS = [('integrase', r'\bintegrase|site-specific recombinase|tyrosine recombinase|serine recombinase|\bxerC|\bxerD'),
       ('transposase/IS', r'transposase|insertion element|\bIS\d|IS[0-9]+ family'),
       ('phage', r'phage|capsid|tail fiber|tail protein|terminase|portal protein|baseplate|holin|prophage'),
       ('conjugation/T4SS', r'conjugal|conjugative|\btra[A-Z]\b|\btrb[A-Z]\b|\bvir[BD]\d|type IV secretion|relaxase|mobilization|\bmob[A-Z]\b'),
       ('restriction-modification', r'restriction|methyltransferase.*(DNA|adenine|cytosine)|DNA methylase|modification methylase'),
       ('CRISPR/defence', r'\bcas\d|CRISPR|abortive infection|\bAbi[A-Z]|toxin-antitoxin|antitoxin|\btoxin\b|\bRelE\b|\bHicB|\bVapC|\bMazF|\bHigB'),
       ('metal resistance', r'mercur|\bmer[A-Z]\b|arsen|\bars[A-Z]\b|copper|\bcop[A-Z]\b|cadmium|\bczc|cobalt-zinc|heavy metal|chromate|tellurite'),
       ('antibiotic resistance', r'beta-lactamase|aminoglycoside|tetracycline|chloramphenicol|multidrug|efflux RND|bleomycin'),
       ('EPS/cell surface', r'glycosyltransferase|polysaccharide|lipopolysaccharide|capsul|O-antigen|exopolysaccharide|wzx|wzy|\bwza|flippase'),
       ('secretion/adhesion', r'type [I1]+ secretion|type VI|hemolysin|RTX|adhesin|pilus|pilin|fimbri'),
       ('transport', r'transporter|permease|ABC transporter|TonB|porin|symporter|antiporter|channel'),
       ('regulation', r'transcriptional regulator|response regulator|sensor histidine kinase|sigma factor|anti-sigma'),
       ('hypothetical', r'^hypothetical protein$|DUF\d|domain-containing protein$')]
def cls(p):
    for k, rx in CLS:
        if re.search(rx, p, re.I): return k
    return 'other enzyme/function'
G['class'] = G['product'].map(cls)
GM = pd.read_pickle(D + 'genemap.pkl').set_index('gene_id')
G = G.join(GM[['KO', 'bin']], on='gene_id')
hc = pd.read_csv(H + 'hc_genes.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
ch = pd.read_csv('/data/dev/rnaseq-exploration/canthyd/canthyd_hits.tsv', sep='\t').sort_values('score', ascending=False).drop_duplicates('gene_id').set_index('gene_id')
G['hc_category'] = G.gene_id.map(hc.category).fillna(''); G['hc_tier'] = G.gene_id.map(hc.tier).fillna('')
G['canthyd'] = G.gene_id.map(ch.hmm).fillna('')
G['hydrocarbon'] = (G.hc_tier == 'confident') | (G.canthyd != '')
G.loc[G.hydrocarbon, 'class'] = 'hydrocarbon degradation'

# other MGE evidence on the island
pl = pd.read_csv(H + 'mge/plasmid_summary.tsv', sep='\t').set_index('seq_name')
pv = pd.read_csv(H + 'mge/provirus.tsv', sep='\t')
ig = pd.read_csv(H + 'mge/integrons.tsv', sep='\t')
ig = ig.groupby(['ID_integron', 'ID_replicon']).agg(start=('pos_beg', 'min'), end=('pos_end', 'max')).reset_index()
def overlap(df, c, a, b, sc='source_seq'):
    x = df[df[sc] == c]
    return bool(((x.start <= b) & (x.end >= a)).any())
GI['plasmid_score'] = GI.contig.map(pl.plasmid_score)
GI['provirus'] = [overlap(pv, r.contig, r.start, r.end) for r in GI.itertuples()]
GI['integron'] = [overlap(ig, r.contig, r.start, r.end, 'ID_replicon') for r in GI.itertuples()]
c2b = pd.read_csv(D + 'contig2bin.tsv', sep='\t', header=None, names=['contig', 'bin']).drop_duplicates('contig').set_index('contig').bin
GI['bin'] = GI.contig.map(c2b).fillna('unbinned')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
def lab(b):
    if b == 'unbinned': return 'unbinned'
    if plastid_label(b): return plastid_label(b)
    t = [x.split('__', 1)[1] for x in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in x and x.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or 'unclassified'
GI['taxon'] = GI.bin.map(lab)
GI['family'] = [next((x[3:] for x in str(gt.get(b.replace('dastool-', ''), '')).split(';') if x.startswith('f__')), '') for b in GI.bin]
cnt = G[G.type == 'cds'].groupby(['island_id', 'class']).size().unstack(fill_value=0)
GI = GI.join(cnt, on='island_id').fillna({k: 0 for k in cnt.columns})
GI['n_cds'] = GI.island_id.map(G[G.type == 'cds'].groupby('island_id').size()).fillna(0).astype(int)
GI['hc_genes'] = GI.island_id.map(G[G.hydrocarbon].groupby('island_id').size()).fillna(0).astype(int)
GI['hc_categories'] = GI.island_id.map(G[G.hydrocarbon].groupby('island_id').apply(lambda x: ', '.join(sorted(set(x.hc_category.replace('', np.nan).dropna()) | set(x.canthyd.replace('', np.nan).dropna()))))).fillna('')
GI['mobility_score'] = (GI.get('integrase', 0) > 0).astype(int) + (GI.get('transposase/IS', 0) > 0).astype(int) + GI.tRNA_edge.astype(int) + \
    (GI.get('phage', 0) >= 3).astype(int) + (GI.get('conjugation/T4SS', 0) >= 2).astype(int) + GI.integron.astype(int) + GI.provirus.astype(int) + (GI.dGC.abs() >= 3).astype(int)

# expression: island CPM (sum of its CDS) per library, corrected metadata
X, GM2 = genes(); S = samples(); m = meta().loc[S]
cds = GM2['type'].values == 'CDS'; tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
row = pd.Series(np.arange(len(GM2)), index=GM2.gene_id.values)
G['row'] = G.gene_id.map(row)
gr = G.dropna(subset=['row'])
gene_cpm = pd.DataFrame(X[gr.row.astype(int).values].toarray() / tot * 1e6, index=gr.gene_id.values, columns=S)
E = gene_cpm.groupby(gr.island_id.values).sum()
E.to_csv(O + 'island_expr.tsv', sep='\t')
oil = m.phase.eq('oiled'); ctl = m.phase.eq('control')
for s in ('2024-25', '2025-26'):
    for med in ('SW',):
        a = m.index[oil & (m.season == s) & (m.medium == med)]; b = m.index[ctl & (m.season == s) & (m.medium == med)]
        GI[f'cpm_oil_{s}'] = GI.island_id.map(E[a].median(1)); GI[f'cpm_ctl_{s}'] = GI.island_id.map(E[b].median(1))
        GI[f'log2_oil_ctl_{s}'] = np.log2((GI[f'cpm_oil_{s}'] + 0.5) / (GI[f'cpm_ctl_{s}'] + 0.5))
        ga = gene_cpm[a].median(1); gb = gene_cpm[b].median(1)
        G[f'cpm_oil_{s}'] = G.gene_id.map(ga); G[f'cpm_ctl_{s}'] = G.gene_id.map(gb)
GI['cpm_max'] = GI.island_id.map(E.max(1))

# IslandPath-DIMOB false positives: highly expressed core-gene clusters (ribosomal proteins, RNA polymerase, EF-Tu/G) have atypical
# codon usage and get called as islands. Flag islands with >= 2 such genes; they track host abundance, not mobile DNA.
CORE = r'ribosomal protein|ribosomal subunit protein|RNA polymerase subunit (alpha|beta)|elongation factor (Tu|G|Ts)|DNA gyrase|translation initiation factor|ATP synthase|chaperonin|\bDnaK\b|Clp protease'
GI['core_genes'] = GI.island_id.map(G[(G.type == 'cds') & G['product'].str.contains(CORE, case=False, regex=True)].groupby('island_id').size()).fillna(0).astype(int)
GI['core_like'] = GI.core_genes >= 2

# host-normalised expression: island share of its host MAG's CDS transcripts, diesel vs control; and in diesel samples relative to
# the island's share of host genes (activity index > 1 = island transcribed above the host's per-gene average)
import scipy.sparse as sp
bins = GM2['bin'].fillna('').values
ub, bi = np.unique(bins, return_inverse=True)
ind = sp.csr_matrix((cds.astype(float), (bi, np.arange(len(bins)))), shape=(len(ub), len(bins)))
BC = pd.DataFrame((ind @ X).toarray() / tot * 1e6, index=ub, columns=S)   # bin CPM of CDS
ncds_bin = pd.Series(np.asarray(ind.sum(1)).ravel(), index=ub)
hb = GI.set_index('island_id').bin
host = BC.reindex(hb.values); host.index = hb.index
share = E.reindex(hb.index) / host.where(host >= 20)          # need >= 20 host CPM for a usable share
for s in ('2024-25', '2025-26'):
    a = m.index[oil & (m.season == s) & (m.medium == 'SW')]; b = m.index[ctl & (m.season == s) & (m.medium == 'SW')]
    sa, sb = share[a].median(1), share[b].median(1)
    GI[f'host_cpm_oil_{s}'] = GI.island_id.map(host[a].median(1)); GI[f'host_cpm_ctl_{s}'] = GI.island_id.map(host[b].median(1))
    GI[f'log2_share_oil_ctl_{s}'] = GI.island_id.map(np.log2((sa + 1e-4) / (sb + 1e-4)).where(sa.notna() & sb.notna()))
    frac_genes = GI.set_index('island_id').n_cds / ncds_bin.reindex(hb.values).values
    GI[f'activity_oil_{s}'] = GI.island_id.map(sa / frac_genes)
GI.to_csv(O + 'islands.tsv', sep='\t', index=False)
G.drop(columns=['mid', 'row']).to_csv(O + 'island_genes.tsv', sep='\t', index=False)

# function enrichment: class shares inside islands vs the same contigs outside islands (CDS only)
F['class'] = F['product'].map(cls)
F = F[F.type == 'cds']; F['in_island'] = F.gene_id.isin(set(G.gene_id))
F.loc[F.gene_id.isin(hc.index[hc.tier == 'confident']) | F.gene_id.isin(ch.index), 'class'] = 'hydrocarbon degradation'
en = F.groupby(['class', 'in_island']).size().unstack(fill_value=0); en.columns = ['outside', 'inside']
en['share_in'] = en.inside / en.inside.sum(); en['share_out'] = en.outside / en.outside.sum(); en['enrich'] = en.share_in / en.share_out
en.sort_values('enrich', ascending=False).to_csv(O + 'function_enrichment.tsv', sep='\t')
pd.set_option('display.width', 250)
print(len(GI), 'islands;', (GI.bin != 'unbinned').sum(), 'in MAGs;', GI.bin.nunique() - 1, 'MAGs')
print(en.round(3).to_string())
print(GI[GI.hc_genes > 0].sort_values('hc_genes', ascending=False)[['island_id', 'contig', 'len', 'taxon', 'n_cds', 'hc_genes', 'hc_categories', 'mobility_score', 'dGC', 'log2_oil_ctl_2024-25', 'log2_oil_ctl_2025-26', 'cpm_max']].head(40).round(2).to_string())
