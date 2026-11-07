"""Cycloclasticus hydrocarbon-gene inventory, operon neighbourhoods and per-gene expression time courses."""
import sys, re, json
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
C = '/data/dev/rnaseq-exploration/cyclo/'
BINS = ['dastool-lorbin_184', 'dastool-metabat_9342', 'dastool-metabat_5432', 'dastool-lorbin_722_sub', 'dastool-comebin_3787_sub']
# category: (regex on Bakta product, KO list)
CAT = [
    ('RHO alpha (PAH/aromatic dioxygenase large)', r'(ring.hydroxylating dioxygenase.*(alpha|large)|PAH dioxygenase.*large|naphthalene 1,2-dioxygenase(?! .*(beta|small))|benzene 1,2-dioxygenase subunit alpha|benzoate/toluate 1,2-dioxygenase alpha|phenylpropionate.*dioxygenase.*alpha|anthranilate 1,2-dioxygenase large|aromatic dioxygenase large|carbazole.*dioxygenase|terephthalate.*alpha|biphenyl dioxygenase.*alpha|Rieske non-heme iron oxygenase|ring.hydroxylating dioxygenase$|aromatic.ring.hydroxylating dioxygenase$)', ['K14579', 'K11943', 'K08689', 'K03268', 'K05549']),
    ('RHO beta (small subunit)', r'(ring.hydroxylating dioxygenase.*(beta|small)|PAH dioxygenase.*small|benzene 1,2-dioxygenase subunit beta|benzoate.*1,2-dioxygenase.*(beta|small)|ethylbenzene dioxygenase beta|aromatic oxygenase small|phenylpropionate.*beta)', ['K14580', 'K11944', 'K15750', 'K16268', 'K05550']),
    ('RHO electron transfer (ferredoxin/reductase)', r'(PAH dioxygenase component ferredoxin|dioxygenase ferredoxin|ring-hydroxylating dioxygenases|naphthalenesulfonate dioxygenase ferredoxin)', ['K14578', 'K14581', 'K18087', 'K18088']),
    ('dihydrodiol dehydrogenase', r'dihydrodiol dehydrogenase|dihydro.*diol dehydrogenase', ['K14582', 'K08690']),
    ('extradiol ring cleavage', r'(biphenyl-2,3-diol 1,2-dioxygenase|dihydroxybiphenyl 1,2-dioxygenase|dihydroxynaphthalene dioxygenase|catechol 2,3-dioxygenase|extradiol|LigB family|dihydroxy-p-cumate|2,3-dihydroxybenzoate 3,4-dioxygenase)', ['K00446', 'K07104', 'K00462', 'K14583', 'K10621']),
    ('intradiol / ortho cleavage', r'(catechol 1,2-dioxygenase|protocatechuate 3,4-dioxygenase|muconate cycloisomerase|muconolactone|maleylacetate|3-oxoadipate)', ['K03381', 'K00448', 'K00449', 'K00217']),
    ('meta-cleavage lower pathway', r'(2-hydroxymuconic semialdehyde|2-hydroxymuconate|hydroxymuconate|oxopent|4-oxalocrotonate|hydratase-aldolase|MetA-pathway|Phenol degradation protein meta|salicylaldehyde dehydrogenase)', ['K10216', 'K10217', 'K01617', 'K18364', 'K14584', 'K14585']),
    ('salicylate / gentisate', r'(salicylate hydroxylase|salicylate 1-monooxygenase|gentisate)', ['K00480', 'K00450']),
    ('phenol/toluene monooxygenase (multicomponent)', r'(phenol hydroxylase|phenol 2-monooxygenase|DmpK|DmpM|DmpN|DmpP|toluene.*monooxygenase|xylene monooxygenase|aromatic/alkene monooxygenase)', ['K16242', 'K16243', 'K16244', 'K16245', 'K16246', 'K16249', 'K15760', 'K15761', 'K15762', 'K15763', 'K15764', 'K15765']),
    ('alkane hydroxylation', r'(alkane 1-monooxygenase|alkane hydroxylase|long-chain alkane monooxygenase|flavin-binding monooxygenase AlmA|cytochrome P450|rubredoxin)', ['K00496', 'K05297']),
    ('Cu-membrane monooxygenase (pmo/amo-like)', r'methane monooxygenase|ammonia monooxygenase', ['K10944', 'K10945', 'K10946']),
    ('chemotaxis', r'chemotaxis|methyl-accepting', []),
]
X, g = genes(); S = samples(); m = meta().loc[S]
co = pd.read_csv(C + 'cyclo_genes_coords.tsv', sep='\t')[['gene_id', 'start', 'end', 'strand']]
cg = pd.read_csv(C + 'cyclo_contigs.tsv', sep='\t', header=None, names=['contig', 'bin'])
G = g.reset_index(drop=True)
G = G[G.contig.isin(cg.contig)].copy()
G['row'] = G.index
G = G.merge(co, on='gene_id', how='left')
G['cbin'] = G.contig.map(cg.set_index('contig')['bin'])
def cat(r):
    # KOs from mag_marine eggNOG/kofam are keyed to an older Bakta run (IDs shifted) -> Bakta product only
    for name, rx, kl in CAT:
        if re.search(rx, r['product'], re.I): return name
    return ''
G['hc'] = G.apply(cat, axis=1)
# expression
cdsrows = np.where(g['type'].values == 'CDS')[0]
tot = np.asarray(X[cdsrows].sum(0)).ravel()
E = X[G.row.values].toarray().astype(float)           # pairs
cpm = E / tot * 1e6
# within-genome TPM (per bin, CDS only)
wtpm = np.zeros_like(E)
for b in G.cbin.unique():
    ix = np.where((G.cbin == b).values & (G['type'] == 'CDS').values)[0]
    rpk = E[ix] / (G.length.values[ix, None] / 1000)
    s = rpk.sum(0); s[s == 0] = 1
    wtpm[ix] = rpk / s * 1e6
G.drop(columns=['row']).to_csv(C + 'cyclo_genes_annot.tsv', sep='\t', index=False)
np.save(C + 'cyclo_cpm.npy', cpm); np.save(C + 'cyclo_wtpm.npy', wtpm); np.save(C + 'cyclo_pairs.npy', E)
print(G[G.hc != ''].groupby(['cbin', 'hc']).size().unstack(0).fillna(0).astype(int).to_string())
