"""Viral transcription per library from CDS gene counts on confident viral contigs (geNomad score >= 0.95, >= 2 hallmarks,
CheckV viral > host genes; CheckV proviruses excluded here and handled separately). CPM = per million CDS read pairs (as hc.py).
Out: out/vrna_cpm.tsv (vOTU-level later), out/vrna_classes.tsv (lytic/lysogeny marker classes x samples), out/vrna_contig_cpm.pkl"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
V = '/data/dev/rnaseq-exploration/viruses/'
vc = pd.read_csv(V + 'out/viral_confident.tsv', sep='\t')
lyt = vc[vc.provirus != 'Yes']
X, G = genes(); S = samples()
cds = G['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
sel = cds & G.contig.isin(set(lyt.seq_name)).values
idx = np.where(sel)[0]
Gv = G.iloc[idx][['gene_id', 'contig', 'product']].reset_index(drop=True)
C = (X[idx].toarray() / tot * 1e6)
pd.to_pickle({'genes': Gv, 'cpm': C.astype(np.float32), 'samples': S}, V + 'out/vrna_gene_cpm.pkl')
CL = [('major capsid', r'major capsid|capsid protein'), ('terminase', r'terminase'), ('portal', r'portal'),
      ('tail', r'tail'), ('lysis (holin/endolysin/spanin)', r'holin|lysozyme|endolysin|spanin|lysis protein|peptidoglycan hydrolase|N-acetylmuramoyl'),
      ('integrase/excisionase', r'integrase|excisionase|recombinase'), ('repressor (cI/Cro)', r'repressor|\bCro\b|\bcI\b'),
      ('DNA replication', r'DNA polymerase|primase|helicase|ribonucleotide reductase|thymidylate'),
      ('auxiliary metabolic (psbA/D, phoH, mazG, cobS)', r'photosystem|\bpsb|phoH|mazG|cobS|high light inducible|PhoH')]
import re
cls = Gv['product'].fillna('').map(lambda p: next((k for k, rx in CL if re.search(rx, p, re.I)), 'other'))
R = pd.DataFrame(C, columns=S).groupby(cls.values).sum()
R.loc['all viral CDS'] = C.sum(0)
R.to_csv(V + 'out/vrna_classes.tsv', sep='\t')
m = meta().loc[S]
pd.set_option('display.width', 250)
x = R.loc['all viral CDS']
g = pd.DataFrame({'cpm': x, 'season': m.season, 'group': m.group, 'medium': m.medium, 'phase': m.phase})
print(g.groupby(['season', 'group', 'medium', 'phase']).cpm.median().round(0).to_string())
print(len(idx), 'viral CDS;', lyt.shape[0], 'contigs')
