"""Assign eukaryotic IBP contigs to a lineage via MarFERReT best hits of the contig's MetaEuk proteins (majority at class/phylum level).
-> ../ibp_euk_taxa.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
G = C_DIR + 'grex/'
f = pd.read_csv(G + 'euk_ibp_marferret.tsv', sep='\t'); f['contig'] = f.protein_id.str.split('|').str[1]
mt = pd.read_csv(G + 'mft_tax.tsv', sep='\t', header=None, names=['mft', 'taxid']).set_index('mft').taxid
nodes = pd.read_csv(G + 'nodes_min.tsv', sep='\t', header=None, names=['id', 'parent', 'rank']).set_index('id')
names = pd.read_csv(G + 'names_sci.tsv', sep='\t', header=None, names=['id', 'name']).set_index('id').name
par, rank = nodes.parent.to_dict(), nodes['rank'].to_dict()
def lin(t):
    out = {}
    while t in par and t != 1:
        r = rank.get(t)
        if r in ('phylum', 'class', 'order', 'genus', 'clade') and r not in out: out[r] = names.get(t, '')
        if names.get(t) in ('Bacillariophyta', 'Haptophyta', 'Dinophyceae', 'Cryptophyceae', 'Chlorophyta', 'Ciliophora', 'Metazoa', 'Fungi', 'Labyrinthulomycetes', 'Stramenopiles'):
            out.setdefault('group', names.get(t))
        t = par[t]
    return out
f['taxid'] = f.best_hit.map(mt)
L = f.taxid.dropna().astype(int).drop_duplicates().map(lambda t: lin(t))
lk = {t: l for t, l in zip(f.taxid.dropna().astype(int).drop_duplicates(), L)}
f['group'] = f.taxid.map(lambda t: lk.get(int(t), {}).get('group', '') if t == t else '')
f['class'] = f.taxid.map(lambda t: lk.get(int(t), {}).get('class', '') if t == t else '')
f['genus'] = f.taxid.map(lambda t: lk.get(int(t), {}).get('genus', '') if t == t else '')
f = f[f.pident >= 40]
maj = f.groupby('contig').agg(n=('group', 'size'), group=('group', lambda x: x.value_counts().index[0] if len(x) else ''),
                              cls=('class', lambda x: x.value_counts().index[0] if len(x) else ''), mf_genus=('genus', lambda x: x.value_counts().index[0] if len(x) else ''))
I = pd.read_csv(C_DIR + 'ibp_context.tsv', sep='\t', keep_default_na=False)
E = I[I.domain.isin(['eukaryote', 'unknown'])].join(maj, on='contig')
E[['gene_id', 'contig', 'domain', 'plen', 'cpm_SI', 'cpm_SW', 'ice_log2', 'n', 'group', 'cls', 'mf_genus']].to_csv(C_DIR + 'ibp_euk_taxa.tsv', sep='\t', index=False)
S = E.fillna({'group': 'no MarFERReT hit', 'cls': ''}).groupby(['group', 'cls']).agg(n=('gene_id', 'size'), cpm_SI=('cpm_SI', 'sum'), cpm_SW=('cpm_SW', 'sum')).sort_values('cpm_SI', ascending=False)
print(S.round(2).head(20).to_string())
print(E.sort_values('cpm_SI', ascending=False).head(10)[['gene_id', 'domain', 'plen', 'cpm_SI', 'cpm_SW', 'group', 'cls', 'mf_genus']].round(2).to_string())
