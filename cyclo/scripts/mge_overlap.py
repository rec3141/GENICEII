"""Flag Cycloclasticus HC genes on plasmid-like contigs, genomic islands, integrons, or next to transposases/integrases."""
import pandas as pd, numpy as np
C = '/data/dev/rnaseq-exploration/cyclo/'
G = pd.read_csv(C + 'cyclo_genes_annot.tsv', sep='\t', keep_default_na=False).sort_values(['contig', 'start']).reset_index(drop=True)
sc = pd.read_csv(C + 'cyclo_genomad_scores.tsv', sep='\t').set_index('seq_name')
mge = pd.read_csv(C + 'cyclo_mge.tsv', sep='\t', header=None, dtype=str, names=list(range(20)))
isl = mge[mge[0].str.contains('islandpath')][[3, 4, 5]].rename(columns={3: 'contig', 4: 's', 5: 'e'}).astype({'s': int, 'e': int})
plas = mge[mge[0].str.contains('plasmid_summary')].set_index(2)
integ = mge[mge[0].str.contains('integrons')][[3, 4, 5, 6, 10, 11]].rename(columns={3: 'id', 4: 'contig', 6: 's', 7: 'e'})
integ = mge[mge[0].str.contains('integrons')].rename(columns={2: 'iid', 3: 'contig', 5: 's', 6: 'e', 10: 'elt', 11: 'ann'})
integ = integ[integ.contig.str.startswith('contig')]
integ[['s', 'e']] = integ[['s', 'e']].astype(int)
ispan = integ.groupby(['contig', 'iid']).agg(s=('s', 'min'), e=('e', 'max'), n=('elt', 'size'), intI=('ann', lambda a: (a == 'intI').any())).reset_index()
mob = G['product'].str.contains(r'transposase|integrase|IS[0-9]+|recombinase|resolvase|insertion element|conjugal|relaxase|MobA|TraG|phage', case=False, regex=True)
G['mob'] = mob
G['plasmid_score'] = G.contig.map(sc['plasmid_score'])
G['virus_score'] = G.contig.map(sc['virus_score'])
G['genomad_plasmid'] = G.contig.isin(plas.index)
G['island'] = [((isl.contig == c) & (isl.s <= e) & (isl.e >= s)).any() for c, s, e in zip(G.contig, G.start, G.end)]
G['integron'] = [((ispan.contig == c) & (ispan.s - 2000 <= e) & (ispan.e + 2000 >= s)).any() for c, s, e in zip(G.contig, G.start, G.end)]
G['ord'] = G.groupby('contig').cumcount()
near = np.zeros(len(G), bool)
for c, d in G.groupby('contig'):
    mo = d.ord[d.mob].values
    if len(mo): near[d.index] = [np.abs(mo - o).min() <= 5 for o in d.ord]
G['mob_within5'] = near
G.to_csv(C + 'cyclo_genes_mge.tsv', sep='\t', index=False)
hc = G[(G.hc != '') & (G.hc != 'chemotaxis')]
print(hc.groupby('cbin')[['genomad_plasmid', 'island', 'integron', 'mob_within5']].sum().join(hc.groupby('cbin').size().rename('n_hc')).to_string())
pd.set_option('display.width', 250, 'display.max_colwidth', 60)
print(hc[hc.genomad_plasmid | hc.island | hc.integron][['cbin', 'contig', 'start', 'product', 'plasmid_score', 'genomad_plasmid', 'island', 'integron', 'mob_within5']].to_string())
print('\nplasmid contigs (geNomad) in Cyclo bins:')
for c, r in plas.iterrows():
    if r[1] == 'unbinned_kaiju_cyclo': continue
    d = G[G.contig == c]
    print(r[1], c, r[3], 'len', 'genes', r[5], 'score', r[7], 'conj', str(r[12])[:80], '| HC genes:', '; '.join(d['product'][(d.hc != '') & (d.hc != 'chemotaxis')].str.slice(0, 40)))
print(ispan[ispan.contig.isin(G.contig[G.cbin != 'unbinned_kaiju_cyclo'])].groupby('contig').agg(n=('n', 'sum'), intI=('intI', 'any')).sort_values('n', ascending=False).head(10))
