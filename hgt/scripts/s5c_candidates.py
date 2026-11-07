"""Strict HGT candidates: >=99% nt identity, different genus AND (harmonized) class or family, low depth correlation."""
import sys; sys.path.insert(0, '.')
from hcdef import *
H = '/data/dev/rnaseq-exploration/hgt/'
a = pd.read_csv(H + 'shared_pairs.tsv', sep='\t')
for c in a.columns:
    if a[c].dtype == object: a[c] = a[c].fillna('')
CLS = {'Flavobacteriia': 'Bacteroidia', 'Sphingobacteriia': 'Bacteroidia', 'Cytophagia': 'Bacteroidia', 'Saprospiria': 'Bacteroidia',
       'Chitinophagia': 'Bacteroidia', 'Betaproteobacteria': 'Gammaproteobacteria', 'Actinomycetes': 'Actinomycetia', 'Actinobacteria': 'Actinomycetia'}
FAM = {'Piscirickettsiaceae': 'Cycloclasticaceae', 'Roseobacteraceae': 'Rhodobacteraceae', 'Paracoccaceae': 'Rhodobacteraceae',
       'Candidatus Pelagibacteraceae': 'Pelagibacteraceae', 'Comamonadaceae': 'Burkholderiaceae', 'Oxalobacteraceae': 'Burkholderiaceae',
       'Phyllobacteriaceae': 'Rhizobiaceae', 'Colwelliaceae': 'Alteromonadaceae', 'Pseudoalteromonadaceae': 'Alteromonadaceae',
       'Shewanellaceae': 'Alteromonadaceae', 'Idiomarinaceae': 'Alteromonadaceae', 'Psychromonadaceae': 'Alteromonadaceae',
       'Oceanospirillaceae': 'Marinomonadaceae', 'Marinobacteraceae': 'Marinobacteraceae', 'Hyphomicrobiaceae': 'Hyphomicrobiaceae'}
for s in 'qt':
    a[f'{s}_cls'] = a[f'{s}_class'].map(lambda x: CLS.get(x, x)); a[f'{s}_fam'] = a[f'{s}_family'].map(lambda x: FAM.get(x, x))
    a[f'{s}_gen'] = a[f'{s}_genus'].str.replace(r'_[A-Z]$', '', regex=True).str.replace('Candidatus ', '')
known = lambda x: (x != '') & ~x.str.contains('unclassified|Incertae', case=False)
gen_diff = known(a.q_gen) & known(a.t_gen) & (a.q_gen != a.t_gen)
a['level'] = np.select([known(a.q_cls) & known(a.t_cls) & (a.q_cls != a.t_cls) & gen_diff,
                        known(a.q_fam) & known(a.t_fam) & (a.q_fam != a.t_fam) & gen_diff, gen_diff], ['cross-class', 'cross-family', 'cross-genus'], '')
a['independent_hosts'] = a.depth_rho < 0.2      # same-MAG contig pairs: median rho 0.60; random pairs: median -0.52
c = a[(a.level != '') & (a.nt_ident >= 99)]
print(c.groupby(['level', 'independent_hosts']).size().unstack(fill_value=0).to_string())
s = c[c.independent_hosts | c.depth_rho.isna()]
s.to_csv(H + 'hgt_candidates.tsv', sep='\t', index=False)
pd.set_option('display.width', 260)
print(s[['level', 'q_category', 'q_genus', 'q_family', 't_genus', 't_family', 'q_tax_source', 't_tax_source', 'nt_ident', 'depth_rho', 'q', 't', 'q_contig', 't_contig']].sort_values(['level', 'q_category']).to_string())
