"""Illumina metagenome k-mer screen -> per-target abundance table + site data (out/mg_kmer.tsv, site key 'mg')."""
import sys, glob, re, json; sys.path.insert(0, '.')
from common import *
S = pd.read_csv('../mgscreen/samples.tsv', sep='\t')
rows = [pd.read_csv(f, sep='\t', header=None, names=['sample', 'target', 'hits', 'tot', 'size', 'abund']) for f in glob.glob('../mgscreen/out/*.tsv')]
A = pd.concat(rows); A['pct'] = 100 * A.hits / A.tot  # % of sampled k-mers that are unique to the target
R = A.pivot(index='sample', columns='target', values='pct')
info = S.set_index('sample')
# label corrections / QC flags (identity check), joined by Genohub id parsed from the fastq path
w = info.fq.str.extract(r'raw/[^/]+/[^/]+/(?:(\d+)-|(plate\d+)-)([A-H]\d\d)_')
info['gid'] = np.where(w[0].notna(), '2599682_plate' + w[0].fillna('') + '_' + w[2], '1343765_' + w[1].fillna('') + '_' + w[2])
C = corrections('metagenome'); info['label_original'] = info.name
fix = info.gid.map(C.sample_name_corrected).fillna(''); info['name'] = info.name.where(fix == '', fix)
info['qc_flag'] = info.gid.map(C.qc_flag).fillna('')
u = info.name.str.upper()
info['season'] = np.where(info.project.str.startswith('1343765'), '2025-26', '2024-25')
info['pool'] = np.where(u.str.contains(r'(^|_)EXP(_|$)'), 'EXP', np.where(u.str.contains(r'(^|_)CTL(_|$)'), 'CTL', 'other'))
info['medium'] = np.where(u.str.contains(r'_SI'), 'SI', np.where(u.str.contains(r'_SW'), 'SW', 'other'))
def day(s):
    m_ = re.search(r'(PRE|DAY)\s*(\d+)', s.upper()); return np.nan if not m_ else (-int(m_[2]) if m_[1] == 'PRE' else int(m_[2]))
info['day'] = info.name.map(day)
R = R.join(info[['name', 'label_original', 'qc_flag', 'gid', 'season', 'pool', 'medium', 'day']])
R.to_csv('../out/mg_kmer.tsv', sep='\t')
T = [c for c in R.columns if c not in ('name', 'label_original', 'qc_flag', 'gid', 'season', 'pool', 'medium', 'day')]
R = R[~R.qc_flag.isin(EXCLUDE_FLAGS)]   # QC-excluded libraries (failed, unidentifiable, contaminated) stay in mg_kmer.tsv but not on the site
LAB = {'metabat_9342': 'Cycloclasticus metabat_9342 (2024-25 bloom)', 'lorbin_184': 'Cycloclasticus lorbin_184 (2025-26 bloom)', 'comebin_3787_sub': 'Cycloclasticus comebin_3787_sub',
       'lorbin_722_sub': 'Cycloclasticus lorbin_722_sub', 'metabat_5432': 'Cycloclasticus metabat_5432', 'lorbin_167': 'Sulfurimonas lorbin_167', 'semibin_088': 'Sulfurimonas semibin_088',
       'lorbin_032': 'Sulfitobacter lorbin_032', 'metabat_1122': 'Sulfitobacter metabat_1122', 'comebin_2596': 'Casp-alpha2 comebin_2596', 'comebin_2576': 'Porticoccaceae 50-400-T64 comebin_2576',
       'lorbin_369': 'Oleispira lorbin_369', 'lorbin_042': 'Yoonia lorbin_042 (control-pool marker)', 'semibin_072': 'UBA5078 semibin_072 (control-pool marker)', 'metabat_8244': 'Cognaticolwellia metabat_8244',
       'plasmid_contig_0009679': 'Sphingorhabdus CYP153/alkB plasmid contig_0009679', 'plasmid_contig_0009104': 'Pseudosulfitobacter alkB plasmid-like contig_0009104'}
out = {'targets': [{'id': t, 'label': LAB.get(t, t)} for t in T],
       'samples': [{'name': x['name'], 'season': x.season, 'pool': x.pool, 'medium': x.medium, 'day': None if pd.isna(x.day) else int(x.day)} for _, x in R.iterrows()],
       'values': {t: [round(float(v), 4) for v in R[t]] for t in T}}
json.dump(out, open('../out/mg_site.json', 'w'), separators=(',', ':'))
pd.set_option('display.width', 250)
key = ['metabat_9342', 'lorbin_184', 'semibin_088', 'lorbin_032', 'comebin_2596', 'plasmid_contig_0009679']
print(R[R.pool != 'other'].groupby(['season', 'pool', pd.cut(R.day, [-10, 0, 5, 12, 25, 70, 100])])[key].median().round(3).dropna(how='all').to_string())
