"""DNA (nanopore) MAG abundance in CMO samples and RNA:DNA activity ratios for matched samples."""
import sys, re; sys.path.insert(0, '.')
from common import *
O = '../out/'
d = pd.read_csv('../dna/bin_depths.tsv', sep='\t', index_col=0); d.columns = [c.replace('.sorted.bam', '') for c in d.columns]
L = d.pop('length')
ns = pd.read_csv('../dna/nanopore_samples.tsv', sep='\t'); ns = ns[ns.run_barcode.isin(d.columns)].copy()
C = corrections('nanopore'); ns['label_original'] = ns.sample_name
fix = ns.run_barcode.map(C.sample_name_corrected).fillna(''); ns['sample_name'] = ns.sample_name.where(fix == '', fix)
ns['season'] = ns.log.map({2024: '2024-25', 2025: '2025-26'}).fillna(ns.log.astype(str).map({'2024': '2024-25', '2025': '2025-26'}))
u = ns.sample_name.str.upper()
ns['group'] = np.where(u.str.contains('_EXP_'), 'EXP', np.where(u.str.contains('_CTL_'), 'CTL', 'other'))
ns['medium'] = np.where(u.str.contains('_SI'), 'SI', np.where(u.str.contains('_SW'), 'SW', 'other'))
def day(s):
    m_ = re.search(r'(PRE|DAY)(\d+)', s.upper())
    return np.nan if not m_ else (-int(m_[2]) if m_[1] == 'PRE' else int(m_[2]))
ns['day'] = ns.sample_name.map(day)
bins = d.drop(index='unbinned')
# relative abundance: share of binned coverage (depth × length) per sample
cov = bins.mul(L.drop('unbinned'), axis=0)
rel = cov / cov.sum() * 100
# collapse resequenced libraries of the same sample (same season+name)
ns['key'] = ns.season + '|' + ns.sample_name
relS = rel[ns.run_barcode].T.groupby(ns.key.values).mean().T
info = ns.drop_duplicates('key').set_index('key')[['sample_name', 'season', 'group', 'medium', 'day']]
relS.to_csv(O + 'dna_rel_abundance.tsv', sep='\t'); info.to_csv(O + 'dna_samples.tsv', sep='\t')
# RNA:DNA for matched samples
m = meta(); tpm = table('bin_tpm').drop(index='unbinned', errors='ignore')
rk = m.season + '|' + m.sample_name.str.replace(r'_SW([ab])$', lambda x: '_SW_' + x[1].upper(), regex=True)
match = [(k, s) for s, k in rk.items() if k in relS.columns]
print('matched RNA/DNA samples:', len(match))
rna_rel = tpm / tpm.sum() * 100
ratio = pd.DataFrame({k: np.log2((rna_rel[s] + 0.01) / (relS[k] + 0.01)) for k, s in match})
ratio.to_csv(O + 'rna_dna_log2ratio.tsv', sep='\t')
pd.Series(dict(match)).to_csv(O + 'rna_dna_matches.tsv', sep='\t', header=False)
key = {'Cycloclasticus': ['dastool-lorbin_184', 'dastool-metabat_9342', 'dastool-metabat_5432'], 'Sulfitobacter lorbin_032': ['dastool-lorbin_032'],
       'Casp-alpha2 comebin_2596': ['dastool-comebin_2596'], 'Porticoccaceae x4': ['dastool-comebin_2576', 'dastool-metabat_1191', 'dastool-comebin_3425', 'dastool-lorbin_288'],
       'Sulfurimonas x2': ['dastool-lorbin_167', 'dastool-semibin_088'], 'Oleispira lorbin_369': ['dastool-lorbin_369']}
t = pd.DataFrame({k: relS.loc[v].sum() for k, v in key.items()}).join(info).sort_values(['season', 'group', 'day'])
pd.set_option('display.width', 250)
print(t[t.group != 'other'].round(2).to_string())
print(t[t.group == 'other'][['sample_name', 'Cycloclasticus', 'Sulfurimonas x2', 'Sulfitobacter lorbin_032']].round(2).to_string())
