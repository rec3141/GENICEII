"""Cross-platform identity table (2024-25, day <= 7 seawater + every corrected/flagged library) -> swapcheck/identity.tsv.
Diesel-pool signature = Sulfitobacter lorbin_032 + Casp-alpha2 comebin_2596 + Sulfurimonas semibin_088; control markers = Yoonia lorbin_042 + UBA5078 semibin_072.
Units differ by platform: nanopore % of binned coverage, Illumina % of read k-mers, RNA % of binned TPM."""
import sys, re; sys.path.insert(0, 'scripts')
import pandas as pd, numpy as np
from common import meta, corrections
DIE = ['lorbin_032', 'comebin_2596', 'semibin_088']; CTL = ['lorbin_042', 'semibin_072']
early = lambda n: bool(re.search(r'CMO_(CTL|EXP)_(PRE0\d|DAY0[0-7])_SW', n or ''))
rows = []
D = pd.read_csv('dna/bin_depths.tsv', sep='\t', index_col=0); L = D.pop('length'); D.columns = [c.replace('.sorted.bam', '') for c in D.columns]
rel = D.mul(L, axis=0).drop(index='unbinned'); rel = 100 * rel / rel.sum()
ns = pd.read_csv('dna/nanopore_samples.tsv', sep='\t'); C = corrections('nanopore')
for x in ns[ns.run_barcode.isin(rel.columns) & (ns.log.astype(str) == '2024')].itertuples():
    if early(x.sample_name) or x.run_barcode in C.index:
        c = C.loc[x.run_barcode] if x.run_barcode in C.index else None
        rows.append(dict(platform='nanopore', library=x.run_barcode, label=x.sample_name, corrected=c.sample_name_corrected if c is not None else '', flag=c.qc_flag if c is not None else '',
                         diesel=rel.loc[['dastool-' + b for b in DIE], x.run_barcode].sum(), control=rel.loc[['dastool-' + b for b in CTL], x.run_barcode].sum()))
I = pd.read_csv('out/mg_kmer.tsv', sep='\t', index_col=0)
for s, x in I[I.season == '2024-25'].iterrows():
    if early(x.label_original) or isinstance(x.qc_flag, str):
        rows.append(dict(platform='Illumina', library=x.gid, label=x.label_original, corrected=x['name'] if x['name'] != x.label_original else '', flag=x.qc_flag if isinstance(x.qc_flag, str) else '',
                         diesel=np.nan if 'lorbin_042' not in I else x[DIE].sum(), control=x[CTL].sum()))
T = pd.read_csv('data/marine_coassembly.bin_tpm.tsv', sep='\t', index_col=0).drop(index='unbinned', errors='ignore'); T = 100 * T / T.sum(); M = meta()
for lib in T.columns:
    if lib in M.index and M.loc[lib, 'season'] == '2024-25':
        lab = M.loc[lib, 'label_original'] if isinstance(M.loc[lib, 'label_original'], str) else M.loc[lib, 'sample_name']
        if early(lab):
            rows.append(dict(platform='RNA', library=lib, label=lab, corrected=M.loc[lib, 'sample_name'] if lab != M.loc[lib, 'sample_name'] else '', flag=M.loc[lib, 'qc_flag'] if isinstance(M.loc[lib, 'qc_flag'], str) else '',
                             diesel=T.loc[['dastool-' + b for b in DIE], lib].sum(), control=T.loc[['dastool-' + b for b in CTL], lib].sum()))
R = pd.DataFrame(rows)
R['day'] = R.label.str.extract(r'_(PRE|DAY)(\d+)').apply(lambda r: -int(r[1]) if r[0] == 'PRE' else int(r[1]), axis=1)
R = R.sort_values(['day', 'label', 'platform', 'library'])
R.to_csv('swapcheck/identity.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(R.round(3).to_string(index=False))
