"""Nanopore DNA: viral share of mapped bases per library (depth x length over confident non-provirus viral contigs / all contigs),
per class, and per-vOTU depth. Sample labels via common.nanopore_samples() (corrected). Out: out/vdna_samples.tsv, out/vdna_votu.pkl"""
import sys, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
V = '/data/dev/rnaseq-exploration/viruses/'
d = pd.read_csv(V + 'grex/viral_depths.tsv.gz', sep='\t', index_col=0)
L = d.pop('len'); d.columns = [c.replace('.sorted.bam', '') for c in d.columns]
tot = pd.read_csv(V + 'grex/depth_totals.tsv', sep='\t', index_col=0).total_bp_depth; tot.index = tot.index.str.replace('.sorted.bam', '')
vc = pd.read_csv(V + 'out/viral_confident.tsv', sep='\t').set_index('region')
bp = d.mul(L, axis=0)
share = bp.div(tot[d.columns], axis=1) * 100
ns = nanopore_samples(); ns = ns[ns.run_barcode.isin(d.columns)].copy()
ns['season'] = ns.log.astype(str).map({'2024': '2024-25', '2025': '2025-26'})
u = ns.sample_name.str.upper()
ns['group'] = np.where(u.str.contains('_EXP_'), 'EXP', np.where(u.str.contains('_CTL_'), 'CTL', 'other'))
ns['medium'] = np.where(u.str.contains('_SI'), 'SI', np.where(u.str.contains('_SW'), 'SW', 'other'))
def day(s):
    m_ = re.search(r'(PRE|DAY)(\d+)', s.upper()); return np.nan if not m_ else (-int(m_[2]) if m_[1] == 'PRE' else int(m_[2]))
ns['day'] = ns.sample_name.map(day)
cls = vc.cls.reindex(share.index).fillna('').replace('', 'unclassified')
by = share.groupby(cls.values).sum().T
by['viral_pct'] = share.sum()
ns = ns.set_index('run_barcode').join(by)
ns['tot_bp_depth'] = tot
ns.to_csv(V + 'out/vdna_samples.tsv', sep='\t')
vot = bp.groupby(vc.votu_rep.reindex(bp.index).fillna(pd.Series(bp.index, index=bp.index)).values).sum()
pd.to_pickle({'bp': vot.astype(np.float32), 'tot': tot}, V + 'out/vdna_votu.pkl')
x = ns[(ns.group != 'other') & ns.day.notna()]
pd.set_option('display.width', 250)
print(x.groupby(['season', 'group', 'medium', pd.cut(x.day, [-10, 0, 5, 12, 25, 90])])[['viral_pct', 'Caudoviricetes', 'Megaviricetes']].median().round(2).dropna(how='all').to_string())
print(len(ns), 'nanopore libraries;', len(x), 'CMO pool libraries')
