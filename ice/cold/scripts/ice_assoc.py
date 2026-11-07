"""Per-MAG ice vs water (and TOP vs BOT) association from RNA (bin TPM) and nanopore DNA, paired by season x pool x day.
-> ../mag_ice_assoc.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
m = meta().loc[samples()]
m['layer'] = m.ice_layer.fillna('')
T = table('bin_tpm'); T = T[[c for c in T.columns if c in m.index]]
T = T.loc[[i for i in T.index if i.startswith('dastool-')]]
C = clr(T.values + 0.0, pseudo=1.0); C = pd.DataFrame(C, index=T.index, columns=T.columns)

def paired(C, info, a_mask, b_mask, keys=('season', 'group', 'day')):
    """median and sign-consistency of mean(a) - mean(b) over blocks having both"""
    d = []
    for k, g in info.groupby(list(keys)):
        a = g.index[a_mask(g)]; b = g.index[b_mask(g)]
        if len(a) and len(b): d.append(C[a].mean(1) - C[b].mean(1))
    D = pd.concat(d, axis=1)
    return D.median(1), (np.sign(D).mean(1)), D.shape[1]

x = m[(m.group != 'other') & m.medium.isin(['SI', 'SW'])]
ice_med, ice_sgn, n1 = paired(C, x, lambda g: g.medium == 'SI', lambda g: g.medium == 'SW')
tb_med, tb_sgn, n2 = paired(C, x[x.layer.isin(['TOP', 'BOT'])], lambda g: g.layer == 'TOP', lambda g: g.layer == 'BOT')
print('RNA blocks ice/water', n1, 'top/bot', n2)
# DNA (nanopore, per-bin depth -> relative abundance)
dd = pd.read_csv('/data/dev/rnaseq-exploration/dna/bin_depths.tsv', sep='\t', index_col=0).drop(columns='length')
dd.columns = [c.replace('.sorted.bam', '') for c in dd.columns]
ns = pd.read_csv('/data/dev/rnaseq-exploration/dna/nanopore_samples.tsv', sep='\t').set_index('run_barcode')
ns = ns.loc[[c for c in dd.columns if c in ns.index]]
u = ns.sample_name.str.upper()
ns['group'] = np.where(u.str.contains('CTL'), 'CTL', np.where(u.str.contains('EXP'), 'EXP', 'other'))
ns['medium'] = np.where(u.str.contains('_SI'), 'SI', np.where(u.str.contains('_SW'), 'SW', 'other'))
ns['layer'] = np.where(u.str.contains('_TOP'), 'TOP', np.where(u.str.contains('_BOT'), 'BOT', ''))
ns['day'] = u.str.extract(r'(?:DAY|PRE)(\d+)')[0].astype(float)
ns['season'] = ns['log'].astype(str)
nsx = ns[(ns.group != 'other') & (ns.medium != 'other') & ns.day.notna()]
R = dd[nsx.index]; R = R / R.sum(0) * 100
CD = pd.DataFrame(clr(R.values, pseudo=1e-3), index=R.index, columns=R.columns)
dice_med, dice_sgn, n3 = paired(CD, nsx, lambda g: g.medium == 'SI', lambda g: g.medium == 'SW')
dtb_med, dtb_sgn, n4 = paired(CD, nsx[nsx.layer != ''], lambda g: g.layer == 'TOP', lambda g: g.layer == 'BOT')
print('DNA blocks ice/water', n3, 'top/bot', n4)
q = pd.read_csv('/data/dev/rnaseq-exploration/data/summary.tsv', sep='\t', index_col=0)
out = pd.DataFrame({'rna_ice': ice_med, 'rna_ice_sign': ice_sgn, 'rna_top': tb_med, 'rna_top_sign': tb_sgn,
                    'dna_ice': dice_med, 'dna_ice_sign': dice_sgn, 'dna_top': dtb_med, 'dna_top_sign': dtb_sgn,
                    'max_tpm': T.max(1), 'mean_dna_pct': R.mean(1)})
out = out.join(q[['SCG_completeness', 'SCG_redundancy', 'size']])
def cls(r):
    ok_r = r.max_tpm >= 50; ok_d = r.mean_dna_pct >= 0.02 if pd.notna(r.mean_dna_pct) else False
    s_ice = (ok_r and r.rna_ice >= 1 and r.rna_ice_sign >= 0.5) or (ok_d and r.dna_ice >= 1 and r.dna_ice_sign >= 0.5)
    s_w = (ok_r and r.rna_ice <= -1 and r.rna_ice_sign <= -0.5) or (ok_d and r.dna_ice <= -1 and r.dna_ice_sign <= -0.5)
    if s_ice and not s_w: return 'ice'
    if s_w and not s_ice: return 'water'
    if ok_r or ok_d: return 'both' if not (s_ice or s_w) else 'conflict'
    return 'rare'
out['habitat'] = out.apply(cls, axis=1)
out.index.name = 'bin'
out.round(4).to_csv('../mag_ice_assoc.tsv', sep='\t')
print(out.habitat.value_counts())
print(out[['rna_ice', 'dna_ice']].corr(method='spearman'))
