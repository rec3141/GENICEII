"""Do ice-associated MAGs have more 'cold-adapted' proteomes or more cold genes? Controls: GC, family fixed effects, within-genus pairs.
Uses the parent's paired ice/water MAG contrast (ice/out/mag_icewater.tsv, delta = median log2 ice/water across 24 blocks).
-> ../mag_cold_table.tsv, ../mag_test.tsv, ../genus_pairs.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
from scipy.stats import spearmanr, wilcoxon
MI = pd.read_csv(C_DIR + 'mag_indices.tsv', sep='\t', index_col=0)
IW = pd.read_csv('/data/dev/rnaseq-exploration/ice/out/mag_icewater.tsv', sep='\t', index_col=0)
TB = pd.read_csv('/data/dev/rnaseq-exploration/ice/out/mag_topbot.tsv', sep='\t', index_col=0)
inv = pd.read_csv(C_DIR + 'cold_inventory_mag.tsv', sep='\t', index_col=0)
tx = pd.read_csv('/data/dev/rnaseq-exploration/data/gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1])
tx['bin'] = np.where(tx.user_genome.str.startswith('dastool-'), tx.user_genome, 'dastool-' + tx.user_genome)
tx = tx.drop_duplicates('bin').set_index('bin').classification.str.split(';', expand=True)
tx.columns = list('dpcofgs'); tx = tx.apply(lambda c: c.str[3:])
D = MI.join(IW[['delta', 'n_up', 'n_down', 'q']].rename(columns={'delta': 'ice_delta', 'q': 'ice_q'})).join(TB[['delta']].rename(columns={'delta': 'top_delta'}))
D = D.join(tx[['p', 'c', 'f', 'g', 's']]).join(inv.fillna(0))
D = D[(D.SCG_completeness >= 50)]
D['genome_mb'] = D.size_bp / 1e6
D['ibp_n'] = D.get('Ice-binding proteins', 0)
D.round(5).to_csv(C_DIR + 'mag_cold_table.tsv', sep='\t')
E = D.dropna(subset=['ice_delta'])
print('MAGs >=50% complete:', len(D), ' with ice/water contrast:', len(E))
IDX = ['arg_lys', 'pro', 'acidic', 'charged', 'gly', 'ser', 'gln_thr', 'aliphatic', 'gravy', 'polar_unch']
CATS = [c for c in inv.columns]
rows = []
def resid_fe(y, d, fe='f'):
    """residual of y after GC (linear) + family fixed effects"""
    Xd = pd.get_dummies(d[fe].fillna('NA'), drop_first=False).astype(float)
    Xm = np.column_stack([d.gc.values, Xd.values])
    b, *_ = np.linalg.lstsq(Xm, y, rcond=None); return y - Xm @ b
for v in IDX + CATS + ['gc', 'genome_mb']:
    y = E[v].astype(float).values
    if v in CATS: y = np.log1p(y / E.genome_mb.values)   # genes per Mb (log1p)
    r0, p0 = spearmanr(y, E.ice_delta)
    ok = E.f.notna() & (E.groupby('f').f.transform('size') >= 3)
    if v not in ('gc',):
        ry = resid_fe(y[ok.values], E[ok]); rx = resid_fe(E.ice_delta.values[ok.values], E[ok])
        r1, p1 = spearmanr(ry, rx)
    else: r1 = p1 = np.nan
    rows.append({'variable': v, 'rho_raw': r0, 'p_raw': p0, 'rho_gc_family': r1, 'p_gc_family': p1, 'n_fe': int(ok.sum())})
R = pd.DataFrame(rows); R['q_gc_family'] = bh(R.p_gc_family.fillna(1))
R.round(4).to_csv(C_DIR + 'mag_test.tsv', sep='\t', index=False)
pd.set_option('display.width', 200)
print(R.round(3).to_string(index=False))
# within-genus pairs: most ice-leaning vs most water-leaning MAG of the same genus (|difference| >= 1 log2)
pr = []
for gname, d in E[E.g.fillna('') != ''].groupby('g'):
    if len(d) < 2: continue
    a, b = d.ice_delta.idxmax(), d.ice_delta.idxmin()
    if d.ice_delta[a] - d.ice_delta[b] < 1: continue
    r = {'genus': gname, 'ice_mag': a, 'water_mag': b, 'ice_delta_ice': d.ice_delta[a], 'ice_delta_water': d.ice_delta[b]}
    for v in IDX + ['gc', 'genome_mb']: r['d_' + v] = d.loc[a, v] - d.loc[b, v]
    for v in CATS: r['d_' + v] = d.loc[a, v] - d.loc[b, v]
    pr.append(r)
P = pd.DataFrame(pr); P.round(4).to_csv(C_DIR + 'genus_pairs.tsv', sep='\t', index=False)
print(len(P), 'genus pairs')
for v in IDX + ['gc', 'genome_mb', 'Ice-binding proteins', 'Cold-shock proteins', 'Compatible solute synthesis', 'Compatible solute transport', 'Fatty-acid desaturases', 'PUFA synthase (pfa)', 'EPS export', 'Carotenoid synthesis']:
    x = P['d_' + v].dropna(); x = x[x != 0]
    p = wilcoxon(x).pvalue if len(x) >= 5 else np.nan
    print(f'{v:32s} median diff (ice-water MAG) {P["d_" + v].median():+.4f}  ice>water {(P["d_" + v] > 0).sum()}/{len(P)}  p={p:.3g}')
