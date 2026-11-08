"""Shared differential-abundance contrasts (used by da.py logic and hull_mags.py). See da.py for the contrast definitions."""
import pandas as pd, numpy as np, sys
from scipy.stats import mannwhitneyu, wilcoxon
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import meta, nanopore_samples, bh, OIL_DAY
import contrasts as C
def dna_meta(cols):
    ns = nanopore_samples(); di = pd.read_csv('/data/dev/rnaseq-exploration/out/dna_samples.tsv', sep='\t')
    di['sample_name'] = di.sample_name.astype(str); dmap = di.drop_duplicates('sample_name').set_index('sample_name'); md = []
    for rb in cols:
        sn = ns.sample_name[ns.run_barcode == rb]; sn = sn.iloc[0] if len(sn) else ''
        r = dmap.loc[sn] if sn in dmap.index else None
        if r is None or pd.isna(r.day) or r.group not in ('CTL', 'EXP'): md.append(dict(lib=rb, season=None, group=None, medium=None, day=np.nan, phase='other')); continue
        od = OIL_DAY['CMO ' + r.season]
        ph = 'pre' if r.day < 0 else ('oiled' if (r.group == 'EXP' and r.day >= od) else ('control' if r.group == 'CTL' else 'pre'))
        md.append(dict(lib=rb, season=r.season, group=r.group, medium=r.medium, day=r.day, phase=ph))
    return pd.DataFrame(md).set_index('lib')
def run(ab_r, ab_d):
    """ab_r: features x RNA libraries (abundance); ab_d: features x DNA libraries. Returns long table with q (BH per data x contrast)."""
    mr = meta().loc[ab_r.columns, ['season', 'group', 'medium', 'phase', 'day']]; md = dna_meta(ab_d.columns)
    res = []
    def add(blob, dt, con, eff, p, n1, n2):
        res.append(dict(blob=blob, data=dt, contrast=con, log2_effect=eff, p=p, n_a=n1, n_b=n2))
    def lg(v):
        nz = v[v > 0]; pc = 0.01 * (nz.median() if len(nz) else 1e-9); return np.log2(v + pc)
    res = []
    def add(blob, dt, con, eff, p, n1, n2):
        res.append(dict(blob=blob, data=dt, contrast=con, log2_effect=eff, p=p, n_a=n1, n_b=n2))
    for dt, A, M in (('RNA', ab_r, mr), ('DNA', ab_d, md)):
        for b in A.index:
            v = lg(A.loc[b])
            for s in ('2024-25', '2025-26'):
                o = M.index[(M.season == s) & (M.phase == 'oiled') & M.medium.isin(['SW', 'SI'])]; c = M.index[(M.season == s) & (M.phase == 'control') & M.medium.isin(['SW', 'SI'])]
                if len(o) >= 3 and len(c) >= 3:
                    add(b, dt, f'diesel_vs_control_{s}', v[o].median() - v[c].median(), mannwhitneyu(v[o], v[c]).pvalue, len(o), len(c))
            s = '2025-26'; sw = (M.season == s) & (M.medium == 'SW')
            pe, pc_, o, c = [M.index[sw & (M.phase == 'pre') & (M.group == g)] for g in ('EXP', 'CTL')] + [M.index[sw & (M.phase == ph)] for ph in ('oiled', 'control')]
            if min(len(pe), len(pc_)) >= 2 and min(len(o), len(c)) >= 3:
                eff = (v[o].median() - v[pe].median()) - (v[c].median() - v[pc_].median())
                add(b, dt, 'baci_2025-26_SW', eff, mannwhitneyu(v[o] - v[pe].median(), v[c] - v[pc_].median()).pvalue, len(o), len(c))
            # ice vs water pairs
            if dt == 'RNA':
                P = [dict(p, a=[x for x in p['a'] if x in A.columns], b=[x for x in p['b'] if x in A.columns]) for p in C.PI]
                P = [p for p in P if p['a'] and p['b']]
                dl = [v[p['a']].mean() - v[p['b']].mean() for p in P]
            else:
                dl = []
                for (s_, g_, d_), h in M[M.medium.isin(['SW', 'SI'])].groupby(['season', 'group', 'day']):
                    if (h.medium == 'SI').any() and (h.medium == 'SW').any(): dl.append(v[h.index[h.medium == 'SI']].mean() - v[h.index[h.medium == 'SW']].mean())
            dl = np.array(dl)
            if len(dl) >= 5: add(b, dt, 'ice_vs_water', float(np.mean(dl)), wilcoxon(dl).pvalue if np.any(dl != 0) else 1.0, int((dl > 0).sum()), len(dl))

    R = pd.DataFrame(res); R['q'] = np.nan
    for (dt, con), g in R.groupby(['data', 'contrast']): R.loc[g.index, 'q'] = bh(g.p.values)
    return R
