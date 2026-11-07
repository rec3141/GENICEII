"""Simple under-ice light model per pool: clear-sky surface PAR (par_model.py) attenuated by the measured snow and ice
(ice/phys/cores.tsv, CMO 'Ice physical' sheets). No roof, cloud or oil.
  T = (1 - albedo) * exp(-k_snow * h_snow - k_ice * h_ice)
  albedo 0.85 with snow (>= 1 cm), 0.50 bare ice (> 2 cm), 0.07 open water; k_snow 15 m-1, k_ice 1.5 m-1 (Perovich 1996 ranges).
Snow and ice are interpolated linearly between core days, held flat before the first and after the last core
(snow before the first record = 0 if the first record is 0). Under-ice PAR is the light reaching the ice bottom / water below.
Out: light/out/underice_daily.tsv, underice_samples.tsv, underice_corr.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr
O = '/data/dev/rnaseq-exploration/light/out/'
A_SNOW, A_ICE, A_WATER, K_SNOW, K_ICE = 0.85, 0.50, 0.07, 15.0, 1.5
D0 = {'2024-25': pd.Timestamp('2024-11-28'), '2025-26': pd.Timestamp('2025-11-20')}   # same day-0 as light_site.py
P = pd.read_csv(O + 'par_daily.tsv', sep='\t', parse_dates=['date'])
P['season'] = np.where((P.date.dt.year == 2024) | ((P.date.dt.year == 2025) & (P.date.dt.month < 6)), '2024-25', '2025-26')
P['day'] = [(d - D0[s]).days for d, s in zip(P.date, P.season)]
C = pd.read_csv('/data/dev/rnaseq-exploration/ice/phys/cores.tsv', sep='\t')

def fill(days, d, v, zero_lead=False):
    ok = v.notna()
    if not ok.any(): return np.full(len(days), np.nan)
    dd, vv = d[ok].values.astype(float), v[ok].values.astype(float)
    out = np.interp(days, dd, vv)                      # np.interp holds end values flat
    if zero_lead and vv[0] == 0: out[days < dd[0]] = 0
    return out

def trans(hi_cm, hs_cm):
    hi, hs = np.nan_to_num(hi_cm) / 100, np.nan_to_num(hs_cm) / 100
    alb = np.where(hs >= 0.01, A_SNOW, np.where(hi > 0.02, A_ICE, A_WATER))
    return (1 - alb) * np.exp(-K_SNOW * hs - K_ICE * hi)

rows = []
for s in ('2024-25', '2025-26'):
    p = P[P.season == s]
    for pool in ('CTL', 'EXP'):
        c = C[(C.season == s) & (C.pool == pool)].sort_values('day')
        hi = fill(p.day.values, c.day, c.thickness_cm); hs = fill(p.day.values, c.day, c.snow_cm, zero_lead=True)
        hs = np.where(np.isnan(hs), 0, hs)
        T = trans(hi, hs)
        rows.append(pd.DataFrame({'season': s, 'pool': pool, 'date': p.date.dt.date.astype(str).values, 'day': p.day.values, 'par_surface': p.par_mol.values,
                                  'ice_cm': hi, 'snow_cm': hs, 'transmittance': T, 'par_under': p.par_mol.values * T}))
U = pd.concat(rows); U.to_csv(O + 'underice_daily.tsv', sep='\t', index=False)

# per sample (seawater and ice in CTL/EXP pools)
L = pd.read_csv(O + 'sample_light.tsv', sep='\t', index_col=0)
key = U.set_index(['season', 'pool', 'day'])
def look(x):
    try: return key.loc[(x.season, x.pool, int(x.day))]
    except Exception: return None
vals = {s: look(x) for s, x in L.iterrows() if x.pool in ('CTL', 'EXP') and pd.notna(x.day)}
L['par_under'] = pd.Series({s: v.par_under for s, v in vals.items() if v is not None})
L['snow_cm'] = pd.Series({s: v.snow_cm for s, v in vals.items() if v is not None})
L['ice_model_cm'] = pd.Series({s: v.ice_cm for s, v in vals.items() if v is not None})
L.to_csv(O + 'underice_samples.tsv', sep='\t')

# within-genome light-gene shares vs under-ice PAR (seawater, within pool and season; and pooled across pools per season)
SH = pd.read_pickle(O + 'genome_share_samples.pkl')
GS = pd.read_csv(O + 'genome_share.tsv', sep='\t')
out = []
for (kind, b), v in SH.items():
    v = v.dropna()
    for s in ('2024-25', '2025-26'):
        for g in ('EXP', 'CTL', 'both'):
            ix = [i for i in v.index if i in L.index and L.season[i] == s and L.medium[i] == 'SW' and (g == 'both' or L.pool[i] == g) and pd.notna(L.par_under[i]) and L.phase[i] != 'other']
            if len(ix) < 6 or v[ix].std() == 0: continue
            out.append(dict(kind=kind, bin=b, season=s, pool=g, n=len(ix), r_under=spearmanr(v[ix], L.par_under[ix])[0], r_surface=spearmanr(v[ix], L.par_mol[ix])[0],
                            r_day=spearmanr(v[ix], L.day[ix])[0]))
R = pd.DataFrame(out).merge(GS[['kind', 'bin', 'taxon', 'host_mean']], on=['kind', 'bin'], how='left'); R.to_csv(O + 'underice_corr.tsv', sep='\t', index=False)
pd.set_option('display.width', 220)
print(U.groupby(['season', 'pool']).agg(T_min=('transmittance', 'min'), T_max=('transmittance', 'max'), under_min=('par_under', 'min'), under_max=('par_under', 'max')).round(3))
print(R.groupby(['kind', 'season', 'pool']).agg(n=('bin', 'size'), r_under=('r_under', 'median'), r_surface=('r_surface', 'median'), r_day=('r_day', 'median')).round(2).to_string())
print(R[R.bin.str.contains('metabat_1122')].round(2).to_string())

# ---- correlations with the extra light values (2026-11-07 rerun) ----
# r_under: Spearman with modelled under-ice PAR; r_T: with ice/snow transmittance alone (no seasonal sun);
# partial: rank correlation after removing a linear trend in rank(day) within each pool, pooled over both pools of a season.
# Under-ice light is not monotonic in time (open water early, a sharp drop when snow settles), so the partial isolates that part.
from scipy.stats import rankdata, pearsonr
def partial(v, x, day, grp):
    rv, rx = np.zeros(len(v)), np.zeros(len(v))
    for g in np.unique(grp):
        ix = grp == g
        if ix.sum() < 4: continue
        rd = rankdata(day[ix]); A = np.c_[np.ones(ix.sum()), rd]
        for src, dst in ((v, rv), (x, rx)):
            y = rankdata(src[ix]); dst[ix] = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    ok = np.isin(grp, [g for g in np.unique(grp) if (grp == g).sum() >= 4])
    if ok.sum() < 8 or rv[ok].std() == 0 or rx[ok].std() == 0: return np.nan, np.nan
    return pearsonr(rv[ok], rx[ok])
L['T'] = pd.Series({s: v.transmittance for s, v in vals.items() if v is not None})
SW = L[(L.medium == 'SW') & L.pool.isin(['CTL', 'EXP']) & L.par_under.notna() & (L.phase != 'other')]
E = pd.read_csv(O + 'cat_cpm.tsv', sep='\t', index_col=0)
crow = []
for name in E.index:
    for s in ('2024-25', '2025-26'):
        x = SW[SW.season == s]; v = np.log2(E.loc[name, x.index].astype(float) + 0.5)
        row = dict(category=name, season=s)
        for g in ('EXP', 'CTL'):
            ix = x.pool == g
            if ix.sum() >= 6 and v[ix].std() > 0:
                row[f'r_under_{g}'] = spearmanr(v[ix], x.par_under[ix])[0]; row[f'r_T_{g}'] = spearmanr(v[ix], x['T'][ix])[0]
        if v.std() > 0:
            pr, pp = partial(v.values, x.par_under.values, x.day.values, x.pool.values); row['partial_r'] = pr; row['partial_p'] = pp
        crow.append(row)
CR = pd.DataFrame(crow); CR['partial_q'] = bh(CR.partial_p.fillna(1)); CR.to_csv(O + 'cat_under_corr.tsv', sep='\t', index=False)
# genome level: add r_T and the day-controlled partial
grow = []
for (kind, b), v in SH.items():
    v = v.dropna()
    for s in ('2024-25', '2025-26'):
        x = SW[(SW.season == s) & SW.index.isin(v.index)]
        if len(x) < 8 or v[x.index].std() == 0: continue
        pr, pp = partial(v[x.index].values, x.par_under.values, x.day.values, x.pool.values)
        grow.append(dict(kind=kind, bin=b, season=s, n=len(x), partial_r=pr, partial_p=pp, r_T=spearmanr(v[x.index], x['T'])[0]))
GP = pd.DataFrame(grow); GP['partial_q'] = bh(GP.partial_p.fillna(1))
GP = GP.merge(GS[['kind', 'bin', 'taxon', 'host_mean']], on=['kind', 'bin'], how='left'); GP.to_csv(O + 'underice_partial.tsv', sep='\t', index=False)
print(CR.round(2).to_string())
print(GP.groupby(['kind', 'season']).agg(n=('bin', 'size'), med_partial=('partial_r', 'median'), n_q05=('partial_q', lambda q: int((q < 0.05).sum())), med_rT=('r_T', 'median')).round(2))
print(GP[GP.partial_q < 0.1].sort_values('partial_q').round(3).to_string())
print(GP[GP.bin.str.contains('metabat_1122')].round(3).to_string())
