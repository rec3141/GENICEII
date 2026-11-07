"""Ice physics from the CMO 'Ice physical' log sheets -> per-core layer summary (T, bulk S, brine volume, brine salinity).
Out: ice/phys/profiles.tsv (tidy), ice/phys/cores.tsv (one row per pool x day)."""
import re, numpy as np, pandas as pd, datetime as dt
D = '/home/grid/Downloads/'; O = '/data/dev/rnaseq-exploration/ice/phys/'

def rng(v):
    """Depth cell -> (lo, hi) cm. Excel turned '4-6' into dates (month-day)."""
    if pd.isna(v): return (np.nan, np.nan)
    if isinstance(v, (dt.datetime, pd.Timestamp)): return (float(v.month), float(v.day))
    s = str(v).strip().replace('>', '')
    m = re.match(r'^(-?[\d.]+)\s*-\s*([\d.]+)$', s)
    if m: return (float(m[1]), float(m[2]))
    try: x = float(s); return (x, x)
    except ValueError: return (np.nan, np.nan)
num = lambda v: pd.to_numeric(pd.Series([v]), errors='coerce')[0]

rows = []
a = pd.read_excel(D + 'CMO 2024 MICB.xlsx', sheet_name='Ice physical')
a.columns = ['id', 'pool', 'day', 'date', 'time', 'loc', 'tdepth', 'T', 'snow', 'sdepth', 'S', 'com']
a[['id', 'pool', 'day']] = a[['id', 'pool', 'day']].ffill()
for _, r in a.iterrows():
    pool = re.search(r'CMO_(EXP|CTL)', str(r.id))[1]  # sheet's Pool column disagrees with the ID once (DAY08 salinity); trust the ID
    day = int(re.sub(r'\D', '', str(r.day)))
    if pd.notna(r['T']) and 'broke' not in str(r.com):  # DAY20 EXP 30 cm: core broke, air-cooled
        lo, hi = rng(r.tdepth)
        rows.append(dict(season='2024-25', pool=pool, day=day, var='T', lo=lo, hi=hi, val=float(r['T']), snow=num(r.snow)))
    if pd.notna(r.S):
        lo, hi = rng(r.sdepth)
        if np.isnan(lo):  # 'depth 1 - bottom' .. 'depth 5 - top' (DAY16 EXP), unknown cm: index 1..5 from bottom
            k = int(re.search(r'depth (\d)', str(r.sdepth))[1]); lo = hi = -k
        rows.append(dict(season='2024-25', pool=pool, day=day, var='S', lo=lo, hi=hi, val=float(r.S), snow=num(r.snow)))
b = pd.read_excel(D + 'CMO 2025 MICB.xlsx', sheet_name='Ice physical')
b.columns = ['id', 'usage', 'pool', 'core', 'day', 'date', 'time', 'loc', 'wT', 'thick', 'kovac', 'freeb', 'trange', 'iceT', 'isT', 'ssT', 'airT', 'snow', 'sTop', 'sBot', 'com']
thick25 = {}
for _, r in b.iterrows():
    day = int(re.sub(r'\D', '', str(r.day))); pool = r.pool
    t = num(r.thick)
    if pd.notna(t): thick25.setdefault((pool, day), []).append(t)
    base = dict(season='2025-26', pool=pool, day=day, snow=num(r.snow) if not isinstance(r.snow, dt.datetime) else np.nan)
    if pd.notna(num(r.iceT)):
        lo, hi = rng(r.trange); rows.append(dict(base, var='T', lo=lo, hi=hi, val=num(r.iceT)))
    if pd.notna(num(r.sTop)):
        lo, hi = rng(r.trange); rows.append(dict(base, var='S', lo=lo, hi=hi, val=num(r.sTop)))
    if pd.notna(num(r.sBot)): rows.append(dict(base, var='Sbot', lo=np.nan, hi=np.nan, val=num(r.sBot)))
    if pd.notna(num(r.wT)): rows.append(dict(base, var='Twater', lo=np.nan, hi=np.nan, val=num(r.wT)))
    if pd.notna(num(r.airT)): rows.append(dict(base, var='Tair', lo=np.nan, hi=np.nan, val=num(r.airT)))
P = pd.DataFrame(rows)
# 2025-26 CTL DAY12 temperatures get colder with depth (-15 -> -19 C at the ice bottom), physically implausible; drop them
P = P[~((P.season == '2025-26') & (P.pool == 'CTL') & (P.day == 12) & (P['var'] == 'T'))]
# air temperature rows in 2024-25 are depth -100
P.loc[(P['var'] == 'T') & (P.lo < -50), 'var'] = 'Tair'
P.loc[(P['var'] == 'T') & (P.lo < 0) & (P.lo > -50), 'var'] = 'Tsnow'
P.to_csv(O + 'profiles.tsv', sep='\t', index=False)

# ice thickness from the RNA/DNA sample logs (2024-25) and the 2025-26 sheet
m = pd.read_csv('/data/dev/genice-sample-metadata/outputs/cmo_rna_metadata.tsv', sep='\t')
th = {}
for _, r in m.iterrows():
    t = pd.to_numeric(re.sub(r'\(.*', '', str(r.ice_thickness_cm)), errors='coerce')
    if pd.notna(t) and r.experiment == 'CMO 2024-25' and r.treatment in ('EXP', 'CTL') and pd.notna(r.day): th[('2024-25', r.treatment, int(r.day))] = t
for (pool, day), v in thick25.items(): th[('2025-26', pool, day)] = float(np.median(v))

def brine(T, S):
    """Brine volume fraction (Frankenstein & Garner 1967) and brine salinity (Assur 1958 / Cox & Weeks fit), valid -0.5..-22.9 C."""
    if not (np.isfinite(T) and np.isfinite(S)) or T >= -0.5: return np.nan, np.nan
    vb = S * (49.185 / abs(T) + 0.532) / 1000
    sb = 1000 * (1 - 54.11 / T) ** -1
    return vb, sb

out = []
for (season, pool, day), g in P.groupby(['season', 'pool', 'day']):
    T = g[g['var'] == 'T'].dropna(subset=['lo']); S = g[g['var'] == 'S'].dropna(subset=['lo'])
    H = th.get((season, pool, day), np.nan)
    if np.isnan(H) and len(T): H = np.nan
    mid = lambda d: (d.lo + d.hi) / 2
    rec = dict(season=season, pool=pool, day=day, thickness_cm=H, snow_cm=g.snow.dropna().max() if g.snow.notna().any() else np.nan,
               T_air=g.loc[g['var'] == 'Tair', 'val'].mean())
    if len(T):
        z = mid(T); rec['T_top'] = T.val[z == z.min()].mean(); rec['T_bot'] = T.val[z == z.max()].mean(); rec['z_Tbot'] = z.max()
        rec['T_mean'] = T.val.mean()
    if len(S):
        if (S.lo < 0).all():  # index-from-bottom layout
            rec['S_top'] = S.val[S.lo == S.lo.min()].mean(); rec['S_bot'] = S.val[S.lo == S.lo.max()].mean()
        else:
            z = mid(S); rec['S_top'] = S.val[z == z.min()].mean(); rec['S_bot'] = S.val[z == z.max()].mean()
        rec['S_mean'] = S.val.mean()
    sb = g.loc[g['var'] == 'Sbot', 'val']
    if len(sb): rec['S_bot'] = sb.mean()
    for L in ('top', 'bot'):
        vb, br = brine(rec.get(f'T_{L}', np.nan), rec.get(f'S_{L}', np.nan))
        rec[f'Vb_{L}'] = vb; rec[f'Sbrine_{L}'] = br
    out.append(rec)
C = pd.DataFrame(out)
C.to_csv(O + 'cores.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(C.round(3).to_string())
