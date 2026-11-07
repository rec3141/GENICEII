"""Low-oxygen expression index for the CMO diesel mesocosms (no dissolved-O2 data exist; do_pct is empty).
Two ways:
 1. Community: log2 ratios of O2-responsive gene groups (CPM of CDS reads): high-affinity terminal oxidases (cbb3, bd) vs low-affinity (aa3, bo3);
    O2-independent vs O2-dependent paralogues (hemN/hemF, nrdDG/nrdAB); plus CPM of anaerobic/microaerobic markers.
 2. Within genome: for each genome carrying both high- and low-affinity oxidases, log2(HA/LA) read ratio per sample, centred on that genome's own mean;
    sample index = weighted mean over genomes (weight sqrt(oxidase reads)). A bloom of a cbb3-heavy genome does not move this index.
    Plus a decomposition of the binned community HA fraction into composition (genome mix at fixed per-genome fractions) and regulation (the rest).
Out: oxygen/out/{group_cpm.tsv, sample_index.tsv, genome_shift.tsv, genome_ha.tsv, contrasts.tsv, covariates.tsv}"""
import sys, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr, wilcoxon, mannwhitneyu
O = '/data/dev/rnaseq-exploration/oxygen/out/'; os.makedirs(O, exist_ok=True)
K = lambda s: s.split()
GROUPS = {  # name: (role, KOs)
 'cbb3 oxidase (ccoNOPQ)': ('HA', K('K00404 K00405 K00406 K00407')),
 'bd oxidase (cydAB, appX)': ('HA', K('K00425 K00426 K22501')),
 'aa3 oxidase (coxABC/ctaCDE)': ('LA', K('K02274 K02275 K02276')),
 'bo3 oxidase (cyoABCD)': ('LA', K('K02297 K02298 K02299 K02300')),
 'hemN (O2-independent coproporphyrinogen oxidase)': ('anaer_paralog', K('K02495')),
 'hemF (O2-dependent coproporphyrinogen oxidase)': ('aer_paralog', K('K00228')),
 'nrdDG (anaerobic class III RNR)': ('anaer_paralog', K('K21636 K04068')),
 'nrdAB (O2-dependent class I RNR)': ('aer_paralog', K('K00525 K00526')),
 'narGHI (membrane nitrate reductase)': ('anaer', K('K00370 K00371 K00374')),
 'napAB (periplasmic nitrate reductase)': ('anaer', K('K02567 K02568')),
 'hoxFUYH (bidirectional hydrogenase)': ('anaer', K('K00436 K18005 K18006 K18007')),
 'H2-evolving NiFe group 4 (hyc, hyf, ech, mbh)': ('anaer', K('K15830 K12136 K12137 K12138 K12139 K12140 K12141 K12142 K12143 K12144 K12145 K14086 K14087 K14088 K14089 K14090 K14091 K18016 K18017 K18023')),
 'pflB + pflA (pyruvate formate-lyase)': ('anaer', K('K00656 K04069')),
 'ldh (lactate dehydrogenase)': ('anaer', K('K00016')),
 'adhE (fermentative alcohol dehydrogenase)': ('anaer', K('K04072')),
 'fnr (CRP/FNR family, includes non-O2 sensors)': ('reg', K('K01420')),
 'arcAB (aerobic respiration control)': ('reg', K('K07773 K07648')),
 'hmp (flavohemoglobin)': ('globin', K('K05916')),
 'glbN (truncated hemoglobin)': ('globin', K('K06886')),
}
X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
kos = G['KO'].fillna('').str.split(',')
allko = set(k for _, ks in GROUPS.values() for k in ks)
idx = np.where(kos.apply(lambda x: bool(allko.intersection(x))).values & cds)[0]
sub = G.iloc[idx].reset_index(drop=True); sub['bin'] = sub['bin'].fillna('')
R = pd.DataFrame(X[idx].toarray().astype(float), columns=S)          # reads
CPM = R / tot * 1e6
skos = sub['KO'].fillna('').str.split(',')
mem = {n: skos.apply(lambda x: bool(set(ks).intersection(x))).values for n, (_, ks) in GROUPS.items()}
E = pd.DataFrame({n: CPM[s].sum(0) for n, s in mem.items()}).T
E.to_csv(O + 'group_cpm.tsv', sep='\t')
role = {n: r for n, (r, _) in GROUPS.items()}
grp = lambda r: [n for n in GROUPS if role[n] == r]
ha = np.any([mem[n] for n in grp('HA')], 0); la = np.any([mem[n] for n in grp('LA')], 0)

I = pd.DataFrame(index=S)
I['comm_HA_LA'] = np.log2((E.loc[grp('HA')].sum() + 0.5) / (E.loc[grp('LA')].sum() + 0.5))
I['comm_cbb3_aa3'] = np.log2((E.loc['cbb3 oxidase (ccoNOPQ)'] + 0.5) / (E.loc['aa3 oxidase (coxABC/ctaCDE)'] + 0.5))
I['comm_bd_aa3'] = np.log2((E.loc['bd oxidase (cydAB, appX)'] + 0.5) / (E.loc['aa3 oxidase (coxABC/ctaCDE)'] + 0.5))
I['comm_hemN_hemF'] = np.log2((E.loc['hemN (O2-independent coproporphyrinogen oxidase)'] + 0.5) / (E.loc['hemF (O2-dependent coproporphyrinogen oxidase)'] + 0.5))
I['comm_nrdDG_nrdAB'] = np.log2((E.loc['nrdDG (anaerobic class III RNR)'] + 0.5) / (E.loc['nrdAB (O2-dependent class I RNR)'] + 0.5))
I['comm_anaer_cpm'] = np.log2(E.loc[grp('anaer')].sum() + 0.5)

# within genome
B = sub.bin.values
HAr = R[ha].groupby(B[ha]).sum(); LAr = R[la].groupby(B[la]).sum()
gen = sorted(g for g in set(HAr.index) & set(LAr.index) - {''} if is_prokaryote(g))   # prokaryote genomes only (bin_domain)
HAr, LAr = HAr.reindex(gen).fillna(0), LAr.reindex(gen).fillna(0)
T = HAr + LAr
elig = [g for g in gen if T.loc[g].sum() >= 300 and HAr.loc[g].sum() >= 30 and LAr.loc[g].sum() >= 30]
L = np.log2((HAr.loc[elig] + 0.5) / (LAr.loc[elig] + 0.5)).where(T.loc[elig] >= 20)
Lc = L.sub(L.mean(1), axis=0)
W = np.sqrt(T.loc[elig]).where(T.loc[elig] >= 20)
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def blab(b):
    if special_label(b): return special_label(b)
    t = [y.split('__', 1)[1] for y in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in y and y.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')
lab = {g: blab(g) for g in gen}
cyc = [g for g in elig if 'Cycloclasticus' in lab[g]]
def wmean(Lc, W, rows):
    a, w = Lc.loc[rows], W.loc[rows]
    return (a * w).sum() / w.where(a.notna()).sum()
I['within_HA_LA'] = wmean(Lc, W, elig)
I['within_HA_LA_noCyclo'] = wmean(Lc, W, [g for g in elig if g not in cyc])
I['within_n_genomes'] = Lc.notna().sum()
# decomposition of binned community HA fraction (log2 odds): composition vs regulation
HAall = R[ha].groupby(B[ha]).sum(); LAall = R[la].groupby(B[la]).sum()
allg = sorted((set(HAall.index) | set(LAall.index)) - {''})
HAall, LAall = HAall.reindex(allg).fillna(0), LAall.reindex(allg).fillna(0)
Tall = HAall + LAall
fbar = (HAall.sum(1) + 0.5) / (Tall.sum(1) + 1)
F = (HAall.sum() + 0.5) / (Tall.sum() + 1)
Fc = (Tall.mul(fbar, axis=0).sum() + 0.5) / (Tall.sum() + 1)
lo = lambda p: np.log2(p / (1 - p))
I['binned_HA_logodds'] = lo(F); I['composition_HA_logodds'] = lo(Fc); I['regulation_HA_logodds'] = lo(F) - lo(Fc)
I['binned_oxidase_share'] = Tall.sum() / (R[ha | la].sum() + 1e-9)
# Cycloclasticus bloom size (TPM of all Cycloclasticus bins)
tpm = table('bin_tpm')
cycbins = [b for b in tpm.index if 'Cycloclasticus' in blab(b)]
I['cyclo_tpm'] = tpm.loc[cycbins, S].sum()
# metadata
for c in ('sample_name', 'season', 'group', 'medium', 'phase', 'day', 'ice_layer', 'water_temp_c', 'ice_thickness_cm', 'roof_open', 'mixers_on'):
    I[c] = m[c]
I.to_csv(O + 'sample_index.tsv', sep='\t')

# per-genome shift: diesel vs control after dosing, seawater, per season (centred log2 HA/LA)
gs = []
for s in ('2024-25', '2025-26'):
    for md in ('SW', 'SI'):
        o = m.index[(m.season == s) & (m.medium == md) & (m.phase == 'oiled')]; c = m.index[(m.season == s) & (m.medium == md) & (m.phase == 'control')]
        for g in elig:
            a, b = L.loc[g, o].dropna(), L.loc[g, c].dropna()
            if len(a) >= 3 and len(b) >= 3:
                gs.append(dict(season=s, medium=md, bin=g, taxon=lab[g], n_oil=len(a), n_ctl=len(b), oil=a.median(), ctl=b.median(), shift=a.median() - b.median(),
                               reads_oil=T.loc[g, o].median(), reads_ctl=T.loc[g, c].median(), p=mannwhitneyu(a, b).pvalue))
GS = pd.DataFrame(gs); GS['q'] = bh(GS.p.values); GS.to_csv(O + 'genome_shift.tsv', sep='\t', index=False)
pd.DataFrame({'bin': elig, 'taxon': [lab[g] for g in elig], 'HA_reads': HAr.loc[elig].sum(1).values, 'LA_reads': LAr.loc[elig].sum(1).values,
              'mean_log2': L.mean(1).values, 'n_samples': L.notna().sum(1).values}).to_csv(O + 'genome_ha.tsv', sep='\t', index=False)
pd.to_pickle({'L': L, 'Lc': Lc, 'W': W, 'lab': lab}, O + 'within.pkl')

# contrasts
sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts'); import contrasts as C
IDX = ['comm_HA_LA', 'comm_cbb3_aa3', 'comm_bd_aa3', 'comm_hemN_hemF', 'comm_nrdDG_nrdAB', 'comm_anaer_cpm', 'within_HA_LA', 'within_HA_LA_noCyclo', 'regulation_HA_logodds', 'composition_HA_logodds']
rows = []
def paired(P, name):
    for k in IDX:
        d = np.array([I.loc[p['a'], k].mean() - I.loc[p['b'], k].mean() for p in P], float); d = d[np.isfinite(d)]
        rows.append(dict(contrast=name, index=k, n=len(d), mean=d.mean(), n_up=int((d > 0).sum()), p=wilcoxon(d).pvalue if len(d) > 5 and np.any(d != 0) else np.nan))
paired(C.PI, 'ice minus water'); paired(C.PT, 'top minus bottom ice')
for s in ('2024-25', '2025-26'):
    for k in IDX:
        x = I[(I.season == s) & (I.medium == 'SW')]
        med = lambda g, ph: x[(x.group == g) & (x.phase == ph)][k].median()
        o, c = x[x.phase == 'oiled'][k].dropna(), x[x.phase == 'control'][k].dropna()
        rows.append(dict(contrast=f'diesel minus control SW {s}', index=k, n=len(o) + len(c), mean=o.median() - c.median(), n_up=None, p=mannwhitneyu(o, c).pvalue if len(o) > 2 and len(c) > 2 else np.nan))
        rows.append(dict(contrast=f'BACI SW {s}', index=k, n=None, mean=(med('EXP', 'oiled') - med('EXP', 'pre')) - (med('CTL', 'control') - med('CTL', 'pre')), n_up=None, p=np.nan))
CT = pd.DataFrame(rows); CT.to_csv(O + 'contrasts.tsv', sep='\t', index=False)

# covariates (Spearman within subsets)
cov = []
cores = pd.read_csv('/data/dev/rnaseq-exploration/ice/phys/cores.tsv', sep='\t')
def core_val(s, pool, day, col):
    c = cores[(cores.season == s) & (cores.pool == pool)].dropna(subset=[col])
    if c.empty or pd.isna(day): return np.nan
    j = (c.day - day).abs().idxmin(); return c.loc[j, col] if abs(c.loc[j, 'day'] - day) <= 2 else np.nan
for col in ('T_top', 'Vb_top', 'Sbrine_top', 'Vb_bot', 'thickness_cm'):
    I['core_' + col] = [core_val(s, g, d, col) for s, g, d in zip(I.season, I.group, I.day)]
I['water_temp_c'] = pd.to_numeric(I.water_temp_c, errors='coerce'); I['ice_thickness_cm'] = pd.to_numeric(I.ice_thickness_cm, errors='coerce')
I['log_cyclo'] = np.log10(I.cyclo_tpm + 1)
subsets = {'all CTL/EXP water': (I.group.isin(['CTL', 'EXP'])) & (I.medium == 'SW'), 'control pool water': (I.group == 'CTL') & (I.medium == 'SW'),
           'diesel pool water': (I.group == 'EXP') & (I.medium == 'SW'), 'all ice': I.medium == 'SI', 'top ice': (I.medium == 'SI') & I.ice_layer.fillna('').str.contains('top', case=False)}
for sn, mask in subsets.items():
    for v in ('log_cyclo', 'water_temp_c', 'ice_thickness_cm', 'core_thickness_cm', 'core_T_top', 'core_Vb_top', 'core_Sbrine_top', 'day'):
        for k in ('comm_HA_LA', 'within_HA_LA', 'within_HA_LA_noCyclo', 'regulation_HA_logodds'):
            x = I.loc[mask, [v, k]].dropna()
            if len(x) >= 6 and x[v].nunique() > 2:
                r, p = spearmanr(x[v], x[k]); cov.append(dict(subset=sn, covariate=v, index=k, n=len(x), rho=r, p=p))
for fld in ('roof_open', 'mixers_on'):
    for k in ('comm_HA_LA', 'within_HA_LA'):
        x = I[I.group.isin(['CTL', 'EXP'])].dropna(subset=[k])
        a, b = x[x[fld] == 'Yes'][k], x[x[fld] == 'No'][k]
        if len(a) >= 3 and len(b) >= 3: cov.append(dict(subset='CTL/EXP all media', covariate=f'{fld} Yes vs No', index=k, n=len(a) + len(b), rho=a.median() - b.median(), p=mannwhitneyu(a, b).pvalue))
CV = pd.DataFrame(cov); CV.to_csv(O + 'covariates.tsv', sep='\t', index=False)
I.to_csv(O + 'sample_index.tsv', sep='\t')
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
print('eligible genomes', len(elig), 'cyclo', [lab[g] + ' ' + g for g in cyc]); print('binned oxidase share median', I.binned_oxidase_share.median().round(2))
print(CT.round(3).to_string())
print(CV[CV.p < 0.1].round(3).to_string())
print(GS.sort_values('q').head(25).round(3).to_string())
mos = I[I.sample_name.str.contains('MOS', na=False)]
print(mos[IDX].round(2).to_string()); print(I[(I.season == '2025-26') & (I.group == 'EXP') & (I.day.between(3, 5))][['sample_name', 'medium'] + IDX].round(2).to_string())
