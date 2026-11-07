"""site/states.json from states/out (states.py)."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
O = '/data/dev/rnaseq-exploration/states/out/'; E0 = '/data/dev/rnaseq-exploration/'
KEY = ['lorbin_184', 'metabat_9342', 'metabat_5432', 'comebin_3787_sub', 'comebin_2576', 'semibin_088', 'vamb_549', 'metabat_5376', 'vamb_4149', 'lorbin_095', 'metabat_9833']
def r(x, n=3):
    try: x = float(x)
    except (TypeError, ValueError): return None
    return round(x, n) if np.isfinite(x) else None
R = pd.read_csv(O + 'mag_states.tsv', sep='\t'); F = pd.read_csv(O + 'functions.tsv', sep='\t')
CP = pd.read_csv(O + 'changepoints.tsv', sep='\t'); RA = pd.read_csv(O + 'rarefy.tsv', sep='\t')
SC = pd.read_pickle(O + 'scores.pkl'); LO = pd.read_pickle(O + 'loadings.pkl'); GR = pd.read_pickle(O + 'growth.pkl')
S = samples(); m = meta().loc[S]
core = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
sidx = {s: i for i, s in enumerate(core.index)}
J = {'samples': [{'id': s, 'name': core.sample_name[s], 'season': core.season[s], 'pool': core.group[s], 'medium': core.medium[s],
                  'layer': core.ice_layer[s] if isinstance(core.ice_layer[s], str) else '', 'day': int(core.day[s]), 'phase': core.phase[s]} for s in core.index]}
R['short'] = R.bin.str.replace('dastool-', '', regex=False); R['key'] = R.short.isin(KEY)
R['rank'] = R['R2_2025-26'].fillna(-1)
J['n_pass'] = int(len(R)); J['min_reads'] = 3000; J['min_samples'] = 20
detail = set(R.sort_values('rank', ascending=False).head(60).bin) | set(R.bin[R.key])
J['mags'] = []
for x in R.sort_values('rank', ascending=False).itertuples():
    J['mags'].append({'bin': x.short, 'taxon': x.taxon, 'n': int(x.n), 'genes': int(x.genes), 'ev': [r(x.pc1), r(x.pc2), r(x.pc3)], 'rd1': r(x.r_depth_pc1, 2), 'rd2': r(x.r_depth_pc2, 2),
        'gmm_k': int(x.gmm_k), 'dbic': r(x.gmm_dbic, 1), 'km_k': int(x.km_k), 'sil': r(x.sil, 2),
        'key': bool(x.key), 'detail': x.bin in detail})
# dashed column names are not attributes: fill from the frame directly
cols = {'R2_25': 'R2_2025-26', 'q_25': 'q_2025-26', 'R2_24': 'R2_2024-25', 'q_24': 'q_2024-25', 'R2ref_25': 'R2ref_2025-26', 'qref_25': 'qref_2025-26',
        'shift_25': 'shift_2025-26', 'shift_24': 'shift_2024-25', 'sep_25': 'sep_2025-26', 'sep_24': 'sep_2024-25',
        'R2ice_25': 'R2ice_2025-26', 'qice_25': 'qice_2025-26', 'R2ice_24': 'R2ice_2024-25', 'qice_24': 'qice_2024-25'}
Rs = R.sort_values('rank', ascending=False).reset_index(drop=True)
for i, mg in enumerate(J['mags']):
    for k, c in cols.items(): mg[k] = r(Rs.loc[i, c], 4 if k.startswith('q') else 3) if c in Rs else None
J['detail'] = {}
for b in detail:
    s = SC[b]; sh = b.replace('dastool-', '')
    f = F[F.bin == b]
    J['detail'][sh] = {'idx': [sidx[x] for x in s['samples']], 'pc': [[r(v, 2) for v in row] for row in s['pc']], 'ev': [r(v) for v in s['ev']],
        'cl': [int(c) for c in s['cluster']], 'depth': [int(v) for v in s['depth']],
        'growth': [r(GR.get(b, {}).get(x), 5) for x in s['samples']],
        'func': [{'cat': z.category, 'n': int(z.n_genes), 'share': r(z.mean_share, 5), 'r1': r(z.r_pc1, 2), 'r2': r(z.r_pc2, 2)} for z in f.itertuples()],
        'load': {str(j): {k: [{'g': g, 'p': p, 'w': r(w, 3)} for g, p, w in LO[b][j][k][:8]] for k in ('neg', 'pos')} for j in (0, 1)},
        'cp': [{'season': z.season, 'pool': z.pool, 'medium': z.medium, 'day': r(z.day, 1), 'eff': r(z.effect, 2)} for z in CP[CP.bin == b].itertuples()]}
RA['short'] = RA.bin.str.replace('dastool-', '', regex=False); rmap = RA.set_index('short').r_pc1
for mg in J['mags']: mg['rr1'] = r(rmap.get(mg['bin']), 2)
J['rarefy'] = [{'bin': z.bin.replace('dastool-', ''), 'depth': int(z.depth), 'r1': r(z.r_pc1, 2), 'r2': r(z.r_pc2, 2), 'R2': r(getattr(z, 'R2_2025_26_rare', np.nan)), 'p': r(getattr(z, 'p_2025_26_rare', np.nan))} for z in RA[RA.short.isin(KEY)].itertuples()]
AS = pd.read_csv(O + 'assoc.tsv', sep='\t'); AS['short'] = AS.bin.str.replace('dastool-', '', regex=False); AS = AS.set_index('short')
TR = pd.read_csv(O + 'trends.tsv', sep='\t')
LAB = {'R2_season': 'winter', 'R2_group': 'pool', 'R2_medium': 'ice vs water'}
for mg in J['mags']:
    a = AS.loc[mg['bin']]; v = {k: a.get(k) for k in LAB if pd.notna(a.get(k))}
    mg['axis'] = LAB[max(v, key=v.get)] if v else None; mg['axis_r2'] = r(max(v.values()), 2) if v else None
    mg['R2_season'] = r(a.get('R2_season'), 2); mg['R2_pool'] = r(a.get('R2_group'), 2); mg['R2_medium'] = r(a.get('R2_medium'), 2)
for b, d in J['detail'].items():
    t = TR[TR.bin == 'dastool-' + b]
    d['trend'] = [{'season': z.season, 'pool': z.pool, 'n': int(z.n), 'r1': r(z.r_pc1_day, 2), 'r2': r(z.r_pc2_day, 2), 'rg': r(z.r_growth_day, 2),
                   'g0': r(z.growth_first, 1), 'g1': r(z.growth_last, 1)} for z in t.itertuples()]
cnt = pd.Series([mg['axis'] for mg in J['mags']]).value_counts()
J['axis_counts'] = {k: int(v) for k, v in cnt.items()}
J['ice_sig'] = {s: [int((R[f'qice_{s}'] < 0.05).sum()), int(R[f'qice_{s}'].notna().sum())] for s in ('2024-25', '2025-26')}
J['depth_flag'] = int((R.r_depth_pc1.abs() > 0.6).sum()); J['thin_flag'] = int((RA.r_pc1 < 0.7).sum())
Fd = F.assign(short=F.bin.str.replace('dastool-', '', regex=False)).set_index(['short', 'category'])
def fr(b, c, k='r_pc1'):
    try: return f"{Fd.loc[(b, c), k]:+.2f}"
    except KeyError: return '–'
def mg(b): return next(x for x in J['mags'] if x['bin'] == b)
def tr(b, se, P):
    t = TR[(TR.bin == 'dastool-' + b) & (TR.season == se) & (TR.pool == P)]
    return t.iloc[0] if len(t) else None
def g(t): return f"{t.growth_first:.0f}→{t.growth_last:.0f}‰"
S88, C184, C9342, P2576, S5376, U4149 = 'semibin_088', 'lorbin_184', 'metabat_9342', 'comebin_2576', 'metabat_5376', 'vamb_4149'
ac = J['axis_counts']
J['reading'] = [
 f"<b>Ice versus water is the biggest switch.</b> For {ac.get('ice vs water', 0)} of {len(J['mags'])} genomes the main axis of within-genome variation separates ice from water, more than winter ({ac.get('winter', 0)}) or pool ({ac.get('pool', 0)}). After dosing, ice and water states differ significantly in {J['ice_sig']['2024-25'][0]} of {J['ice_sig']['2024-25'][1]} testable genomes in 2024-25 and {J['ice_sig']['2025-26'][0]} of {J['ice_sig']['2025-26'][1]} in 2025-26. The same population behaves differently in brine than in the water below.",
 f"<b>Sulfurimonas differs between pools, but not because of diesel.</b> semibin_088 has one of the strongest pool separations (R² {mg(S88)['R2_25']:.2f}), yet the diesel-pool water moved less from its own pre-dosing state than the control did (shift {mg(S88)['shift_25']:.1f}×). The pools carry different strains of this genome (see Sulfur &amp; hydrogen), which this method cannot tell apart from physiology. Along PC1 it trades ribosomal proteins (r {fr(S88, 'ribosomal proteins (growth)')}) against motility and chemotaxis (r {fr(S88, 'motility / chemotaxis')}), and it slides from the growth end to the motile end through the 2025-26 season in both pools (growth index {g(tr(S88, '2025-26', 'EXP'))} diesel, {g(tr(S88, '2025-26', 'CTL'))} control).",
 f"<b>Cycloclasticus lorbin_184 winds down its hydrocarbon genes as its bloom ages.</b> Its PC1 rises with time in the 2025-26 diesel pool (r {tr(C184, '2025-26', 'EXP').r_pc1_day:+.2f} with day). Moving up PC1, hydrocarbon-degradation genes (r {fr(C184, 'hydrocarbon degradation')}) and storage (r {fr(C184, 'storage (PHA, polyP, glycogen)')}) lose share and chaperones and proteases gain it (r {fr(C184, 'chaperones / proteases')}). Its growth index falls ({g(tr(C184, '2025-26', 'EXP'))}). This is the pattern of a population running down its substrate and coping with protein stress, not of a fresh bloom.",
 f"<b>The 2024-25 Cycloclasticus (metabat_9342) moves between growth and scavenging.</b> Its PC1 sets ribosomal proteins (r {fr(C9342, 'ribosomal proteins (growth)')}) and cold-shock proteins (r {fr(C9342, 'cold shock')}) against motility (r {fr(C9342, 'motility / chemotaxis')}), transporters (r {fr(C9342, 'transporters')}), storage (r {fr(C9342, 'storage (PHA, polyP, glycogen)')}) and nitrogen acquisition (r {fr(C9342, 'N acquisition')}): a fast-growth state versus a foraging, nutrient-limited state.",
 f"<b>50-400-T64 is still in growth mode at the end.</b> Its growth index rises slightly over the 2025-26 diesel season ({g(tr(P2576, '2025-26', 'EXP'))}), consistent with its late bloom (see the 50-400-T64 page). Its main axis pairs chaperones (r {fr(P2576, 'chaperones / proteases')}), oxidative-stress (r {fr(P2576, 'oxidative stress (ROS)')}) and cold-shock genes against transporters (r {fr(P2576, 'transporters')}) and nitrogen acquisition (r {fr(P2576, 'N acquisition')}).",
 f"<b>Background heterotrophs slow down in the diesel pool.</b> Two abundant genomes that do not eat diesel lose growth-index share in the diesel pool but not in the control in 2025-26: Sulfitobacter metabat_5376 ({g(tr(S5376, '2025-26', 'EXP'))} diesel, {g(tr(S5376, '2025-26', 'CTL'))} control; shift {mg(S5376)['shift_25']:.1f}×, R² {mg(S5376)['R2_25']:.2f}, {('q < 0.001' if mg(S5376)['q_25'] < 0.001 else 'q = ' + format(mg(S5376)['q_25'], '.2g'))}) and UBA9145 vamb_4149 ({g(tr(U4149, '2025-26', 'EXP'))} diesel, {g(tr(U4149, '2025-26', 'CTL'))} control). In UBA9145, phage-shock (psp) genes run opposite to PC1 (r {fr(U4149, 'phage shock / envelope stress')}), the same envelope-stress response seen on the Metals &amp; toxicity page.",
 f"<b>What this cannot do.</b> Genomes that are abundant in only one pool (the blooming Cycloclasticus and 50-400-T64) cannot be compared between pools, because the control has too few reads on them. A MAG can hide several strains, and different strains in different pools or winters look like different states. Nitrosopumilus had too few reads for a profile. One pool per treatment, so no test separates diesel from other pool differences.",
]
json.dump(J, open(E0 + 'site/states.json', 'w'), separators=(',', ':'), allow_nan=False)
print(len(J['mags']), 'mags,', len(J['detail']), 'with detail,', os.path.getsize(E0 + 'site/states.json') // 1000, 'kB')
