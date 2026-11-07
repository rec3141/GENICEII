"""site/light.json from light/out (light.py, par_model.py): categories with hosts, time courses, modelled clear-sky PAR (daily and per
sample), within-genome shares for AAP/rhodopsin/plastid genomes, plastid photosynthesis share of plastid transcripts,
metabat_1122 (Sulfitobacter) AAP share vs its PCA state and N acquisition, AAP vs rhodopsin partitioning among heterotrophs."""
import sys, json, os, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr, mannwhitneyu
O = '/data/dev/rnaseq-exploration/light/out/'; E0 = '/data/dev/rnaseq-exploration/'
def r(x, n=2):
    try: x = float(x)
    except (TypeError, ValueError): return None
    return round(x, n) if np.isfinite(x) else None
def text(J):
    C = {c['name'].split(' (')[0]: c for c in J['categories']}
    ps2, rc, lh, pr_, blu, pps, phr = C['photosystem II'], C['reaction centre'], C['light-harvesting LH1'], C['proton-pumping rhodopsins'], C['BLUF blue-light photoreceptors'], C['PpsR/AppA photosynthesis regulators'], C['photolyases and cryptochromes']
    Pd = pd.DataFrame(J['par_daily']); sm = pd.DataFrame(J['samples']); pl = pd.Series(J['plastid']['share']); tp = pd.Series(J['plastid']['tpm'])
    post = lambda se, pool, md: (sm.season == se) & (sm.pool == pool) & (sm.medium == md) & (sm.day >= (5.5 if se == '2024-25' else 1))
    pre = lambda se, pool: (sm.season == se) & (sm.pool == pool) & (sm.medium == 'SW') & (sm.day < 0)
    q = lambda v: f"{v:,.0f}"
    pe, pc_ = pl[post('2025-26', 'EXP', 'SW')].median(), pl[post('2025-26', 'CTL', 'SW')].median()
    te, tc, tpe, tpc = tp[post('2025-26', 'EXP', 'SW')].median(), tp[post('2025-26', 'CTL', 'SW')].median(), tp[pre('2025-26', 'EXP')].median(), tp[pre('2025-26', 'CTL')].median()
    G_ = pd.DataFrame(J['genomes'])
    def wbaci(k):
        g = G_[(G_.kind == k)].dropna(subset=['oil_2025-26', 'ctl_2025-26', 'pre_exp_2025-26', 'pre_ctl_2025-26'])
        g = g[(g[['oil_2025-26', 'ctl_2025-26', 'pre_exp_2025-26', 'pre_ctl_2025-26']] > 0).all(1)]
        b = np.log2(g['oil_2025-26'] / g['pre_exp_2025-26']) - np.log2(g['ctl_2025-26'] / g['pre_ctl_2025-26']); return b.median(), int((b > 0).sum()), len(b)
    ab, au, an = wbaci('AAP'); rb, ru, rn = wbaci('rhodopsin')
    M = J['m1122']; part = pd.DataFrame(J['partition'])
    rcA = part[(part.cat == 'reaction centre') & (part.klass == 'Alphaproteobacteria')].share.sum(); rcG = part[(part.cat == 'reaction centre') & (part.klass == 'Gammaproteobacteria')].share.sum()
    prA = part[(part.cat == 'proton-pumping rhodopsins') & (part.klass == 'Alphaproteobacteria')].share.sum(); prAc = part[(part.cat == 'proton-pumping rhodopsins') & (part.klass.isin(['Actinomycetes', 'Acidimicrobiia']))].share.sum()
    prF = part[(part.cat == 'proton-pumping rhodopsins') & (part.klass.isin(['Flavobacteriia', 'Bacteroidia']))].share.sum()
    dl = Pd.dl.min(); parmin, parmax = Pd.par.min(), Pd[(Pd.season == '2025-26') & (Pd.day <= 90)].par.max()
    f = lambda v: f"{v:+.1f}".replace('-', '−')
    g2 = lambda v: f"{v:+.2f}".replace('-', '−')
    T = {}
    T['kicker'] = ("Even under ice in a Churchill winter, light is a resource. Algae photosynthesise with chlorophyll in their plastids. Many heterotrophic bacteria top up their energy with light as well, either with bacteriochlorophyll reaction centres (aerobic anoxygenic photosynthesis, AAP, typical of roseobacters) or with rhodopsin proton pumps. "
        "This page follows 17 light-gene groups through both winters against a modelled light curve. No light was measured, so the model gives the clear-sky maximum above the ice. Expression is CPM of protein-coding reads, with hosts from MAGs or, for unbinned contigs, Kaiju; claims about regulation use each genome’s own transcripts. One pool per treatment, so all contrasts are descriptive.")
    T['stats'] = [[f"{f(ps2['lfc25'])}", f"log₂ diesel ÷ control, photosystem II transcripts in 2025-26 seawater (BACI {f(ps2['baci25'])}): the algae are lost"],
        [f"{pe:.0f} vs {pc_:.0f}‰", "photosynthesis-gene share of plastid transcripts, diesel vs control 2025-26: the plastids that remain run photosynthesis as before"],
        [f"{f(rc['lfc25'])}", f"log₂ diesel ÷ control, AAP reaction-centre transcripts 2025-26 (BACI {f(rc['baci25'])}); within genomes BACI {f(ab)}: composition, not induction"],
        [f"{blu['top_up']}/{blu['top_n']}", "top-ice vs bottom-ice pairs with more BLUF blue-light photoreceptor transcripts at the top"]]
    T['sky'] = (f"Churchill (58.8° N) never has polar night: on the winter solstice the sun is up for {dl:.1f} hours and reaches about 8° above the horizon. The model gives clear-sky daily PAR from {parmin:.1f} mol photons m⁻² d⁻¹ at the solstice to about {parmax:.0f} by mid-February, a roughly {parmax / parmin:.0f}-fold rise over each experiment. "
        "Both experiments start near the darkest point, so every sample series also runs along a rising light curve. The coloured curves estimate how much of it reached the water under each pool’s snow and ice.")
    ur = J['underice_range']
    T['underice'] = (f"Under the ice, a simple model (no roof, cloud or oil) attenuates the clear-sky light by each pool’s measured snow and ice: transmitted fraction = (1 − albedo) × exp(−15 m⁻¹ × snow − 1.5 m⁻¹ × ice), with albedo 0.85 for snow, 0.50 for bare ice and 0.07 for open water. "
        f"In 2025-26 both pools start nearly open (up to {ur['2025-26|EXP'][1]:.0%} of surface light transmitted in the diesel pool, {ur['2025-26|CTL'][1]:.0%} in the control), then fall to about 1% once snow settles on 30–80 cm of ice after about day 20, so the water darkens even as daylight returns. "
        f"In 2024-25 the cores start at about 20 cm of ice with a little snow, so only {ur['2024-25|CTL'][0]:.0%}–{ur['2024-25|EXP'][1]:.0%} gets through all season. "
        "Snow and ice are interpolated between core days and held flat before the first core. Within-genome light-gene shares show no consistent link to modelled under-ice light (median r near zero for AAP, rhodopsin and plastid photosynthesis), apart from metabat_1122 below. Vertical lines mark the logged roof state on 2025-26 sampling days (solid = open, dotted = closed); it was not logged in 2024-25.")
    T['algae'] = (f"Plastid photosynthesis transcripts collapse in the 2025-26 diesel pool (photosystem II {f(ps2['lfc25'])} log₂, BACI {f(ps2['baci25'])}). That is the algal loss shown on the Eukaryotes page: plastid MAGs fall to {q(te)} TPM in diesel-pool water against {q(tc)} in the control (both pools {q(tpc)}–{q(tpe)} before dosing). "
        f"Within the plastids, though, the share of transcripts that are photosynthesis genes does not change: {pe:.0f}‰ in the diesel pool vs {pc_:.0f}‰ in the control after dosing. Across 16 plastid and cyanobacterial genomes, the within-genome BACI is 0.0. The surviving algae are not visibly shaded or switched off, at the level of plastid transcription. Photosynthesis genes are lower in top ice than bottom ice ({ps2['top_up']} of {ps2['top_n']} pairs higher at the top), as expected for bottom-ice algae.")
    T['genomes'] = (f"Community AAP transcripts rise with diesel because AAP roseobacters, mainly Sulfitobacter, are more abundant in the diesel pool. Within each genome the AAP share does not rise: median within-genome BACI {g2(ab)} ({au} of {an} genomes up). Rhodopsins are flat as well ({g2(rb)}, {ru} of {rn} up). "
        "So diesel did not make the heterotrophs lean harder on light. Correlations with modelled PAR inside each pool are near zero (median r −0.1 to +0.1), so no seasonal switch-on shows up either. Most AAP and rhodopsin genomes express less of their light genes in ice than in the water below.")
    T['partition'] = (f"AAP is almost entirely alphaproteobacterial ({rcA:.0%} of reaction-centre transcripts; roseobacters such as Sulfitobacter, Octadecabacter, Planktomarina and RS62), with gammaproteobacteria at {rcG:.0%}. "
        f"Proton-pumping rhodopsins are spread more widely: alphaproteobacteria {prA:.0%} (mostly SAR11-type), actinobacteria {prAc:.0%}, Bacteroidota {prF:.0%}. Rhodopsins carry about ten times as many transcripts as AAP reaction centres overall.")
    T['m1122'] = (f"On the Physiological states page this Sulfitobacter’s water state loads PufH next to the ammonium transporter. The AAP share confirms it: AAP genes are {M['sw_med']:.1f}‰ of its transcripts in water against {M['si_med']:.1f}‰ in ice (Mann–Whitney p = {M['p_sw_si']:.4f}). "
        f"Within seawater, AAP share runs with the nitrogen-hungry end of PC1 (Spearman r {g2(M['r_pc1_sw'])} with PC1) and with its own nitrogen-acquisition share (r {g2(M['r_nacq_sw'])}). Against the under-ice light model it runs the other way: in 2025-26 water the AAP share falls as modelled under-ice light rises (r {g2(M['r_under']['2025-26|EXP'])} in the diesel pool, {g2(M['r_under']['2025-26|CTL'])} in the control) and climbs as snow and ice thicken through the season (r with day {g2(M['r_day']['2025-26|EXP'])}). That fits roseobacter AAP biology, where bacteriochlorophyll synthesis is repressed by light and made in the dark. Season, darkening and nitrogen hunger still move together, so these cannot be separated.")
    T['reading'] = [
        f"<b>Observed: diesel removed algae, not their photosynthesis.</b> Plastid photosystem transcripts fall about {2 ** -ps2['lfc25']:.0f}× in 2025-26 diesel seawater, in step with plastid abundance. The surviving plastids keep the same photosynthesis share ({pe:.0f} vs {pc_:.0f}‰). Shading by oil, if any, does not show in plastid transcription.",
        f"<b>Observed: light-harvesting bacteria did not switch on with diesel.</b> AAP transcripts rise ({f(rc['lfc25'])} log₂) only because roseobacters bloom; within genomes AAP and rhodopsin shares are flat (BACI {g2(ab)} and {g2(rb)}). Light is not obviously used as a stress subsidy here.",
        f"<b>Observed: the top of the ice senses light.</b> BLUF blue-light photoreceptors ({blu['top_up']}/{blu['top_n']} pairs), PpsR/AppA (the light- and redox-sensing repressor of AAP genes; {pps['top_up']}/{pps['top_n']}) and photolyases ({phr['top_up']}/{phr['top_n']}) are higher in top ice than bottom ice. This fits more light, and more UV damage to repair, near the ice surface. Composition differs between layers, so this is not a within-genome result.",
        "<b>Observed, in one genome: more AAP as the water darkens.</b> Sulfitobacter metabat_1122 expresses more AAP genes in water than in ice, most in its nitrogen-hungry state, and more as snow and ice cut the modelled under-ice light through 2025-26. Making light-harvesting machinery in dim light fits roseobacter AAP regulation; light energy could also stretch budgets while nitrogen limits growth. Season, darkening and nitrogen hunger change together.",
        f"<b>Model caveat.</b> Surface PAR is a clear-sky maximum above the ice (Haurwitz irradiance; PAR = 0.45 × GHI; 4.57 µmol J⁻¹). It rises steadily through both experiments, so it is confounded with experiment day; correlations with it can only be called consistent with seasonal light return. The under-ice curve adds measured snow and ice but uses textbook attenuation values and ignores the roof, cloud and oil; under-ice PAR sensors in each pool would make light a usable covariate.",
    ]
    return T
summ = pd.read_csv(O + 'cat_summary.tsv', sep='\t', index_col=0); E = pd.read_csv(O + 'cat_cpm.tsv', sep='\t', index_col=0)
G = pd.read_csv(O + 'genes.tsv', sep='\t', keep_default_na=False); M = pd.read_pickle(O + 'gene_cpm.pkl'); mem = pd.read_pickle(O + 'members.pkl')
GS = pd.read_csv(O + 'genome_share.tsv', sep='\t'); SH = pd.read_pickle(O + 'genome_share_samples.pkl')
L = pd.read_csv(O + 'sample_light.tsv', sep='\t', index_col=0); P = pd.read_csv(O + 'par_daily.tsv', sep='\t', parse_dates=['date'])
S = list(E.columns); m = meta().loc[S]
oil = list(m.index[m.phase == 'oiled']); ctl = list(m.index[m.phase == 'control'])
J = {'groups': list(dict.fromkeys(summ.group)), 'categories': [], 'timecourse': {}}
for c in summ.index:
    x = summ.loc[c]
    if x.n_genes == 0: continue
    idx = mem[c]; sub = G.iloc[idx]; mm = M.iloc[idx].values
    hosts = {}
    for w, cols in (('diesel', oil), ('control', ctl), ('all', S)):
        ci = [S.index(s) for s in cols]; v = pd.Series(mm[:, ci].sum(1)).groupby(sub.host.values).sum(); t = v.sum()
        hosts[w] = [{'host': h, 'share': r(s / t, 3)} for h, s in v.sort_values(ascending=False).head(5).items()] if t > 0 else []
    z = lambda *v: all((not np.isfinite(t)) or t < 0.5 for t in v)
    J['categories'].append({'name': c, 'group': x.group, 'n_genes': int(x.n_genes), 'mean': r(x.mean_cpm, 1),
        'lfc24': None if z(x['oil_2024-25_SW'], x['ctl_2024-25_SW']) else r(x['log2_2024-25_SW']),
        'lfc25': None if z(x['oil_2025-26_SW'], x['ctl_2025-26_SW']) else r(x['log2_2025-26_SW']),
        'baci25': None if z(x['oil_2025-26_SW'], x['ctl_2025-26_SW'], x['pre_exp_2025-26'], x['pre_ctl_2025-26']) else r(x['baci_2025-26']),
        'oil24': r(x['oil_2024-25_SW'], 1), 'ctl24': r(x['ctl_2024-25_SW'], 1), 'oil25': r(x['oil_2025-26_SW'], 1), 'ctl25': r(x['ctl_2025-26_SW'], 1),
        'pre_exp25': r(x['pre_exp_2025-26'], 1), 'pre_ctl25': r(x['pre_ctl_2025-26'], 1),
        'ice': r(x.ice_log2), 'ice_up': int(x.ice_up), 'ice_n': int(x.ice_n), 'top': r(x.top_log2), 'top_up': int(x.top_up), 'top_n': int(x.top_n),
        'rpar': {f'{s}|{g}': r(x[f'rpar_{s}_{g}']) for s in ('2024-25', '2025-26') for g in ('EXP', 'CTL')}, 'hosts': hosts})
    J['timecourse'][c] = [r(v, 2) for v in E.loc[c, S]]
J['samples'] = [{'id': s, 'name': m.sample_name[s], 'season': m.season[s], 'pool': m.group[s], 'medium': m.medium[s],
                 'layer': m.ice_layer[s] if isinstance(m.ice_layer[s], str) else '', 'day': None if pd.isna(m.day[s]) else int(m.day[s]),
                 'date': str(L.date[s]), 'par': r(L.par_mol[s], 2), 'dl': r(L.day_length_h[s], 2),
                 'roof': L.roof_open[s] if isinstance(L.roof_open[s], str) else ''} for s in S]
# modelled daily PAR, with experiment day per season (day 0 = 2024-11-28 and 2025-11-20; the dominant collection_date - day offset)
D0 = {'2024-25': pd.Timestamp('2024-11-28'), '2025-26': pd.Timestamp('2025-11-20')}
J['par_daily'] = []
for _, x in P.iterrows():
    se = '2024-25' if x.date.year == 2024 or (x.date.year == 2025 and x.date.month < 6) else '2025-26'
    J['par_daily'].append({'season': se, 'date': str(x.date.date()), 'day': int((x.date - D0[se]).days), 'par': r(x.par_mol, 2), 'dl': r(x.day_length_h, 2), 'toa': r(x.toa_MJ, 2)})
# within-genome shares
gs = []
KEYS = ['oil_2025-26', 'ctl_2025-26', 'pre_exp_2025-26', 'pre_ctl_2025-26', 'oil_2024-25', 'ctl_2024-25', 'ice_log2',
        'rpar_2025-26_EXP', 'rpar_2025-26_CTL', 'rpar_2024-25_EXP', 'rpar_2024-25_CTL']
for _, x in GS.iterrows():
    if x.n_ok < 10: continue
    if x.kind == 'AAP' and re.search(r'Synechococc|Cyanobac|Plastid|Prochloro', str(x.taxon)): continue   # chlorophyll-synthesis homologs, not AAP
    v = SH[(x.kind, x.bin)]
    g = {'kind': x.kind, 'bin': x.bin.replace('dastool-', ''), 'taxon': x.taxon, 'klass': x.klass if isinstance(x.klass, str) else '',
         'n': int(x.n_ok), 'host': r(x.host_mean, 0), 'share': r(x.share_med, 2), 'ice_up': int(x.ice_up), 'ice_n': int(x.ice_n),
         'series': [r(v.get(s), 3) for s in S]}
    g.update({k: r(x[k], 2) for k in KEYS})
    gs.append(g)
J['genomes'] = gs
# plastid photosynthesis share of plastid transcripts (photosynthesis genes on plastid MAGs / all plastid-MAG CDS transcripts)
BIN = pd.read_pickle(E0 + 'nitrogen/out/bin_cds_cpm.pkl')
pl = [b for b in BIN.index if plastid_label(b)]
phot_cats = [c for c in summ.index if summ.loc[c, 'group'] == 'Oxygenic photosynthesis']
pidx = np.unique(np.concatenate([mem[c] for c in phot_cats]))
pmask = G.iloc[pidx].plastid.astype(str).eq('True').values
psum = M.iloc[pidx[pmask]].sum(0); ptot = BIN.loc[pl, S].sum(0)
J['plastid'] = {'share': [r(a / b * 1e3, 1) if b > 50 else None for a, b in zip(psum[S], ptot[S])], 'tpm': [r(v, 1) for v in ptot[S]], 'n_bins': len(pl)}
# metabat_1122
b = 'dastool-metabat_1122'; st = json.load(open(E0 + 'site/states.json')); d = st['detail']['metabat_1122']
pc1 = {st['samples'][i]['id']: d['pc'][k][0] for k, i in enumerate(d['idx'])}
aap = SH.get(('AAP', b)); ng = pd.read_csv(E0 + 'nitrogen/out/genes.tsv', sep='\t', keep_default_na=False); nm = pd.read_pickle(E0 + 'nitrogen/out/gene_cpm.pkl')
nmem = pd.read_pickle(E0 + 'nitrogen/out/members.pkl')
NST = ['amtB (ammonium transporter)', 'glnB/glnK (PII nitrogen regulators)', 'ntrBC (Ntr two-component system)', 'glnA (glutamine synthetase)', 'urtABCDE (urea ABC transporter)', 'nitrate/nitrite transporters (nrtABC, narK, nasA-type)']
ni = np.unique(np.concatenate([nmem[c] for c in NST])); ngb = ng.iloc[ni]; ids = ngb.gene_id[ngb.bin == b]
nshare = nm.loc[ids, S].sum() / BIN.loc[b, S].replace(0, np.nan) * 1e3
rows = []
for s in S:
    if aap is None or pd.isna(aap.get(s)) or s not in pc1: continue
    rows.append({'id': s, 'aap': r(aap[s], 3), 'pc1': r(pc1[s], 2), 'nacq': r(nshare[s], 3)})
R = pd.DataFrame(rows).set_index('id')
mk = m.loc[R.index]
def sp(a, bb): return r(spearmanr(a, bb)[0], 2) if len(a) > 5 else None
J['m1122'] = {'points': rows,
    'r_pc1_all': sp(R.aap, R.pc1), 'r_nacq_all': sp(R.aap, R.nacq),
    'r_pc1_sw': sp(R.aap[mk.medium == 'SW'], R.pc1[mk.medium == 'SW']), 'r_nacq_sw': sp(R.aap[mk.medium == 'SW'], R.nacq[mk.medium == 'SW']),
    'r_par_sw25': sp(R.aap[(mk.medium == 'SW') & (mk.season == '2025-26')], L.par_mol[R.index][(mk.medium == 'SW') & (mk.season == '2025-26')]),
    'sw_med': r(R.aap[mk.medium == 'SW'].median(), 2), 'si_med': r(R.aap[mk.medium == 'SI'].median(), 2),
    'p_sw_si': r(mannwhitneyu(R.aap[mk.medium == 'SW'], R.aap[mk.medium == 'SI'])[1], 4) if (mk.medium == 'SI').sum() > 3 else None}
# partitioning among heterotrophs: AAP vs rhodopsin pumps, by class (transcripts, all samples, CPM sums)
part = []
for c in ('reaction centre (pufL, pufM, pufC, puhA)', 'proton-pumping rhodopsins (proteo-, xantho-, actino-)'):
    idx = mem[c]; sub = G.iloc[idx]; v = M.iloc[idx].sum(1).groupby(sub.klass.replace('', 'unclassified').values).sum()
    for k_, val in v.sort_values(ascending=False).head(8).items(): part.append({'cat': c.split(' (')[0], 'klass': k_, 'cpm_sum': r(val, 0), 'share': r(val / v.sum(), 3)})
J['partition'] = part
# under-ice light model (underice.py) and logged roof state (2025-26 sampling days only)
U = pd.read_csv(O + 'underice_daily.tsv', sep='\t'); UC = pd.read_csv(O + 'underice_corr.tsv', sep='\t')
J['underice'] = [{'season': x.season, 'pool': x.pool, 'day': int(x.day), 'under': r(x.par_under, 3), 'ice': r(x.ice_cm, 1), 'snow': r(x.snow_cm, 1), 'T': r(x.transmittance, 3)} for x in U.itertuples()]
J['underice_range'] = {f'{a}|{b}': [r(g.transmittance.min(), 3), r(g.transmittance.max(), 3), r(g.par_under.min(), 3), r(g.par_under.max(), 2)] for (a, b), g in U.groupby(['season', 'pool'])}
J['underice_corr'] = {f'{k}|{a}|{b}': r(g.r_under.median(), 2) for (k, a, b), g in UC.groupby(['kind', 'season', 'pool'])}
mm = m[m.roof_open.isin(['Yes', 'No']) & m.group.isin(['CTL', 'EXP']) & m.day.notna()]
J['roof'] = [{'season': a, 'pool': b, 'day': int(c), 'open': d == 'Yes'} for (a, b, c, d) in sorted(set(zip(mm.season, mm.group, mm.day, mm.roof_open)))]
x1122 = UC[UC.bin == 'dastool-metabat_1122'].set_index(['season', 'pool'])
J['m1122']['r_under'] = {f'{a}|{b}': r(x1122.loc[(a, b), 'r_under'], 2) for a, b in x1122.index}
J['m1122']['r_day'] = {f'{a}|{b}': r(x1122.loc[(a, b), 'r_day'], 2) for a, b in x1122.index}
J['text'] = text(J)
json.dump(J, open(E0 + 'site/light.json', 'w'), separators=(',', ':'), allow_nan=False)
print(len(J['categories']), 'cats;', len(gs), 'genomes;', os.path.getsize(E0 + 'site/light.json') // 1000, 'kB')
print(json.dumps({k: v for k, v in J['m1122'].items() if k != 'points'}))
print(pd.DataFrame(part).to_string())
