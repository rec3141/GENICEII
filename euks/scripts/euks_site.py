"""site/euks.json from euks/out (euks.py) + plastid bin TPM."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr
O = '/data/dev/rnaseq-exploration/euks/out/'; E0 = '/data/dev/rnaseq-exploration/'
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
GS = pd.read_csv(O + 'group_share.tsv', sep='\t', index_col=0); GEN = pd.read_csv(O + 'genus_share.tsv', sep='\t', index_col=0)
gsum = pd.read_csv(O + 'group_summary.tsv', sep='\t', index_col=0); nsum = pd.read_csv(O + 'genus_summary.tsv', sep='\t', index_col=0)
st = pd.read_csv(O + 'sample_stats.tsv', sep='\t', index_col=0); ssum = pd.read_csv(O + 'stat_summary.tsv', sep='\t', index_col=0)
VL = pd.read_csv(O + 'virus_links.tsv', sep='\t')
S = list(GS.columns); m = meta().loc[S]
def row(x, name, extra=None):
    d = {'name': name, 'mean': r(x['mean'], 2)}
    for s, k in (('2024-25', '24'), ('2025-26', '25')):
        for f in ('oil', 'ctl', 'pre_exp', 'pre_ctl'): d[f'{f}{k}'] = r(x[f'{f}_{s}'], 2)
        d[f'lfc{k}'] = r(x[f'lfc_{s}']); d[f'baci{k}'] = r(x[f'baci_{s}']); d[f'ice_lfc{k}'] = r(x[f'lfc_ice_{s}'])
    d.update({'ice': r(x.ice), 'ice_up': int(x.ice_up), 'ice_n': int(x.ice_n), 'top': r(x.top), 'top_up': int(x.top_up), 'top_n': int(x.top_n)})
    if extra: d.update(extra)
    return d
J = {'samples': [{'id': s, 'name': m.sample_name[s], 'season': m.season[s], 'pool': m.group[s], 'medium': m.medium[s],
                  'day': None if pd.isna(m.day[s]) else int(m.day[s])} for s in S]}
J['groups'] = [row(gsum.loc[n], n, {'role': gsum.role[n]}) for n in gsum.index]
J['group_tc'] = {n: [r(v, 3) for v in GS.loc[n, S]] for n in GS.index}
top = nsum[nsum['mean'] >= 0.3].index
J['genera'] = [row(nsum.loc[n], n.split(' [')[0], {'group': n.split(' [')[1].rstrip(']')}) for n in top]
J['genus_tc'] = {n.split(' [')[0]: [r(v, 3) for v in GEN.loc[n, S]] for n in top}
LAB = {'euk_rrna_pct': 'Eukaryotic share of rRNA reads (%)', 'plastid_mrna_pct': 'Plastid-encoded share of mRNA (Kaiju “cyanobacterial”, %)',
       'plastid_mag_tpm': 'Plastid genomes (MAG TPM)', 'giant_virus_cpm': 'Giant-virus transcripts (Megaviricetes, CPM)', 'virophage_cpm': 'Virophage contig_0049668 (CPM)',
       'grazer_pct': 'Grazers, % of eukaryotic 18S', 'phototroph_pct': 'Strict phototrophs, % of eukaryotic 18S'}
J['stats'] = [row(ssum.loc[k], LAB[k], {'key': k}) for k in ssum.index]
J['stat_tc'] = {k: [r(v, 3) for v in st.loc[S, k]] for k in LAB}
# plastid bins
tpm = table('bin_tpm'); core = m[m.group.isin(['CTL', 'EXP']) & (m.medium == 'SW')]
pl = []
for b in tpm.index:
    l = special_label(b) if is_plastid(b) else ''
    if not l: continue
    d = {'bin': b.replace('dastool-', ''), 'label': l}
    for s, k in (('2024-25', '24'), ('2025-26', '25')):
        for f, sel in (('oil', core.phase == 'oiled'), ('ctl', core.phase == 'control'), ('pre_exp', (core.phase == 'pre') & (core.group == 'EXP')), ('pre_ctl', (core.phase == 'pre') & (core.group == 'CTL'))):
            d[f'{f}{k}'] = r(tpm.loc[b, core.index[(core.season == s) & sel]].median(), 0)
    d['mean'] = r(tpm.loc[b, S].mean(), 0)
    if d['mean'] >= 50: pl.append(d)
J['plastids'] = sorted(pl, key=lambda d: -d['mean'])
# eukaryotic MAGs: hull-curated EE## (euks/mags: survey -> compact -> final -> expr -> finalize -> kernel -> da -> register_blobs give
# the E## source groups; hulls drawn by eye on the contig map, hulls.tsv -> hull_mags.py -> hull_expr.py), copied into data/ by
# scripts/bin_domains.py. Expression = reads on the MAG's contigs per million contig-assigned reads (any read on the contig, so
# intron-containing genes count), from euks/mags/out/hull_mags_expr_samples.tsv; DA from hull_mags_da.tsv (dacore).
HM = pd.read_csv(D + 'hull_mags.tsv', sep='\t', keep_default_na=False)
EG = HM[HM.mag.str.startswith('EE')]
GX = hull_mag_expression()       # EE/P x libraries (common.py)
DA = hull_mag_da()               # per-MAG RNA/DNA differential abundance
BD = pd.read_csv(D + 'bin_domain.tsv', sep='\t', index_col=0)
dbin = {g: [b.replace('dastool-', '') for b in BD.index[BD.get('hull_mag', pd.Series(dtype=str)).reindex(BD.index) == g]] for g in EG.mag}
core_all = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
def fnum(v, n=2):
    return r(float(v), n) if str(v) not in ('', 'nan') else None
eg = []
for x in EG.itertuples():
    d = {'id': x.mag, 'label': hull_mag_label(x.mag), 'status': x.status, 'network': str(x.network_ok) == 'True', 'mbp': r(x.bp / 1e6, 1),
         'contigs': int(x.contigs), 'gc': r(x.gc, 1), 'tiara': fnum(x.tiara_euk), 'who': fnum(x.who_euk), 'rms': r(x.tsne_rms, 1), 'hint': x.taxonomy_hint, 'tax': x.taxonomy, 'tax_rank': x.taxonomy_rank, 'tax_conf': x.taxonomy_confidence,
         'src': x.source_group, 'gained': r(x.bp_gained / 1e6, 1), 'lost': fnum(float(x.bp_lost) / 1e6 if str(x.bp_lost) not in ('', 'nan') else '', 1),
         'caution': x.caution, 'dastool': dbin.get(x.mag, [])}
    a = DA[DA.group == x.mag]
    if len(a): d['da'] = {f'{y.data}|{y.contrast}': [r(y.log2_effect, 2), r(y.q, 4)] for y in a.itertuples()}
    if x.mag in GX.index:
        v = GX.loc[x.mag].reindex(S)
        for s_, k in (('2024-25', '24'), ('2025-26', '25')):
            for f, sel in (('oil', core.phase == 'oiled'), ('ctl', core.phase == 'control'), ('pre_exp', (core.phase == 'pre') & (core.group == 'EXP')), ('pre_ctl', (core.phase == 'pre') & (core.group == 'CTL'))):
                d[f'{f}{k}'] = r(v[core.index[(core.season == s_) & sel]].median(), 0)
        d['median'] = r(v[core_all.index].median(), 0)
        d['expr'] = [r(t, 1) for t in v]
    eg.append(d)
order = {'network': 0, '<4 Mbp': 1, 'little CMO RNA': 2, 'doubtful': 3}
J['euk_groups'] = sorted(eg, key=lambda d: (order.get(d['status'], 3), -(d.get('median') or 0)))
J['da_contrasts'] = [['RNA', 'diesel_vs_control_2024-25', 'RNA diesel ÷ control 24-25'], ['RNA', 'diesel_vs_control_2025-26', 'RNA diesel ÷ control 25-26'],
                     ['RNA', 'baci_2025-26_SW', 'RNA BACI 25-26'], ['RNA', 'ice_vs_water', 'RNA ice ÷ water'],
                     ['DNA', 'diesel_vs_control_2024-25', 'DNA diesel ÷ control 24-25'], ['DNA', 'ice_vs_water', 'DNA ice ÷ water']]
J['euk_sources'] = {'groups': int(pd.read_csv(D + 'euk_groups.tsv', sep='\t').shape[0]), 'split': sorted(EG.mag[EG.source_group == 'E21'])}
# contig-map figure with the hull outlines (euks/mags/scripts/plot_hulls.py) -> site/img, palette-compressed for the web
from PIL import Image
os.makedirs(E0 + 'site/img', exist_ok=True)
im = Image.open(E0 + 'euks/mags/out/euk_hull_mags_tsne.png').convert('RGB')
im.quantize(colors=128, method=Image.Quantize.MEDIANCUT).save(E0 + 'site/img/euk_hull_mags_tsne.png', optimize=True)
if os.path.getsize(E0 + 'site/img/euk_hull_mags_tsne.png') > 1e6:
    im.resize((im.width * 3 // 4, im.height * 3 // 4), Image.LANCZOS).quantize(colors=128, method=Image.Quantize.MEDIANCUT).save(E0 + 'site/img/euk_hull_mags_tsne.png', optimize=True)
print('figure', im.size, os.path.getsize(E0 + 'site/img/euk_hull_mags_tsne.png') // 1000, 'kB')
J['euk_missed'] = {'tiara_euk_5kb_mbp': 883, 'in_groups_mbp': round(EG.bp.sum() / 1e6)}   # EE MAGs (all contigs >= 1 kb)
# virus links: top 6 positive per season x virus (q < 0.05)
vl = []
for (s, v), x in VL.groupby(['season', 'virus']):
    x = x[(x.q < 0.05) & (x.rho > 0)].sort_values('rho', ascending=False).head(8)
    vl += [{'season': s, 'virus': v, 'taxon': t.split(' [')[0], 'group': t.split(' [')[1].rstrip(']') if '[' in t else t, 'level': lv, 'rho': r(rh), 'q': r(q, 4)}
           for t, lv, rh, q in zip(x.taxon, x.level, x.rho, x.q)]
J['virus_links'] = vl
# grazers vs the Cycloclasticus bloom (2025-26 diesel pool samples, SW+SI)
gz = []
ix = m.index[(m.season == '2025-26') & (m.group == 'EXP') & m.medium.isin(['SW', 'SI'])]
for n in GS.index:
    rho, p = spearmanr(GS.loc[n, ix], st.loc[ix, 'cycloclasticus_tpm'])
    gz.append({'group': n, 'role': gsum.role[n], 'rho': r(rho), 'p': r(p, 4)})
J['bloom_links'] = sorted(gz, key=lambda d: -(d['rho'] or -9))
json.dump(J, open(E0 + 'site/euks.json', 'w'), separators=(',', ':'), allow_nan=False)
print(os.path.getsize(E0 + 'site/euks.json') // 1000, 'kB;', len(J['genera']), 'genera;', len(J['plastids']), 'plastids;', len(J['euk_groups']), 'EE MAGs')
print(pd.DataFrame(J['bloom_links']).to_string()); print(pd.DataFrame(J['plastids'])[['bin', 'label', 'mean', 'pre_exp25', 'pre_ctl25', 'oil25', 'ctl25', 'oil24', 'ctl24']].to_string())
