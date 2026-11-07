"""Light use: oxygenic photosynthesis (plastids, cyanobacteria), aerobic anoxygenic photosynthesis (AAP: puf/bch), rhodopsins,
light sensing/protection. Per-sample CPM of CDS reads by gene group; hosts (MAG or Kaiju); diesel vs control, BACI, ice vs water,
top vs bottom; within-genome shares (group CPM / genome CDS CPM) for AAP and rhodopsin genomes; plastid photosynthesis share of
plastid transcripts; modelled clear-sky PAR per sample (light/out/par_daily.tsv) as a seasonal covariate.
Out: light/out/{cat_cpm.tsv, cat_summary.tsv, genes.tsv, gene_cpm.pkl, members.pkl, genome_share.tsv, genome_share_samples.pkl, sample_light.tsv}"""
import sys, os, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
import contrasts as C
from scipy.stats import spearmanr, mannwhitneyu
O = '/data/dev/rnaseq-exploration/light/out/'; E0 = '/data/dev/rnaseq-exploration/'
K = lambda s: s.split()
HELIO = r'[Hh]eliorhodopsin'
PUMP = r'(?<![Hh]elio)rhodopsin|Rhodopsin|Xanthorodopsin'
CATS = [  # group, name, KOs, product regex, host class restriction: None | 'phot' (plastid/cyano only) | 'nonphot' (exclude plastid/cyano)
 ('Oxygenic photosynthesis', 'photosystem II (psbA, psbD)', K('K02703 K02706'), r'^Photosystem II (protein D1|D2 protein)|^Photosystem Q\(B\)', 'phot'),
 ('Oxygenic photosynthesis', 'photosystem I (psaA, psaB)', K('K02689 K02690'), r'^Photosystem I P700 chlorophyll a apoprotein', 'phot'),
 ('Oxygenic photosynthesis', 'cytochrome b6f (petB, petD)', K('K02635 K02637'), None, 'phot'),
 ('Oxygenic photosynthesis', 'phycobilisome (cpc, apc, cpe)', K('K02284 K02285 K02092 K02093 K05376 K05377'), r'[Pp]hycobili|[Pp]hycocyanin|[Pp]hycoerythrin|[Aa]llophycocyanin', 'phot'),
 ('Oxygenic photosynthesis', 'RuBisCO large subunit (rbcL), phototroph hosts', K('K01601'), None, 'phot'),
 ('Aerobic anoxygenic photosynthesis', 'reaction centre (pufL, pufM, pufC, puhA)', K('K08928 K08929 K13991 K13994'), r'[Pp]hotosynthetic reaction center|reaction center subunit [HLM]', 'nonphot'),
 ('Aerobic anoxygenic photosynthesis', 'light-harvesting LH1 (pufA, pufB)', K('K08926 K08927'), r'[Ll]ight-harvesting protein B-?875|LH1', 'nonphot'),
 ('Aerobic anoxygenic photosynthesis', 'bacteriochlorophyll synthesis (bchXYZ, bchC/F/G/P, acsF)', K('K11333 K11334 K11335 K11337 K11336 K04040 K10960 K04035'), r'[Bb]acteriochlorophyll', 'nonphot'),
 ('Aerobic anoxygenic photosynthesis', 'spheroidene carotenoids (crtA, crtC, crtD, crtF)', K('K09845 K09844 K09843 K09846'), None, 'nonphot'),
 ('Rhodopsins', 'proton-pumping rhodopsins (proteo-, xantho-, actino-)', K('K04641'), PUMP, None),
 ('Rhodopsins', 'heliorhodopsins (function unknown)', [], HELIO, None),
 ('Rhodopsins', 'retinal synthesis (blh)', K('K21817'), r'beta-carotene 15,15|\bBlh\b', None),
 ('Light sensing and protection', 'photolyases and cryptochromes (phr)', K('K01669'), r'[Pp]hoto-?lyase|[Cc]ryptochrome', None),
 ('Light sensing and protection', 'BLUF blue-light photoreceptors', [], r'BLUF', None),
 ('Light sensing and protection', 'bacteriophytochromes', K('K11354'), r'[Pp]hytochrome', None),
 ('Light sensing and protection', 'carotenoid synthesis (crtB, crtI, crtY, crtZ)', K('K02291 K10027 K06443 K02294'), None, None),
 ('Light sensing and protection', 'PpsR/AppA photosynthesis regulators', K('K21023'), r'\bPpsR\b|\bAppA\b', 'nonphot'),
]
AAP = [c[1] for c in CATS if c[0] == 'Aerobic anoxygenic photosynthesis']
PUMPS = 'proton-pumping rhodopsins (proteo-, xantho-, actino-)'

X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
kos = G['KO'].fillna('').str.split(','); prod = G['product'].fillna('')
allko = set(k for c in CATS for k in c[2]); rxall = '|'.join(f'(?:{c[3]})' for c in CATS if c[3])
sel = (kos.apply(lambda x: bool(allko.intersection(x))).values | prod.str.contains(rxall, regex=True).values) & cds
idx = np.where(sel)[0]
sub = G.iloc[idx].reset_index(drop=True); sub['bin'] = sub['bin'].fillna('')
M = pd.DataFrame(X[idx].toarray() / tot * 1e6, columns=S)
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def lin(b): return str(gt.get(str(b).replace('dastool-', ''), ''))
def rank(b, r):
    x = [y[3:] for y in lin(b).split(';') if y.startswith(r + '__')]; return x[0] if x else ''
def blab(b):
    if not b: return 'unbinned'
    if special_label(b): return special_label(b)
    t = [y.split('__', 1)[1] for y in lin(b).split(';') if '__' in y and y.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')
need = set(sub.contig[sub.bin == '']); kj = {}
with open(D + 'kaiju_contigs.tsv') as fh:
    next(fh)
    for line in fh:
        c_, st, tid, l = line.rstrip('\n').split('\t')
        if c_ in need: kj[c_] = [x.strip() for x in l.split(';')]
def kgen(c_):
    l = kj.get(c_)
    if not l: return 'unbinned, unclassified'
    g = l[5] if len(l) > 5 and l[5] not in ('', 'NA') else (l[4] if len(l) > 4 else l[-1])
    return f'unbinned ({g}, Kaiju)'
sub['phylum'] = [rank(b, 'p') if b else (kj.get(c_, ['', ''])[1] if len(kj.get(c_, [])) > 1 else '') for b, c_ in zip(sub.bin, sub.contig)]
sub['klass'] = ['Eukaryota (Tiara)' if b and bin_domain(b) == 'eukaryote' else 'Plastid/organelle' if b and is_plastid(b) else rank(b, 'c') if b else (kj.get(c_, ['', '', ''])[2] if len(kj.get(c_, [])) > 2 else '') for b, c_ in zip(sub.bin, sub.contig)]
sub['plastid'] = [is_plastid(b) for b in sub.bin]
bid = sub.bin.str.replace('dastool-', '', regex=False)
sub['taxon'] = [blab(b) for b in sub.bin]
sub['host'] = np.where(sub.bin == '', sub.contig.map(kgen), np.where(sub.taxon == bid, sub.taxon, sub.taxon + ' ' + bid))
# phototroph hosts: plastid MAGs, Cyanobacteriota MAGs (GTDB puts plastids there too), unbinned with Kaiju Cyanobacteria (mostly plastid-derived)
phot = np.array(sub.plastid.values, dtype=bool) | sub.phylum.str.contains('Cyanobacteri', regex=True).values
sub['phot'] = phot
skos = sub['KO'].fillna('').str.split(','); sprod = sub['product'].fillna('')
members, rows = {}, {}
for grp, name, ks, rx, rest in CATS:
    s_ = skos.apply(lambda x: bool(set(ks).intersection(x))).values if ks else np.zeros(len(sub), bool)
    s_ = np.array(s_, dtype=bool, copy=True)
    if rx: s_ = s_ | sprod.str.contains(rx, regex=True).values
    if name == PUMPS: s_ &= ~sprod.str.contains(HELIO, regex=True).values
    if rest == 'phot': s_ &= phot
    if rest == 'nonphot': s_ &= ~phot
    members[name] = np.where(s_)[0]; rows[name] = M[s_].sum(0)
E = pd.DataFrame(rows).T
E.to_csv(O + 'cat_cpm.tsv', sep='\t'); pd.to_pickle(members, O + 'members.pkl')
sub[['gene_id', 'contig', 'product', 'KO', 'bin', 'phylum', 'klass', 'taxon', 'host', 'plastid', 'phot']].to_csv(O + 'genes.tsv', sep='\t', index=False)
Mi = M.copy(); Mi.index = sub.gene_id; Mi.to_pickle(O + 'gene_cpm.pkl')

# modelled light per sample (collection date)
P = pd.read_csv(O + 'par_daily.tsv', sep='\t', parse_dates=['date']).set_index('date')
cd = pd.to_datetime(m.collection_date, errors='coerce').dt.normalize()
L = pd.DataFrame({'date': cd.dt.date, 'par_mol': cd.map(P.par_mol), 'day_length_h': cd.map(P.day_length_h), 'roof_open': m.roof_open,
                  'ice_cm': pd.to_numeric(m.ice_thickness_cm, errors='coerce'), 'season': m.season, 'pool': m.group, 'medium': m.medium,
                  'layer': m.ice_layer, 'day': m.day, 'phase': m.phase}, index=S)
L.to_csv(O + 'sample_light.tsv', sep='\t')

core = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
names = [c[1] for c in CATS]; med = lambda cols: E[cols].median(1)
summ = pd.DataFrame({'group': [c[0] for c in CATS], 'n_genes': [len(members[n]) for n in names]}, index=names)
for s in ('2024-25', '2025-26'):
    for md in ('SW', 'SI'):
        o = core.index[(core.season == s) & (core.medium == md) & (core.phase == 'oiled')]; c = core.index[(core.season == s) & (core.medium == md) & (core.phase == 'control')]
        summ[f'oil_{s}_{md}'] = med(o); summ[f'ctl_{s}_{md}'] = med(c); summ[f'log2_{s}_{md}'] = np.log2((med(o) + .5) / (med(c) + .5))
    for g, k in (('EXP', 'exp'), ('CTL', 'ctl')):
        p = core.index[(core.season == s) & (core.medium == 'SW') & (core.group == g) & (core.phase == 'pre')]; summ[f'pre_{k}_{s}'] = med(p)
    summ[f'baci_{s}'] = np.log2((summ[f'oil_{s}_SW'] + .5) / (summ[f'pre_exp_{s}'] + .5)) - np.log2((summ[f'ctl_{s}_SW'] + .5) / (summ[f'pre_ctl_{s}'] + .5))
Lg = np.log2(E + 0.5)
r1, _ = C.contrast(Lg, C.PI); r2, _ = C.contrast(Lg, C.PT)
summ['ice_log2'] = r1.delta; summ['ice_up'] = r1.n_up; summ['ice_n'] = r1.n; summ['top_log2'] = r2.delta; summ['top_up'] = r2.n_up; summ['top_n'] = r2.n
summ['mean_cpm'] = E.mean(1)
# seasonal light covariate: Spearman of category CPM with modelled PAR, seawater, within each pool x season after dosing (and all SW)
for s in ('2024-25', '2025-26'):
    for g in ('EXP', 'CTL'):
        ix = core.index[(core.season == s) & (core.group == g) & (core.medium == 'SW')]
        summ[f'rpar_{s}_{g}'] = [spearmanr(E.loc[n, ix], L.par_mol[ix])[0] if E.loc[n, ix].std() > 0 else np.nan for n in names]
summ.to_csv(O + 'cat_summary.tsv', sep='\t')

# within-genome shares for light genomes (AAP, rhodopsin pumps, plastids): share of the genome's CDS transcripts (per mille)
BIN = pd.read_pickle(E0 + 'nitrogen/out/bin_cds_cpm.pkl')
gsets = {'AAP': np.unique(np.concatenate([members[n] for n in AAP])), 'rhodopsin': members[PUMPS],
         'photosynthesis (plastid/cyano)': np.unique(np.concatenate([members[n] for n in names[:5]]))}
share_rows, share_samples = [], {}
for kind, ix_ in gsets.items():
    sb = sub.iloc[ix_]; V = M.iloc[ix_].groupby(sb.bin.values).sum(); V = V[V.index != '']
    for b in V.index:
        if b not in BIN.index: continue
        if kind != 'photosynthesis (plastid/cyano)' and not is_prokaryote(b): continue   # AAP/rhodopsin hosts: prokaryotes only
        h = BIN.loc[b, S]; ok = h >= 200                           # genome CPM >= 200 (~ >= 300 CDS reads at median depth)
        if ok.sum() < 10: continue
        sh = (V.loc[b, S] / h.replace(0, np.nan) * 1e3).where(ok)
        share_samples[(kind, b)] = sh
        row = dict(kind=kind, bin=b, taxon=blab(b), phylum=rank(b, 'p'), klass=('Eukaryota (Tiara)' if bin_domain(b) == 'eukaryote' else 'Plastid/organelle' if is_plastid(b) else rank(b, 'c')), n_ok=int(ok.sum()), host_mean=float(h.mean()), share_med=float(sh.median()))
        for s in ('2024-25', '2025-26'):
            for k, sel_ in (('oil', (core.phase == 'oiled')), ('ctl', (core.phase == 'control')), ('pre_exp', (core.phase == 'pre') & (core.group == 'EXP')), ('pre_ctl', (core.phase == 'pre') & (core.group == 'CTL'))):
                cols = core.index[(core.season == s) & (core.medium == 'SW') & sel_]; v = sh[cols].dropna()
                row[f'{k}_{s}'] = float(v.median()) if len(v) >= 2 else np.nan; row[f'n_{k}_{s}'] = len(v)
            for g in ('EXP', 'CTL'):
                cols = core.index[(core.season == s) & (core.group == g) & (core.medium == 'SW')]; v = sh[cols].dropna()
                row[f'rpar_{s}_{g}'] = spearmanr(v, L.par_mol[v.index])[0] if len(v) >= 6 and v.std() > 0 else np.nan
        # ice vs water within genome: paired where both have the genome
        d = []
        for pr_ in C.PI:                                           # a = ice libraries, b = paired water libraries
            va = sh[pr_['a']].dropna(); vb = sh[pr_['b']].dropna()
            if len(va) and len(vb) and va.mean() > 0 and vb.mean() > 0: d.append(np.log2(va.mean() / vb.mean()))
        row['ice_log2'] = float(np.median(d)) if len(d) >= 4 else np.nan; row['ice_n'] = len(d); row['ice_up'] = int(np.sum(np.array(d) > 0))
        share_rows.append(row)
GS = pd.DataFrame(share_rows); GS.to_csv(O + 'genome_share.tsv', sep='\t', index=False); pd.to_pickle(share_samples, O + 'genome_share_samples.pkl')
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(summ[['n_genes', 'mean_cpm', 'log2_2024-25_SW', 'log2_2025-26_SW', 'baci_2025-26', 'ice_log2', 'ice_up', 'ice_n', 'top_log2', 'top_up', 'rpar_2025-26_EXP', 'rpar_2025-26_CTL']].round(2).to_string())
print(GS.sort_values('host_mean', ascending=False).groupby('kind').head(12)[['kind', 'bin', 'taxon', 'n_ok', 'host_mean', 'share_med', 'oil_2025-26', 'ctl_2025-26', 'pre_exp_2025-26', 'pre_ctl_2025-26', 'ice_log2', 'ice_up', 'ice_n', 'rpar_2025-26_EXP', 'rpar_2025-26_CTL']].round(2).to_string())
