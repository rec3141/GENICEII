"""Nitrogen cycling genes: per-sample CPM (of CDS read pairs), hosts, diesel vs control, 2025-26 seawater BACI, ice vs water, top vs bottom.
Gene sets from the coordinate-corrected KOs (genemap.pkl). KOs shared between processes are split by host taxonomy:
narGH (K00370/1) = nitrite oxidoreductase nxrAB in Nitrospinota/Nitrospirota/Nitrococcus; amoABC (K10944-6) = ammonia monooxygenase only in nitrifier hosts, other copper monooxygenases (CuMMO) kept apart.
Host-normalised N-limitation response: for each genome, share of its transcripts in N-acquisition genes, diesel vs control (controls for blooms).
Out: nitrogen/out/{cat_cpm.tsv, cat_summary.tsv, genes.tsv, gene_cpm.pkl, members.pkl, host_nstress.tsv, bin_cds_cpm.pkl}"""
import sys, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
import contrasts as C
O = '/data/dev/rnaseq-exploration/nitrogen/out/'; os.makedirs(O, exist_ok=True)
K = lambda s: s.split()
NITRIFIER_PHYLA = r'Nitrospinota|Nitrospirota'
NITRIFIER_GENERA = r'Nitrococcus|Nitrospina|Nitrospira|Nitrobacter|Nitrotoga|Nitrolancea'
AOB_GENERA = r'Nitrosomonas|Nitrosospira|Nitrosococcus|Nitrosopumilus|Nitrosopelagicus|Nitrosarchaeum|Nitrosomarinus|Nitrososphaera|Nitrosocosmicus'
AMO_HOSTS = r'Nitroso|Nitrosopumil|Nitrososphaer|Nitrosomonad|Nitrosococc|Thermoproteota|Thaumarchaeota|Nitrospira'  # comammox Nitrospira too
CATS = [  # (group, name, KOs)
 ('N2 fixation', 'nifHDK (nitrogenase)', K('K02588 K02586 K02591')),
 ('Nitrification', 'amoABC (ammonia monooxygenase, nitrifier hosts)', K('K10944 K10945 K10946')),
 ('Nitrification', 'CuMMO in other hosts (pmo/amo-like; alkane or methane)', K('K10944 K10945 K10946')),
 ('Nitrification', 'hao (hydroxylamine dehydrogenase)', K('K10535')),
 ('Nitrification', 'nxrAB (nitrite oxidoreductase)', K('K00370 K00371')),
 ('Denitrification', 'narGHI (membrane nitrate reductase)', K('K00370 K00371 K00374')),
 ('Denitrification', 'napAB (periplasmic nitrate reductase)', K('K02567 K02568')),
 ('Denitrification', 'nirK (Cu nitrite reductase)', K('K00368')),
 ('Denitrification', 'nirS (cd1 nitrite reductase)', K('K15864')),
 ('Denitrification', 'norBC (nitric oxide reductase)', K('K04561 K02305')),
 ('Denitrification', 'nosZ (nitrous oxide reductase)', K('K00376')),
 ('DNRA', 'nrfAH (cytochrome c nitrite reductase)', K('K03385 K15876')),
 ('DNRA', 'nirBD (NADH nitrite reductase)', K('K00362 K00363')),
 ('Anammox', 'hzsABC + hdh (hydrazine synthase/dehydrogenase)', K('K20932 K20933 K20934 K20935')),
 ('Nitrate assimilation', 'nasA/narB (assimilatory nitrate reductase)', K('K00372 K00367')),
 ('Nitrate assimilation', 'nirA (ferredoxin nitrite reductase)', K('K00366')),
 ('Nitrate assimilation', 'nitrate/nitrite transporters (nrtABC, narK, nasA-type)', K('K15576 K15577 K15578 K15579 K02575 K02574')),
 ('Ammonium and N regulation', 'amtB (ammonium transporter)', K('K03320')),
 ('Ammonium and N regulation', 'glnB/glnK (PII nitrogen regulators)', K('K04751 K04752')),
 ('Ammonium and N regulation', 'ntrBC (Ntr two-component system)', K('K07708 K07712')),
 ('Ammonium and N regulation', 'glnA (glutamine synthetase)', K('K01915')),
 ('Ammonium and N regulation', 'gltBD (glutamate synthase, GOGAT)', K('K00265 K00266 K00284')),
 ('Ammonium and N regulation', 'gdhA (NADP glutamate dehydrogenase)', K('K00262')),
 ('Organic N', 'ureABC (urease)', K('K01430 K01429 K01428 K14048')),
 ('Organic N', 'urtABCDE (urea ABC transporter)', K('K11959 K11960 K11961 K11962 K11963')),
 ('Organic N', 'cynS (cyanase)', K('K01725')),
 ('Organic N', 'amino acid/peptide uptake (livK, dppA, oppA)', K('K01999 K12368 K15580')),
 ('NO detoxification', 'hmp (flavohemoglobin NO dioxygenase)', K('K05916')),
 ('NO detoxification', 'norVW (flavorubredoxin)', K('K12264 K12265')),
 ('NO detoxification', 'hcp (hybrid cluster protein)', K('K05601')),
]
NSTRESS = ['amtB (ammonium transporter)', 'glnB/glnK (PII nitrogen regulators)', 'ntrBC (Ntr two-component system)', 'glnA (glutamine synthetase)',
           'urtABCDE (urea ABC transporter)', 'nitrate/nitrite transporters (nrtABC, narK, nasA-type)']

X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'
cidx = np.where(cds)[0]
tot = np.asarray(X[cidx].sum(0)).ravel()
kos = G['KO'].fillna('').str.split(',')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def lin(b): return str(gt.get(str(b).replace('dastool-', ''), ''))
def rank(b, r):
    x = [y[3:] for y in lin(b).split(';') if y.startswith(r + '__')]
    return x[0] if x else ''
def blab(b):
    if not isinstance(b, str) or not b: return 'unbinned'
    if plastid_label(b): return plastid_label(b)
    t = [y.split('__', 1)[1] for y in lin(b).split(';') if '__' in y and y.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')

# per-bin total CDS CPM (host normaliser)
bins = G['bin'].fillna('').values
if os.path.exists(O + 'bin_cds_cpm.pkl'): BIN = pd.read_pickle(O + 'bin_cds_cpm.pkl')
else:
    import scipy.sparse as sp
    bc = pd.Series(bins[cidx]); codes, ub = pd.factorize(bc)
    A = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(ub), len(codes)))
    BIN = pd.DataFrame((A @ X[cidx]).toarray() / tot * 1e6, index=ub, columns=S); BIN = BIN[BIN.index != '']
    BIN.to_pickle(O + 'bin_cds_cpm.pkl')

allko = set(k for _, _, ks in CATS for k in ks)
idx = np.where(kos.apply(lambda x: bool(allko.intersection(x))).values & cds)[0]
sub = G.iloc[idx].reset_index(drop=True)
M = pd.DataFrame(X[idx].toarray() / tot * 1e6, columns=S)
sub['bin'] = sub['bin'].fillna('')
sub['lineage'] = [lin(b) if b else '' for b in sub['bin']]
sub['phylum'] = [rank(b, 'p') if b else '' for b in sub['bin']]
sub['taxon'] = [blab(b) for b in sub['bin']]
# Kaiju lineage for unbinned contigs (for nxr/amo split)
need = set(sub.contig[sub.bin == ''])
kj = {}
with open(D + 'kaiju_contigs.tsv') as fh:
    next(fh)
    for line in fh:
        c, st, tid, l = line.rstrip('\n').split('\t')
        if c in need: kj[c] = l
sub['kaiju'] = [kj.get(c, '') if b == '' else '' for c, b in zip(sub.contig, sub.bin)]
fulllin = sub.lineage + ';' + sub.kaiju
is_nob = fulllin.str.contains(NITRIFIER_PHYLA + '|' + NITRIFIER_GENERA, regex=True).values
is_amo = fulllin.str.contains(AMO_HOSTS, regex=True).values
skos = sub['KO'].fillna('').str.split(',')
members, rows, hosts = {}, {}, {}
for grp, name, ks in CATS:
    sel = skos.apply(lambda x: bool(set(ks).intersection(x))).values.copy()
    if name.startswith('nxrAB'): sel &= is_nob
    if name.startswith('narGHI'): sel &= ~is_nob
    if name.startswith('amoABC'): sel &= is_amo
    if name.startswith('CuMMO'): sel &= ~is_amo
    members[name] = np.where(sel)[0]; rows[name] = M[sel].sum(0)
E = pd.DataFrame(rows).T
E.to_csv(O + 'cat_cpm.tsv', sep='\t')
pd.to_pickle(members, O + 'members.pkl')
sub[['gene_id', 'contig', 'product', 'KO', 'bin', 'phylum', 'taxon', 'lineage', 'kaiju']].to_csv(O + 'genes.tsv', sep='\t', index=False)
M.index = sub.gene_id; M.to_pickle(O + 'gene_cpm.pkl')

core = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
med = lambda cols: E[cols].median(1)
names = [n for _, n, _ in CATS]
summ = pd.DataFrame({'group': [g for g, *_ in CATS], 'n_genes': [len(members[n]) for n in names]}, index=names)
for s in ('2024-25', '2025-26'):
    for md in ('SW', 'SI'):
        o = core.index[(core.season == s) & (core.medium == md) & (core.phase == 'oiled')]
        c = core.index[(core.season == s) & (core.medium == md) & (core.phase == 'control')]
        if len(o) and len(c):
            summ[f'oil_{s}_{md}'] = med(o); summ[f'ctl_{s}_{md}'] = med(c)
            summ[f'log2_{s}_{md}'] = np.log2((med(o) + 0.5) / (med(c) + 0.5))
    for g, k in (('EXP', 'exp'), ('CTL', 'ctl')):
        p = core.index[(core.season == s) & (core.medium == 'SW') & (core.group == g) & (core.phase == 'pre')]
        summ[f'pre_{k}_{s}'] = med(p) if len(p) else np.nan
for s in ('2024-25', '2025-26'):
    summ[f'baci_{s}'] = np.log2((summ[f'oil_{s}_SW'] + .5) / (summ[f'pre_exp_{s}'] + .5)) - np.log2((summ[f'ctl_{s}_SW'] + .5) / (summ[f'pre_ctl_{s}'] + .5))
L = np.log2(E + 0.5)
r1, _ = C.contrast(L, C.PI); r2, _ = C.contrast(L, C.PT)
summ['ice_log2'] = r1.delta; summ['ice_up'] = r1.n_up; summ['ice_n'] = r1.n; summ['ice_q'] = r1.get('q', np.nan)
summ['top_log2'] = r2.delta; summ['top_up'] = r2.n_up; summ['top_n'] = r2.n; summ['top_q'] = r2.get('q', np.nan)
summ['mean_cpm'] = E.mean(1)
summ.to_csv(O + 'cat_summary.tsv', sep='\t')

# host-normalised N-acquisition share per genome (seawater, diesel vs control after dosing, per season)
ns_idx = np.unique(np.concatenate([members[n] for n in NSTRESS]))
NS = M.iloc[ns_idx].groupby(sub.bin.values[ns_idx]).sum()
NS = NS[NS.index != '']
hr = []
for s in ('2024-25', '2025-26'):
    o = list(core.index[(core.season == s) & (core.medium == 'SW') & (core.phase == 'oiled')])
    c = list(core.index[(core.season == s) & (core.medium == 'SW') & (core.phase == 'control')])
    for b in NS.index:
        ho, hc = BIN.loc[b, o], BIN.loc[b, c]
        if ho.median() < 100 and hc.median() < 100: continue
        so = (NS.loc[b, o] / ho.replace(0, np.nan)).median(); sc = (NS.loc[b, c] / hc.replace(0, np.nan)).median()
        hr.append(dict(season=s, bin=b, taxon=blab(b), phylum=rank(b, 'p'), host_oil=ho.median(), host_ctl=hc.median(),
                       share_oil=so, share_ctl=sc, log2_share=np.log2((so + 1e-4) / (sc + 1e-4))))
H = pd.DataFrame(hr); H.to_csv(O + 'host_nstress.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
cols = ['n_genes', 'mean_cpm', 'log2_2024-25_SW', 'log2_2025-26_SW', 'baci_2024-25', 'baci_2025-26', 'ice_log2', 'ice_up', 'ice_n', 'top_log2', 'top_up']
print(summ[cols].round(2).to_string())
print(H.sort_values('host_oil', ascending=False).head(30).round(3).to_string())
