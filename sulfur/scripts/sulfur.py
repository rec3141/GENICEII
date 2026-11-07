"""Sulfur, hydrogen, carbon-fixation and electron-acceptor genes: per-sample CPM (of CDS read pairs), host genomes, diesel vs control,
before-after-control-impact (2025-26 seawater), ice vs water. Gene sets from the coordinate-corrected KOs (genemap.pkl), plus Bakta products
for dsr accessory genes that have no KO. dsr/apr/sat direction is assigned by host taxonomy (Desulfobacterota etc. = reductive).
Out: sulfur/out/{cat_cpm.tsv, cat_host_cpm.pkl, cat_summary.tsv, cat_hosts.tsv, dsr_genes.tsv, sulfurimonas_inventory.tsv}"""
import sys, re, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
import contrasts as C
O = '/data/dev/rnaseq-exploration/sulfur/out/'; os.makedirs(O, exist_ok=True)

K = lambda s: s.split()
CATS = [  # (group, name, KOs, product regex or None)
 ('Sulfur oxidation', 'sqr (sulfide:quinone oxidoreductase)', K('K17218'), None),
 ('Sulfur oxidation', 'fccAB (flavocytochrome c sulfide dehydrogenase)', K('K17229 K17230'), None),
 ('Sulfur oxidation', 'soxAX', K('K17222 K17223'), None),
 ('Sulfur oxidation', 'soxYZ', K('K17226 K17227'), None),
 ('Sulfur oxidation', 'soxB', K('K17224'), None),
 ('Sulfur oxidation', 'soxCD (sulfane dehydrogenase)', K('K17225 K22622'), None),
 ('Sulfur oxidation', 'soeABC (cytoplasmic sulfite oxidation)', K('K21307 K21308 K21309'), None),
 ('Sulfur oxidation', 'sorA (periplasmic sulfite dehydrogenase)', K('K05301'), None),
 ('Sulfur oxidation', 'tsdA (thiosulfate dehydrogenase)', K('K19713'), None),
 ('Sulfur oxidation', 'sdo/ETHE1 (persulfide dioxygenase)', K('K17725'), None),
 ('Dissimilatory sulfur (direction by host)', 'dsrAB (dissimilatory sulfite reductase)', K('K11180 K11181'), None),
 ('Dissimilatory sulfur (direction by host)', 'aprAB (APS reductase)', K('K00394 K00395'), None),
 ('Dissimilatory sulfur (direction by host)', 'qmoABC (reductive-pathway electron transfer)', K('K16885 K16886 K16887'), None),
 ('Dissimilatory sulfur (direction by host)', 'dsrD (marker of reductive Dsr)', [], r'\bdsrD\b|DsrD'),
 ('Dissimilatory sulfur (direction by host)', 'dsrEFH (marker of oxidative rDsr)', [], r'\bdsr[EFH]\b|DsrE|DsrF|DsrH'),
 ('Sulfur reduction', 'ttrABC (tetrathionate reductase)', K('K08357 K08358 K08359'), None),
 ('Sulfur reduction', 'phsA/psrA (thiosulfate/polysulfide reductase)', K('K08352'), None),
 ('DMSP and methanethiol', 'dmdA (DMSP demethylase)', K('K17486'), None),
 ('DMSP and methanethiol', 'dmdBCD (MMPA-CoA pathway)', K('K20034 K20035 K20036'), None),
 ('DMSP and methanethiol', 'DMSP lyases (dddL/P/W/Y/X/D)', K('K16953 K28067 K28068 K28069 K28072 K28073'), None),
 ('DMSP and methanethiol', 'mddA (methanethiol S-methyltransferase)', K('K21310'), None),
 ('Assimilatory sulfate reduction (background)', 'cysDN + cysH + cysIJ', K('K00957 K00956 K00390 K00381 K00380'), None),
 ('Hydrogen', 'NiFe group 1 respiratory uptake (hyaAB/hybOC/hupSL)', K('K06281 K06282 K03620'), None),
 ('Hydrogen', 'NiFe group 1b quinone-reactive (hydAB, Campylobacterota)', K('K05922 K05927'), None),
 ('Hydrogen', 'NiFe group 2b sensory (hupUV)', K('K23548 K23549'), None),
 ('Hydrogen', 'NiFe group 3 bidirectional (hoxFUYH)', K('K00436 K18005 K18006 K18007 K05586 K05587 K05588'), None),
 ('Hydrogen', 'NiFe group 3 F420 / methyl-viologen (frh, mvh)', K('K00440 K00441 K00443 K14126 K14128'), None),
 ('Hydrogen', 'NiFe group 4 H2-evolving (hyc, hyf, ech, mbh)', K('K15830 K12136 K12137 K12138 K12139 K12140 K12141 K12142 K12143 K12144 K12145 K14086 K14087 K14088 K14089 K14090 K14091 K18016 K18017 K18023'), None),
 ('Hydrogen', '[FeFe] hydrogenases (hydA, hnd)', K('K00532 K00533 K00534 K06441 K17997 K17998 K17999 K18330 K18331 K18332 K17992 K25123'), None),
 ('Hydrogen', 'NiFe maturation (hypA–F)', K('K04651 K04652 K04653 K04654 K04655 K04656'), None),
 ('Hydrogen', 'nitrogenase (nifHDK, H2 by-product)', K('K02588 K02586 K02591'), None),
 ('Alternative donors', 'aerobic CO dehydrogenase (coxLMS)', K('K03520 K03519 K03518'), None),
 ('Alternative donors', 'formate dehydrogenase (fdoGH/fdhF)', K('K00123 K00124 K22015'), None),
 ('Carbon fixation', 'rTCA: ATP-citrate lyase (aclAB)', K('K15230 K15231'), None),
 ('Carbon fixation', 'rTCA: pyruvate:ferredoxin oxidoreductase (porABCD)', K('K00169 K00170 K00171 K00172'), None),
 ('Carbon fixation', 'rTCA: 2-oxoglutarate:ferredoxin oxidoreductase (korAB)', K('K00174 K00175'), None),
 ('Carbon fixation', 'Calvin: RuBisCO (rbcL) + prk', K('K01601 K00855'), None),
 ('Electron acceptors', 'nitrate: periplasmic napAB', K('K02567 K02568'), None),
 ('Electron acceptors', 'nitrate: membrane narGH', K('K00370 K00371'), None),
 ('Electron acceptors', 'nitrite to NO: nirS / nirK', K('K15864 K00368'), None),
 ('Electron acceptors', 'NO reductase norBC', K('K04561 K02305'), None),
 ('Electron acceptors', 'N2O reductase nosZ', K('K00376'), None),
 ('Electron acceptors', 'O2, high affinity: cbb3 oxidase (ccoNO)', K('K00404 K00405'), None),
 ('Electron acceptors', 'O2, high affinity: bd oxidase (cydAB)', K('K00425 K00426'), None),
 ('Electron acceptors', 'O2, low affinity: aa3 oxidase (coxA/ctaD)', K('K02274'), None),
]
REDUCERS = r'Desulfobacterota|Desulfobacterota_[A-Z]|Thermodesulfobacteriota|Desulfuromonadota|Nitrospirota|Bacillota'

X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
kos = G['KO'].fillna('').str.split(',')
prod = G['product'].fillna('')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def rank(b, r):
    t = str(gt.get(str(b).replace('dastool-', ''), ''))
    x = [y[3:] for y in t.split(';') if y.startswith(r + '__')]
    return x[0] if x else ''
def blab(b):
    if not isinstance(b, str) or not b: return 'unbinned'
    if special_label(b): return special_label(b)
    t = [y.split('__', 1)[1] for y in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in y and y.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')

allko = set(k for _, _, ks, _ in CATS for k in ks)
selko = kos.apply(lambda x: bool(allko.intersection(x))).values
selrx = prod.str.contains(r'\bdsrD\b|DsrD|\bdsr[EFH]\b|DsrE|DsrF|DsrH|issimilatory sulfite reductase', regex=True).values
idx = np.where((selko | selrx) & cds)[0]
sub = G.iloc[idx].reset_index(drop=True)
M = pd.DataFrame(X[idx].toarray() / tot * 1e6, columns=S)
sub['bin'] = sub['bin'].fillna('')
sub['phylum'] = [rank(b, 'p') if b else '' for b in sub['bin']]
sub['taxon'] = [blab(b) for b in sub['bin']]
skos = sub['KO'].fillna('').str.split(',')
rows, hosts, members = {}, {}, {}
for grp, name, ks, rx in CATS:
    sel = skos.apply(lambda x: bool(set(ks).intersection(x))).values if ks else np.zeros(len(sub), bool)
    if rx: sel |= sub['product'].str.contains(rx, regex=True).values
    members[name] = np.where(sel)[0]
    rows[name] = M[sel].sum(0)
    hosts[name] = M[sel].groupby(sub['bin'][sel].values).sum()
E = pd.DataFrame(rows).T          # categories x samples
E.to_csv(O + 'cat_cpm.tsv', sep='\t')
pd.to_pickle({'hosts': hosts, 'members': members}, O + 'cat_host_cpm.pkl')
sub.drop(columns=[c for c in sub.columns if c not in ('gene_id', 'contig', 'product', 'KO', 'bin', 'phylum', 'taxon', 'length')]).to_csv(O + 'genes.tsv', sep='\t', index=False)
M.index = sub.gene_id; M.to_pickle(O + 'gene_cpm.pkl')

# summaries
core = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
def med(cols): return E[cols].median(1)
summ = pd.DataFrame({'group': [g for g, *_ in CATS], 'n_genes': [len(members[n]) for _, n, *_ in CATS]}, index=[n for _, n, *_ in CATS])
for s in ('2024-25', '2025-26'):
    for md in ('SW', 'SI'):
        o = core.index[(core.season == s) & (core.medium == md) & (core.phase == 'oiled')]
        c = core.index[(core.season == s) & (core.medium == md) & (core.phase == 'control')]
        if len(o) and len(c):
            summ[f'oil_{s}_{md}'] = med(o); summ[f'ctl_{s}_{md}'] = med(c)
            summ[f'log2_{s}_{md}'] = np.log2((med(o) + 0.5) / (med(c) + 0.5))
    pre_e = core.index[(core.season == s) & (core.medium == 'SW') & (core.group == 'EXP') & (core.phase == 'pre')]
    pre_c = core.index[(core.season == s) & (core.medium == 'SW') & (core.group == 'CTL') & (core.phase == 'pre')]
    summ[f'pre_exp_{s}'] = med(pre_e) if len(pre_e) else np.nan; summ[f'pre_ctl_{s}'] = med(pre_c) if len(pre_c) else np.nan
# BACI 2025-26 seawater: log2(oil post / oil pre) - log2(ctl post / ctl pre), medians, pseudocount 0.5
s = '2025-26'
if 'oil_2025-26_SW' in summ:
    summ['baci_2025-26'] = (np.log2((summ['oil_2025-26_SW'] + .5) / (summ['pre_exp_2025-26'] + .5)) - np.log2((summ['ctl_2025-26_SW'] + .5) / (summ['pre_ctl_2025-26'] + .5)))
# ice vs water, top vs bottom (paired)
L = np.log2(E + 0.5)
r1, _ = C.contrast(L, C.PI); r2, _ = C.contrast(L, C.PT)
summ['ice_log2'] = r1.delta; summ['ice_up'] = r1.n_up; summ['ice_n'] = r1.n
summ['top_log2'] = r2.delta; summ['top_up'] = r2.n_up; summ['top_n'] = r2.n
summ['mean_cpm'] = E.mean(1)
summ.to_csv(O + 'cat_summary.tsv', sep='\t')

# hosts: share of each category's transcripts by genome, all samples and diesel-pool-after-dosing
hrows = []
oil = m.index[m.phase == 'oiled']; ctl = m.index[m.phase == 'control']
for name, H in hosts.items():
    for where, cols in (('all', S), ('diesel', list(oil)), ('control', list(ctl))):
        tot_ = H[cols].sum().sum()
        if tot_ <= 0: continue
        sh = (H[cols].sum(1) / tot_).sort_values(ascending=False)
        for b, v in sh.head(6).items():
            hrows.append(dict(category=name, where=where, bin=b, taxon=blab(b), phylum=rank(b, 'p') if b else '', share=v))
pd.DataFrame(hrows).to_csv(O + 'cat_hosts.tsv', sep='\t', index=False)

# dsr/apr/sat genes with host phylum -> reductive vs oxidative call
d = sub.iloc[np.concatenate([members[n] for n in ('dsrAB (dissimilatory sulfite reductase)', 'aprAB (APS reductase)', 'qmoABC (reductive-pathway electron transfer)', 'dsrD (marker of reductive Dsr)', 'dsrEFH (marker of oxidative rDsr)')])].copy()
d['mean_cpm'] = M.iloc[d.index].mean(1).values
d['direction'] = np.where(d.phylum.str.contains(REDUCERS, regex=True), 'reductive (host taxonomy)', np.where(d.phylum == '', 'unbinned', 'oxidative (host taxonomy)'))
d.sort_values('mean_cpm', ascending=False).to_csv(O + 'dsr_genes.tsv', sep='\t', index=False)

pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
cols = ['n_genes', 'mean_cpm', 'log2_2024-25_SW', 'log2_2025-26_SW', 'baci_2025-26', 'ice_log2', 'ice_up', 'top_log2', 'top_up']
print(summ[cols].round(2).to_string())
print(d.groupby(['direction']).mean_cpm.agg(['size', 'sum']).round(1))
