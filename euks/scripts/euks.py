"""Eukaryotes from 18S rRNA reads (SILVA, marine_coassembly.taxonomy_rrna.tsv; genus-rank rows partition each domain).
Kaiju (prokaryote + virus DB) cannot see eukaryotic mRNA; plastid transcription is tracked by Kaiju Cyanobacteriota mRNA (plastid-encoded
genes hit cyanobacteria) and plastid MAG TPM. Groups by SILVA lineage; likely contaminant metazoa (mammals, fish, insects) dropped.
Out: euks/out/{group_share.tsv, genus_share.tsv, sample_stats.tsv, group_summary.tsv, genus_summary.tsv, virus_links.tsv}"""
import sys, re, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
import contrasts as C
from scipy.stats import spearmanr
O = '/data/dev/rnaseq-exploration/euks/out/'; E0 = '/data/dev/rnaseq-exploration/'
S = samples(); m = meta().loc[S]
t = pd.read_csv(D + 'marine_coassembly.taxonomy_rrna.tsv', sep='\t')
dom = t[t['rank'] == 'domain'].set_index('lineage')[S]
g = t[(t['rank'] == 'genus') & t.lineage.str.startswith('Eukaryota')].set_index('lineage')[S]
CONTAM = r'Mammalia|Euteleostomi|Pterygota|Actinopteri|Aves|Streptophyta;Embryophyta|Phragmoplastophyta|Hexapoda'
g = g[~g.index.str.contains(CONTAM)]
GROUPS = [  # (group, trophic role, regex on lineage) first match wins
 ('Diatoms', 'phototroph', r'Diatomea'),
 ('Haptophytes', 'phototroph', r'Haptophyta'),
 ('Green algae (prasinophytes etc.)', 'phototroph', r'Chloroplastida|prasinophyte|Archaeplastida'),
 ('Cryptophytes', 'phototroph', r'Cryptomonadales|Cryptophyceae;(?!Kathablepharidae)'),
 ('Chrysophytes and other ochrophytes', 'mixotroph', r'Ochrophyta'),
 ('Dinoflagellates', 'mixotroph', r'Dinophyceae|Dinoflagellata;(?!MALV|Syndiniales)'),
 ('Syndiniales (parasitic dinoflagellates)', 'parasite', r'MALV|Syndiniales'),
 ('Ciliates', 'grazer', r'Ciliophora'),
 ('Cercozoa', 'grazer', r'Cercozoa'),
 ('Other Rhizaria', 'grazer', r'Rhizaria'),
 ('Choanoflagellates', 'grazer', r'Diaphanoeca|Helgoeca|Bicosta|Savillea|Stephanoeca|Acanthoeca|Choanoflag|Monosiga|Salpingoeca|Codosiga|Pleurasiga|Calliacantha|Parvicorbicula'),
 ('Katablepharids', 'grazer', r'Kathablepharidae'),
 ('Bicosoecids (Cafeteria etc.)', 'grazer', r'Bicosoecida|Cafeteria'),
 ('MAST and other heterotrophic stramenopiles', 'grazer', r'MAST|Developayella|Stramenopiles;(?!Ochrophyta|Labyrinthul|Hyphochytri|Peronospor|Oomycota)'),
 ('Telonemids', 'grazer', r'Telonema'),
 ('Kinetoplastids and diplonemids', 'grazer', r'Kinetoplastea|Diplonemea'),
 ('Euglenids', 'mixotroph', r'Euglenida|Euglenozoa'),
 ('Labyrinthulids and hyphochytrids', 'osmotroph', r'Labyrinthul|Hyphochytri|Peronospor|Oomycota'),
 ('Fungi', 'osmotroph', r'Nucletmycea'),
 ('Metazoa (hydrozoans, copepods, …)', 'animal', r'Hydrozoa|Leptothecata|Metazoa|Arthropoda|Copepoda|Cnidaria|Annelida|Mollusca|Nematoda|Rotifera|Holozoa;(?!unclassified)'),
 ('Other and unclassified eukaryotes', 'other', r'.'),
]
def grp(l):
    for n, r, rx in GROUPS:
        if re.search(rx, l): return n
lab = pd.Series({l: grp(l) for l in g.index})
role = dict((n, r) for n, r, _ in GROUPS)
etot = g.sum(0)
GS = g.groupby(lab.values).sum() / etot * 100          # % of eukaryotic 18S reads
GS = GS.loc[[n for n, *_ in GROUPS if n in GS.index]]
GS.to_csv(O + 'group_share.tsv', sep='\t')
gs = g.div(etot, axis=1) * 100
keep = gs.max(1) >= 2
short = lambda l: ';'.join([x for x in l.split(';')[3:] if x not in ('unclassified', 'Incertae Sedis')][-2:]) or l
GEN = gs[keep]; GEN.index = [f'{short(l)} [{lab[l]}]' for l in GEN.index]; GEN = GEN.groupby(level=0).sum()
GEN.to_csv(O + 'genus_share.tsv', sep='\t')

# per-sample stats
k = pd.read_csv(D + 'marine_coassembly.taxonomy_mrna_kaiju.tsv', sep='\t')
kd = k[k['rank'] == 'domain'].set_index('lineage')[S]
kp = k[k['rank'] == 'phylum'].set_index('lineage')[S]
ktot = kd.loc[['Bacteria', 'Archaea', 'Viruses']].sum(0)
tpm = table('bin_tpm')
pl = [b for b in tpm.index if is_plastid(b)]
vs = pd.read_csv(E0 + 'viruses/out/vrna_samples.tsv', sep='\t', index_col=0)
vc = pd.read_pickle(E0 + 'viruses/out/votu_rna_cpm.pkl')
VP = 'contig_0049668'
st = pd.DataFrame({'euk_rrna_pct': dom.loc['Eukaryota'] / dom.sum(0) * 100, 'euk_reads': etot,
                   'plastid_mrna_pct': kp.loc['Bacteria;Cyanobacteriota'] / ktot * 100, 'plastid_mag_tpm': tpm.loc[pl, S].sum(0),
                   'giant_virus_cpm': vs.reindex(S)['Megaviricetes'], 'virophage_class_cpm': vs.reindex(S)['Maveriviricetes'],
                   'virophage_cpm': vc.loc[VP, S] if VP in vc.index else np.nan})
for r in ('phototroph', 'mixotroph', 'grazer', 'parasite', 'osmotroph'):
    st[f'{r}_pct'] = GS.loc[[n for n in GS.index if role[n] == r]].sum(0)
st['grazer_to_phototroph'] = np.log2((st.grazer_pct + .1) / (st.phototroph_pct + .1))
cy = tpm.loc[[b for b in tpm.index if 'Cycloclasticus' in str(kaiju_bin_label(b)) or b in ('dastool-metabat_9342', 'dastool-lorbin_184', 'dastool-metabat_5432')], S].sum(0)
st['cycloclasticus_tpm'] = cy
st.to_csv(O + 'sample_stats.tsv', sep='\t')

# summaries: groups and stats rows; seawater diesel vs control, BACI per season; ice/water and top/bottom paired
core = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
def summarise(X, pc):
    out = pd.DataFrame(index=X.index)
    for s in ('2024-25', '2025-26'):
        sel = lambda ph, g=None, md='SW': core.index[(core.season == s) & (core.medium == md) & (core.phase == ph) & ((core.group == g) if g else True)]
        o, c, pe, pcn = X[sel('oiled')].median(1), X[sel('control')].median(1), X[sel('pre', 'EXP')].median(1), X[sel('pre', 'CTL')].median(1)
        out[f'oil_{s}'], out[f'ctl_{s}'], out[f'pre_exp_{s}'], out[f'pre_ctl_{s}'] = o, c, pe, pcn
        out[f'lfc_{s}'] = np.log2((o + pc) / (c + pc))
        out[f'baci_{s}'] = np.log2((o + pc) / (pe + pc)) - np.log2((c + pc) / (pcn + pc))
        ice_o, ice_c = X[sel('oiled', md='SI')].median(1), X[sel('control', md='SI')].median(1)
        out[f'lfc_ice_{s}'] = np.log2((ice_o + pc) / (ice_c + pc))
    L = np.log2(X + pc)
    r1, _ = C.contrast(L, C.PI); r2, _ = C.contrast(L, C.PT)
    out['ice'], out['ice_up'], out['ice_n'], out['ice_q'] = r1.delta, r1.n_up, r1.n, r1.q
    out['top'], out['top_up'], out['top_n'], out['top_q'] = r2.delta, r2.n_up, r2.n, r2.q
    out['mean'] = X.mean(1)
    return out
gsum = summarise(GS, 0.1); gsum['role'] = [role[n] for n in gsum.index]
gsum.to_csv(O + 'group_summary.tsv', sep='\t')
nsum = summarise(GEN, 0.1); nsum.to_csv(O + 'genus_summary.tsv', sep='\t')
ssum = summarise(st[['euk_rrna_pct', 'plastid_mrna_pct', 'plastid_mag_tpm', 'giant_virus_cpm', 'virophage_cpm', 'grazer_pct', 'phototroph_pct']].T.astype(float), 0.1)
ssum.to_csv(O + 'stat_summary.tsv', sep='\t')

# virus links: Spearman across pool samples (CTL/EXP, SW+SI) per season, giant viruses / virophage vs each genus share
vl = []
for s in ('2024-25', '2025-26'):
    ix = core.index[core.season == s]
    for v in ('giant_virus_cpm', 'virophage_cpm'):
        y = st.loc[ix, v]
        if y.notna().sum() < 8 or y.std() == 0: continue
        for n in list(GEN.index) + list(GS.index):
            x = (GEN if n in GEN.index else GS).loc[n, ix]
            if x.std() == 0: continue
            rho, p = spearmanr(x, y)
            vl.append(dict(season=s, virus=v, taxon=n, level='genus' if n in GEN.index else 'group', rho=rho, p=p, n=len(ix)))
VL = pd.DataFrame(vl); VL['q'] = bh(VL.p.fillna(1).values); VL.to_csv(O + 'virus_links.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(gsum[['role', 'mean', 'oil_2024-25', 'ctl_2024-25', 'lfc_2024-25', 'baci_2024-25', 'oil_2025-26', 'ctl_2025-26', 'lfc_2025-26', 'baci_2025-26', 'ice', 'ice_up', 'top', 'top_up']].round(2).to_string())
print(ssum[['oil_2024-25', 'ctl_2024-25', 'pre_exp_2024-25', 'pre_ctl_2024-25', 'oil_2025-26', 'ctl_2025-26', 'pre_exp_2025-26', 'pre_ctl_2025-26', 'baci_2025-26', 'ice', 'ice_up', 'top', 'top_up']].round(2).to_string())
for (s, v), x in VL.groupby(['season', 'virus']):
    print(s, v); print(x.sort_values('rho', ascending=False).head(8)[['taxon', 'rho', 'q']].round(3).to_string(index=False))
