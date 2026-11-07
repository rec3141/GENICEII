"""Porticoccaceae 50-400-T64 (dastool-comebin_2576): per-gene expression (CPM of CDS reads) in diesel vs control seawater, within-genome share,
functional toolkit (product/KO regex + CANT-HYD + hgt HC catalogue), N/S/H2 genes, islands/MGEs.
Out: porticocc/out/{genes_expr.tsv, toolkit.tsv, summary.json}"""
import sys, os, json, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
B = '/data/dev/rnaseq-exploration/'; O = B + 'porticocc/out/'
BIN = 'dastool-comebin_2576'
X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
idx = np.where((G.bin.values == BIN) & cds)[0]
sub = G.iloc[idx].reset_index(drop=True)
M = pd.DataFrame(X[idx].toarray() / tot * 1e6, columns=S)
host = M.sum(0)
core = m[m.group.isin(['CTL', 'EXP']) & (m.medium == 'SW')]
cols = {}
for s in ('2024-25', '2025-26'):
    for k, sel in (('oil', core.phase == 'oiled'), ('ctl', core.phase == 'control'), ('pre_exp', (core.phase == 'pre') & (core.group == 'EXP')), ('pre_ctl', (core.phase == 'pre') & (core.group == 'CTL'))):
        cols[f'{k}_{s}'] = list(core.index[(core.season == s) & sel])
for k, c in cols.items():
    sub['cpm_' + k] = M[c].median(1).values
    sub['share_' + k] = (M[c].div(host[c].replace(0, np.nan), axis=1)).median(1).values * 1e3   # per mille of the genome's transcripts
    sub['pshare_' + k] = M[c].sum(1).values / max(host[c].sum(), 1e-9) * 1e3                       # pooled over samples (robust when the genome is rare)
sub['pshare_bg_2025-26'] = M[cols['ctl_2025-26'] + cols['pre_exp_2025-26'] + cols['pre_ctl_2025-26']].sum(1).values / host[cols['ctl_2025-26'] + cols['pre_exp_2025-26'] + cols['pre_ctl_2025-26']].sum() * 1e3
sub['bg_reads_cpm_sum'] = M[cols['ctl_2025-26'] + cols['pre_exp_2025-26'] + cols['pre_ctl_2025-26']].sum(1).values
sub['log2_share_2025-26'] = np.log2((sub['share_oil_2025-26'] + 0.01) / (sub['share_ctl_2025-26'] + 0.01))
ch = pd.read_csv(B + 'canthyd/canthyd_hits.tsv', sep='\t').sort_values('score', ascending=False).drop_duplicates('gene_id').set_index('gene_id')
sub['canthyd'] = sub.gene_id.map(ch.hmm).fillna('')
hc = pd.read_csv(B + 'hgt/hc_genes.tsv', sep='\t', keep_default_na=False).set_index('gene_id')
sub['hc_category'] = sub.gene_id.map(hc.category).fillna('')
RX = [  # (toolkit class, regex on product, KOs)
 ('CuMMO (pmo/amo-like)', r'monooxygenase subunit [ABC]|AmoB|PmoA', 'K10944 K10945 K10946 K28504'),
 ('alkane hydroxylases: AlkB, cytochrome P450', r'[Aa]lkane (1-)?monooxygenase|alkB|[Cc]ytochrome P450', 'K00496 K00493'),
 ('Baeyer-Villiger monooxygenases (ketone to ester)', r'[Cc]yclohexanone monooxygenase|[Aa]cetone monooxygenase|Baeyer', 'K03379 K18371 K18091 K21730'),
 ('cycloalkanol / secondary-alcohol dehydrogenases', r'[Cc]yclopentanol|[Cc]yclohexanol|secondary alcohol', 'K17754'),
 ('other monooxygenases', r'[Mm]onooxygenase|[Hh]ydroxylase', ''),
 ('ring-cleaving / aromatic dioxygenases', r'[Dd]ioxygenase', ''),
 ('alcohol dehydrogenases (incl. PQQ)', r'[Aa]lcohol dehydrogenase|PQQ|[Mm]ethanol dehydrogenase|quinoprotein', 'K00114 K00119 K14028 K23995 K00001 K13953 K13954'),
 ('aldehyde dehydrogenases', r'[Aa]ldehyde dehydrogenase', 'K00128 K00138 K00140 K14519'),
 ('acyl-CoA / beta-oxidation', r'[Aa]cyl-CoA (dehydrogenase|oxidase|synthetase|ligase)|[Ee]noyl-CoA|3-hydroxyacyl-CoA|[Tt]hiolase|acetyl-CoA acyltransferase', ''),
 ('propionate / methylcitrate (odd-chain)', r'[Mm]ethylcitrate|prpB|prpC|prpD|[Pp]ropionyl-CoA carboxylase|methylmalonyl', 'K01720 K01659 K03417 K01965 K01966 K01847'),
 ('glyoxylate shunt', r'[Ii]socitrate lyase|[Mm]alate synthase', 'K01637 K01638'),
 ('copper handling', r'[Cc]opper|Cu\(I\)|CopZ|CopA|CusA|multicopper|[Cc]opper oxidase|CopC|CopD', ''),
 ('PHA storage', r'[Pp]oly\(3-hydroxy|PhaC|[Pp]hasin|PhaR', ''),
 ('rhodopsin / light', r'[Rr]hodopsin|[Pp]hotolyase|[Bb]acteriochlorophyll|puf[LM]', ''),
 ('nitrogen acquisition', r'[Aa]mmonium transporter|[Nn]itrogen regulatory protein P-II|[Gg]lutamine synthetase|[Uu]rea ABC|[Uu]rease|[Nn]itrate|[Nn]itrite', ''),
 ('sulfur / hydrogen', r'[Ss]ulfide|[Tt]hiosulfate|SoxB|\b[Hh]ydrogenase|[Ss]ulfite', ''),
 ('terminal oxidases', r'cbb3|[Cc]ytochrome bd|[Cc]ytochrome c oxidase|[Cc]ytochrome-c oxidase', ''),
 ('motility / chemotaxis', r'[Ff]lagell|[Cc]hemotaxis|MCP|methyl-accepting', ''),
]
kos = sub.KO.fillna('').str.split(',')
sub['toolkit'] = ''
for name, rx, ks in RX:
    sel = sub['product'].str.contains(rx, regex=True) | kos.apply(lambda x: bool(set(ks.split()).intersection(x)))
    sub.loc[sel & (sub.toolkit == ''), 'toolkit'] = name
sub.to_csv(O + 'genes_expr.tsv', sep='\t', index=False)
T = []
for name, *_ in RX:
    x = sub[sub.toolkit == name]
    T.append(dict(toolkit=name, n=len(x), oil25=x['cpm_oil_2025-26'].sum(), ctl25=x['cpm_ctl_2025-26'].sum(),
                  share_oil25=x['share_oil_2025-26'].sum(), share_ctl25=x['share_ctl_2025-26'].sum(),
                  oil24=x['cpm_oil_2024-25'].sum(), ctl24=x['cpm_ctl_2024-25'].sum(),
                  pshare_oil25=x['pshare_oil_2025-26'].sum(), pshare_bg25=x['pshare_bg_2025-26'].sum(), pshare_oil24=x['pshare_oil_2024-25'].sum()))
T = pd.DataFrame(T); T['log2_share25'] = np.log2((T.share_oil25 + .01) / (T.share_ctl25 + .01)); T['log2_pshare25'] = np.log2((T.pshare_oil25 + .05) / (T.pshare_bg25 + .05)); T.to_csv(O + 'toolkit.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
print(T.round(2).to_string())
print('\nTop expressed genes in the 2025-26 diesel pool (per mille of genome transcripts):')
print(sub.sort_values('share_oil_2025-26', ascending=False).head(40)[['gene_id', 'contig', 'gene', 'product', 'KO', 'toolkit', 'cpm_oil_2025-26', 'share_oil_2025-26', 'share_ctl_2025-26', 'log2_share_2025-26']].round(2).to_string())
print('\nMost induced (share, 2025-26, genes with >=20 CPM in diesel):')
print(sub[sub['cpm_oil_2025-26'] >= 20].sort_values('log2_share_2025-26', ascending=False).head(30)[['gene_id', 'product', 'toolkit', 'cpm_oil_2025-26', 'share_oil_2025-26', 'share_ctl_2025-26', 'log2_share_2025-26']].round(2).to_string())
print('\nCANT-HYD / HC catalogue hits:'); print(sub[(sub.canthyd != '') | (sub.hc_category != '')][['gene_id', 'product', 'canthyd', 'hc_category', 'cpm_oil_2025-26', 'cpm_ctl_2025-26']].round(1).to_string())
pd.to_pickle({'M': M.set_index(sub.gene_id), 'host': host}, O + 'gene_cpm.pkl')
print('\nhost CPM medians:', {k: round(float(host[c].median()), 0) for k, c in cols.items()})
