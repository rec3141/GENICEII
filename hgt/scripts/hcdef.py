"""Hydrocarbon-degradation gene categories from Bakta (extra) products; KO remap / CANT-HYD added in s1_select.py."""
import re, sys
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
import numpy as np, pandas as pd

# category, pathway, step, product regex (lower-case Bakta product), KOs (for the coordinate-remapped KO table)
CAT = [
    ('alkB', 'alkane', 'activation', r'^alkane (1-)?monooxygenase( \d)?$|^alkane hydroxylase|^alkane 1-monooxygenase alkb', ['K00496']),
    ('alkG/alkT (rubredoxin system)', 'alkane', 'accessory', r'^rubredoxin( \d)?$|^rubredoxin-nad\(\+\) reductase|^rubredoxin--nad', ['K05297']),
    ('CYP153', 'alkane', 'activation', r'^$a', []),
    ('ladA/almA', 'alkane', 'activation', r'long-chain alkane monooxygenase|^flavin-binding monooxygenase( alma)?$', []),
    ('C1-C4 alkane MO', 'short alkane', 'activation', r'^(soluble )?methane monooxygenase(?!/ammonia)|propane (2-)?monooxygenase|butane monooxygenase', ['K10944', 'K16157', 'K18223', 'K18224']),
    ('RHO dioxygenase alpha/beta', 'aromatic (PAH/BTEX)', 'activation',
     r'ring[- ]hydroxylating dioxygenase (subunit )?(alpha|beta|large|small)|ring[- ]hydroxylating dioxygenase$|^naphthalene 1,2-dioxygenase( subunit (alpha|beta)| large| small|$)|^biphenyl (2,3-)?dioxygenase( subunit (alpha|beta)| large| small)?$|^benzene 1,2-dioxygenase subunit|^toluene (1,2-)?dioxygenase|^benzoate 1,2-dioxygenase (large|small|subunit)|anthranilate 1,2-dioxygenase (large|small)|phthalate (4,5|3,4)-dioxygenase$|^pah dioxygenase',
     ['K14579', 'K14580', 'K11943', 'K11944', 'K08689', 'K15750', 'K03268', 'K16268', 'K05549', 'K05550', 'K18068', 'K16319', 'K16320']),
    ('RHO ferredoxin/reductase', 'aromatic (PAH/BTEX)', 'accessory', r'naphthalene 1,2-dioxygenase.*(ferredoxin|reductase)|biphenyl dioxygenase (ferredoxin|reductase)|benzoate 1,2-dioxygenase (electron|reductase)', ['K14578', 'K14581', 'K18087', 'K18088', 'K05784']),
    ('BMM toluene/phenol monooxygenase', 'aromatic (BTEX/phenol)', 'activation',
     r'toluene[- ](4-|3-|2-)?monooxygenase|toluene methyl-monooxygenase|^phenol (2-)?monooxygenase|^phenol hydroxylase( p\d| component)?|^multicomponent phenol hydroxylase',
     ['K15760', 'K15761', 'K15762', 'K15763', 'K15764', 'K15765', 'K16242', 'K16243', 'K16244', 'K16245', 'K16246', 'K16249', 'K15757', 'K15758']),
    ('extradiol dioxygenase', 'meta cleavage', 'ring cleavage',
     r'^catechol[- ]2,3-dioxygenase( \d)?$|biphenyl-2,3-diol 1,2-dioxygenase|2,3-dihydroxybiphenyl 1,2-dioxygenase|1,2-dihydroxynaphthalene dioxygenase|protocatechuate 4,5-dioxygenase|3,4-dihydroxyphenylacetate 2,3-dioxygenase|^3-methylcatechol 2,3-dioxygenase|2,3-dihydroxy-p-cumate',
     ['K00446', 'K07104', 'K00462', 'K14583', 'K04100', 'K04101', 'K16270', 'K10621']),
    ('intradiol dioxygenase', 'ortho cleavage', 'ring cleavage', r'^catechol 1,2-dioxygenase|^protocatechuate 3,4-dioxygenase|^chlorocatechol 1,2-dioxygenase', ['K03381', 'K00448', 'K00449', 'K15253']),
    ('gentisate 1,2-dioxygenase', 'gentisate', 'ring cleavage', r'^gentisate 1,2-dioxygenase', ['K00450']),
    ('salicylate hydroxylase', 'PAH lower', 'downstream', r'^salicylate (1-)?(hydroxylase|monooxygenase)', ['K00480']),
    ('PAH lower (nahB/D/E/F)', 'PAH lower', 'downstream', r'^(cis-)?(1,2-)?dihydro.*diol dehydrogenase|naphthalene dihydrodiol|trans-o-hydroxybenzylidenepyruvate|salicylaldehyde dehydrogenase', ['K14582']),
    ('maleylacetate reductase', 'ortho/chloroaromatic', 'downstream', r'maleylacetate reductase', ['K00217']),
    ('ortho lower (cat/pca)', 'ortho cleavage', 'downstream',
     r'^muconate cycloisomerase|^muconolactone (delta-)?isomerase|^3-oxoadipate enol-lactonase|^3-oxoadipate coa-transferase|^4-carboxymuconolactone decarboxylase( pcac)?$|^3-carboxy-cis,cis-muconate cycloisomerase|^beta-ketoadipate',
     ['K01856', 'K03464', 'K01055', 'K01031', 'K01032', 'K01857', 'K01607']),
    ('meta lower (xyl/dmp)', 'meta cleavage', 'downstream',
     r'2-hydroxymuconic semialdehyde|2-hydroxymuconate[- ]semialdehyde|2-hydroxymuconate-6-semialdehyde|2-oxopent-4-enoate hydratase|2-keto-4-pentenoate hydratase|^4-hydroxy-2-oxovalerate aldolase|4-oxalocrotonate (decarboxylase|tautomerase)$|acetaldehyde dehydrogenase \(acylating\)',
     ['K10216', 'K10217', 'K01617', 'K02554', 'K01666', 'K04073']),
    ('paa phenylacetate', 'phenylacetate', 'activation',
     r'phenylacetyl-coa epoxidase|^phenylacetate--coa ligase( paak)?$|^phenylacetic acid degradation protein|oxepin-coa hydrolase|1,2-epoxyphenylacetyl-coa isomerase',
     ['K02609', 'K02610', 'K02611', 'K02612', 'K02613', 'K02618', 'K01912']),
    ('anaerobic (bss/bcr/bam)', 'anaerobic aromatics', 'activation', r'benzylsuccinate synthase|^benzoyl-coa reductase|alkylsuccinate synthase|^bam[b-i] ', ['K07540', 'K29040', 'K29041', 'K04112', 'K04113', 'K04114', 'K04115', 'K19515', 'K19516']),
]
PATHWAY = {c: p for c, p, _, _, _ in CAT}
STEP = {c: s for c, _, s, _, _ in CAT}


def classify(g, ko=None):
    """-> (category Series, evidence Series) for each CDS in genemap g. ko: Series gene_id -> 'Kxxxx,...' (remapped)."""
    cds = g['type'].values == 'CDS'
    prod = g['product'].str.lower().values
    reg = np.array(['regulator' in p or 'regulatory' in p for p in prod])
    cat = np.array([''] * len(g), dtype=object); ev = cat.copy()
    for c, _, _, rx, _ in CAT:
        r = re.compile(rx)
        hit = np.array([bool(r.search(p)) for p in prod]) & cds & (cat == '') & ~reg
        cat[hit] = c; ev[hit] = 'product'
    if ko is not None:
        ko2cat = {k: c for c, _, _, _, ks in CAT for k in ks}
        kk = g['gene_id'].map(ko).fillna('').values
        for i in np.where(cds & (kk != ''))[0]:
            cs = [ko2cat[k] for k in kk[i].split(',') if k in ko2cat]
            if not cs: continue
            if cat[i] == '': cat[i] = cs[0]; ev[i] = 'KO'
            elif cat[i] in cs: ev[i] = 'product+KO'
    return pd.Series(cat, index=g.index), pd.Series(ev, index=g.index)

CANTHYD = {'AlkB': 'alkB', 'CYP153': 'CYP153', 'NdoB': 'RHO dioxygenase alpha/beta', 'NdoC': 'RHO dioxygenase alpha/beta',
           'non_NdoB_type': 'RHO dioxygenase alpha/beta', 'MAH_alpha': 'RHO dioxygenase alpha/beta', 'MAH_beta': 'RHO dioxygenase alpha/beta',
           'TmoA_BmoA': 'BMM toluene/phenol monooxygenase', 'TmoE': 'BMM toluene/phenol monooxygenase', 'TomA3': 'BMM toluene/phenol monooxygenase',
           'TomA4': 'BMM toluene/phenol monooxygenase', 'DmpO': 'BMM toluene/phenol monooxygenase', 'LadA_alpha': 'ladA/almA', 'LadA_beta': 'ladA/almA',
           'AlmA_GroupI': 'ladA/almA', 'AlmA_GroupIII': 'ladA/almA', 'PrmA': 'C1-C4 alkane MO', 'CmdA': 'RHO dioxygenase alpha/beta'}
