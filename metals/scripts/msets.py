"""Metal and toxicity gene sets (corrected KOs + Bakta product regexes) -> metals/out/members.tsv (gene, category, group, evidence).
A gene joins the first category it matches, in the order listed. Product regexes catch genes the KO transfer missed."""
import sys, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
O = '/data/dev/rnaseq-exploration/metals/out/'
# (category, group, KOs, product regex or None)
SETS = [
 ('Mercury (mer)', 'toxic metal', 'K00520 K08363 K08364 K08365 K19057 K19058 K19059', r'mercur|organomercurial|alkylmercury'),
 ('Arsenic (ars, aio)', 'toxic metal', 'K00537 K03741 K18701 K03325 K03893 K11811 K03892 K23988 K07755 K08355 K08356', r'\barsenic|\barsenate|\barsenite|arsenical'),
 ('Copper (cop, cue, pco, cus sensors)', 'toxic metal', 'K17686 K01533 K14588 K07156 K07233 K07245 K07213 K21600 K11923 K19591 K07644 K07665', r'copper resistance|copper-(translocating|exporting)|multicopper oxidase|\bCopZ\b|\bCueO\b'),
 ('Heavy-metal RND efflux (czc, cus, cnr, sil)', 'toxic metal', 'K15725 K15726 K15727 K07787 K07796 K07798 K07810 K11326', r'cobalt-zinc-cadmium|heavy metal efflux|\bCzc[ABC]\b|\bCus[ABCF]\b|silver efflux'),
 ('Zn/Cd/Pb export (zntA, cadA, czcD) + sensors', 'toxic metal', 'K01534 K16264 K13638 K21885 K21903 K22043 K22298', r'cadmium|lead.*ATPase|zinc.*cadmium|cation diffusion facilitator|\bCzcD\b|\bZntA\b|\bCadA\b'),
 ('Chromate (chrA)', 'toxic metal', 'K07240', r'chromate'),
 ('Tellurite (ter, teh)', 'toxic metal', 'K03304 K16868 K05791 K05792 K05793 K05794 K05795', r'tellur'),
 ('Nickel/cobalt export (rcnA, nreB)', 'toxic metal', 'K07785 K08970 K23240 K21886', r'nickel/cobalt efflux|\bRcnA\b|nickel resistance'),
 ('Metallothionein (smtA, bmtA)', 'toxic metal', 'K21904', r'metallothionein'),
 ('Iron storage and Fur regulation', 'micronutrient', 'K02217 K03594 K03711 K09826', r'\bferritin\b|bacterioferritin'),
 ('Iron uptake (siderophores, TonB, Feo)', 'micronutrient', 'K03832 K03561 K03559 K04759 K16090 K16089 K02014 K19611', r'siderophore|ferrichrome|ferric (enterobactin|aerobactin|hydroxamate|citrate)|enterobactin|vibriobactin|ferrous iron transport'),
 ('Mn/Zn uptake (mnt, sit, znu, Zur)', 'micronutrient', 'K03322 K11601 K11602 K11603 K19973 K19975 K19976 K11604 K11707 K11708 K11709 K11710 K09815 K09816 K09817 K02076 K09823 K12950', r'manganese transport|zinc (uptake|transport)|\bZnuA\b'),
 ('Ni/Co uptake (nik, cbi, nixA)', 'micronutrient', 'K15584 K15585 K02006 K02007 K02008 K02009 K07241', r'nickel (uptake|transport|ABC)|cobalt (uptake|transport)'),
 ('Solvent/multidrug RND efflux (acrAB-tolC, mexAB-oprM, ttg, srp)', 'solvent stress', 'K03585 K18138 K18139 K12340 K18131 K18135', r'solvent efflux|toluene efflux|\bTtg[ABC]\b|\bSrp[ABC]\b'),
 ('Membrane cis-trans isomerase (cti)', 'solvent stress', '', r'fatty acid cis.?trans.?isomerase|\bCti\b'),
 ('Phage-shock membrane stress (psp)', 'solvent stress', 'K03969 K03970 K03971 K03973 K03974', r'phage shock protein'),
 ('Oxidative stress (kat, sod, ahp, dps, osmC)', 'oxidative stress', 'K03781 K03782 K04564 K04565 K04047 K24119 K24126 K03387 K04063 K09825 K13639 K00432 K11065', None),
 ('Protein-folding stress (dnaK, groEL, ibpA, clpB, htpG)', 'proteotoxic stress', 'K04043 K04077 K04080 K04081 K03695 K04079', None),
]
if __name__ == '__main__':
    G = pd.read_pickle(D + 'genemap.pkl')
    cds = (G['type'] == 'CDS').values
    kos = G.KO.fillna('').str.split(',')
    prod = G['product'].fillna('')
    taken = np.zeros(len(G), bool); rows = []
    # peptidyl-prolyl cis-trans isomerases must not match cti; GET3/arsA-like eukaryotic tail-anchoring ATPases stay in (arsA KO is shared)
    excl = prod.str.contains(r'peptidyl-prolyl', case=False).values
    for cat, grp, k, rx in SETS:
        ks = set(k.split())
        bk = kos.apply(lambda x: bool(ks.intersection(x))).values if ks else np.zeros(len(G), bool)
        bp = prod.str.contains(rx, case=False, regex=True).values if rx else np.zeros(len(G), bool)
        sel = (bk | bp) & cds & ~taken & ~excl
        taken |= sel
        idx = np.where(sel)[0]
        rows.append(pd.DataFrame({'row': idx, 'gene_id': G.gene_id.values[idx], 'contig': G.contig.values[idx], 'bin': G.bin.fillna('').values[idx],
                                  'category': cat, 'group': grp, 'evidence': np.where(bk[idx] & bp[idx], 'KO+product', np.where(bk[idx], 'KO', 'product')),
                                  'KO': G.KO.fillna('').values[idx], 'product': prod.values[idx]}))
    M = pd.concat(rows); M.to_csv(O + 'members.tsv', sep='\t', index=False)
    print(M.groupby(['group', 'category', 'evidence']).size().unstack(fill_value=0).to_string())
    print('contigs:', M.contig.nunique())
    # product sanity: top products per category
    for c, d in M.groupby('category'): print(c, '|', d['product'].value_counts().head(4).to_dict())
