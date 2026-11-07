"""Biodegradation pathway catalogue: per-MAG completeness (KEGG module logic on eggNOG KOs + product rules)
and expression (CPM of CDS pairs) per sample; writes ../out/pathways/*."""
import sys, os, re, glob, json; sys.path.insert(0, '.')
from common import *
O = '../out/pathways/'; os.makedirs(O, exist_ok=True)
MOD = {}
for f in sorted(glob.glob(D + 'modules/*.txt')):
    t = open(f).read()
    MOD[os.path.basename(f)[:-4]] = (re.search(r'^NAME\s+(.*)$', t, re.M)[1], re.search(r'^DEFINITION\s+(.*)$', t, re.M)[1])
CUSTOM = {  # id: (name, definition, family)
    'C_alkB': ('Alkane hydroxylase AlkB (+rubredoxin reductase)', 'K00496 -K05297', 'alkane'),
    'C_pca': ('Protocatechuate ortho-cleavage (beta-ketoadipate)', 'K00448+K00449 K01857 K01607 (K01055,K14727) K01031+K01032 K07823', 'lower aromatic'),
    'C_gen': ('Gentisate pathway', 'K00450 K01800 K16163', 'lower aromatic'),
    'C_hmg': ('Homogentisate pathway', 'K00451 K01800 K01555', 'lower aromatic'),
    'C_mar': ('Maleylacetate reductase (hydroxyquinol / chlorocatechol branch)', 'K00217', 'lower aromatic'),
}
PROD = {  # product-defined single-step functions
    'P_cyp153': ('Cytochrome P450 alkane hydroxylase (CYP153)', r'(?=.*alkane)(?=.*(p450|cytochrome))', 'alkane'),
    'P_lada': ('Long-chain alkane monooxygenase LadA / AlmA', r'long-chain alkane monooxygenase|\blada\b|\balma\b', 'alkane'),
    'P_rhd': ('Aromatic ring-hydroxylating dioxygenase (any)', r'ring-hydroxylating dioxygenase|aromatic.*dioxygenase.*alpha|rieske.*dioxygenase', 'upper aromatic'),
    'P_ass': ('Alkylsuccinate / benzylsuccinate synthase (anaerobic)', r'alkylsuccinate synthase|benzylsuccinate synthase|\bassa\b|\bbssa\b', 'anaerobic'),
}
FAM = {'M00174': 'C1', 'M00418': 'anaerobic', 'M00541': 'anaerobic', 'M00540': 'anaerobic', 'M00568': 'lower aromatic', 'M00569': 'lower aromatic',
       'M00533': 'lower aromatic', 'M00637': 'lower aromatic', 'M00638': 'lower aromatic', 'M00878': 'lower aromatic', 'M00623': 'lower aromatic',
       'M00636': 'lower aromatic', 'M00624': 'lower aromatic'}

def tokenize(s): return re.findall(r'K\d{5}|[(),+\- ]', s)
def parse(tok):
    """KEGG definition -> nested ('and'|'or'|'cx', [..]) tree; '-' marks optional components."""
    def expr(i, stop):
        steps, cur, alts, opt = [], [], [], False
        while i < len(tok) and tok[i] not in stop:
            t = tok[i]
            if t == ' ':
                if cur: alts.append(('cx', cur) if len(cur) > 1 else cur[0]); cur = []
                if alts: steps.append(('or', alts) if len(alts) > 1 else alts[0]); alts = []
                i += 1; continue
            if t == ',':
                if cur: alts.append(('cx', cur) if len(cur) > 1 else cur[0]); cur = []
                i += 1; continue
            if t == '+': i += 1; continue
            if t == '-': opt = True; i += 1; continue
            if t == '(':
                node, i = expr(i + 1, (')',)); i += 1
            else:
                node = ('k', t); i += 1
            cur.append(('opt', node) if opt else node); opt = False
        if cur: alts.append(('cx', cur) if len(cur) > 1 else cur[0])
        if alts: steps.append(('or', alts) if len(alts) > 1 else alts[0])
        return (('and', steps) if len(steps) > 1 else steps[0]), i
    return expr(0, ())[0]
def score(n, have):
    t, v = n
    if t == 'k': return 1.0 if v in have else 0.0
    if t == 'opt': return None
    if t == 'or': return max(s for s in (score(x, have) for x in v) if s is not None)
    sc = [s for s in (score(x, have) for x in v) if s is not None]
    return sum(sc) / len(sc) if sc else 1.0
def kos_of(n):
    t, v = n
    if t == 'k': return {v}
    if t == 'opt': return kos_of(v)
    return set().union(*[kos_of(x) for x in v])
def steps(n): return n[1] if n[0] == 'and' else [n]

DEFS = {**{k: (v[0], v[1], FAM.get(k, 'upper aromatic')) for k, v in MOD.items()}, **CUSTOM}
TREES = {k: parse(tokenize(d)) for k, (_, d, _) in DEFS.items()}

X, g = genes(); S = samples(); m = meta().loc[S]
cds = (g['type'] == 'CDS').values
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
ko = g['KO'].fillna('')
prod = g['product'].fillna('').str.lower()
# gene index per KO (only KOs we need)
need = set().union(*[kos_of(t) for t in TREES.values()])
gi = {}
for i, ks in zip(np.where(cds & ko.ne('').values)[0], ko.values[cds & ko.ne('').values]):
    for k in ks.split(','):
        if k in need: gi.setdefault(k, []).append(i)
for pid, (name, rx, fam) in PROD.items():
    gi[pid] = list(np.where(cds & prod.str.contains(rx, regex=True).values)[0])
    DEFS[pid] = (name, pid, fam); TREES[pid] = ('k', pid)
bins = g['bin'].values
allbins = sorted(set(bins) - {''})
# per pathway: per MAG completeness, gene list, expression per sample
comp, rows, expr_mag, expr_comm = {}, [], {}, {}
for pid, tree in TREES.items():
    ks = kos_of(tree)
    idx = sorted(set(i for k in ks for i in gi.get(k, [])))
    if not idx: comp[pid] = {}; expr_comm[pid] = np.zeros(len(S)); continue
    idx = np.array(idx)
    have_by_bin = {}
    for k in ks:
        for i in gi.get(k, []):
            have_by_bin.setdefault(bins[i] or 'unbinned', set()).add(k)
    comp[pid] = {b: round(score(tree, h), 3) for b, h in have_by_bin.items()}
    M = X[idx].toarray() / tot * 1e6
    expr_comm[pid] = M.sum(0)
    for b in set(bins[idx]):
        sel = bins[idx] == b
        expr_mag[(pid, b or 'unbinned')] = M[sel].sum(0)
    for i in idx:
        rows.append((pid, g.gene_id.values[i], g.contig.values[i], bins[i] or 'unbinned', ko.values[i], g['product'].values[i], int(g.length.values[i])))
C = pd.DataFrame(comp).fillna(0)   # bins × pathways
C.to_csv(O + 'completeness.tsv', sep='\t')
pd.DataFrame(rows, columns=['pathway', 'gene_id', 'contig', 'bin', 'KO', 'product', 'length']).to_csv(O + 'pathway_genes.tsv', sep='\t', index=False)
E = pd.DataFrame(expr_mag, index=S).T; E.index.names = ['pathway', 'bin']
E.to_pickle(O + 'expr_mag.pkl')
EC = pd.DataFrame(expr_comm, index=S).T; EC.to_csv(O + 'expr_community.tsv', sep='\t')
pd.DataFrame({k: {'name': v[0], 'definition': v[1], 'family': v[2]} for k, v in DEFS.items()}).T.to_csv(O + 'defs.tsv', sep='\t')
print('pathways', len(TREES), 'MAG×pathway with genes', len(E))
print((C.drop(index=['unbinned'], errors='ignore') >= 0.75).sum().sort_values(ascending=False).rename(lambda k: f'{k} {DEFS[k][0][:55]}').to_string())
