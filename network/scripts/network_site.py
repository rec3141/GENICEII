"""Site data for the Network page (Food web). Run after nodes.py, export.py, the four methods, consensus.py and layout.py.
site/network_{dna,rna}.json, one per layer: nodes (t-SNE position, domain, kind, status, role, label) and the FULL edge sets of all
four methods as column arrays (a, b node indices; s sign; m / x bit masks of the methods agreeing with / opposing that sign;
fw FlashWeave |weight| x100; sps / ses SPRING / SPIEC-EASI stability x100; spb / seb their |beta|
x1000; rho residual-Spearman |rho| x100; q -log10 Spearman q x10). The score columns are sparse: one value per edge that
the method reports, in edge order (the page expands them); 'strong' is recomputed on the page from these scores. Spearman edges that no other method reports are kept only at
|rho| >= SP_FLOOR. Also FlashWeave's meta-variable ('driver') edges, method counts and overlaps.
site/network.json: methods, permissivity steps (STEPS, filtered client-side), edges per step, the known-pair check, top partners per
eukaryote and per giant virus / virophage (consensus edges). out/permissivity_steps.tsv: known pairs recovered at each step."""
import os, sys, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
import common

BIT = {m: 1 << i for i, m in enumerate(METHODS)}
SITE = os.path.join(REPO, 'site', 'network.json')
def r(x, n=3):
    try: x = float(x)
    except Exception: return None
    return round(x, n) if np.isfinite(x) else None

# trophic role of eukaryotic MAGs from their taxonomy call (data/ee_taxonomy.tsv: call, protein lineage, 18S), first match wins
ROLES = [(r'Choano|Acanthoecid', 'grazer', 'choanoflagellate'), (r'Ciliophora|Oligohymeno|Spirotrich', 'grazer', 'ciliate'),
         (r'Calanidae|Copepoda|Arthropoda|Metazoa', 'grazer', 'copepod or other animal'), (r'Chrysophyc|Chromulin|Ochromonad', 'mixotroph', 'chrysophyte'),
         (r'Bacillario|Diatom|Minidiscus|Attheya|Chaetoc|Thalassios|Navicul|Staurosira', 'alga', 'diatom'),
         (r'Micromonas|Bathycoccus|Picochlorum|Chlorophyt|Mamiell', 'alga', 'green alga'), (r'Hapto|Phaeocystis', 'alga', 'haptophyte'),
         (r'Laminariales|Phaeophyc', 'alga', 'brown alga'), (r'Cercozoa', 'protist', 'cercozoan'), (r'Ascomycota|Fungi', 'fungus', 'fungus')]
NA = nodes_all()
ET = pd.read_csv(common.D + 'ee_taxonomy.tsv', sep='\t', keep_default_na=False).set_index('mag')
def role(mid):
    # MAGs that PR #18 calls 'doubtful, not a genome' get no trophic role (their taxonomy describes a mixture)
    if str(NA.genome_status.get(mid, '') if 'genome_status' in NA else '').startswith('doubtful'): return '', ''
    if mid not in ET.index: return '', ''
    x = ET.loc[mid]
    full = ' '.join(str(x.get(c, '')) for c in ('call', 'protein_phylum', 'protein_class', 'protein_order', 'protein_family', 'protein_genus', 'rrna_18S'))
    for s in (str(x.get('call', '')), full):                   # the taxonomy call first, then the supporting lineages
        for rx, ro, grp in ROLES:
            if re.search(rx, s): return ro, grp + ('?' if x.get('confidence') == 'low' else '')
    return 'protist', ''

def nlabel(i, x):
    if x.domain == 'prokaryote': return f'{x.label} ({i})'
    if x.domain == 'organelle' and x.kind == 'plastid bin': return f'{x.label} ({i})'
    return x.label

# ---------------------------------------------------------------- permissivity steps (one ordinal slider on the page)
# Every step filters the stored results of the four methods; nothing is recomputed. A method SUPPORTS an edge at a step when it
# reports the edge with the edge's sign and its score passes the step's level:
#   'strong': SPRING and SPIEC-EASI stability (share of StARS subsamples with the edge) >= STRONG_STAB; residual Spearman
#             |rho| >= STRONG_RHO and q < STRONG_Q; FlashWeave: reported (its own conditional test at alpha 0.01 is the strict part);
#   'any':    reported (Spearman: q < 0.05).
# A clause {k, lvl, design} passes an edge when at least k of the methods switched on support it, no switched-on method reports
# the opposite sign, and (design clauses) one supporting method accounts for the design (FlashWeave with meta-variables, or the
# residual Spearman). With fewer than k methods switched on, k becomes the number switched on. Each step is the union of its own
# clause and every stricter step's, so the steps are nested: every link of a step is also in all more permissive steps.
# Step 'cons' equals the consensus of consensus.py.
STRONG_STAB, STRONG_RHO, STRONG_Q = 0.9, 0.6, 1e-3
SP_FLOOR = 0.4          # Spearman-only edges with |rho| below this are not exported (RNA has 69,711 Spearman edges at q < 0.05)
_C = [(dict(k=4, lvl='strong', design=False), 'All four methods, strong scores'),
      (dict(k=3, lvl='strong', design=False), 'Three methods with strong scores'),
      (dict(k=3, lvl='any', design=False), 'Three methods'),
      (dict(k=2, lvl='strong', design=True), 'Three methods, or two with strong scores (one design-aware)'),
      (dict(k=2, lvl='any', design=True), 'Two methods, one design-aware (consensus)'),
      (dict(k=2, lvl='any', design=False), 'Any two methods'),
      (dict(k=1, lvl='strong', design=False), 'Any two methods, or one with a strong score'),
      (dict(k=1, lvl='any', design=False), 'Any one method')]
STEPS = [dict(label=l, clauses=[c for c, _ in _C[:i + 1]], **({'cons': True} if 'consensus' in l else {})) for i, (_, l) in enumerate(_C)]
DEFAULT_STEP = 4

def read_methods(layer, N):
    meta_cols = list(pd.read_csv(os.path.join(CACHE, f'{layer}_meta.tsv'), sep='\t', index_col=0).columns)
    F = pd.read_csv(os.path.join(OUT, f'flashweave_{layer}.edgelist'), sep='\t', comment='#', header=None, names=['a', 'b', 'weight'])
    F = F[~(F.a.isin(meta_cols) | F.b.isin(meta_cols))].copy(); F['sign'] = np.sign(F.weight).astype(int)
    F[['a', 'b']] = [edge_key(a, b) for a, b in zip(F.a, F.b)]
    M = {'flashweave': F}
    for m in ('spring', 'spieceasi', 'spearman'):
        M[m] = pd.read_csv(os.path.join(OUT, f'edges_{m}_{layer}.tsv'), sep='\t')
        M[m][['a', 'b']] = [edge_key(a, b) for a, b in zip(M[m].a, M[m].b)]     # same key order for R and Python output
    T = None
    for m, E in M.items():
        cols = {'sign': f'{m}_sign', 'weight': f'{m}_w'}
        if 'stability' in E: cols['stability'] = f'{m}_stab'
        if m == 'spearman': cols['q'] = 'spearman_q'
        e = E.set_index(['a', 'b'])[list(cols)].rename(columns=cols)
        T = e if T is None else T.join(e, how='outer')
    for m in METHODS: T[f'{m}_sign'] = T[f'{m}_sign'].fillna(0).astype(int)
    S = T[[f'{m}_sign' for m in METHODS]]
    pos, neg = (S > 0).sum(1), (S < 0).sum(1)
    T['sign'] = np.where(pos >= neg, 1, -1)
    T['mask'] = sum(np.where(S[f'{m}_sign'].to_numpy() == T.sign.to_numpy(), BIT[m], 0) for m in METHODS)
    T['anti'] = sum(np.where(S[f'{m}_sign'].to_numpy() == -T.sign.to_numpy(), BIT[m], 0) for m in METHODS)
    only_sp = T['mask'].eq(BIT['spearman']) & T.anti.eq(0)
    n_drop = int((only_sp & (T.spearman_w.abs() < SP_FLOOR)).sum())
    T = T[~(only_sp & (T.spearman_w.abs() < SP_FLOOR))]
    return T.reset_index(), n_drop

def strong_mask(T):
    # on the exported (rounded) values, exactly as the page recomputes it
    st = np.full(len(T), BIT['flashweave']) * (T.flashweave_sign != 0)
    st = st + BIT['spring'] * (np.round(T.spring_stab.fillna(0) * 100) >= STRONG_STAB * 100) + BIT['spieceasi'] * (np.round(T.spieceasi_stab.fillna(0) * 100) >= STRONG_STAB * 100)
    st = st + BIT['spearman'] * ((np.round(T.spearman_w.abs().fillna(0) * 100) >= STRONG_RHO * 100) & (q10(T.spearman_q.to_numpy(float)) > -np.log10(STRONG_Q) * 10))
    return st.astype(int).to_numpy()

def popcount(v): return np.array([bin(int(x)).count('1') for x in v])
DESIGN_BITS = BIT['flashweave'] | BIT['spearman']
def passes(T, strong, step, bits=sum(BIT.values())):
    n_on = bin(bits).count('1'); m = T['mask'].to_numpy() & bits; clean = (T.anti.to_numpy() & bits) == 0; out = np.zeros(len(T), bool)
    for c in step['clauses']:
        sup = m & strong if c['lvl'] == 'strong' else m
        ok = (popcount(sup) >= min(c['k'], n_on)) & clean
        if c['design']: ok &= (sup & DESIGN_BITS) > 0
        out |= ok
    return out

J = {'methods': [{'key': m, 'label': METHOD_LABEL[m], 'bit': BIT[m], 'design_aware': m in ('flashweave', 'spearman')} for m in METHODS],
     'params': {'prev': PREV, 'min_ppm': MIN_PPM, 'prev_forced': PREV_FORCED, 'strong_stab': STRONG_STAB, 'strong_rho': STRONG_RHO,
                'strong_q': STRONG_Q, 'sp_floor': SP_FLOOR}, 'steps': STEPS, 'default_step': DEFAULT_STEP, 'layers': {}}
DRV = pd.read_csv(os.path.join(OUT, 'drivers.tsv'), sep='\t')
SUM = pd.read_csv(os.path.join(OUT, 'method_summary.tsv'), sep='\t'); OVL = pd.read_csv(os.path.join(OUT, 'method_overlap.tsv'), sep='\t')
KNOWN = [('contig_0000564', 'EE07'), ('contig_0000623', 'EE07'), ('contig_0049668', 'EE07'), ('contig_0000362', 'EE04'), ('contig_0000610', 'EE22'), ('P01', 'EE05')]
REC = []
q10 = lambda q: np.where(np.isfinite(q), np.minimum(np.round(-np.log10(np.maximum(q, 1e-300)) * 10), 3000), 0).astype(int)
for layer in LAYERS:
    X, lib, N = load_layer(layer)
    pos = pd.read_csv(os.path.join(OUT, f'layout_{layer}.tsv'), sep='\t', index_col=0)
    share = X.mul(1.0 if layer == 'rna' else DNA_UNIT).div(lib['size'], axis=0) * 1e6
    ids = list(N.index); ix = {k: i for i, k in enumerate(ids)}
    nodes = []
    for i, x in N.iterrows():
        ro, grp = role(i) if x.domain == 'eukaryote' else ('', '')
        st = []
        if str(x.get('doubtful', '')) == 'True': st.append('doubtful')
        if str(x.get('multi', '')) == 'True': st.append('multi-genome')
        if layer == 'rna' and x.domain in ('eukaryote', 'organelle') and str(x.get('rna_thin', '')) == 'True': st.append('thin RNA profile')
        nodes.append({'id': i, 'label': nlabel(i, x), 'domain': x.domain, 'kind': x.kind, 'status': st, 'ee_status': (x.genome_status if isinstance(x.get('genome_status'), str) and x.genome_status else x.status) if x.domain == 'eukaryote' else '',
                      'caution': x.get('caution', '') if isinstance(x.get('caution', ''), str) else '', 'role': ro, 'group': grp,
                      'x': r(pos.x[i], 4), 'y': r(pos.y[i], 4), 'ab': r(share[i].median(), 2), 'mean': r(share[i].mean(), 2), 'prev': r((X[i] > 0).mean(), 2)})
    T, n_drop = read_methods(layer, N)
    strong = strong_mask(T)
    i100 = lambda v, f: np.round(v.fillna(0).abs() * f).astype(int).tolist()
    def sp(m, v, f):           # sparse: one value per edge that method reports (either sign), in edge order
        has = T[f'{m}_sign'].to_numpy() != 0
        return np.round(v[has].fillna(0).abs() * f).astype(int).tolist()
    edges = {'a': [ix[a] for a in T.a], 'b': [ix[b] for b in T.b], 's': T.sign.astype(int).tolist(), 'm': T['mask'].astype(int).tolist(),
             'x': T.anti.astype(int).tolist(),
             'fw': sp('flashweave', T.flashweave_w, 100), 'sps': sp('spring', T.spring_stab, 100), 'spb': sp('spring', T.spring_w, 1000),
             'ses': sp('spieceasi', T.spieceasi_stab, 100), 'seb': sp('spieceasi', T.spieceasi_w, 1000), 'rho': sp('spearman', T.spearman_w, 100),
             'q': q10(T.spearman_q.to_numpy(float))[T.spearman_sign.to_numpy() != 0].tolist()}
    counts = [int(passes(T, strong, s).sum()) for s in STEPS]
    key = {edge_key(a, b): k for k, (a, b) in enumerate(zip(T.a, T.b))}
    for si, s in enumerate(STEPS):
        ok = passes(T, strong, s)
        for v, h in KNOWN:
            k = key.get(edge_key(v, h)); REC.append(dict(layer=layer, step=si + 1, label=s['label'], edges=counts[si], pair=f'{v}-{h}', found=bool(k is not None and ok[k])))
    dr = DRV[(DRV.layer == layer) & DRV.node.isin(ids)]
    L = {'n_libs': int(len(lib)), 'nodes': nodes, 'edges': edges, 'step_edges': counts, 'sp_dropped': n_drop,
         'drivers': [[ix[d.node], d.driver, r(d.weight, 3)] for d in dr.itertuples()],
         'summary': SUM[SUM.layer == layer].drop(columns='layer').replace({np.nan: None}).to_dict('records'),
         'overlap': OVL[OVL.layer == layer].drop(columns='layer').replace({np.nan: None}).to_dict('records'),
         'by_domain': N.domain.value_counts().to_dict()}
    json.dump(L, open(os.path.join(REPO, 'site', f'network_{layer}.json'), 'w'), separators=(',', ':'), allow_nan=False)
    J['layers'][layer] = {k: L[k] for k in ('n_libs', 'step_edges', 'sp_dropped', 'summary', 'overlap', 'by_domain')} | {'n_nodes': len(nodes)}
    print(layer, len(nodes), 'nodes;', len(T), 'edges exported (', n_drop, 'Spearman-only edges below |rho|', SP_FLOOR, 'dropped); per step:', counts)
REC = pd.DataFrame(REC); REC.to_csv(os.path.join(OUT, 'permissivity_steps.tsv'), sep='\t', index=False)
print(REC.pivot_table(index=['layer', 'step', 'edges'], columns='pair', values='found', aggfunc='first').astype(int).to_string())
# known pairs
K = pd.read_csv(os.path.join(OUT, 'known_pairs.tsv'), sep='\t', keep_default_na=False)
J['known'] = [{**{k: (bool(v) if k in ('in_layer', 'consensus') else v) for k, v in x.items() if k in ('layer', 'a', 'b', 'pair', 'in_layer', 'support', 'consensus', 'eukaryote_partners')},
               'methods': [m for m in METHODS if int(x[m]) != 0], 'signs': {m: int(x[m]) for m in METHODS}} for x in K.to_dict('records')]
# top partners per eukaryote and per eukaryotic virus: consensus edges, both layers
C = pd.read_csv(os.path.join(OUT, 'edges_consensus.tsv'), sep='\t')
lab = {i: nlabel(i, x) for i, x in NA.iterrows()}
focus = [i for i, x in NA.iterrows() if x.domain == 'eukaryote' or x.kind in ('giant virus', 'virophage')]
P = []
for f in focus:
    c = C[(C.a == f) | (C.b == f)].copy()
    if c.empty: continue
    c['other'] = np.where(c.a == f, c.b, c.a)
    rows = []
    for o, g in c.groupby('other'):
        sg = int(np.sign(g.sign.sum())) or int(g.sign.iloc[0])
        rows.append({'id': o, 'label': lab.get(o, o), 'domain': NA.domain.get(o, ''), 'sign': sg, 'n_layers': int(g.layer.nunique()),
                     'support': int(g.support.max()), 'layers': {l: gg.methods.iloc[0].split(',') for l, gg in g.groupby('layer')}})
    rows.sort(key=lambda d: (-d['n_layers'], -d['support'], d['sign'] < 0, d['label']))
    P.append({'id': f, 'label': lab[f], 'domain': NA.domain[f], 'kind': NA.kind[f], 'role': role(f)[1] if NA.domain[f] == 'eukaryote' else '', 'n': len(rows), 'partners': rows[:8]})
J['partners'] = P
json.dump(J, open(SITE, 'w'), separators=(',', ':'), allow_nan=False)
print(len(P), 'focus nodes with partners;', {f: os.path.getsize(os.path.join(REPO, 'site', f)) // 1000 for f in ('network.json', 'network_dna.json', 'network_rna.json')}, 'kB')
