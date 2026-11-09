"""Site data for the Network page (Food web). Run after nodes.py, export.py, the five methods, consensus.py, layout.py and
forcelayout.py.
site/network_{layer}.json, one per layer: nodes (domain, kind, status, role, label), positions of every layout (t-SNE and one force
layout per method plus the consensus; -100..100, 1 decimal) and the edge sets of all five methods as column arrays (a, b node
indices; f flags, see below; fw FlashWeave |weight| x100; sps / ses SPRING /
SPIEC-EASI stability x100; rho residual-Spearman |rho| x100; pr Presence Pearson |r| x100; f packs, with M = number of methods,
bit 0 = positive sign, bits 1..M the mask of methods agreeing with that sign, the next M bits the mask of methods with the opposite
sign, the next M bits the methods whose score is strong there (netcommon.strong_bits). The score columns are sparse: one value per edge that the
method reports, in edge order (the page expands them). Edges that only the
residual Spearman reports are kept at |rho| >= SP_FLOOR (netcommon.edge_table). Also FlashWeave's meta-variable ('driver') edges, method counts and overlaps.
site/network.json: methods, layouts, layers, permissivity steps (STEPS, filtered client-side), edges per step, the known-pair check,
top partners per eukaryote and per giant virus / virophage (consensus edges). out/permissivity_steps.tsv: known pairs recovered at
each step."""
import os, sys, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netcommon import *
sys.path.insert(0, os.path.join(REPO, 'scripts'))
import common

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

# genera the page always labels: Cycloclasticus, the diesel-bloom genus of the site (its own page, 'cyclo')
FEATURED = ('Cycloclasticus',)
def genus(x):
    g = x.get('genus', '')
    if isinstance(g, str) and g: return re.sub(r'_[A-Z]$', '', g)
    m = re.match(r'([A-Z][a-z]+)', str(x.label)); return m[1] if m and x.label != 'unclassified' else ''

# possible same genome: two prokaryotic MAGs of the same GTDB genus (or, both without a GTDB genus, the same majority Kaiju
# genus) whose RNA profiles agree almost perfectly over the libraries where both are detected (Spearman rho >= SAME_RHO, at least
# SAME_N such libraries). Their links are flagged (the page draws them dashed and can hide them); nothing else changes.
SAME_RHO, SAME_N = 0.95, 10
def same_genome_pairs():
    from scipy.stats import spearmanr
    X, lib, N = load_layer('rna'); P = N[N.domain == 'prokaryote']
    def gkey(i):
        g = P.genus.get(i, '')
        if isinstance(g, str) and g: return 'g:' + re.sub(r'_[A-Z]$', '', g)
        k = common.kaiju_bin_label(i); return 'k:' + k.split(' (')[0] if k else ''
    G = pd.Series({i: gkey(i) for i in P.index}); G = G[G != '']
    out = {}
    for g, ids in G.groupby(G):
        ids = list(ids.index)
        for a in range(len(ids)):
            for b in range(a + 1, len(ids)):
                u, v = X[ids[a]].to_numpy(), X[ids[b]].to_numpy(); m = (u > 0) & (v > 0)
                if m.sum() < SAME_N: continue
                rho = spearmanr(u[m], v[m]).correlation
                if rho >= SAME_RHO: out[edge_key(ids[a], ids[b])] = (g[2:], round(float(rho), 3), int(m.sum()))
    return out
SAME = same_genome_pairs()
print(len(SAME), 'possible same-genome pairs (RNA rho >=', SAME_RHO, '):', sorted(SAME.items(), key=lambda kv: -kv[1][1])[:12])

def nlabel(i, x):
    if x.domain == 'prokaryote': return f'{x.label} ({i})'
    if x.domain == 'organelle' and x.kind == 'plastid bin': return f'{x.label} ({i})'
    return x.label

# ---------------------------------------------------------------- permissivity steps (one ordinal slider on the page)
# Every step filters the stored results of the five methods; nothing is recomputed. A method SUPPORTS an edge at a step when it
# reports the edge with the edge's sign and its score passes the step's level:
#   'strong': netcommon.strong_bits (SPRING / SPIEC-EASI stability >= 0.9; residual Spearman |rho| >= 0.6 and q < 0.001;
#             Presence Pearson |r| >= 0.9, q < 0.001, co-present in >= 30 libraries; FlashWeave: every edge, its own conditional
#             test at alpha 0.01 being the strict part);
#   'any':    reported (Spearman and Presence Pearson: q < 0.05).
# A clause {k, lvl, design} passes an edge when at least k of the methods switched on support it, no switched-on method reports
# the opposite sign, and (design clauses) one supporting method accounts for the design (FlashWeave with meta-variables, or the
# residual Spearman). With fewer than k methods switched on, k becomes the number switched on. Each step is the union of its own
# clause and every stricter step's, so the steps are nested: every link of a step is also in all more permissive steps.
# With five methods, three or even four supporters can all be design-blind (SPRING, SPIEC-EASI, Presence Pearson), so every
# clause down to the consensus carries the design condition; otherwise the union would put design-blind links into the consensus.
# Step 'cons' equals the consensus of consensus.py.
_C = [(dict(k=4, lvl='strong', design=True), 'Four or five methods, strong scores'),
      (dict(k=3, lvl='strong', design=True), 'Three methods with strong scores, one design-aware'),
      (dict(k=3, lvl='any', design=True), 'Three methods, one design-aware'),
      (dict(k=2, lvl='strong', design=True), 'Three methods, or two with strong scores (one design-aware)'),
      (dict(k=2, lvl='any', design=True), 'Two methods, one design-aware (consensus)'),
      (dict(k=2, lvl='any', design=False), 'Any two methods'),
      (dict(k=1, lvl='strong', design=False), 'Any two methods, or one with a strong score'),
      (dict(k=1, lvl='any', design=False), 'Any one method')]
STEPS = [dict(label=l, clauses=[c for c, _ in _C[:i + 1]], **({'cons': True} if 'consensus' in l else {})) for i, (_, l) in enumerate(_C)]
DEFAULT_STEP = 4

def popcount(v): return np.array([bin(int(x)).count('1') for x in v])
DESIGN_BITS = sum(BIT[m] for m in DESIGN_AWARE)
def passes(T, strong, step, bits=sum(BIT.values())):
    # the sign is the one the switched-on methods give; any disagreement among them drops the edge (as nwPass on the page)
    n_on = bin(bits).count('1'); ag = T['mask'].to_numpy() & bits; op = T.anti.to_numpy() & bits
    on = np.where(op == 0, ag, np.where(ag == 0, op, 0)); out = np.zeros(len(T), bool)
    for c in step['clauses']:
        sup = on & strong if c['lvl'] == 'strong' else on
        ok = popcount(sup) >= max(1, min(c['k'], n_on))
        if c['design']: ok &= (sup & DESIGN_BITS) > 0
        out |= ok
    return out

LAYOUTS = [('tsne', 't-SNE')] + [(m, METHOD_LABEL[m]) for m in METHODS] + [('consensus', 'Consensus')]   # forcelayout.py
# single-method line styling on the page: the score column of each method and its class limits (5 classes, same in every layer)
SCORE = {'flashweave': ('fw', 'FlashWeave |weight|', [0.3, 0.4, 0.5, 0.7]), 'spring': ('sps', 'SPRING stability', [0.4, 0.6, 0.8, 0.9]),
         'spieceasi': ('ses', 'SPIEC-EASI stability', [0.4, 0.6, 0.8, 0.9]), 'spearman': ('rho', 'residual Spearman |ρ|', [0.3, 0.4, 0.5, 0.6]),
         'pearson': ('pr', 'Presence Pearson |r|', [0.87, 0.9, 0.93, 0.96])}   # its edges start at |r| 0.85
J = {'methods': [{'key': m, 'label': METHOD_LABEL[m], 'bit': BIT[m], 'design_aware': m in DESIGN_AWARE, 'abbr': METHOD_ABBR[m],
                  'col': SCORE[m][0], 'score': SCORE[m][1], 'bins': SCORE[m][2]} for m in METHODS],
     'params': {'prev': PREV, 'min_ppm': MIN_PPM, 'prev_forced': PREV_FORCED, 'strong_stab': STRONG_STAB, 'strong_rho': STRONG_RHO,
                'strong_q': STRONG_Q, 'strong_pp_r': STRONG_PP_R, 'strong_pp_n': STRONG_PP_N, 'sp_floor': SP_FLOOR,
                'pp_n_min': PP_N_MIN, 'pp_q': PP_Q, 'pp_r': PP_R, 'pp_n': PP_N},
     'steps': STEPS, 'default_step': DEFAULT_STEP, 'layouts': [{'key': k, 'label': l} for k, l in LAYOUTS],
     'layer_list': [{'key': l, **{k: LAYER_CFG[l][k] for k in ('label', 'long', 'what')}} for l in sorted(LAYERS, key=lambda l: LAYER_CFG[l]['page'])],
     'layers': {}}
DRV = pd.read_csv(os.path.join(OUT, 'drivers.tsv'), sep='\t')
SUM = pd.read_csv(os.path.join(OUT, 'method_summary.tsv'), sep='\t'); OVL = pd.read_csv(os.path.join(OUT, 'method_overlap.tsv'), sep='\t')
KNOWN = [('contig_0000564', 'EE07'), ('contig_0000623', 'EE07'), ('contig_0049668', 'EE07'), ('contig_0000362', 'EE04'), ('contig_0000610', 'EE22'), ('P01', 'EE05')]
REC, SAMEC = [], []

def scaled(x, y):
    """One layout to -100..100 on its larger axis (aspect kept), 1 decimal, as [x..., y...] lists."""
    x, y = np.asarray(x, float), np.asarray(y, float); x, y = x - (x.max() + x.min()) / 2, y - (y.max() + y.min()) / 2
    f = 100 / max(np.abs(x).max(), np.abs(y).max(), 1e-9)
    return {'x': np.round(x * f, 1).tolist(), 'y': np.round(y * f, 1).tolist()}

MED_WORD = {'SW': 'seawater', 'SI': 'sea ice', 'other': ''}
ICE_RANK = {'TOP': (1, 'top'), 'HIGH': (1, 'top'), 'BOT': (3, 'bottom'), 'LOW': (3, 'bottom')}
def sample_order(lib):
    """Sample-picker order and labels: season > pool (control, diesel, then the rest) > medium (seawater, sea ice top, other ice,
    sea ice bottom) > day; samples outside the two pools (source sites, enrichments, blanks) at the end, by name."""
    rows = []
    for k, x in lib.iterrows():
        g = {'CTL': 0, 'EXP': 1}.get(x.group, 2); lay = str(x.get('layer', '') if isinstance(x.get('layer', ''), str) else '').upper()
        ir, iw = ICE_RANK.get(lay, (2, lay.lower())) if x.medium == 'SI' else ((0, '') if x.medium == 'SW' else (4, ''))
        d = pd.to_numeric(x.day, errors='coerce'); pool = g < 2 and x.medium in ('SW', 'SI')
        lab = ' · '.join(t for t in (x.season, ['control', 'diesel'][g] if g < 2 else '', (MED_WORD.get(x.medium, '') + (' ' + iw if iw else '')).strip(),
                                       f'day {int(d)}'.replace('-', '−') if pool and np.isfinite(d) else '') if t)
        if not pool: lab = f'{x.season} · {x["name"]}'
        rows.append((not pool, x.season, g, ir, d if np.isfinite(d) else 999, str(x['name']), k, lab))
    rows.sort()
    return [(r[6], r[7], str(r[5])) for r in rows]

def sig2(v):
    if v <= 0 or not np.isfinite(v): return 0
    v = float(f'{v:.2g}'); return int(v) if v >= 10 else v

for layer in LAYERS:
    X, lib, N = load_layer(layer)
    # per-sample node values for the sample picker (lazily loaded file): RNA reads per kb per million (the robust profile),
    # nanopore DNA parts per million of the library's aligned bases, Illumina DNA per million mapped read pairs; 2 significant digits
    if layer == 'rna': V = X.div(N.bp / 1e3, axis=1).div(lib['size'] / 1e6, axis=0); unit = 'reads per kb per million'
    else: V = rel_abund(X, lib, layer) * 1e6; unit = 'per million mapped read pairs' if layer == 'mgdna' else 'per million aligned bases'
    so = sample_order(lib)
    smp = {'unit': unit, 'libs': [{'id': k, 'label': l, 'name': nm} for k, l, nm in so],
           'v': [[sig2(v) for v in V.loc[k, N.index]] for k, _, _ in so]}
    json.dump(smp, open(os.path.join(REPO, 'site', f'network_{layer}_samples.json'), 'w'), separators=(',', ':'), allow_nan=False)
    pos = pd.read_csv(os.path.join(OUT, f'layout_{layer}.tsv'), sep='\t', index_col=0).loc[N.index]
    fl = pd.read_csv(os.path.join(OUT, f'forcelayout_{layer}.tsv'), sep='\t', index_col=0).loc[N.index]
    layouts = {'tsne': scaled(pos.x, pos.y)} | {k: scaled(fl[f'{k}_x'], fl[f'{k}_y']) for k, _ in LAYOUTS[1:]}
    share = rel_abund(X, lib, layer) * 1e6
    ids = list(N.index); ix = {k: i for i, k in enumerate(ids)}
    nodes = []
    for i, x in N.iterrows():
        ro, grp = role(i) if x.domain == 'eukaryote' else ('', '')
        st = []
        if str(x.get('doubtful', '')) == 'True': st.append('doubtful')
        if str(x.get('multi', '')) == 'True': st.append('multi-genome')
        if str(x.get('fragment', '')) == 'True': st.append('fragment')
        if str(x.get('contam_10_20', '')) == 'True': st.append('contamination 10–20%')
        if layer == 'rna' and x.domain in ('eukaryote', 'organelle') and str(x.get('rna_thin', '')) == 'True': st.append('thin RNA profile')
        nodes.append({'id': i, 'label': nlabel(i, x), 'domain': x.domain, 'kind': x.kind, 'status': st, 'ee_status': (x.genome_status if isinstance(x.get('genome_status'), str) and x.genome_status else x.status) if x.domain == 'eukaryote' else '',
                      'caution': x.get('caution', '') if isinstance(x.get('caution', ''), str) else '', 'role': ro, 'group': grp,
                      **({'genus': genus(x)} if x.domain == 'prokaryote' else {}), **({'feat': True} if x.domain == 'prokaryote' and genus(x) in FEATURED else {}),
                      'ab': r(share[i].median(), 2), 'mean': r(share[i].mean(), 2), 'prev': r((X[i] > 0).mean(), 2)})
    T, n_drop = edge_table(layer)
    strong = strong_bits(T)
    def sp(m, v, f=100):       # sparse: one value per edge that method reports (either sign), in edge order
        has = T[f'{m}_sign'].to_numpy() != 0
        return np.round(np.abs(np.nan_to_num(np.asarray(v, float)[has])) * f).astype(int).tolist()
    nb = len(METHODS)
    flags = (T.sign.to_numpy() > 0).astype(int) | (T['mask'].to_numpy() << 1) | (T.anti.to_numpy() << (1 + nb)) | (strong << (1 + 2 * nb))
    edges = {'a': [ix[a] for a in T.a], 'b': [ix[b] for b in T.b], 'f': flags.astype(int).tolist(),
             'fw': sp('flashweave', T.flashweave_w), 'sps': sp('spring', T.spring_stab), 'ses': sp('spieceasi', T.spieceasi_stab),
             'rho': sp('spearman', T.spearman_w), 'pr': sp('pearson', T.pearson_w)}
    edges['same'] = [k for k, (a, b) in enumerate(zip(T.a, T.b)) if (a, b) in SAME]       # rows flagged 'possible same genome'
    counts = [int(passes(T, strong, s).sum()) for s in STEPS]
    key = {edge_key(a, b): k for k, (a, b) in enumerate(zip(T.a, T.b))}
    cons = set(map(tuple, pd.read_csv(os.path.join(OUT, 'edges_consensus.tsv'), sep='\t').query('layer == @layer')[['a', 'b']].to_numpy()))
    for si, s in enumerate(STEPS):
        ok = passes(T, strong, s)
        if s.get('cons'):     # the consensus step must reproduce consensus.py exactly
            got = {(a, b) for a, b, o in zip(T.a, T.b, ok) if o}
            assert got == cons, (layer, len(got), len(cons), len(got - cons), len(cons - got))
        for v, h in KNOWN:
            k = key.get(edge_key(v, h)); REC.append(dict(layer=layer, step=si + 1, label=s['label'], edges=counts[si], pair=f'{v}-{h}', found=bool(k is not None and ok[k])))
    dr = DRV[(DRV.layer == layer) & DRV.node.isin(ids)]
    okc = passes(T, strong, STEPS[DEFAULT_STEP])
    SAMEC.append(dict(layer=layer, pairs_in_layer=sum(a in ix and b in ix for a, b in SAME), edges=len(edges['same']),
                      consensus_edges=int(sum(okc[k] for k in edges['same']))))
    L = {'n_libs': int(len(lib)), 'nodes': nodes, 'layouts': layouts, 'edges': edges, 'step_edges': counts, 'dropped': n_drop,
         'drivers': [[ix[d.node], d.driver, r(d.weight, 3)] for d in dr.itertuples()],
         'summary': SUM[SUM.layer == layer].drop(columns='layer').replace({np.nan: None}).to_dict('records'),
         'overlap': OVL[OVL.layer == layer].drop(columns='layer').replace({np.nan: None}).to_dict('records'),
         'by_domain': N.domain.value_counts().to_dict()}
    json.dump(L, open(os.path.join(REPO, 'site', f'network_{layer}.json'), 'w'), separators=(',', ':'), allow_nan=False)
    J['layers'][layer] = {k: L[k] for k in ('n_libs', 'step_edges', 'dropped', 'summary', 'overlap', 'by_domain')} | {'n_nodes': len(nodes)}
    print(layer, len(nodes), 'nodes;', len(T), 'edges exported (single-method edges dropped below the floors:', n_drop, '); per step:', counts)
print('possible same genome:', SAMEC)
pd.DataFrame([dict(a=a, b=b, genus=g, rho=r, n=n) for (a, b), (g, r, n) in SAME.items()]).to_csv(os.path.join(OUT, 'same_genome_pairs.tsv'), sep='\t', index=False)
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
print(len(P), 'focus nodes with partners;', {f: os.path.getsize(os.path.join(REPO, 'site', f)) // 1000 for f in ['network.json'] + [f'network_{l}{x}.json' for l in LAYERS for x in ('', '_samples')]}, 'kB')
