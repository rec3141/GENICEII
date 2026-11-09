"""site/network.json for the Network page (Food web): per layer the nodes (t-SNE position, domain, kind, status, role, label),
every edge any conditional method reports plus every edge with at least two methods (with the agreeing methods as a bit mask),
FlashWeave's meta-variable ('driver') edges, method counts and overlaps; the known-pair check; top partners per eukaryote and per
giant virus / virophage. Run after nodes.py, export.py, the four methods, consensus.py and layout.py."""
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

J = {'methods': [{'key': m, 'label': METHOD_LABEL[m], 'bit': BIT[m], 'design_aware': m in ('flashweave', 'spearman')} for m in METHODS],
     'params': {'prev': PREV, 'min_ppm': MIN_PPM, 'prev_forced': PREV_FORCED}, 'layers': {}}
ALL = pd.read_csv(os.path.join(OUT, 'edges_all.tsv'), sep='\t')
DRV = pd.read_csv(os.path.join(OUT, 'drivers.tsv'), sep='\t')
SUM = pd.read_csv(os.path.join(OUT, 'method_summary.tsv'), sep='\t'); OVL = pd.read_csv(os.path.join(OUT, 'method_overlap.tsv'), sep='\t')
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
        nodes.append({'id': i, 'label': nlabel(i, x), 'domain': x.domain, 'kind': x.kind, 'status': st, 'ee_status': (x.genome_status if isinstance(x.get('genome_status'), str) and x.genome_status else x.status) if x.domain == 'eukaryote' else '',
                      'caution': x.get('caution', '') if isinstance(x.get('caution', ''), str) else '', 'role': ro, 'group': grp,
                      'x': r(pos.x[i], 4), 'y': r(pos.y[i], 4), 'ab': r(share[i].median(), 2), 'mean': r(share[i].mean(), 2), 'prev': r((X[i] > 0).mean(), 2)})
    E = ALL[ALL.layer == layer]
    edges = []
    for e in E.itertuples():
        mask = sum(BIT[m] for m in METHODS if getattr(e, f'{m}_sign') == e.sign)
        anti = sum(BIT[m] for m in METHODS if getattr(e, f'{m}_sign') == -e.sign)
        edges.append([ix[e.a], ix[e.b], int(e.sign), mask, anti, r(e.spring_w, 3), r(e.flashweave_w, 3), r(e.spearman_w, 2)])
    dr = DRV[(DRV.layer == layer) & DRV.node.isin(ids)]
    J['layers'][layer] = {'n_libs': int(len(lib)), 'nodes': nodes, 'edges': edges,
                          'drivers': [[ix[d.node], d.driver, r(d.weight, 3)] for d in dr.itertuples()],
                          'summary': SUM[SUM.layer == layer].drop(columns='layer').replace({np.nan: None}).to_dict('records'),
                          'overlap': OVL[OVL.layer == layer].drop(columns='layer').replace({np.nan: None}).to_dict('records'),
                          'by_domain': N.domain.value_counts().to_dict()}
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
print({l: (len(v['nodes']), len(v['edges']), len(v['drivers'])) for l, v in J['layers'].items()}, len(P), 'focus nodes with partners;', os.path.getsize(SITE) // 1000, 'kB')
