"""Pathway section data: MAG × pathway completeness/expression heatmap, community time courses, host attribution per gene family."""
import sys, json; sys.path.insert(0, '.')
from common import *
P = '../out/pathways/'
m = meta(); S = samples(); m = m.loc[S]
defs = pd.read_csv(P + 'defs.tsv', sep='\t', index_col=0)
C = pd.read_csv(P + 'completeness.tsv', sep='\t', index_col=0)
E = pd.read_pickle(P + 'expr_mag.pkl')
EC = pd.read_csv(P + 'expr_community.tsv', sep='\t', index_col=0)
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
kj = pd.read_csv(D + 'kaiju_contigs.tsv', sep='\t', usecols=['contig_id', 'lineage']).set_index('contig_id')['lineage']
def lab(b):
    if b == 'unbinned': return 'unbinned contigs'
    t = gt.get(b.replace('dastool-', ''), '')
    d = dict(x.split('__', 1) for x in t.split(';') if '__' in x)
    for k, w in (('g', ''), ('f', ' (f)'), ('o', ' (o)'), ('c', ' (c)')):
        if d.get(k) and not (k != 'g' and kaiju_bin_label(b)): return d[k] + w
    return kaiju_bin_label(b) or 'unclassified'
def phylum(b):
    t = gt.get(b.replace('dastool-', ''), ''); d = dict(x.split('__', 1) for x in t.split(';') if '__' in x)
    return (d.get('c') or d.get('p') or '')
oil = (m.phase == 'oiled') & m.medium.isin(['SW', 'SI']); ctl = (m.phase == 'control') & m.medium.isin(['SW', 'SI'])
paths = [p for p in defs.index if p in EC.index and EC.loc[p].max() > 0]
# MAG ranking: total expression of degradation pathways in oiled samples, among MAGs with any pathway >= 0.5 complete
deg = [p for p in paths if defs.loc[p, 'family'] != 'C1']
Ebin = E.reset_index()
tot_oil = Ebin[Ebin.pathway.isin(deg)].groupby('bin')[list(m.index[oil])].sum().mean(1)
cand = [b for b in tot_oil.sort_values(ascending=False).index if b in C.index and (C.loc[b, deg] >= 0.5).any()]
rows = cand[:48]
heat = []
for b in rows:
    for p in paths:
        c = float(C.loc[b, p]) if p in C.columns and b in C.index else 0.0
        if (p, b) in E.index:
            v = E.loc[(p, b)]
            mo, mc = v[oil].mean(), v[ctl].mean()
            heat.append({'bin': b.replace('dastool-', ''), 'p': p, 'c': round(c, 2), 'oil': round(float(mo), 2), 'ctl': round(float(mc), 2),
                         'lfc': round(float(np.log2((mo + 0.5) / (mc + 0.5))), 2)})
        elif c > 0:
            heat.append({'bin': b.replace('dastool-', ''), 'p': p, 'c': round(c, 2), 'oil': 0, 'ctl': 0, 'lfc': 0})
out = {'defs': {p: {'name': defs.loc[p, 'name'], 'family': defs.loc[p, 'family'], 'def': defs.loc[p, 'definition']} for p in paths},
       'rows': [{'bin': b.replace('dastool-', ''), 'label': lab(b), 'class': phylum(b)} for b in rows],
       'heat': heat,
       'community': {p: [round(float(v), 2) for v in EC.loc[p, S]] for p in paths}}
# gene-family host attribution
hc = pd.read_csv('../out/hc_cpm.tsv', sep='\t', index_col=0)
mem = pd.read_pickle('../out/hc_members.pkl')
X, g = genes()
# CANT-HYD (curated hydrocarbon HMMs, trusted cutoffs) as an independent family set
ch = pd.read_csv('../canthyd/canthyd_hits.tsv', sep='\t')
gid = pd.Series(np.arange(len(g)), index=g.gene_id.values)
cds = (g['type'] == 'CDS').values; tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
for h, d in ch.groupby('hmm'):
    idx = gid.reindex(d.gene_id).dropna().astype(int).values
    if len(idx) == 0: continue
    k = 'CANT-HYD ' + h; mem[k] = idx
    hc[k] = np.asarray(X[idx].sum(0)).ravel() / tot * 1e6
fam = {}
col = {s: i for i, s in enumerate(S)}
for k, idx in mem.items():
    if len(idx) == 0: continue
    M = X[idx].toarray()
    res = {}
    for se in ['2024-25', '2025-26']:
        for grp, mask in (('oiled', oil), ('control', ctl)):
            ss = [col[s] for s in m.index[mask & (m.season == se)]]
            v = M[:, ss].sum(1); t = v.sum()
            if t == 0: continue
            host = pd.Series(v, index=[(g.bin.values[i] or 'unbinned') if g.bin.values[i] else 'kaiju:' + str(kj.get(g.contig.values[i], '')).split(';')[-2].strip() if ';' in str(kj.get(g.contig.values[i], '')) else 'unbinned' for i in idx])
            host = host.groupby(level=0).sum().sort_values(ascending=False)
            top = [{'host': h.replace('dastool-', ''), 'label': (lab(h) if not h.startswith('kaiju:') else h[6:] + ' (unbinned, Kaiju)'), 'pct': round(100 * x / t, 1)} for h, x in host.head(5).items()]
            res[f'{se} {grp}'] = {'pairs': int(t), 'top': top}
    fam[k] = {'cpm': [round(float(v), 2) for v in hc.loc[S, k]], 'hosts': res, 'n_genes': int(len(idx))}
# BACI (2025-26 seawater, medians): (EXP post / EXP pre) / (CTL post / CTL pre)
sw = m.medium == 'SW'; s26 = m.season == '2025-26'
grpI = {k: m.index[s26 & sw & (m.group == g) & (m.day >= 1 if post else m.day < 0)] for k, (g, post) in
        {'ep': ('EXP', True), 'eb': ('EXP', False), 'cp': ('CTL', True), 'cb': ('CTL', False)}.items()}
for k in fam:
    v = hc[k]
    med = {g: float(v.reindex(ix).median()) + 0.5 for g, ix in grpI.items()}
    fam[k]['baci26'] = round((med['ep'] / med['eb']) / (med['cp'] / med['cb']), 2)
out['families'] = fam
# potential vs realised across all MAGs
pr = []
Cb = C.drop(index=['unbinned'], errors='ignore')
for p in paths:
    if p not in Cb.columns: continue
    enc = Cb.index[Cb[p] >= 0.75]
    ex = []
    for b in enc:
        if (p, b) in E.index: ex.append((b, float(E.loc[(p, b)][oil].mean())))
    ex = sorted(ex, key=lambda t: -t[1]); tot = sum(v for _, v in ex) or 1
    top = ex[0] if ex else None
    pr.append({'p': p, 'encoded': int(len(enc)), 'expressed': int(sum(v > 1 for _, v in ex)),
               'top_share': round(100 * top[1] / tot, 1) if top else 0, 'top': lab(top[0]) + ' · ' + top[0].replace('dastool-', '') if top else '',
               'classes': pd.Series([phylum(b) for b in enc]).value_counts().head(4).to_dict()})
out['potential'] = pr
json.dump(out, open(P + 'site_pathways.json', 'w'), separators=(',', ':'))
json.dump(rows, open(P + 'key_bins.json', 'w'))
print(len(rows), len(heat), len(fam))
print([ (r['bin'], r['label']) for r in out['rows'][:20]])
