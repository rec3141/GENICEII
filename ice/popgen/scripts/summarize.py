"""Collect abundance + pop/*.json + gene FST -> popgen.json (site), popgen_summary.tsv, fst_outlier_genes.tsv."""
import json, glob, os, numpy as np, pandas as pd
D = '/data/dev/rnaseq-exploration/ice/popgen'; E = '/data/dev/rnaseq-exploration'
S = pd.read_csv(f'{D}/samples_map.tsv', sep='\t', index_col=0)
import sys as _s; _s.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); from common import excluded_names
S = S[~S.index.isin(excluded_names('metagenome'))]   # QC-excluded metagenomes (metadata qc_flag)
P = pd.read_csv(f'{D}/abund_pct.tsv', sep='\t', index_col=0); S = S.loc[[s for s in S.index if s in P.columns]]; P = P[S.index]
Cv = pd.read_csv(f'{D}/abund_cov.tsv', sep='\t', index_col=0)[S.index]
tax = pd.read_csv(f'{E}/data/gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome').classification
def T(b):
    t = tax.get(b, tax.get(b.replace('dastool-', ''), ''))
    r = {k: '' for k in 'pcofgs'}
    for x in str(t).split(';'):
        if len(x) > 3 and x[1:3] == '__': r[x[0]] = x[3:]
    return r
def label(b):
    t = T(b); best = t['s'] or t['g'] or t['f'] or t['o'] or t['c'] or t['p'] or 'unclassified'
    return f"{best} {b.replace('dastool-', '')}"
# paired log2 ratios (ice layer vs water of the same pool-day), pseudocount = 0.001 % of reads
S['pd'] = S.season + '|' + S.pool + '|' + S.day.astype(int).astype(str)
def paired(l1, l2):
    out = []
    for pd_, sub in S.groupby('pd'):
        a = sub[sub.layer.isin(l1)].index; b = sub[sub.layer.isin(l2)].index
        if len(a) and len(b): out.append((pd_, np.log2(P[a].mean(1) + 1e-3) - np.log2(P[b].mean(1) + 1e-3)))
    return pd.DataFrame({k: v for k, v in out})
ICE = ['TOP', 'BOT', 'MID', 'SI', 'PANCAKE']
IW = paired(ICE, ['SW']); TB = paired(['TOP'], ['BOT'])
res = pd.DataFrame({'label': [label(b) for b in P.index], 'mean_pct': P.mean(1), 'max_pct': P.max(1),
                    'ice_pct': P.loc[:, S.layer.isin(ICE).values].mean(1), 'water_pct': P.loc[:, (S.layer == 'SW').values].mean(1),
                    'log2_ice_water': IW.median(1), 'n_pairs_up': (IW > 0).sum(1), 'n_pairs': IW.shape[1],
                    'log2_top_bot': TB.median(1), 'n_tb_up': (TB > 0).sum(1), 'n_tb': TB.shape[1]}, index=P.index)
for k in 'cofg': res[k] = [T(b)[k] for b in P.index]
res.sort_values('mean_pct', ascending=False).to_csv(f'{D}/abundance_icewater_dna.tsv', sep='\t')
# siblings: congeneric MAG pairs (skani 80-95 % ANI)
sk = pd.read_csv(f'{D}/ref/skani_pairs.tsv', sep='\t')
sk = sk[(sk.ANI < 95) & (sk.ANI >= 80)].copy(); sk['A'] = 'dastool-' + sk.a; sk['B'] = 'dastool-' + sk.b
sk = sk[sk.A.isin(res.index) & sk.B.isin(res.index)]
sk = sk[(res.mean_pct.reindex(sk.A).values > .05) & (res.mean_pct.reindex(sk.B).values > .05)]
sk['dA'] = res.log2_ice_water.reindex(sk.A).values; sk['dB'] = res.log2_ice_water.reindex(sk.B).values
sk['contrast'] = (sk.dA - sk.dB).abs()
sk.sort_values('contrast', ascending=False).to_csv(f'{D}/sibling_pairs.tsv', sep='\t', index=False)
# popgen per bin
gm = pd.read_pickle(f'{E}/data/genemap.pkl').set_index('gene_id')[['product', 'KO', 'gene']].rename(columns={'gene': 'symbol'})
rows, pg, outl = [], [], []
for f in sorted(glob.glob(f'{D}/pop/*.json')):
    j = json.load(open(f)); b = 'dastool-' + j['bin']; r = {'bin': b, 'label': label(b)}
    for g, v in j['groups'].items(): r[f'pi|{g}'] = v['pi']; r[f'snvkb|{g}'] = v['snv_per_kb']; r[f'cov_bp|{g}'] = v['covered_bp']
    for c, v in j['contrasts'].items(): r[f'fst|{c}'] = v['fst']; r[f'fixed|{c}'] = v['fixed']; r[f'sites|{c}'] = v['sites']
    for c, v in j['pnps'].items(): r[f'pnps|{c}'] = v['pNpS']
    pr = pd.DataFrame(j['pairs'])
    if len(pr):
        for k, sub in pr.groupby(['season', 'kind']): r[f'pairmed|{k[0]}|{k[1]}'] = sub.fst.median()
    rows.append(r)
    pg.append({'bin': b, 'label': label(b), 'groups': j['groups'], 'contrasts': j['contrasts'], 'pnps': j['pnps'], 'pairs': j['pairs']})
    gf = f'{D}/pop/{j["bin"]}.genes.tsv'
    if os.path.exists(gf):
        G = pd.read_csv(gf, sep='\t'); G = G[G.sites >= 5]
        for se, sub in G.groupby('season'):
            if len(sub) < 50: continue
            gw = sub.num.sum() / sub.den.sum()
            mu, sd = sub.fst.mean(), sub.fst.std()
            o = sub[(sub.fst > max(gw * 3, mu + 3 * sd)) & (sub.fst > .1)].copy()
            o['bin'] = b; o['label'] = label(b); o['genome_fst'] = gw; o['z'] = (o.fst - mu) / sd
            outl.append(o)
R = pd.DataFrame(rows); R.to_csv(f'{D}/popgen_summary.tsv', sep='\t', index=False)
O = pd.concat(outl) if outl else pd.DataFrame(columns=['gene'])
if len(O):
    O = O.join(gm, on='gene'); O.sort_values(['bin', 'fst'], ascending=[True, False]).to_csv(f'{D}/fst_outlier_genes.tsv', sep='\t', index=False)
# site json
top = res[res.mean_pct > .2].sort_values('mean_pct', ascending=False).index.tolist()
keep = list(dict.fromkeys(top + list(sk.A.head(40)) + list(sk.B.head(40)) + list(R.bin)))
site = {'samples': [{'id': s, 'name': r['name'], 'season': r.season, 'pool': r.pool, 'day': int(r.day), 'layer': r.layer} for s, r in S.iterrows()],
        'abundance': {b: {'label': res.label[b], 'family': res.f[b], 'genus': res.g[b], 'pct': [round(float(v), 4) for v in P.loc[b]],
                          'log2_ice_water': None if pd.isna(res.log2_ice_water[b]) else round(float(res.log2_ice_water[b]), 2),
                          'n_pairs_up': int(res.n_pairs_up[b]), 'n_pairs': int(res.n_pairs[b]),
                          'log2_top_bot': None if pd.isna(res.log2_top_bot[b]) else round(float(res.log2_top_bot[b]), 2)} for b in keep},
        'siblings': sk.head(40)[['A', 'B', 'ANI', 'dA', 'dB']].round(2).to_dict('records'),
        'popgen': pg,
        'outlier_genes': [] if not len(O) else O[['bin', 'label', 'season', 'gene', 'sites', 'fst', 'genome_fst', 'product', 'KO']].round(3).fillna('').to_dict('records')}
json.dump(site, open(f'{D}/popgen.json', 'w'), separators=(',', ':'), default=lambda x: None if x != x else x)
print(len(keep), 'bins in site json;', len(pg), 'popgen bins;', len(O), 'outlier genes')
