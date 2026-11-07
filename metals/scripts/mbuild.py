"""site/metals.json for the Metals & toxicity page, from metals/out."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/metals/scripts')
from common import *
from msets import SETS
R = '/data/dev/rnaseq-exploration/'; O = R + 'metals/out/'
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
def lab(b):
    if not isinstance(b, str) or not b: return 'unbinned contigs'
    if special_label(b): return special_label(b)
    t = [x.split('__', 1)[1] for x in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in x and x.split('__', 1)[1]]
    return ((t[-1] if t else '') or kaiju_bin_label(b) or 'unclassified') + ' ' + b.replace('dastool-', '')
M = pd.read_csv(O + 'members.tsv', sep='\t', keep_default_na=False, low_memory=False)
E = pd.read_csv(O + 'cat_cpm.tsv', sep='\t', index_col=0)
DD = pd.read_csv(O + 'diesel.tsv', sep='\t')
IW = pd.read_csv(O + 'icewater.tsv', sep='\t', index_col=0)
pairs = {k: pd.read_csv(O + f'{k}_pairs.tsv', sep='\t', index_col=0) for k in ('icewater', 'topbot')}
Hs = pd.read_csv(O + 'hosts.tsv', sep='\t', index_col=[0, 1], keep_default_na=False)
MS = pd.read_csv(O + 'mobility_summary.tsv', sep='\t', index_col=[0, 1])
DG = pd.read_csv(O + 'degraders.tsv', sep='\t')
cats = []
for cat, grp, _, _ in SETS:
    c = {'name': cat, 'group': grp, 'n_all': int((M.category == cat).sum()), 'n_conf': int(((M.category == cat) & (M.confident.astype(str) == 'True')).sum())}
    c['diesel'] = [{'season': x.season, 'medium': x.medium, 'oil': r(x.oil, 1), 'ctl': r(x.ctl, 1), 'lfc': r(x.log2_oil_ctl), 'baci': r(x.get('log2_baci')),
                    'pre_exp': r(x.get('pre_exp'), 1), 'pre_ctl': r(x.get('pre_ctl'), 1), 'n_oil': int(x.n_oil), 'n_ctl': int(x.n_ctl), 'n_pre_exp': int(x.n_pre_exp)}
                   for _, x in DD[DD.category == cat].iterrows()]
    x = IW.loc[cat]
    c['ice'] = {'iw': r(x.icewater_delta), 'iw_up': int(x.icewater_n_up), 'iw_n': int(x.icewater_n), 'iw_q': r(x.icewater_q, 4),
                'tb': r(x.topbot_delta), 'tb_up': int(x.topbot_n_up), 'tb_n': int(x.topbot_n), 'tb_q': r(x.topbot_q, 4),
                'cpm': {k: r(x[f'cpm_{k}'], 1) for k in ('SW', 'TOP', 'BOT')},
                'iw_pairs': [r(v) for v in pairs['icewater'].loc[cat]], 'tb_pairs': [r(v) for v in pairs['topbot'].loc[cat]]}
    if cat in Hs.index.get_level_values(0):
        h = Hs.loc[cat]
        hosts = {}
        for col in h.columns:
            tot = h[col].sum(); top = h[col].sort_values(ascending=False).head(6)
            hosts[col] = {'total': r(tot, 1), 'top': [{'host': lab(b), 'cpm': r(v, 1), 'share': r(v / tot if tot else 0, 3)} for b, v in top.items() if v > 0]}
        c['hosts'] = hosts
    if (grp, cat) in MS.index:
        mm = MS.loc[(grp, cat)]
        c['mob'] = {k: r(mm[k], 2) for k in ('plasmid_strict', 'genomic_island', 'integron', 'IS_within_10kb', 'mge_borne', 'mobile_assoc', 'hc_within_10kb')}
    cats.append(c)
ref = MS.loc[('reference', 'Hydrocarbon genes (confident, hgt/)')]
J = {'samples': list(E.columns), 'cats': cats, 'series': {c: [r(v, 2) for v in E.loc[c]] for c in E.index},
     'hc_ref': {k: r(ref[k], 2) for k in ('plasmid_strict', 'genomic_island', 'integron', 'IS_within_10kb', 'mge_borne', 'mobile_assoc')},
     'pair_cols': {k: v.columns.tolist() for k, v in pairs.items()}}
# within-host expression of the oil degraders and other diesel-pool responders
deg = []
for (b, labl, cat), d in DG.groupby(['bin', 'label', 'category']):
    get = lambda s, p, ph: (lambda x: r(x['median'].iloc[0], 0) if len(x) else None)(d[(d.season == s) & (d.pool == p) & (d.phase == ph)])
    deg.append({'bin': b.replace('dastool-', ''), 'label': labl, 'category': cat, 'genes': int(d.genes.iloc[0]),
                'v': {f'{s}|{p}|{ph}': get(s, p, ph) for s in ('2024-25', '2025-26') for p, ph in (('EXP', 'pre'), ('EXP', 'oiled'), ('CTL', 'pre'), ('CTL', 'control'))}})
J['degraders'] = deg
# plasmids carrying both toxic-metal and hydrocarbon genes (co-selection candidates)
Mo = pd.read_csv(O + 'mobility.tsv', sep='\t', keep_default_na=False, low_memory=False)
Mo = Mo[(Mo.confident.astype(str) == 'True') & (Mo.group == 'toxic metal')]
hc = pd.read_csv(R + 'hgt/hc_context.tsv', sep='\t', keep_default_na=False); hc = hc[hc.tier == 'confident']
p = pd.read_csv(R + 'hgt/mge/plasmid_summary.tsv', sep='\t').set_index('seq_name')
strict = set(p[(p.plasmid_score >= 0.9) & ((p.n_hallmarks >= 1) | p.conjugation_genes.notna())].index)
both = sorted(set(Mo.contig) & set(hc.contig) & strict, key=lambda c: -p.loc[c, 'length'])
J['coplasmids'] = [{'contig': c, 'len': int(p.loc[c, 'length']), 'score': r(p.loc[c, 'plasmid_score']), 'conj': isinstance(p.loc[c, 'conjugation_genes'], str),
                    'host': lab(Mo[Mo.contig == c].bin.iloc[0]), 'metal': sorted(set(Mo[Mo.contig == c].category)), 'hc': sorted(set(hc[hc.contig == c].category))} for c in both]
J['plasmid_counts'] = {'metal': len(set(Mo.contig) & strict), 'hc': len(set(hc.contig) & strict), 'both': len(both)}
json.dump(J, open(R + 'site/metals.json', 'w'), separators=(',', ':'), allow_nan=False)
print(len(cats), 'categories;', len(deg), 'degrader rows;', len(both), 'co-plasmids;', os.path.getsize(R + 'site/metals.json') // 1000, 'kB')
