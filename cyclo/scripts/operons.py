"""Group Cycloclasticus HC genes into neighbourhoods and summarise their expression."""
import sys, json
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
C = '/data/dev/rnaseq-exploration/cyclo/'
G = pd.read_csv(C + 'cyclo_genes_annot.tsv', sep='\t', keep_default_na=False)
cpm = np.load(C + 'cyclo_cpm.npy'); wt = np.load(C + 'cyclo_wtpm.npy'); pr = np.load(C + 'cyclo_pairs.npy')
S = samples(); m = meta().loc[S]
core = m.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI']).values
def idx(**k):
    s = core.copy()
    for a, v in k.items(): s &= (m[a] == v).values
    return s
G = G.sort_values(['contig', 'start']).reset_index()      # 'index' = row into expression arrays
G['ord'] = G.groupby('contig').cumcount()
hc = G[(G.hc != '') & (G.hc != 'chemotaxis') & (G.cbin != 'unbinned_kaiju_cyclo')]
clusters = []
for c, d in hc.groupby('contig'):
    o = d.ord.values; grp = [[o[0]]]
    for x in o[1:]:
        (grp[-1].append(x) if x - grp[-1][-1] <= 4 else grp.append([x]))
    for gr in grp:
        if len(gr) >= 3: clusters.append((c, min(gr) - 2, max(gr) + 2))
def summ(r):
    out = {}
    for se in ['2024-25', '2025-26']:
        o = idx(season=se, group='EXP') & (m.phase == 'oiled').values
        ct = idx(season=se, group='CTL') & (m.phase == 'control').values
        p = idx(season=se, group='EXP') & (m.phase == 'pre').values
        out[se] = {'cpm_oiled': round(float(cpm[r, o].mean()), 2), 'cpm_ctl': round(float(cpm[r, ct].mean()), 2),
                   'wtpm_oiled': round(float(wt[r, o].mean()), 1), 'wtpm_pre': round(float(wt[r, p].mean()), 1),
                   'wtpm_ctl': round(float(wt[r, ct].mean()), 1), 'pairs_oiled': int(pr[r, o].sum())}
    return out
ops = []
for i, (c, a, b) in enumerate(clusters):
    d = G[(G.contig == c) & (G.ord >= a) & (G.ord <= b)]
    genes_ = [{'id': x.gene_id, 'start': int(x.start), 'end': int(x.end), 'strand': x.strand, 'product': x['product'],
               'gene': x.gene, 'ko': x.KO, 'cat': x.hc, 'length': int(x.length), 'expr': summ(x['index'])} for _, x in d.iterrows()]
    ops.append({'id': f'op{i+1:02d}', 'bin': d.cbin.iloc[0], 'contig': c, 'start': int(d.start.min()), 'end': int(d.end.max()),
                'n_hc': int((d.hc != '').sum()), 'cats': sorted(set(d.hc) - {''}), 'genes': genes_,
                'pairs_oiled_total': int(sum(gg['expr']['2024-25']['pairs_oiled'] + gg['expr']['2025-26']['pairs_oiled'] for gg in genes_))})
json.dump(ops, open(C + 'operons.json', 'w'))
for o in sorted(ops, key=lambda o: -o['pairs_oiled_total']):
    print(o['id'], o['bin'].replace('dastool-', ''), o['contig'], o['start'], o['end'], 'nHC', o['n_hc'], 'pairs', o['pairs_oiled_total'])
    for gg in o['genes']:
        e = gg['expr']
        print(f"   {gg['strand']} {gg['length']:5d} {gg['product'][:58]:58s} {gg['cat'][:22]:22s} 24/25 wTPM pre {e['2024-25']['wtpm_pre']:7.0f} oil {e['2024-25']['wtpm_oiled']:7.0f} | 25/26 pre {e['2025-26']['wtpm_pre']:7.0f} oil {e['2025-26']['wtpm_oiled']:7.0f} ctl {e['2025-26']['wtpm_ctl']:6.0f}")
