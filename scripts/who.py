"""Which contigs/MAGs carry the expression of a hydrocarbon gene category, per season×pool."""
import sys; sys.path.insert(0, '.')
from common import *
X, g = genes(); S = samples(); m = meta().loc[S]
mem = pd.read_pickle('../out/hc_members.pkl')
kj = pd.read_csv(D + 'kaiju_contigs.tsv', sep='\t', usecols=['contig_id', 'lineage']).set_index('contig_id')['lineage']
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
def short(l, n=3):
    return ';'.join([x.strip() for x in str(l).split(';') if x.strip() not in ('NA', '')][-n:])
grp = {k: m.index[(m.season + ' ' + m.design + ' ' + m.phase) == k] for k in
       ['2024-25 EXP_SW oiled', '2024-25 CTL_SW control', '2024-25 EXP_SI oiled', '2024-25 CTL_SI control',
        '2025-26 EXP_SW oiled', '2025-26 CTL_SW control', '2025-26 EXP_SI oiled', '2025-26 CTL_SI control']}
col = {s: i for i, s in enumerate(S)}
for cat in sys.argv[1:]:
    idx = mem[cat]; M = X[idx].toarray()
    print(f'\n=== {cat}')
    for k, ss in grp.items():
        v = M[:, [col[s] for s in ss]].sum(1); t = v.sum()
        if not t: continue
        o = np.argsort(-v)[:4]
        desc = '; '.join(f"{100*v[i]/t:.0f}% {g.contig.values[idx[i]]}[{g.bin.values[idx[i]].replace('dastool-','') or '-'}|{short(gt.get(g.bin.values[idx[i]].replace('dastool-',''), ''),2) if g.bin.values[idx[i]] else short(kj.get(g.contig.values[idx[i]], ''),2)}]" for i in o)
        print(f'  {k:24s} n={len(ss):2d} pairs={t:7d}  {desc}')
