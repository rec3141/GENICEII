"""Extract CuMMO operon proteins (C, A, B) for comebin_2576 and Cycloclasticus metabat_9342 by translating ORFs of the gene lengths in genemap.
The local genemap has no coordinates, so ORFs are matched by exact length (incl. stop) near the anchor gene. Out: porticocc/out/cummo_query.faa"""
import pandas as pd, re
B = '/data/dev/rnaseq-exploration/'
G = pd.read_pickle(B + 'data/genemap.pkl')
CODON = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
         'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
def rc(s): return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
def tr(s): return ''.join(CODON.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))
def fasta(p):
    d, k = {}, None
    for l in open(p):
        if l.startswith('>'): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l.strip())
    return {k: ''.join(v).upper() for k, v in d.items()}
def orfs(seq):
    """all ORFs: (strand, start, end, nt_len_incl_stop, protein) for every start codon ATG/GTG/TTG -> next in-frame stop"""
    out = []
    for strand, s in (('+', seq), ('-', rc(seq))):
        for f in range(3):
            prot = tr(s[f:])
            for m in re.finditer(r'[^*]*\*', prot):
                seg = m.group(); st = m.start()
                for j, aa in enumerate(seg):
                    codon = s[f + 3 * (st + j): f + 3 * (st + j) + 3]
                    if codon in ('ATG', 'GTG', 'TTG'):
                        out.append((strand, f + 3 * (st + j), 3 * (len(seg) - j), 'M' + seg[j + 1:-1]))
    return out
Q = {}
for bin_, contig, ids, lab in (('comebin_2576', 'contig_0004457', ['JPFHOD_001079467', 'JPFHOD_001079468', 'JPFHOD_001079469'], 'Porticoccaceae_50-400-T64'),
                              ('metabat_9342', 'contig_0001426', ['JPFHOD_000613387', 'JPFHOD_000613386', 'JPFHOD_000613385'], 'Cycloclasticus_metabat_9342')):
    fa = B + (f'mgscreen/ref/dastool-{bin_}.fa' if bin_ == 'comebin_2576' else f'cyclo/bins/dastool-{bin_}.fa')
    seq = fasta(fa)[contig]; O = orfs(seq)
    for gid, sub in zip(ids, 'CAB'):
        L = int(G.length[G.gene_id == gid].iloc[0])
        hits = [o for o in O if o[2] == L]
        Q[f'{lab}_{sub}|{gid}'] = hits
        print(gid, sub, L, len(hits), [(h[0], h[1]) for h in hits][:5])
# disambiguate: keep, per gene, the hit closest to the other operon members (same strand, within 5 kb)
out = open(B + 'porticocc/out/cummo_query.faa', 'w')
for lab in ('Porticoccaceae_50-400-T64', 'Cycloclasticus_metabat_9342'):
    keys = [k for k in Q if k.startswith(lab)]
    best = None
    import itertools
    for combo in itertools.product(*[Q[k] for k in keys]):
        if len({c[0] for c in combo}) > 1: continue
        span = max(c[1] for c in combo) - min(c[1] for c in combo)
        if best is None or span < best[0]: best = (span, combo)
    print(lab, 'operon span', best[0])
    for k, c in zip(keys, best[1]): out.write(f'>{k}\n{c[3]}\n')
out.close()
