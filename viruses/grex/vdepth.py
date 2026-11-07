# Stream mag_marine/mapping/depths.txt once: per-sample total depth*length over all contigs, and the depth rows of the viral contigs.
import sys, gzip
ids = set(l.strip() for l in open('viral_ids.txt'))
out = gzip.open('viral_depths.tsv.gz', 'wt'); tot = None
with open(sys.argv[1]) as f:
    hdr = f.readline().rstrip('\n').split('\t')
    cols = [i for i, h in enumerate(hdr) if h.endswith('.bam')]
    out.write('\t'.join(['contig', 'len'] + [hdr[i] for i in cols]) + '\n')
    tot = [0.0] * len(cols); totlen = 0
    for line in f:
        p = line.rstrip('\n').split('\t'); L = float(p[1]); totlen += L
        d = [float(p[i]) for i in cols]
        for k, x in enumerate(d): tot[k] += x * L
        if p[0] in ids: out.write('\t'.join([p[0], p[1]] + [p[i] for i in cols]) + '\n')
out.close()
with open('depth_totals.tsv', 'w') as o:
    o.write('sample\ttotal_bp_depth\n')
    for i, k in zip(cols, tot): o.write(f'{hdr[i]}\t{k}\n')
print('assembly length', totlen)
