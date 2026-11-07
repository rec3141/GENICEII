"""GC content per bin fasta: bin<TAB>gc<TAB>size."""
import sys, glob, os
for f in sorted(glob.glob(sys.argv[1] + '/*.fa')):
    gc = n = 0
    for line in open(f, 'rb'):
        if line[:1] == b'>': continue
        s = line.strip().upper(); gc += s.count(b'G') + s.count(b'C'); n += len(s) - s.count(b'N')
    print(os.path.basename(f)[:-3], round(gc / max(n, 1), 4), n, sep='\t')
