"""Per-protein amino-acid counts for a (large) faa -> ids.txt + counts.npy (uint16, n x 20, order ACDEFGHIKLMNPQRSTVWY)."""
import sys, numpy as np
AA = b'ACDEFGHIKLMNPQRSTVWY'
LUT = np.full(256, 20, np.int64)
for i, a in enumerate(AA): LUT[a] = i
faa, out = sys.argv[1], sys.argv[2]
ids, rows, seq = [], [], []
def flush():
    s = np.frombuffer(b''.join(seq), np.uint8)
    rows.append(np.bincount(LUT[s], minlength=21)[:20].astype(np.uint16))
with open(faa, 'rb') as f:
    for line in f:
        if line[:1] == b'>':
            if ids: flush()
            ids.append(line[1:].split()[0].decode()); seq = []
        else:
            seq.append(line.strip())
    flush()
np.save(out + '/counts.npy', np.vstack(rows))
open(out + '/ids.txt', 'w').write('\n'.join(ids) + '\n')
print(len(ids))
