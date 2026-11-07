"""K-mer screen of Illumina metagenomes against key MAGs/plasmids (no mapper available locally).
Target-unique canonical 25-mers; for each sample, count k-mer hits in the first N R1 reads.
Output: hits per million sampled k-mers per Mbp of target unique k-mers (a coverage-like abundance)."""
import sys, os, glob, gzip, subprocess, numpy as np
K = 25
LUT = np.full(256, 4, np.uint8)
for i, ch in enumerate(b'ACGT'): LUT[ch] = i; LUT[ch + 32] = i
def kmers(seq_bytes):
    """canonical 2-bit 25-mers of a byte string (N-containing windows dropped)."""
    c = LUT[np.frombuffer(seq_bytes, np.uint8)]
    n = len(c) - K + 1
    if n <= 0: return np.empty(0, np.uint64)
    bad = np.concatenate([[0], np.cumsum(c == 4)])
    ok = (bad[K:] - bad[:-K]) == 0
    c64 = np.where(c == 4, 0, c).astype(np.uint64)
    f = np.zeros(n, np.uint64); r = np.zeros(n, np.uint64)
    for j in range(K):
        f = (f << np.uint64(2)) | c64[j:j + n]
        r = r | ((np.uint64(3) - c64[j:j + n]) << np.uint64(2 * j))
    return np.minimum(f, r)[ok]
def fasta_seqs(path):
    seq = []
    for line in open(path, 'rb'):
        if line.startswith(b'>'):
            if seq: yield b''.join(seq); seq = []
        else: seq.append(line.strip())
    if seq: yield b''.join(seq)

def build(refdir, out):
    targets = {}
    for f in sorted(glob.glob(refdir + '/dastool-*.fa')):
        targets[os.path.basename(f)[8:-3]] = [f]
    for name, s in zip(['plasmid_contig_0009104', 'plasmid_contig_0009679'], [None, None]): pass
    ks, lab, names = [], [], []
    for t, (name, files) in enumerate(targets.items()):
        u = np.unique(np.concatenate([kmers(s + b'N' if False else s) for f in files for s in fasta_seqs(f)]))
        ks.append(u); lab.append(np.full(len(u), t, np.int16)); names.append(name)
    # plasmid contigs: each its own target, removed from the MAG sets
    for s, name in zip(fasta_seqs(refdir + '/plasmids.fa'), [l[1:].strip().decode() for l in open(refdir + '/plasmids.fa', 'rb') if l.startswith(b'>')]):
        u = np.unique(kmers(s)); ks.append(u); lab.append(np.full(len(u), len(names), np.int16)); names.append('plasmid_' + name)
    allk = np.concatenate(ks); alll = np.concatenate(lab)
    o = np.argsort(allk, kind='stable'); allk, alll = allk[o], alll[o]
    dup = np.zeros(len(allk), bool); same = allk[1:] == allk[:-1]; dup[1:] |= same; dup[:-1] |= same
    allk, alll = allk[~dup], alll[~dup]
    np.savez(out, k=allk, l=alll, names=np.array(names))
    print({n: int((alll == i).sum()) for i, n in enumerate(names)})

def screen(db, fq, nreads):
    d = np.load(db); K_, L_, names = d['k'], d['l'], list(d['names'])
    p = subprocess.Popen(['bash', '-c', f'zcat {fq} 2>/dev/null | head -n {4 * nreads}'], stdout=subprocess.PIPE)
    lines = p.stdout.read().split(b'\n')[1::4]
    hits = np.zeros(len(names), np.int64); tot = 0
    for b in range(0, len(lines), 50000):
        km = kmers(b'N'.join(lines[b:b + 50000]))
        tot += len(km)
        i = np.searchsorted(K_, km); i[i >= len(K_)] = 0
        m = K_[i] == km
        hits += np.bincount(L_[i[m]], minlength=len(names))
    size = np.bincount(L_, minlength=len(names))
    return names, hits, tot, size

if __name__ == '__main__':
    if sys.argv[1] == 'build': build(sys.argv[2], sys.argv[3])
    else:
        db, sample, fq, n, out = sys.argv[2:7]
        names, hits, tot, size = screen(db, fq, int(n))
        with open(out, 'w') as fh:
            for nm, h, s in zip(names, hits, size):
                fh.write(f'{sample}\t{nm}\t{h}\t{tot}\t{s}\t{h / tot * 1e6 / (s / 1e6):.4f}\n')
