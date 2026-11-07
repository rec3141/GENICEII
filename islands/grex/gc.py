# GC of each island and of its contig outside the island (assembly.fasta streamed once; only island contigs kept)
import sys
isl = {}
for l in open("genomic_islands.tsv"):
    f = l.split()
    if f[0] == "island_id": continue
    isl.setdefault(f[1], []).append((f[0], int(f[2]), int(f[3])))
out = open("island_gc.tsv", "w"); out.write("island_id\tcontig\tcontig_len\tisland_gc\tcontig_gc_outside\n")
name, seq = None, []
def flush():
    if name in isl:
        s = "".join(seq).upper(); tot = s.count("G") + s.count("C")
        for i, a, b in isl[name]:
            w = s[a - 1:b]; g = w.count("G") + w.count("C")
            out.write(f"{i}\t{name}\t{len(s)}\t{g / max(len(w), 1):.4f}\t{(tot - g) / max(len(s) - len(w), 1):.4f}\n")
for l in open(sys.argv[1]):
    if l.startswith(">"):
        flush(); name = l[1:].split()[0]; seq = []
    elif name in isl: seq.append(l.strip())
flush()
