"""pileup_bin.py BIN : bcftools mpileup of all BAMs over one MAG -> pile/BIN.npz
Stores, for candidate polymorphic sites (pooled minor allele >= 4 reads and >= 2%), per-sample ACGT counts,
plus per-sample and per-group counts of sites covered at >= 5/10/20 reads (denominators for pi / SNV density)."""
import sys, subprocess, numpy as np, pandas as pd, os
D = '/data/dev/rnaseq-exploration/ice/popgen'; BCF = os.path.expanduser('~/apps/bcftools-1.9/bcftools')
b = sys.argv[1]
S = pd.read_csv(f'{D}/samples_map.tsv', sep='\t', index_col=0)
S = S[[os.path.exists(f'{D}/bam/{s}.bam.bai') for s in S.index]]
cb = pd.read_csv(f'{D}/ref/contig_bin_len.tsv', sep='\t')
cb = cb[cb.bin == 'dastool-' + b]
os.makedirs(f'{D}/pile', exist_ok=True)
bed = f'{D}/pile/{b}.bed'
cb.assign(s=0)[['contig', 's', 'len']].to_csv(bed, sep='\t', header=False, index=False)
bl = f'{D}/pile/{b}.bams'; open(bl, 'w').write('\n'.join(f'{D}/bam/{s}.bam' for s in S.index) + '\n')
cmd = f"{BCF} mpileup -f {D}/ref/allbins.fa -R {bed} -b {bl} -a AD -q 20 -Q 20 -B -I -d 100000 -Ou 2>/dev/null | {BCF} query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n'"
p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
n = len(S); IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
groups = {}
for i, (s, r) in enumerate(S.iterrows()):
    lay = 'SW' if r.layer == 'SW' else 'ICE'
    groups.setdefault(f'{r.season}|{lay}', []).append(i)
    if r.layer in ('TOP', 'BOT'): groups.setdefault(f'{r.season}|{r.layer}', []).append(i)
gnames = list(groups); gidx = [np.array(groups[g]) for g in gnames]
TH = np.array([5, 10, 20])
cov_s = np.zeros((n, 3), np.int64); cov_g = np.zeros((len(gnames), 3), np.int64)
contig, pos, ref, cnt = [], [], [], []
for line in p.stdout:
    f = line.rstrip('\n').split('\t')
    if f[3] == '<*>':   # fast path: no alternate allele seen
        dep = np.array([0 if x == '.' else int(x.split(',', 1)[0]) for x in f[4:]])
        cov_s += dep[:, None] >= TH
        gd = np.array([dep[g].sum() for g in gidx]); cov_g += gd[:, None] >= TH
        continue
    al = [f[2]] + f[3].split(',')
    a = np.zeros((n, 4), np.int32)
    for j in range(n):
        ad = f[4 + j]
        if ad == '.': continue
        for k, v in enumerate(ad.split(',')):
            if v != '0' and al[k] in IDX: a[j, IDX[al[k]]] += int(v)
    dep = a.sum(1)
    cov_s += dep[:, None] >= TH
    gd = np.array([dep[g].sum() for g in gidx]); cov_g += gd[:, None] >= TH
    if True:
        tot = a.sum(0); srt = np.sort(tot)
        if srt[-2] >= 4 and srt[-2] >= 0.02 * tot.sum():
            contig.append(f[0]); pos.append(int(f[1])); ref.append(IDX.get(f[2].upper(), 4)); cnt.append(a.astype(np.uint16))
np.savez_compressed(f'{D}/pile/{b}.npz', samples=np.array(S.index), contig=np.array(contig), pos=np.array(pos), ref=np.array(ref, np.int8),
                    counts=np.array(cnt, np.uint16).reshape(-1, n, 4), cov_s=cov_s, cov_g=cov_g, gnames=np.array(gnames), th=TH)
print(b, len(pos), 'sites', cov_g.tolist())
