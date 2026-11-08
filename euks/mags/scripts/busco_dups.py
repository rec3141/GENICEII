"""Identity between duplicated BUSCO copies per MAG (2026-11-08): for each Duplicated BUSCO in
out/busco/run_euk_miniprot/<MAG>.fa/run_*/full_table.tsv, blastn (-task blastn) copy 1 against the other copies; best HSP identity.
Haplotypes of one diploid genome match at >= 98%; strains ~90-97%; different species < 90%.
Out: out/busco/dup_identity.tsv (mag, n_dup, median_nt_id, n_ge98). Needs samtools (BUSCO env) and blastn."""
import glob, re, subprocess
import numpy as np, pandas as pd
O = '/data/dev/rnaseq-exploration/euks/mags/out/busco/'
ST = '/data/scratch/grid-envs/conda/busco/bin/samtools'
rows = []
for f in sorted(glob.glob(O + 'run_euk_miniprot/*.fa/run_*/full_table.tsv')):
    mag = re.search(r'(EE\d+)\.fa', f).group(1); fa = f'{O}mags/{mag}.fa'
    D = pd.DataFrame([r for r in (l.rstrip('\n').split('\t')[:6] for l in open(f) if not l.startswith('#'))
                      if len(r) >= 6 and r[1] == 'Duplicated'], columns=['busco', 'status', 'seq', 'start', 'end', 'strand'])
    ids = []
    for b, g in D.groupby('busco'):
        reg = [f'{x.seq}:{min(int(x.start), int(x.end))}-{max(int(x.start), int(x.end))}' for x in g.itertuples()]
        s = [subprocess.run([ST, 'faidx', fa, r], capture_output=True, text=True).stdout for r in reg]
        open(O + 'dupchk/a.fa', 'w').write(s[0]); open(O + 'dupchk/b.fa', 'w').write(''.join(s[1:]))
        o = subprocess.run(['blastn', '-task', 'blastn', '-query', O + 'dupchk/a.fa', '-subject', O + 'dupchk/b.fa',
                            '-outfmt', '6 pident', '-max_hsps', '1'], capture_output=True, text=True).stdout.split()
        ids.append(max(map(float, o)) if o else np.nan)
    rows.append(dict(mag=mag, n_dup=len(ids), median_nt_id=round(float(np.nanmedian(ids)), 1) if ids else np.nan,
                     n_ge98=sum(x >= 98 for x in ids if x == x)))
pd.DataFrame(rows).set_index('mag').to_csv(O + 'dup_identity.tsv', sep='\t')
print(pd.DataFrame(rows).to_string(index=False))
