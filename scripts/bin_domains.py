"""Per-bin domain call from contig-level Tiara classifications (rule-based, rerunnable when bins change).

Input : /data/dev/fastsne/data/microscape/contigs_microscape.tsv.gz (microscape contig explorer table: id, length, bin, tiara, ...)
Rule  : length-weighted fraction of each bin's contigs per Tiara label;
        'eukaryote'  if eukarya   >= 0.5
        'organelle'  if organelle >= 0.5  (Tiara 'organelle', 'plastid' or 'mitochondrion')
        'prokaryote' otherwise.
Output: data/bin_domain.tsv (bin, bp, n_contigs, frac_*, call, rule)"""
import pandas as pd, numpy as np
SRC = '/data/dev/fastsne/data/microscape/contigs_microscape.tsv.gz'
OUT = '/data/dev/rnaseq-exploration/data/bin_domain.tsv'
EUK, ORG = 0.5, 0.5
T = pd.read_csv(SRC, sep='\t', usecols=['id', 'length', 'bin', 'tiara'], low_memory=False)
T = T[T.bin.notna() & (T.bin != '')].copy()
t = T.tiara.fillna('NA')
T['lab'] = np.select([t.eq('eukarya'), t.isin(['organelle', 'plastid', 'mitochondrion']), t.eq('bacteria'), t.eq('archaea'), t.eq('prokarya')],
                     ['eukarya', 'organelle', 'bacteria', 'archaea', 'prokarya'], 'unknown')
bp = T.groupby('bin').length.sum()
F = T.pivot_table(index='bin', columns='lab', values='length', aggfunc='sum', fill_value=0).div(bp, axis=0)
for c in ('eukarya', 'organelle', 'bacteria', 'archaea', 'prokarya', 'unknown'):
    if c not in F: F[c] = 0.0
D = pd.DataFrame({'bp': bp, 'n_contigs': T.groupby('bin').size()})
for c in ('eukarya', 'organelle', 'bacteria', 'archaea', 'prokarya', 'unknown'): D['frac_' + c] = F[c].round(4)
D['call'] = np.select([D.frac_eukarya >= EUK, D.frac_organelle >= ORG], ['eukaryote', 'organelle'], 'prokaryote')
D['rule'] = f'Tiara length-weighted: eukarya>={EUK} eukaryote; organelle>={ORG} organelle; else prokaryote'
D.sort_index().to_csv(OUT, sep='\t', index_label='bin')
print(D.call.value_counts().to_dict())
print('eukaryote:', list(D.index[D.call == 'eukaryote'])); print('organelle:', len(D.index[D.call == 'organelle']), list(D.index[D.call == 'organelle']))
