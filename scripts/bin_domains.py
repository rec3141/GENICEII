"""Per-bin domain call from contig-level Tiara classifications (rule-based, rerunnable when bins change).

Input : /data/dev/fastsne/data/microscape/contigs_microscape.tsv.gz (microscape contig explorer table: id, length, bin, tiara, ...)
Rule  : length-weighted fraction of each bin's contigs per Tiara label;
        'eukaryote'  if eukarya   >= 0.5
        'organelle'  if organelle >= 0.5  (Tiara 'organelle', 'plastid' or 'mitochondrion')
        'prokaryote' otherwise.
Output: data/bin_domain.tsv (bin, bp, n_contigs, frac_*, call, rule, euk_group)

Eukaryotic genome groups (second layer). DAS Tool keeps almost no eukaryotic bins, so eukaryotic nuclear genomes are
consolidated separately from the raw binner outputs by euks/mags/scripts, run in this order:
  survey.py (candidate bins = config LIST at its top) -> compact.py (t-SNE density-gap split) -> final.py (footprint merge,
  stable E## ids against out/v1) -> expr.py (contig-read expression) -> finalize.py (wiring tables + figure).
This script copies their wiring tables into data/:
  data/euk_groups.tsv          one row per E## group (status network / little CMO RNA / doubtful, network_ok, purity, taxonomy hint)
  data/contig_euk_group.tsv.gz contig -> euk_group (+ Tiara call); one group per contig (contigs shared by two overlapping
                               groups go to the larger group, as in euks/mags finalize)
and adds bin_domain.tsv:euk_group = the E## group holding >= 50% of a eukaryote bin's bp (e.g. comebin_3937 -> E05).
scripts/common.py (euk_group, euk_group_info, euk_group_label, contig_host_label, special_label) reads these."""
import os
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
# eukaryotic genome groups (euks/mags) -> data/, and the group of each eukaryote bin
EM = '/data/dev/rnaseq-exploration/euks/mags/out/'; DD = '/data/dev/rnaseq-exploration/data/'
EG = pd.read_csv(EM + 'euk_groups.tsv', sep='\t')
EG[['group', 'status', 'network_ok', 'bp', 'contigs', 'gc', 'tiara_euk', 'who_euk', 'tsne_rms', 'taxonomy_hint', 'rrna18S', 'expr_ppm_median',
    'baci_2025', 'ice_vs_water_log2', 'note', 'parts']].to_csv(DD + 'euk_groups.tsv', sep='\t', index=False)
CE = pd.read_csv(EM + 'contig_euk_group.tsv', sep='\t')
# euks/mags/out/contig_euk_group.tsv is deduplicated upstream (a contig shared by two overlapping groups goes to the larger
# group); guard: if duplicates ever reappear, apply the same rule here
print('duplicate contigs in contig_euk_group.tsv:', int(CE.contig.duplicated().sum()))
CE = CE.assign(_bp=CE.euk_group.map(EG.set_index('group').bp)).sort_values(['contig', '_bp', 'euk_group'], ascending=[True, False, True]).drop_duplicates('contig')
CE[['contig', 'euk_group', 'tiara']].to_csv(DD + 'contig_euk_group.tsv.gz', sep='\t', index=False, compression='gzip')
eg = T.merge(CE[['contig', 'euk_group']], left_on='id', right_on='contig', how='inner')
share = eg.groupby(['bin', 'euk_group']).length.sum().div(bp, level='bin').reset_index(name='share')
share = share.sort_values('share', ascending=False).drop_duplicates('bin')
D['euk_group'] = share[share.share >= 0.5].set_index('bin').euk_group.reindex(D.index)
D.loc[D.call != 'eukaryote', 'euk_group'] = np.nan
D.sort_index().to_csv(OUT, sep='\t', index_label='bin')
print('eukaryote bins -> groups:', D.euk_group.dropna().to_dict(), '| groups:', len(EG), '| contigs in groups:', len(CE))
print(D.call.value_counts().to_dict())
print('eukaryote:', list(D.index[D.call == 'eukaryote'])); print('organelle:', len(D.index[D.call == 'organelle']), list(D.index[D.call == 'organelle']))
