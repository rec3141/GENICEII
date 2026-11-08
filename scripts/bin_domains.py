"""Per-bin domain call from contig-level Tiara classifications (rule-based, rerunnable when bins change).

Input : /data/dev/fastsne/data/microscape/contigs_microscape.tsv.gz (microscape contig explorer table: id, length, bin, tiara, ...)
Rule  : length-weighted fraction of each bin's contigs per Tiara label;
        'eukaryote'  if eukarya   >= 0.5
        'organelle'  if organelle >= 0.5  (Tiara 'organelle', 'plastid' or 'mitochondrion')
        'prokaryote' otherwise.
Output: data/bin_domain.tsv (bin, bp, n_contigs, frac_*, call, rule, euk_group, hull_mag)

Eukaryotic genome groups (second layer). DAS Tool keeps almost no eukaryotic bins, so eukaryotic nuclear genomes are
consolidated separately from the raw binner outputs by euks/mags/scripts, run in this order:
  survey.py (candidate bins = config LIST at its top) -> compact.py (t-SNE density-gap split) -> final.py (footprint merge,
  split rule, stable ids) -> expr.py (contig-read expression) -> finalize.py (wiring tables) -> kernel.py (expression/coverage
  kernel blobs) -> da.py (per-blob RNA/DNA differential abundance) -> register_blobs.py (blobs -> new E##/V##) -> plot.py.
  Ids come from one registry (euks/mags/out/id_registry.tsv); groups >= 50% geNomad-viral get V## ids and are kept out of
  the eukaryotic tables.
This script copies their wiring tables into data/:
  data/euk_groups.tsv          one row per E## group (status network / little CMO RNA / doubtful, network_ok, purity, taxonomy hint)
  data/contig_euk_group.tsv.gz contig -> euk_group (+ Tiara call); one group per contig (contigs shared by two overlapping
                               groups go to the larger group, as in euks/mags finalize)
  data/virus_groups.tsv        one row per V## virus group (from the same pipeline)
  data/contig_virus_group.tsv.gz contig -> V##
and adds bin_domain.tsv:euk_group = the E## group holding >= 50% of a eukaryote bin's bp (e.g. comebin_3937 -> E05).
Hull-curated MAGs (the site's eukaryotic and viral genome labels; E##/V## are kept as provenance):
  euks/mags/scripts hull_panels.py -> hulls.tsv (polygons picked by eye; curated input) -> hull_mags.py -> hull_expr.py -> plot_hulls.py
  data/hull_mags.tsv           one row per EE## (and VV## passing the viral genome-group rule below; none at present)
  data/contig_hull_mag.tsv.gz  contig -> EE##/VV## (+ source group, Tiara call)
  data/ee_taxonomy.tsv         per-EE taxonomy (euks/mags/scripts/identify.py); hull_mags.tsv gets taxonomy, _rank, _confidence
  data/organelle_groups.tsv    P## plastid / M## mitochondrial groups (P01 = the former V02, a plastid geNomad called viral)
  data/contig_organelle_group.tsv.gz contig -> P##/M##
  bin_domain.tsv:hull_mag      the EE## holding >= 50% of a eukaryote bin's bp
scripts/common.py (hull_mag*, organelle_group*, euk_group*, virus_group*, contig_host_label, special_label) reads these."""
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
VG = pd.read_csv(EM + 'virus_groups.tsv', sep='\t')
CV = pd.read_csv(EM + 'contig_virus_group.tsv', sep='\t').drop_duplicates('contig')
VG.to_csv(DD + 'virus_groups.tsv', sep='\t', index=False)       # provenance only; genome-group rule applied to the hull MAGs below
CV[['contig', 'vid', 'tiara']].to_csv(DD + 'contig_virus_group.tsv.gz', sep='\t', index=False, compression='gzip')
print('virus groups:', list(VG.vid), '| contigs:', len(CV), '| contigs in both a virus and a euk group:', len(set(CV.contig) & set(CE.contig)))
# hull-curated MAGs (EE## eukaryotic, VV## viral; euks/mags/scripts/hull_mags.py on the hand-drawn polygons in euks/mags/hulls.tsv):
# the site's eukaryotic/viral genome labels. E##/V## above stay in data/ as provenance (each EE/VV records its source group).
HM = pd.read_csv(EM + 'hull_mags.tsv', sep='\t', keep_default_na=False)
CH = pd.read_csv(EM + 'contig_hull_mag.tsv', sep='\t')
# Viral genome-group rule (VV## and V##): >= 50% of the group's bp must be geNomad Nucleocytoviricota or another giant-virus /
# eukaryotic-virus class (GIANT_EUK_VIRUS). Tailed-phage (Caudoviricetes) clusters are many phages that cluster on composition,
# not one genome: they get no genome-group label and their contigs keep the per-contig vOTU handling. 2026-11-08: V01/VV01
# (2,165 contigs, 20.9 Mbp, 97% Caudoviricetes, depth 1-405x) fails the rule, so no VV## exist; VV01 is retired in the registry.
GIANT_EUK_VIRUS = {'Megaviricetes', 'Maveriviricetes', 'Pokkesviricetes', 'Polintoviricetes', 'Alsuviricetes', 'Pisoniviricetes', 'Revtraviricetes'}
GV = pd.read_csv(SRC, sep='\t', usecols=['id', 'length', 'genomad_phylum', 'genomad_class'], low_memory=False).set_index('id')
def giant_share(ids):
    g = GV.reindex(ids).dropna(subset=['length']); ok = g.genomad_phylum.eq('Nucleocytoviricota') | g.genomad_class.isin(GIANT_EUK_VIRUS)
    return float(g.length[ok].sum() / g.length.sum()) if len(g) else 0.0
vv_drop = []
for v in HM.mag[HM.mag.str.startswith('VV')]:
    sh = giant_share(CH.contig[CH.mag == v])
    print(f'{v}: giant/eukaryotic-virus share {sh:.1%}', '-> genome group' if sh >= 0.5 else '-> not a genome group (phage cluster)')
    if sh < 0.5: vv_drop.append(v)
HM = HM[~HM.mag.isin(vv_drop)].copy(); CH = CH[~CH.mag.isin(vv_drop)]
EGi = EG.set_index('group')
num = lambda c: pd.to_numeric(HM[c], errors='coerce')
nh = HM.groupby('source_group').mag.transform('size')                       # E21 -> several hulls (EE21, EE32-EE36)
src_bp = HM.source_group.map(EGi.bp).astype(float).where(nh == 1)
src_ex = HM.source_group.map(EGi.expr_ppm_median).astype(float)
# status, same thresholds as the E## groups: doubtful if Tiara-eukarya purity < 0.9; then under 4 Mbp; then under 10 reads per million
HM['status'] = np.select([HM.mag.str.startswith('VV'), num('tiara_euk') < 0.9, num('bp') < 4e6, num('expr_ppm_median') < 10],
                         ['virus', 'doubtful', '<4 Mbp', 'little CMO RNA'], 'network')
# cautions (rules): expression jumped >= 10x over the source group (to >= 10 per million) -> the hull probably took in an expressed
# neighbour; the hull added >= 4 Mbp and >= 50% of the source group's size; the hull was drawn on a sparse or small core
CAU = []
for i, x in HM.iterrows():
    c = []
    e, se = num('expr_ppm_median')[i], src_ex[i]
    if e >= 10 and np.isfinite(se) and e >= 10 * max(se, 0.1): c.append('uncertain: hull may include an expressed neighbour')
    g, sb = num('bp_gained')[i], src_bp[i]
    if np.isfinite(sb) and g >= 4e6 and g >= 0.5 * sb: c.append(f'hull added {g / 1e6:.1f} Mbp to the {sb / 1e6:.1f} Mbp source group')
    if str(x.note).startswith('UNSURE'): c.append('hull drawn on a sparse or small core')
    CAU.append('; '.join(c))
HM['caution'] = CAU
HM['taxonomy_hint'] = HM.source_group.map(EGi.taxonomy_hint).fillna('').where(nh == 1, '')
# taxonomy of each EE (euks/mags/scripts/identify.py: MarFERReT proteins + SILVA 18S/28S, rule-based confidence high/medium/low)
TX = pd.read_csv(EM + 'ee_taxonomy.tsv', sep='\t', keep_default_na=False)
TX.to_csv(DD + 'ee_taxonomy.tsv', sep='\t', index=False)
_c = TX.set_index('mag').call.str.extract(r'^(?P<name>[^?(]+?)\??\s*\((?P<rank>[^)]+)\)$')
HM['taxonomy'] = HM.mag.map(_c.name).fillna(''); HM['taxonomy_rank'] = HM.mag.map(_c['rank']).fillna('')
HM['taxonomy_confidence'] = HM.mag.map(TX.set_index('mag').confidence).fillna('')
HM[['mag', 'hull', 'source_group', 'status', 'network_ok', 'caution', 'contigs', 'bp', 'gc', 'n50', 'tiara_euk', 'who_euk', 'virus_bp', 'tsne_rms',
    'bp_gained', 'bp_lost', 'excluded_prok_bp', 'expr_ppm_median', 'taxonomy_hint', 'taxonomy', 'taxonomy_rank', 'taxonomy_confidence', 'note']].to_csv(DD + 'hull_mags.tsv', sep='\t', index=False)
assert not CH.contig.duplicated().any()
CH[['contig', 'mag', 'source_group', 'in_source', 'tiara']].to_csv(DD + 'contig_hull_mag.tsv.gz', sep='\t', index=False, compression='gzip')
eh = T.merge(CH[['contig', 'mag']], left_on='id', right_on='contig', how='inner')
hs = eh.groupby(['bin', 'mag']).length.sum().div(bp, level='bin').reset_index(name='share').sort_values('share', ascending=False).drop_duplicates('bin')
D['hull_mag'] = hs[hs.share >= 0.5].set_index('bin').mag.reindex(D.index)
D.loc[D.call != 'eukaryote', 'hull_mag'] = np.nan
# organelle groups (P## plastid / M## mitochondrion; >= 50% Tiara organelle, register_blobs.py). P01 = the former V02.
OG = pd.read_csv(EM + 'organelle_groups.tsv', sep='\t')
OG[['oid', 'kind', 'note', 'contigs', 'bp', 'gc', 'n50', 'tiara_organelle', 'virus_bp', 'plastid_markers', 'mito_markers', 'dastool_bins', 'status', 'expr_ppm_median']].to_csv(DD + 'organelle_groups.tsv', sep='\t', index=False)
CO = pd.read_csv(EM + 'contig_organelle_group.tsv', sep='\t').drop_duplicates('contig')
CO[['contig', 'oid', 'tiara']].to_csv(DD + 'contig_organelle_group.tsv.gz', sep='\t', index=False, compression='gzip')
print('hull MAGs:', len(HM), dict(HM.status.value_counts()), '| network:', int((HM.network_ok.astype(str) == 'True').sum()), '| contigs:', len(CH),
      '| organelle groups:', list(OG.oid), '| contigs in both a hull MAG and an organelle group:', len(set(CH.contig) & set(CO.contig)))
print(HM.loc[HM.caution != '', ['mag', 'caution']].to_string(index=False))
print('eukaryote bins -> hull MAGs:', D.hull_mag.dropna().to_dict())
VG['giant_virus_bp'] = [round(giant_share(CV.contig[CV.vid == v]), 3) for v in VG.vid]
VG['genome_group'] = VG.giant_virus_bp >= 0.5                        # same viral genome-group rule as the VV## above
VG.to_csv(DD + 'virus_groups.tsv', sep='\t', index=False)
print('V groups (provenance):', VG[['vid', 'giant_virus_bp', 'genome_group']].to_dict('records'))
D.sort_index().to_csv(OUT, sep='\t', index_label='bin')
print('eukaryote bins -> groups:', D.euk_group.dropna().to_dict(), '| groups:', len(EG), '| contigs in groups:', len(CE))
print(D.call.value_counts().to_dict())
print('eukaryote:', list(D.index[D.call == 'eukaryote'])); print('organelle:', len(D.index[D.call == 'organelle']), list(D.index[D.call == 'organelle']))
