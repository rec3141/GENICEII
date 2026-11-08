"""site/lineage.json: full lineage per bin for chart subtitles.
GTDB-Tk classification where it exists (empty ranks dropped). Where GTDB has no call (refined *_sub bins, eukaryotes), the
length-weighted Kaiju majority at each rank from the contig table, tagged '(Kaiju)'. Bins called eukaryote/organelle by the
Tiara rule (data/bin_domain.tsv, scripts/bin_domains.py) carry that call first, plus the plastid 16S host where known.
Also: hull-curated eukaryotic/viral MAGs (EE##/VV##), organelle groups (P##), and their provenance groups (E##/V##)."""
import os, json, re
import pandas as pd, numpy as np
D = '/data/dev/rnaseq-exploration/data/'; OUT = '/data/dev/rnaseq-exploration/site/lineage.json'
RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus']
T = pd.read_csv('/data/dev/fastsne/data/microscape/contigs_microscape.tsv.gz', sep='\t', low_memory=False,
                usecols=['bin', 'length'] + [f'kaiju_{r}' for r in RANKS])
T = T[T.bin.notna() & (T.bin != '')]
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
dom = pd.read_csv(D + 'bin_domain.tsv', sep='\t').set_index('bin') if os.path.exists(D + 'bin_domain.tsv') else None
pl = pd.read_csv(D + 'plastid_bins.tsv', sep='\t').set_index('bin')

def kaiju(g):
    out = []
    for r in RANKS:
        s = g.groupby(g[f'kaiju_{r}'].fillna('NA')).length.sum()
        s = s[~s.index.isin(['NA', '', 'Unclassified'])]
        if s.empty: break
        top, share = s.idxmax(), s.max() / g.length.sum()
        if share < 0.5: break                       # stop at the first rank without a length majority
        out.append(f'{top}')
    return out

L = {}
for b, g in T.groupby('bin'):
    short = b.replace('dastool-', '')
    parts = [p.split('__', 1)[1] for p in str(gt.get(short, '')).split(';') if '__' in p and p.split('__', 1)[1]]
    src = 'GTDB'
    if not parts or len(parts) <= 1:
        k = kaiju(g)
        if len(k) > len(parts): parts, src = k, 'Kaiju'
    lin = ' › '.join(parts) + (f' ({src})' if parts else '')
    pre = ''
    if dom is not None and b in dom.index and dom.loc[b, 'call'] != 'prokaryote':
        c = dom.loc[b, 'call']
        host = pl.host.get(b) if b in pl.index and isinstance(pl.host.get(b), str) and pl.host.get(b) and not pl.host.get(b).startswith(('uncultured', 'environmental')) else ''
        pre = ('Eukaryote' if c == 'eukaryote' else 'Plastid/organelle') + ' (Tiara)' + (f', plastid 16S: {host}' if host else '')
        if src == 'GTDB' and c != 'prokaryote': lin = f'GTDB placement: {lin}' if lin else ''
    L[short] = ' · '.join(x for x in (pre, lin) if x) or 'no classification'
# eukaryotic genome groups (E##, euks/mags via scripts/bin_domains.py) and the DAS Tool eukaryote bins that map onto them
EG = pd.read_csv(D + 'euk_groups.tsv', sep='\t', keep_default_na=False)
for x in EG.itertuples():
    bits = [f'Eukaryote (Tiara, {x.tiara_euk:.0%} of length)', f'genome group {x.group}, {x.bp / 1e6:.1f} Mbp, GC {x.gc:.1f}%',
            f'18S: {x.taxonomy_hint}' if x.taxonomy_hint else 'no confident 18S hint', {'network': 'clean, expressed', 'little CMO RNA': 'clean; most of its reads come from non-CMO samples of the coassembly', '<4 Mbp': 'clean, under 4 Mbp', 'low RNA': 'clean, barely expressed', 'doubtful': 'doubtful (mixed or impure)'}.get(x.status, x.status)]
    L[x.group] = ' · '.join(bits)
VGT = pd.read_csv(D + 'virus_groups.tsv', sep='\t', keep_default_na=False)
for x in VGT.itertuples():
    vb = f', {float(x.virus_bp):.0%} viral by geNomad' if str(x.virus_bp) else ''
    L[x.vid] = f'Virus group {x.vid} · {x.bp / 1e6:.1f} Mbp of contigs, GC {x.gc:.1f}%{vb} · from the eukaryotic-region binning ({x.note})'
# hull-curated MAGs (EE## eukaryotic, VV## viral) and organelle groups (P##): the labels the pages use; E##/V## above are provenance
STL = {'network': 'clean, expressed', 'little CMO RNA': 'clean; most of its reads come from non-CMO samples of the coassembly', '<4 Mbp': 'clean, under 4 Mbp', 'doubtful': 'doubtful (mixed or impure)'}
HMT = pd.read_csv(D + 'hull_mags.tsv', sep='\t', keep_default_na=False)
for x in HMT.itertuples():
    if x.mag.startswith('VV'):
        L[x.mag] = f'Virus MAG {x.mag} · {x.bp / 1e6:.1f} Mbp of contigs, GC {x.gc:.1f}%, {float(x.virus_bp):.0%} viral by geNomad · hull drawn on the contig map around group {x.source_group}'
        continue
    bits = [f'Eukaryote (Tiara, {float(x.tiara_euk):.0%} of length)', f'MAG {x.mag}, {x.bp / 1e6:.1f} Mbp, GC {x.gc:.1f}%, from group {x.source_group}',
            f'18S: {x.taxonomy_hint}' if x.taxonomy_hint else 'no confident 18S hint', STL.get(x.status, x.status)]
    if x.caution: bits.append(x.caution)
    L[x.mag] = ' · '.join(bits)
OGT = pd.read_csv(D + 'organelle_groups.tsv', sep='\t', keep_default_na=False)
for x in OGT.itertuples():
    L[x.oid] = (f'{"Plastid" if x.kind == "plastid" else "Mitochondrial"} group {x.oid} (Tiara organelle {float(x.tiara_organelle):.0%} of length) · '
                f'{x.bp / 1e6:.1f} Mbp, GC {x.gc:.1f}% · several algal plastid genomes together' + (' · formerly listed as virus V02' if 'V02' in x.note else ''))
if dom is not None and 'hull_mag' in dom.columns:
    for b, g in dom.hull_mag.dropna().items():
        s_ = b.replace('dastool-', '')
        if s_ in L: L[s_] = L[s_].replace('Eukaryote (Tiara)', f'Eukaryote (Tiara), part of MAG {g}', 1)
json.dump(L, open(OUT, 'w'), separators=(',', ':'), ensure_ascii=False)
print(len(L), 'bins;', os.path.getsize(OUT) // 1000, 'kB')
for b in ('EE02', 'EE05', 'EE01', 'EE12', 'EE32', 'VV01', 'P01', 'E02', 'V01', 'comebin_4318_sub', 'comebin_2576', 'lorbin_184', 'comebin_4527', 'comebin_3937', 'metabat_8153', 'semibin_088'): print(b, '→', L.get(b))
