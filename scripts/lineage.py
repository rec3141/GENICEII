"""site/lineage.json: full lineage per bin for chart subtitles.
GTDB-Tk classification where it exists (empty ranks dropped). Where GTDB has no call (refined *_sub bins, eukaryotes), the
length-weighted Kaiju majority at each rank from the contig table, tagged '(Kaiju)'. Bins called eukaryote/organelle by the
Tiara rule (data/bin_domain.tsv, scripts/bin_domains.py) carry that call first, plus the plastid 16S host where known."""
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
json.dump(L, open(OUT, 'w'), separators=(',', ':'), ensure_ascii=False)
print(len(L), 'bins;', os.path.getsize(OUT) // 1000, 'kB')
for b in ('comebin_4318_sub', 'comebin_2576', 'lorbin_184', 'comebin_4527', 'comebin_3937', 'metabat_8153', 'semibin_088'): print(b, '→', L.get(b))
