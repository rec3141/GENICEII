"""BUSCO table for the current EE MAG set (run after merge_mags.py / hull_mags.py; replaces the hand-made table of 2026-11-08).
Parses every short_summary.specific.*.txt under out/busco/ (buscolib.index) and, per current MAG (out/contig_hull_mag.tsv, retired
ids left out), takes the newest run whose input is exactly the MAG's contigs. A MAG with no such run (new or changed membership)
is run now (BUSCO 6.1.0, eukaryota_odb12, miniprot; cached in out/busco/runs/). With --no-run, the newest run under the MAG's id is
used instead and flagged matches_members = False.
Out: out/busco/busco_eukaryota_odb12.tsv (mag, Mbp, C, S, D, F, M, n, matches_members, summary, fasta)"""
import sys, os, re, glob
import pandas as pd
sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
import buscolib as B, eukshared as E
O = E.O
RET = {v for k, v in E.load_ids().items() if k.startswith('retired:')}
CM = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t'); CM = CM[CM.mag.str.startswith('EE') & ~CM.mag.isin(RET)]
norun = '--no-run' in sys.argv
rows = []
for m, g in CM.groupby('mag'):
    ids = sorted(g.contig); c = B.cached(ids)
    if c: r = dict(B.parse_summary(c[0]), summary=c[0], fasta=c[1], matches_members=True)
    elif not norun: r = dict(B.busco(m, ids), matches_members=True)
    else:
        f = sorted(glob.glob(B.O + f'**/*{m}*/short_summary.specific.*.txt', recursive=True), key=os.path.getmtime)
        r = dict(B.parse_summary(f[-1]), summary=f[-1], fasta=B.input_of(f[-1]), matches_members=False) if f else dict(matches_members=False)
    rows.append(dict(mag=m, Mbp=round(g.length.sum() / 1e6, 2), **{k: r.get(k) for k in ('C', 'S', 'D', 'F', 'M', 'n', 'matches_members', 'summary', 'fasta')}))
T = pd.DataFrame(rows)
T['summary'] = T.summary.str.replace(B.O, '', regex=False); T['fasta'] = T.fasta.str.replace(B.O, '', regex=False)
T.to_csv(B.O + 'busco_eukaryota_odb12.tsv', sep='\t', index=False)
pd.set_option('display.width', 220); pd.set_option('display.max_colwidth', 70); print(T.to_string(index=False))
