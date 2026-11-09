"""Merge rule for the hull MAGs (EE##), run after hull_mags.py (user-approved 2026-11-08).
Two hull MAGs are one genome cut in pieces (e.g. by t-SNE) when they rise and fall together, do not share sequence, and their union
is a better genome than either part. A candidate pair (A, B) is merged only if ALL of these hold:
  1. RNA rho >= 0.95: Spearman of the outlier-resistant per-library RNA profiles (robustexpr.robust_profile: per-library 10%
     trimmed mean over the mean-scaled profiles of contigs detected in >= 25% of the 128 libraries, the hull_expr.py metric since 2026-11-08). The summed
     profile (eukshared.rna_profile) let a few contigs decide: EE03+EE31 had summed rho 0.53 but robust rho 0.996;
  2. DNA rho >= 0.8: Spearman of the summed CMO nanopore depth x length shares (eukshared.dna_profile);
  3. skani aligned fraction < 10% in both directions (skani dist --min-af 0; no output = ANI screen < 80% = AF 0): the parts are
     not two copies (strains/haplotypes) of the same sequence;
  4. BUSCO (6.1.0, eukaryota_odb12) on the union: duplication < 5% AND completeness higher than the best part.
Transitive: candidate pairs are taken in decreasing RNA rho; when A or B already is a merged group, rules 1-4 are re-tested on the
groups (summed profiles of the groups, skani between the group FASTAs, BUSCO on the whole union).
The merged MAG keeps the lowest EE number; the others are recorded in out/id_registry.tsv as 'retired:EEnn (merged into EEmm,
<date>)' (dates of merges that persist are kept; a merge that no longer passes is removed, so its EE id returns). hull_mags.py
applies those keys; this script reruns hull_mags.py when the set of merges changes.
Near-miss BUSCO: unions are also BUSCO'd for pairs from the same source group (e.g. the E21 pieces EE21/EE32-EE36) and for pairs
with RNA rho >= 0.9 that fail one rule only, so the report shows what a merge would have given. BUSCO and skani results are cached
(buscolib: out/busco/runs/, keyed by membership), so reruns are cheap.
Pre-merge membership: contig_hull_mag.tsv `hull` + registry 'hull:<hull id>'.
Out: out/mag_merge_tests.tsv (every EE pair: rho, AF, BUSCO of parts and union where run, decision), registry keys."""
import sys, re, os, itertools, subprocess, datetime
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
import eukshared as E, buscolib as B, robustexpr as RX
O = E.O
RNA_MIN, DNA_MIN, AF_MAX, DUP_MAX = 0.95, 0.8, 10.0, 5.0
NEAR_RNA = 0.9
SKANI = '/data/dev/danaSeq/mag_analysis/conda-envs/dana-mag-derep/bin/skani'
TODAY = str(datetime.date.today())

IDS = E.load_ids()
CM = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t'); CM = CM[CM.mag.str.startswith('EE')]
CM['raw'] = [IDS['hull:' + h] for h in CM.hull]
SRC = CM.groupby('raw').source_group.first().to_dict()
MEM = {m: sorted(g.contig) for m, g in CM.groupby('raw')}
key = lambda grp: '+'.join(sorted(grp))
def members(grp): return sorted(set().union(*(MEM[m] for m in grp)))
_prof = {}
def prof(grp):
    k = key(grp)
    if k not in _prof: ids = members(grp); _prof[k] = (RX.robust_profile(ids), E.dna_profile(ids))
    return _prof[k]
def rhos(a, b):
    (ra, da), (rb, db) = prof(a), prof(b); return E.spearman(ra, rb), E.spearman(da, db)
def skani(a, b):
    fa, fb = B.fasta_for(key(a), members(a)), B.fasta_for(key(b), members(b))
    out = subprocess.run([SKANI, 'dist', '-q', fa, '-r', fb, '--min-af', '0', '-t', '4'], capture_output=True, text=True).stdout.strip().split('\n')
    if len(out) < 2: return 0.0, 0.0, np.nan          # below skani's ~80% ANI screen: nothing aligned
    f = out[1].split('\t'); return float(f[4]), float(f[3]), float(f[2])     # AF of a (query), AF of b (ref), ANI
def busco(grp): return B.busco(key(grp).replace('+', '_'), members(grp), threads=5)
fmt = lambda r: f"{r['C']:.1f} [{r['D']:.1f}]" if r else ''

# 1) every pair: profiles; skani for pairs worth a look
rows = []
for a, b in itertools.combinations(sorted(MEM), 2):
    r, d = rhos([a], [b])
    rows.append(dict(a=a, b=b, same_source=SRC[a] == SRC[b], rna_rho=round(r, 3), dna_rho=round(d, 3)))
P = pd.DataFrame(rows)
look = (P.rna_rho >= 0.8) | P.same_source
for i in P.index[look]:
    af_a, af_b, ani = skani([P.a[i]], [P.b[i]]); P.loc[i, ['af_a', 'af_b', 'ani']] = af_a, af_b, ani
P['pass_rna'] = P.rna_rho >= RNA_MIN; P['pass_dna'] = P.dna_rho >= DNA_MIN
P['pass_af'] = (P.af_a < AF_MAX) & (P.af_b < AF_MAX)
P['n_fail'] = (~P.pass_rna).astype(int) + (~P.pass_dna).astype(int) + (~P.pass_af).astype(int)
cand = P[(P.n_fail == 0)].sort_values('rna_rho', ascending=False)
near = P[(P.n_fail > 0) & look & (P.same_source | ((P.rna_rho >= NEAR_RNA) & (P.n_fail == 1)))]
# 2) BUSCO of parts and unions needed (parallel; cached)
need = {key([m]): [m] for i in pd.concat([cand, near]).index for m in (P.a[i], P.b[i])}
need.update({key([P.a[i], P.b[i]]): [P.a[i], P.b[i]] for i in pd.concat([cand, near]).index})
B.index(); B._load()      # build once before the worker threads
with ThreadPoolExecutor(3) as ex: BR = dict(zip(need, ex.map(busco, need.values())))
for i in pd.concat([cand, near]).index:
    a, b = P.a[i], P.b[i]; u = BR[key([a, b])]
    P.loc[i, ['busco_a', 'busco_b', 'busco_union']] = fmt(BR[a]), fmt(BR[b]), fmt(u)
    P.loc[i, 'pass_busco'] = bool(u['D'] < DUP_MAX and u['C'] > max(BR[a]['C'], BR[b]['C']))
# 3) transitive merge with re-tests
grp = {m: [m] for m in MEM}; tests = []
for i in cand.index:
    a, b = P.a[i], P.b[i]; A, Bg = grp[a], grp[b]
    if A is Bg: P.loc[i, 'decision'] = 'already in one group'; continue
    if len(A) == 1 and len(Bg) == 1:
        ok = bool(P.pass_busco[i]); why = 'merge' if ok else 'BUSCO: ' + ('duplication >= 5%' if BR[key([a, b])]['D'] >= DUP_MAX else 'union not more complete than the best part')
    else:      # re-test the growing group
        r, d = rhos(A, Bg); af_a, af_b, ani = skani(A, Bg); ua, ub, uu = busco(A), busco(Bg), busco(A + Bg)
        ok = r >= RNA_MIN and d >= DNA_MIN and af_a < AF_MAX and af_b < AF_MAX and uu['D'] < DUP_MAX and uu['C'] > max(ua['C'], ub['C'])
        why = 'merge (group re-test)' if ok else 'group re-test failed'
        tests.append(dict(a=key(A), b=key(Bg), same_source=False, rna_rho=round(r, 3), dna_rho=round(d, 3), af_a=af_a, af_b=af_b, ani=ani,
                          busco_a=fmt(ua), busco_b=fmt(ub), busco_union=fmt(uu), decision=why))
    P.loc[i, 'decision'] = why
    if ok:
        new = sorted(A + Bg)
        for m in new: grp[m] = new
def reason(x):
    if isinstance(x.decision, str): return x.decision
    f = [n for n, ok in (('RNA rho < 0.95', x.pass_rna), ('DNA rho < 0.8', x.pass_dna)) if not ok]
    if x.af_a == x.af_a:
        if not x.pass_af: f.append('skani AF >= 10%')
    else: f.append('(skani not run)')
    return 'no: ' + ', '.join(f)
P['decision'] = P.apply(reason, axis=1)
P = pd.concat([P, pd.DataFrame(tests)], ignore_index=True)
cols = ['a', 'b', 'same_source', 'rna_rho', 'dna_rho', 'af_a', 'af_b', 'ani', 'busco_a', 'busco_b', 'busco_union', 'decision']
P = P[cols].sort_values(['rna_rho'], ascending=False); P.to_csv(O + 'mag_merge_tests.tsv', sep='\t', index=False)
# 4) registry
merges = {m: min(g, key=lambda x: int(x[2:])) for m, g in grp.items() if len(g) > 1}
merges = {m: k for m, k in merges.items() if m != k}
old = {}
for k in list(IDS):
    mt = E.MERGE_RX.match(k)
    if mt: old[mt.group(1)] = (mt.group(2), re.search(r'(\d{4}-\d\d-\d\d)\)$', k).group(1)); IDS.pop(k)
for m, k in sorted(merges.items()):
    date = old[m][1] if m in old and old[m][0] == k else TODAY
    IDS[f'retired:{m} (merged into {k}, {date})'] = m
E.save_ids(IDS)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
show = P[(P.rna_rho >= 0.8) | P.same_source.astype(bool) | P.decision.str.startswith('merge')]
print(show.to_string(index=False)); print('merges:', merges or 'none')
if {m: k for m, (k, _) in old.items()} != merges:
    print('merge set changed -> rerunning hull_mags.py')
    subprocess.run([sys.executable, '/data/dev/rnaseq-exploration/euks/mags/scripts/hull_mags.py'], check=True, stdout=open(O + 'run_hull.log', 'w'))
