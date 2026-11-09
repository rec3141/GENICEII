"""Shared BUSCO helpers for the EE MAGs (BUSCO 6.1.0, eukaryota_odb12, miniprot; installed locally with the user's approval).
A run is keyed by the md5 of the MAG's sorted contig ids (`members_md5`), so a MAG or merge candidate whose membership has not
changed is never re-run: runs live in out/busco/runs/<label>_<md5[:10]>/ with the input contig list in <run>/members.txt.
FASTAs are cut from out/busco/ee_contigs.fa.gz + ee_contigs_add1.fa.gz (together they hold every contig of the EE hulls; a
contig missing from both raises an error rather than silently shrinking a MAG).
Runs made before this cache (out/busco/run_*/...) are adopted when the input FASTA named in their summary still exists, is older
than the summary, and holds exactly the requested contigs (`index()`).
parse_summary(path) -> dict(C, S, D, F, M, n) from a short_summary.specific.*.txt."""
import os, re, gzip, glob, hashlib, subprocess
O = '/data/dev/rnaseq-exploration/euks/mags/out/busco/'
SRC = [O + 'ee_contigs.fa.gz', O + 'ee_contigs_add1.fa.gz']
BUSCO_BIN = '/data/scratch/grid-envs/conda/busco/bin'
LINEAGE = O + 'lineages/lineages/eukaryota_odb12'
RUNS = O + 'runs/'
_SEQ = {}

def members_md5(ids):
    return hashlib.md5('\n'.join(sorted(map(str, ids))).encode()).hexdigest()

def _load():
    if not _SEQ:
        for f in SRC:
            name, buf = None, []
            for line in gzip.open(f, 'rt'):
                if line.startswith('>'):
                    if name: _SEQ[name] = ''.join(buf)
                    name, buf = line[1:].split()[0], []
                else: buf.append(line.strip())
            if name: _SEQ[name] = ''.join(buf)
    return _SEQ

def write_fasta(ids, path):
    S = _load(); miss = [c for c in ids if c not in S]
    if miss: raise KeyError(f'{len(miss)} contigs not in {SRC} (e.g. {miss[:3]}); extend the contig FASTAs first')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as fh:
        for c in sorted(ids): fh.write(f'>{c}\n{S[c]}\n')

def parse_summary(path):
    t = open(path).read()
    m = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)', t)
    return dict(zip(['C', 'S', 'D', 'F', 'M', 'n'], [float(x) for x in m.groups()]))

def run_dir(label, ids):
    return f'{RUNS}{label}_{members_md5(ids)[:10]}'

def summary_of(rd):
    f = glob.glob(rd + '/short_summary.specific.*.txt'); return f[0] if f else None

def input_of(summary):
    m = re.search(r'for file (\S+)', open(summary).read()); return m.group(1) if m else None

_IDX = {}
def index():
    """{members md5: [(mtime, summary, input fasta), ...]} over every short_summary.specific.*.txt under out/busco/."""
    if not _IDX:
        for f in glob.glob(O + '**/short_summary.specific.*.txt', recursive=True):
            rd = os.path.dirname(f); fa = input_of(f)
            if os.path.exists(rd + '/members.txt'): md5 = members_md5(open(rd + '/members.txt').read().split())
            elif fa and os.path.exists(fa) and os.path.getmtime(fa) < os.path.getmtime(f):
                md5 = members_md5([l[1:].split()[0] for l in open(fa) if l.startswith('>')])
            else: continue
            _IDX.setdefault(md5, []).append((os.path.getmtime(f), f, fa))
        for v in _IDX.values(): v.sort(reverse=True)       # newest first
    return _IDX

def cached(ids):
    """Newest summary whose input is exactly `ids`, or None: (summary, input fasta)."""
    h = index().get(members_md5(ids)); return h[0][1:] if h else None

def busco(label, ids, threads=8):
    """BUSCO on the contig set `ids` (cached by membership). Returns dict(C, S, D, F, M, n, run, summary, fasta)."""
    c = cached(ids)
    if c: return dict(parse_summary(c[0]), run=os.path.dirname(c[0]), summary=c[0], fasta=c[1])
    rd = run_dir(label, ids); s = summary_of(rd)
    if s is None:
        fa = f'{O}fasta/{os.path.basename(rd)}.fa'; write_fasta(ids, fa); os.makedirs(RUNS, exist_ok=True)
        env = dict(os.environ, PATH=BUSCO_BIN + ':' + os.environ['PATH'])
        r = subprocess.run(['timeout', '7200', BUSCO_BIN + '/python', BUSCO_BIN + '/busco', '-i', fa, '-m', 'genome', '-l', LINEAGE,
                            '--offline', '-c', str(threads), '-o', os.path.basename(rd), '--out_path', RUNS, '-f'],
                           env=env, capture_output=True, text=True, cwd=O)
        s = summary_of(rd)
        if s is None: raise RuntimeError(f'BUSCO failed for {label}: {r.stdout[-2000:]} {r.stderr[-2000:]}')
        open(rd + '/members.txt', 'w').write('\n'.join(sorted(ids)) + '\n')
        if _IDX: _IDX.setdefault(members_md5(ids), []).insert(0, (os.path.getmtime(s), s, input_of(s)))
        for junk in glob.glob(rd + '/run_*/busco_sequences') + glob.glob(rd + '/run_*/miniprot_output'):   # keep the cache small
            subprocess.run(['rm', '-rf', junk])
    return dict(parse_summary(s), run=rd, summary=s, fasta=input_of(s))

def fasta_for(label, ids):
    """A FASTA of exactly `ids` (reuses a cached run's input or out/busco/fasta/<label>_<md5>.fa)."""
    c = cached(ids)
    if c and c[1] and os.path.exists(c[1]): return c[1]
    fa = f'{O}fasta/{os.path.basename(run_dir(label, ids))}.fa'
    if not os.path.exists(fa): write_fasta(ids, fa)
    return fa
