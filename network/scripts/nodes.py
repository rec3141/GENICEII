"""Network nodes and their abundance in the two layers (DNA, RNA). Rule-based; rerun whenever bins, EE MAGs or vOTUs change.

Nodes (all read from files at run time, no hard-coded MAG ids):
  - prokaryote: DAS Tool MAGs that common.is_prokaryote() and not common.is_plastid(), of high or medium quality by CheckM2
    (MIMAG-style: completeness >= 50%, contamination < 10%; network/data/checkm2_dastool.tsv, copied from Grex
    mag_marine/binning/checkm2/quality_report.tsv, the dastool-* rows).
  - eukaryote: every EE## hull MAG in data/hull_mags.tsv, whatever its status. Flags: doubtful (status 'doubtful' or genome_status
    'doubtful, not a genome') and multi-genome (genome_status 'multi-genome ...').
  - organelle: P##/M## organelle groups (data/organelle_groups.tsv) and DAS Tool bins that common.is_plastid().
  - virus: the vOTUs of the Viruses page host table (viruses/out/votu_hosts.tsv: the most abundant and the most diesel-responsive
    vOTUs), plus every vOTU of a eukaryotic-virus class (geNomad Megaviricetes = giant viruses, Maveriviricetes = virophages)
    with a mean of at least EUKV_MIN_CPM viral-CDS CPM over the 128 RNA libraries (viruses/out/votu_rna.tsv).
Contigs are assigned to one node only, in the order virus > organelle group > EE MAG > DAS Tool bin (a giant-virus contig
binned into a plastid bin counts for the virus only).

Abundance:
  - RNA: read pairs on the node's contigs (data/marine_coassembly.contig_counts.mtx.gz), the 128 analysed libraries
    (common.samples()); library size = all contig-assigned read pairs. The same metric as the EE expression tables.
  - DNA: nanopore coverage (depth x length) on the node's contigs (euks/mags/dna/depths.txt.gz, jgi depth table of the
    mag_marine mapping) in CMO pool libraries (CTL/EXP, seawater or sea ice, with a day; blanks, source sites and enrichments are
    left out); resequenced runs of the same sample are summed. Library size = coverage summed over all contigs.
Out (network/cache, untracked): rna_counts.tsv, dna_cov.tsv (node x library), rna_lib.tsv, dna_lib.tsv (library size + design);
network/out/nodes_all.tsv (every candidate node with labels, flags and per-layer prevalence)."""
import sys, os, re, gzip
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
import scipy.io, scipy.sparse as sp

NET = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
CACHE, OUT = os.path.join(NET, 'cache'), os.path.join(NET, 'out')
os.makedirs(CACHE, exist_ok=True); os.makedirs(OUT, exist_ok=True)
E0 = '/data/dev/rnaseq-exploration/'
EUKV_MIN_CPM = 5
EUK_VIRUS = {'Megaviricetes': 'giant virus', 'Maveriviricetes': 'virophage'}

# ---------------------------------------------------------------- candidate nodes
ck = pd.read_csv(os.path.join(NET, 'data', 'checkm2_dastool.tsv'), sep='\t').set_index('Name')
summ = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def gtdb_name(b):
    t = [x.split('__', 1)[1] for x in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in x and x.split('__', 1)[1]]
    return t[-1] if t else ''
def gtdb_rank(b, r):
    for x in str(gt.get(b.replace('dastool-', ''), '')).split(';'):
        if x.startswith(r + '__') and x[3:]: return x[3:]
    return ''

N = []
for b in summ.index:
    s = b.replace('dastool-', '')
    if is_plastid(b):
        N.append(dict(id=s, domain='organelle', kind='plastid bin', label=plastid_label(b) or 'Plastid/organelle (Tiara)', bp=int(summ.loc[b, 'size']),
                      completeness=ck.Completeness.get(b, np.nan), contamination=ck.Contamination.get(b, np.nan), status='', member=b))
    elif is_prokaryote(b):
        c, x = ck.Completeness.get(b, np.nan), ck.Contamination.get(b, np.nan)
        if not (c >= 50 and x < 10): continue
        nm = gtdb_name(b) or kaiju_bin_label(b) or 'unclassified'
        N.append(dict(id=s, domain='prokaryote', kind='HQ' if (c >= 90 and x < 5) else 'MQ', label=nm, bp=int(summ.loc[b, 'size']), completeness=c, contamination=x,
                      status='', member=b, phylum=gtdb_rank(b, 'p'), cls=gtdb_rank(b, 'c'), order=gtdb_rank(b, 'o'), family=gtdb_rank(b, 'f'), genus=gtdb_rank(b, 'g')))
HM = pd.read_csv(D + 'hull_mags.tsv', sep='\t', keep_default_na=False)
for x in HM[HM.mag.str.startswith('EE')].itertuples():
    gs = str(getattr(x, 'genome_status', '') or '')           # PR #18: 'doubtful, not a genome', 'multi-genome (>= 2 strains)'
    N.append(dict(id=x.mag, domain='eukaryote', kind='EE MAG', label=hull_mag_label(x.mag), bp=int(x.bp), status=x.status, genome_status=gs,
                  doubtful=x.status == 'doubtful' or gs.startswith('doubtful'), multi=gs.startswith('multi-genome'), caution=x.caution,
                  taxonomy=x.taxonomy if x.taxonomy_confidence in ('high', 'medium') else '', taxonomy_guess=x.taxonomy, taxonomy_rank=x.taxonomy_rank,
                  taxonomy_confidence=x.taxonomy_confidence, member=x.mag))
OG = pd.read_csv(D + 'organelle_groups.tsv', sep='\t', keep_default_na=False)
for x in OG.itertuples():
    # an organelle genome is 0.1-0.4 Mbp, so a group over 1 Mbp holds several (P01: 5.0 Mbp, 486 plastid marker genes)
    N.append(dict(id=x.oid, domain='organelle', kind='organelle group', label=organelle_group_label(x.oid), bp=int(float(x.bp)),
                  status='multi-genome' if float(x.bp) > 1e6 else '', multi=float(x.bp) > 1e6, member=x.oid))
VR = pd.read_csv(E0 + 'viruses/out/votu_rna.tsv', sep='\t', index_col=0)
VH = pd.read_csv(E0 + 'viruses/out/votu_hosts.tsv', sep='\t', index_col=0)
vsel = list(VH.index) + [v for v in VR.index[VR.cls.isin(list(EUK_VIRUS)) & (VR.mean_cpm >= EUKV_MIN_CPM)] if v not in VH.index]
for v in vsel:
    c = VR.cls.get(v, '') if isinstance(VR.cls.get(v, ''), str) else ''
    fam = VR.fam.get(v, '') if isinstance(VR.fam.get(v, ''), str) else ''
    kind = EUK_VIRUS.get(c, 'phage' if c in ('Caudoviricetes', 'Malgrandaviricetes', 'Faserviricetes') else 'virus')
    N.append(dict(id=v, domain='virus', kind=kind, label=f'{kind[0].upper() + kind[1:]} {v}' + (f' ({fam})' if fam else (f' ({c})' if c and kind in ('phage', 'virus') else '')),
                  bp=int(VR.length.get(v, 0)), status='', member=v, vclass=c, vfamily=fam, checkv=VR.checkv_quality.get(v, ''), on_viruses_page=v in VH.index,
                  binned=VR.bin.get(v, '') if isinstance(VR.bin.get(v, ''), str) else ''))
N = pd.DataFrame(N).set_index('id')
print(N.groupby(['domain', 'kind']).size().to_string())

# ---------------------------------------------------------------- contig -> node (one node per contig)
vc = pd.read_csv(E0 + 'viruses/out/viral_confident.tsv', sep='\t', usecols=['seq_name', 'region', 'provirus', 'votu_rep'])
vc = vc[(vc.provirus == 'No') & (vc.region == vc.seq_name) & vc.votu_rep.isin(N.index[N.domain == 'virus'])]
asg = {}
[asg.setdefault(c, v) for c, v in zip(vc.seq_name, vc.votu_rep)]
og = pd.read_csv(D + 'contig_organelle_group.tsv.gz', sep='\t')
[asg.setdefault(c, o) for c, o in zip(og.contig, og.oid) if o in N.index]
hm = pd.read_csv(D + 'contig_hull_mag.tsv.gz', sep='\t')
[asg.setdefault(c, m) for c, m in zip(hm.contig, hm.mag) if m in N.index]
c2b = pd.read_csv(D + 'contig2bin.tsv', sep='\t', header=None, names=['contig', 'bin'])
c2b['node'] = c2b.bin.str.replace('dastool-', '')
[asg.setdefault(c, b) for c, b in zip(c2b.contig, c2b.node) if b in N.index]
A = pd.Series(asg, name='node')
N['contigs'] = A.value_counts().reindex(N.index).fillna(0).astype(int)
print('contigs assigned:', len(A), '; nodes without contigs:', list(N.index[N.contigs == 0]))

# ---------------------------------------------------------------- RNA
ids = np.sort(pd.read_pickle(E0 + 'euks/mags/out/contigs_full.pkl').id.values)
ALL, S = all_samples(), samples(); kc = [ALL.index(x) for x in S]
C = scipy.io.mmread(D + 'marine_coassembly.contig_counts.mtx.gz').tocsr()[:, kc]
assert C.shape[0] == len(ids)
pos = pd.Series(np.arange(len(ids)), index=ids)
codes = pd.Categorical(A.values, categories=N.index).codes
M = sp.csr_matrix((np.ones(len(A)), (codes, pos[A.index].values)), shape=(len(N), len(ids)))
R = pd.DataFrame((M @ C).toarray(), index=N.index, columns=S)
m = meta().loc[S]
light = pd.read_csv(E0 + 'light/out/underice_samples.tsv', sep='\t', index_col=0)
rl = pd.DataFrame({'size': np.asarray(C.sum(0)).ravel(), 'season': m.season, 'group': m.group, 'medium': m.medium, 'day': m.day,
                   'layer': m.ice_layer.fillna(''), 'temp': pd.to_numeric(m.water_temp_c, errors='coerce'), 'par_under': light.par_under.reindex(S),
                   'name': m.sample_name}, index=S)
R.to_csv(os.path.join(CACHE, 'rna_counts.tsv'), sep='\t'); rl.to_csv(os.path.join(CACHE, 'rna_lib.tsv'), sep='\t')

# ---------------------------------------------------------------- DNA
ns = nanopore_samples()
ns['season'] = ns.log.astype(str).map({'2024': '2024-25', '2025': '2025-26'})
u = ns.sample_name.str.upper()
ns['group'] = np.where(u.str.contains('_EXP_'), 'EXP', np.where(u.str.contains('_CTL_'), 'CTL', 'other'))
ns['medium'] = np.where(u.str.contains('_SI'), 'SI', np.where(u.str.contains('_SW'), 'SW', 'other'))
def day(s):
    x = re.search(r'(PRE|DAY)(\d+)', s.upper()); return np.nan if not x else (-int(x[2]) if x[1] == 'PRE' else int(x[2]))
ns['day'] = ns.sample_name.map(day)
ns = ns[(ns.group != 'other') & (ns.medium != 'other') & ns.day.notna() & ~ns.run_barcode.isin(excluded('nanopore'))]
f = E0 + 'euks/mags/dna/depths.txt.gz'
hdr = gzip.open(f, 'rt').readline().rstrip('\n').split('\t')
runs = [c for c in hdr[3:] if not c.endswith('-var') and c.replace('.sorted.bam', '') in set(ns.run_barcode)]
ix = [hdr.index(c) for c in runs]; runs = [c.replace('.sorted.bam', '') for c in runs]
node_of = A.to_dict(); nidx = {n: i for i, n in enumerate(N.index)}
cov = np.zeros((len(N), len(runs))); tot = np.zeros(len(runs))
row_of = {c: nidx[n] for c, n in node_of.items()}
for ch in pd.read_csv(f, sep='\t', usecols=[0, 1] + ix, chunksize=200000):
    v = ch.iloc[:, 2:].to_numpy(float) * ch.iloc[:, 1].to_numpy(float)[:, None]
    tot += v.sum(0)
    k = ch.iloc[:, 0].map(row_of)
    ok = k.notna().to_numpy()
    np.add.at(cov, k[ok].astype(int).to_numpy(), v[ok])
ns = ns.set_index('run_barcode').loc[runs]
ns['key'] = ns.season + '|' + ns.sample_name
Dc = pd.DataFrame(cov, index=N.index, columns=runs).T.groupby(ns.key.values).sum().T
dt = pd.Series(tot, index=runs).groupby(ns.key.values).sum()
li = ns.drop_duplicates('key').set_index('key')
pud = pd.read_csv(E0 + 'light/out/underice_daily.tsv', sep='\t')
par = {(a, b, c): v for a, b, c, v in zip(pud.season, pud.pool, pud.day, pud.par_under)}
# DNA libraries carry no temperature; take the water temperature of the RNA sample of the same winter, pool and day if any
tmp = rl.dropna(subset=['temp']).groupby(['season', 'group', 'day']).temp.mean().to_dict()
dl = pd.DataFrame({'size': dt, 'season': li.season, 'group': li.group, 'medium': li.medium, 'day': li.day,
                   'layer': ['TOP' if '_TOP' in s.upper() else 'BOT' if '_BOT' in s.upper() else '' for s in li.sample_name],
                   'temp': [tmp.get((a, b, c), np.nan) for a, b, c in zip(li.season, li.group, li.day)],
                   'par_under': [par.get((a, b, c), np.nan) for a, b, c in zip(li.season, li.group, li.day)],
                   'name': li.sample_name, 'runs': ns.groupby('key').size()}).loc[Dc.columns]
Dc.to_csv(os.path.join(CACHE, 'dna_cov.tsv'), sep='\t'); dl.to_csv(os.path.join(CACHE, 'dna_lib.tsv'), sep='\t')

# ---------------------------------------------------------------- prevalence per layer
N['rna_prev'] = (R > 0).mean(1).round(3); N['rna_share_ppm'] = (R.div(rl['size'], axis=1).mean(1) * 1e6).round(2)
N['dna_prev'] = (Dc > 0).mean(1).round(3); N['dna_share_ppm'] = (Dc.div(dl['size'], axis=1).mean(1) * 1e6).round(2)
N.to_csv(os.path.join(OUT, 'nodes_all.tsv'), sep='\t')
print('RNA libraries', R.shape[1], '; DNA samples', Dc.shape[1], '(runs', len(runs), ')')
print(N.groupby('domain')[['rna_prev', 'dna_prev', 'rna_share_ppm', 'dna_share_ppm']].describe().T.round(3).to_string())
