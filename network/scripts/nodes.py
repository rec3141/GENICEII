"""Network nodes and their abundance in the three layers (nanopore DNA, Illumina DNA, RNA). Rule-based; rerun whenever bins, EE MAGs or vOTUs change.

Nodes (all read from files at run time, no hard-coded MAG ids):
  - prokaryote: DAS Tool MAGs that common.is_prokaryote() and not common.is_plastid(), of high or medium quality by CheckM2
    (network/data/checkm2_dastool.tsv, copied from Grex mag_marine/binning/checkm2/quality_report.tsv, the dastool-* rows):
    completeness >= MIN_COMPL (15%) and contamination < MAX_CONTAM (20%) (since 2026-11-09, second change; before: >= 30% / < 20%,
    and first >= 50% / < 10%, MIMAG medium quality). It brings in fragments and near-complete bins with co-binned strains, among
    them the Cycloclasticus bloom genomes lorbin_184 (94.7% / 10.2%) and metabat_5432 (100% / 12.8%), and the fragments
    lorbin_722_sub (Cycloclasticus, 23%), lorbin_783 (Moritella, 18%) and metabat_8071_sub (Sulfitobacter, 18%). Flags: fragment (completeness < 50%) and contam_10_20 (contamination >= 10%).
  - eukaryote: every EE## hull MAG in data/hull_mags.tsv, whatever its status. Flags: doubtful (status 'doubtful' or genome_status
    'doubtful, not a genome') and multi-genome (genome_status 'multi-genome ...').
  - organelle: P##/M## organelle groups (data/organelle_groups.tsv) and DAS Tool bins that common.is_plastid().
  Pool expression rule (both layers, before the prevalence rule): prokaryotic MAGs need bin TPM (data/marine_coassembly.bin_tpm.tsv)
    of at least POOL_TPM_MEAN on average over the CMO pool libraries OR at least POOL_TPM_MAX in one of them; EE MAGs, organelle groups
    and plastid bins need their robust RNA profile detected in at least POOL_DET of the pool libraries. vOTUs keep their own rule.
    Contigs of a dropped node are not handed to another node.
  - virus: the vOTUs of the Viruses page host table (viruses/out/votu_hosts.tsv: the most abundant and the most diesel-responsive
    vOTUs), plus every vOTU of a eukaryotic-virus class (geNomad Megaviricetes = giant viruses, Maveriviricetes = virophages)
    with a mean of at least EUKV_MIN_CPM viral-CDS CPM over the 128 RNA libraries (viruses/out/votu_rna.tsv).
Contigs are assigned to one node only, in the order virus > organelle group > EE MAG > DAS Tool bin (a giant-virus contig
binned into a plastid bin counts for the virus only).

Abundance:
  - RNA: the outlier-resistant profile of euks/mags/scripts/robustexpr.py (2026-11-08), applied to every node so the layer is
    consistent: per contig reads per kb per million contig-assigned reads (data/marine_coassembly.contig_counts.mtx.gz, the 128
    analysed libraries, common.samples()); contigs detected in >= 25% of libraries, each scaled to mean 1; per-library 10% trimmed
    mean x the node's median contig level. A few rRNA, conserved-gene or neighbour contigs therefore cannot set a MAG's profile
    (EE04: 94% of its reads on one rRNA contig). The methods need counts, so the profile is turned back into read-pair
    pseudo-counts: profile x node kb x library size / 1e6. Nodes with < 5 such contigs (every vOTU, small bins) fall back to the
    length-weighted mean, which gives back their summed read pairs exactly. Library size = all contig-assigned read pairs. The old
    summed counts are kept in cache/rna_counts_summed.tsv; nodes_all.tsv records rna_contigs_used, rna_thin (< 20 contigs) and
    rna_rho_robust_summed.
  - DNA: nanopore coverage (depth x length) on the node's contigs (euks/mags/dna/depths.txt.gz, jgi depth table of the
    mag_marine mapping) in CMO pool libraries (CTL/EXP, seawater or sea ice, with a day; blanks, source sites and enrichments are
    left out); resequenced runs of the same sample are summed. Library size = coverage summed over all contigs.
  - Illumina DNA (layer 'mgdna', added 2026-11-09): read pairs on the node's contigs in the Illumina metagenomes of the CMO pools
    (common.illumina_dna_samples / illumina_dna_contig_matrix: lanes summed per library; QC-flagged libraries and libraries with
    < 1 M pairs on contigs left out; CTL/EXP, seawater or sea ice, with a day: 85 libraries). Library size = all contig-assigned read
    pairs. With fixed-length reads this is the nanopore metric (depth x length share), so both DNA layers use one metric; it is
    summed over the node's contigs, as in the nanopore layer (DNA is not spiky the way RNA is; see euks/mags/scripts/mgdna.py).
Out (network/cache, untracked): rna_counts.tsv, dna_cov.tsv, mgdna_counts.tsv (node x library), rna_lib.tsv, dna_lib.tsv,
mgdna_lib.tsv (library size + design);
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
MIN_COMPL, MAX_CONTAM = 15, 20          # prokaryotic bins (CheckM2); flags: fragment < 50% complete, contamination >= 10%
# pool expression (2026-11-09): which genomes to show at all, before the per-layer prevalence rule (netcommon), in both layers.
# Pool libraries: common.meta() group CTL or EXP and medium SW or SI (both winters, before and after dosing; site, enrichment,
# MOS and other libraries do not count).
POOL_TPM_MEAN, POOL_TPM_MAX = 10, 100     # prokaryotic MAGs: bin TPM mean over pool libraries >= 10, OR >= 100 in one (short bloomers)
POOL_DET = 0.10                           # EE MAGs, organelle groups and plastid bins: robust RNA profile > 0 in >= 10% of pool libraries
EUK_VIRUS = {'Megaviricetes': 'giant virus', 'Maveriviricetes': 'virophage'}

# ---------------------------------------------------------------- candidate nodes
ck = pd.read_csv(os.path.join(NET, 'data', 'checkm2_dastool.tsv'), sep='\t').set_index('Name')
summ = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def gtdb_cls(b):
    # the DAS Tool bin's own row ('dastool-<bin>') first; the plain binner name only as a fallback (before 2026-11-09 only the
    # plain name was read, so *_sub bins, which exist only as dastool-* rows, came out 'unclassified')
    s = 'dastool-' + b.replace('dastool-', '')
    return str(gt.get(s, gt.get(b.replace('dastool-', ''), '')))
def gtdb_name(b):
    t = [x.split('__', 1)[1] for x in gtdb_cls(b).split(';') if '__' in x and x.split('__', 1)[1]]
    return t[-1] if t else ''
def gtdb_rank(b, r):
    for x in gtdb_cls(b).split(';'):
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
        if not (c >= MIN_COMPL and x < MAX_CONTAM): continue
        nm = gtdb_name(b) or kaiju_bin_label(b) or 'unclassified'
        N.append(dict(id=s, domain='prokaryote', kind='HQ' if (c >= 90 and x < 5) else 'MQ', label=nm, bp=int(summ.loc[b, 'size']), completeness=c, contamination=x,
                      fragment=bool(c < 50), contam_10_20=bool(x >= 10),
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
Rsum = pd.DataFrame((M @ C).toarray(), index=N.index, columns=S)
# Illumina DNA (metagenomes): read pairs on the node's contigs in the CMO pool libraries (common.illumina_dna_samples: CTL/EXP,
# seawater or sea ice, with a day; QC-flagged and < 1 M pair libraries already left out). Each library is its own sample.
MGS = illumina_dna_samples(); MGS = MGS[MGS.pool]
Cm, mlibs = illumina_dna_contig_matrix(); Cm = Cm[:, [mlibs.index(x) for x in MGS.index]]
assert Cm.shape[0] == len(ids)
Mg = pd.DataFrame((M @ Cm).toarray(), index=N.index, columns=MGS.index)
mg_size = np.asarray(Cm.sum(0)).ravel()
sys.path.insert(0, E0 + 'euks/mags/scripts'); import robustexpr as RX
from scipy.stats import spearmanr
RP, RI = RX.robust_profiles(A)                                       # node x library, reads per kb per million
size = np.asarray(C.sum(0)).ravel()
R = RP.reindex(N.index).fillna(0).mul(RI.kb.reindex(N.index).fillna(0), axis=0) * size / 1e6
N['rna_contigs_used'] = RI.n_used.reindex(N.index).fillna(0).astype(int); N['rna_method'] = RI.method.reindex(N.index).fillna('')
N['rna_thin'] = RI.thin.reindex(N.index).fillna(False).astype(bool) & (N.rna_method == 'robust')
N['rna_rho_robust_summed'] = [round(spearmanr(R.loc[n], Rsum.loc[n]).correlation, 3) if Rsum.loc[n].sum() > 0 else np.nan for n in N.index]
# ---------------------------------------------------------------- pool expression rule
mt = meta(); pool = [x for x in S if x in mt.index and mt.group[x] in ('CTL', 'EXP') and mt.medium[x] in ('SW', 'SI')]
tpm = pd.read_csv(D + 'marine_coassembly.bin_tpm.tsv', sep='\t', index_col=0)
tpm.index = tpm.index.str.replace('dastool-', '')
tp = tpm.reindex(N.index)[pool]
N['pool_tpm_mean'] = tp.mean(axis=1).round(2); N['pool_tpm_max'] = tp.max(axis=1).round(1)
N['pool_det'] = (R[pool] > 0).mean(axis=1).round(3)
prok = N.domain == 'prokaryote'; other = N.domain.isin(['eukaryote', 'organelle'])
N['pool_expressed'] = np.where(prok, (N.pool_tpm_mean >= POOL_TPM_MEAN) | (N.pool_tpm_max >= POOL_TPM_MAX),
                               np.where(other, N.pool_det >= POOL_DET, True))
drop = N.index[~N.pool_expressed.astype(bool)]
print(f'pool expression rule ({len(pool)} pool libraries): drops {len(drop)} nodes:', N.loc[drop].domain.value_counts().to_dict(), list(drop[N.loc[drop].domain != 'prokaryote']))
N.loc[drop].to_csv(os.path.join(OUT, 'nodes_dropped_pool_expression.tsv'), sep='\t')
N = N.drop(index=drop); R = R.drop(index=drop); Rsum = Rsum.drop(index=drop); Mg = Mg.drop(index=drop)

print('robust vs summed RNA, Spearman by domain:'); print(N.groupby('domain').rna_rho_robust_summed.describe().round(3).to_string())
print('EE / organelle-group nodes:'); print(N.loc[N.kind.isin(['EE MAG', 'organelle group']), ['rna_contigs_used', 'rna_thin', 'rna_rho_robust_summed']].to_string())
Rsum.to_csv(os.path.join(CACHE, 'rna_counts_summed.tsv'), sep='\t')
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
row_of = {c: nidx[n] for c, n in node_of.items() if n in nidx}
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

# ---------------------------------------------------------------- Illumina DNA library table
ml = pd.DataFrame({'size': mg_size, 'season': MGS.season, 'group': MGS.group, 'medium': MGS.medium, 'day': MGS.day,
                   'layer': MGS.layer, 'temp': [tmp.get((a, b, c), np.nan) for a, b, c in zip(MGS.season, MGS.group, MGS.day)],
                   'par_under': [par.get((a, b, c), np.nan) for a, b, c in zip(MGS.season, MGS.group, MGS.day)],
                   'name': MGS.sample_name}, index=MGS.index)
Mg.to_csv(os.path.join(CACHE, 'mgdna_counts.tsv'), sep='\t'); ml.to_csv(os.path.join(CACHE, 'mgdna_lib.tsv'), sep='\t')

# ---------------------------------------------------------------- prevalence per layer
N['rna_prev'] = (R > 0).mean(1).round(3); N['rna_share_ppm'] = (R.div(rl['size'], axis=1).mean(1) * 1e6).round(2)
N['dna_prev'] = (Dc > 0).mean(1).round(3); N['dna_share_ppm'] = (Dc.div(dl['size'], axis=1).mean(1) * 1e6).round(2)
N['mgdna_prev'] = (Mg > 0).mean(1).round(3); N['mgdna_share_ppm'] = (Mg.div(ml['size'], axis=1).mean(1) * 1e6).round(2)
N.to_csv(os.path.join(OUT, 'nodes_all.tsv'), sep='\t')
print('RNA libraries', R.shape[1], '; DNA samples', Dc.shape[1], '(runs', len(runs), '); Illumina DNA libraries', Mg.shape[1])
print(N.groupby('domain')[['rna_prev', 'dna_prev', 'mgdna_prev', 'rna_share_ppm', 'dna_share_ppm', 'mgdna_share_ppm']].describe().T.round(3).to_string())
