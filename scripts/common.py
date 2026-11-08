"""Shared loaders for the CMO RNA-seq exploration (marine coassembly expression + sample metadata)."""
import re
import numpy as np, pandas as pd

D = '/data/dev/rnaseq-exploration/data/'
META = '/data/dev/genice-sample-metadata/outputs/cmo_rna_metadata.tsv'
# oil addition: 2024-25 diesel on DAY05 16:30 (after the DAY05 11:00 samples); 2025-26 136 L oil on DAY01
OIL_DAY = {'CMO 2024-25': 5.5, 'CMO 2025-26': 1}
CORR = '/data/dev/genice-sample-metadata/sources/sample_corrections.tsv'


# Libraries whose metadata qc_flag marks them for exclusion (sources/sample_corrections.tsv in the metadata archive):
# failed Illumina libraries (low_yield), unidentifiable pool (pool_mismatch), and sampling/handling contamination
# (contamination: the day-18 2024-25 control RNA libraries and their metagenome). 'swap' rows are relabelled, not excluded.
EXCLUDE_FLAGS = ('low_yield', 'pool_mismatch', 'contamination')


def excluded(dataset='rna'):
    """Library ids (RNA sample ids, Illumina genohub ids, nanopore run_barcodes) excluded from all analyses."""
    c = corrections(dataset)
    return set(c.index[c.qc_flag.isin(EXCLUDE_FLAGS)])


def excluded_names(dataset='metagenome'):
    """Sample names (after relabelling) of QC-excluded libraries, for tables keyed by sample name (e.g. ice/popgen)."""
    c = corrections(dataset); c = c[c.qc_flag.isin(EXCLUDE_FLAGS)]
    return set(np.where(c.sample_name_corrected != '', c.sample_name_corrected, c.label_original))


def corrections(dataset):
    """Label corrections / QC flags (rna|metagenome|nanopore) from the cross-platform identity check, indexed by library id."""
    c = pd.read_csv(CORR, sep='\t', comment='#', dtype=str).fillna('')
    return c[c.dataset == dataset].set_index('id')


def meta():
    m = pd.read_csv(META, sep='\t', dtype=str).set_index('sample')
    m = m[~m.index.isin(excluded('rna'))]          # QC-excluded libraries never reach an analysis
    m['season'] = m['experiment'].str[-7:]
    tp = m['timepoint'].fillna('')
    m['day'] = [(-int(t[3:]) if t.startswith('PRE') else int(t[3:])) if t else np.nan for t in tp]
    med = m['env_medium'].fillna('')
    m['medium'] = np.select([med.eq('seawater'), med.eq('sea ice')], ['SW', 'SI'], 'other')
    trt = m['treatment'].fillna('')
    m['group'] = np.where(trt.isin(['CTL', 'EXP']), trt, 'other')
    m.loc[m['sample_name'].str.contains('SF-|sea foam', na=False), 'group'] = 'other'
    m['oiled'] = [(g == 'EXP') and d >= OIL_DAY[e] for g, d, e in zip(m['group'], m['day'], m['experiment'])]
    m['phase'] = np.select([m['day'] < 0, m['oiled'], (m['group'] == 'EXP')], ['pre', 'oiled', 'pre'], 'unoiled')
    m.loc[m['group'] == 'CTL', 'phase'] = np.where(m.loc[m['group'] == 'CTL', 'day'] < 0, 'pre', 'control')
    m.loc[m['group'] == 'other', 'phase'] = 'other'
    m['design'] = m['group'] + '_' + m['medium']
    return m


def table(name, index_cols=1):
    t = pd.read_csv(D + f'marine_coassembly.{name}.tsv', sep='\t', index_col=list(range(index_cols)))
    ex = excluded('rna')                            # QC-excluded libraries, whether samples are columns or rows
    return t.drop(columns=[c for c in t.columns if c in ex]).drop(index=[i for i in t.index if i in ex])


def clr(x, pseudo=0.5):
    l = np.log(x + pseudo)
    return l - l.mean(axis=0)


def all_samples():
    """Every RNA library in the count matrices, in column order (including QC-excluded ones)."""
    return pd.read_csv(D + 'marine_coassembly.samples.tsv', sep='\t')['sample'].tolist()


def genes():
    """Gene x sample counts restricted to samples() (columns in the same order), and the gene table."""
    import scipy.sparse as sp
    X = sp.load_npz(D + 'gene_counts.npz').tocsc()
    A = all_samples(); keep = [i for i, s in enumerate(A) if s not in excluded('rna')]
    if X.shape[1] != len(A): raise ValueError(f'gene_counts has {X.shape[1]} columns, samples.tsv {len(A)}')
    return X[:, keep].tocsr(), pd.read_pickle(D + 'genemap.pkl')


def samples():
    """RNA libraries used in analyses: all_samples() minus QC-excluded libraries."""
    ex = excluded('rna')
    return [s for s in all_samples() if s not in ex]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


_KJ = None
def kaiju_bin_label(b):
    """Majority Kaiju genus of a bin's contigs (used when GTDB-Tk has no call, e.g. *_sub bins); '' if no clear majority."""
    global _KJ
    if _KJ is None:
        _KJ = pd.read_csv(D + 'bin_kaiju_genus.tsv', sep='\t', index_col=0)['kaiju_genus']
    v = _KJ.get(b if b.startswith('dastool-') else 'dastool-' + b, '')
    if not isinstance(v, str) or not v: return ''
    pct = float(v.split('(')[-1].split('%')[0])
    return v.split(' (')[0] + ' (Kaiju)' if pct >= 50 else ''


_PL = None
def plastid_label(b):
    """'Plastid (host)' for bins that are algal chloroplast genomes (GTDB Gloeomargaritales, or <600 kb with a chloroplast 16S
    in SILVA); '' otherwise. Host = SILVA best hit of the plastid 16S where one was recovered (data/plastid_bins.tsv)."""
    global _PL
    if _PL is None:
        _PL = pd.read_csv(D + 'plastid_bins.tsv', sep='\t').set_index('bin')
    b = b if b.startswith('dastool-') else 'dastool-' + b
    if b not in _PL.index: return ''
    x = _PL.loc[b]
    if not (x.gtdb_order or (x['size'] < 600000 and x.n_plastid_ssu > 0)): return ''
    h = x.host if isinstance(x.host, str) and x.host and not x.host.startswith(('uncultured', 'environmental')) else 'unidentified alga'
    return f'Plastid ({h})'


def nanopore_samples():
    """dna/nanopore_samples.tsv with sample_name corrected by the cross-platform identity check (original kept in label_original)."""
    ns = pd.read_csv('/data/dev/rnaseq-exploration/dna/nanopore_samples.tsv', sep='\t')
    fix = ns.run_barcode.map(corrections('nanopore').sample_name_corrected).fillna('')
    ns['label_original'] = ns.sample_name
    ns['sample_name'] = ns.sample_name.where(fix == '', fix)
    return ns


_BD = None
def bin_domain(b):
    """'eukaryote' | 'organelle' | 'prokaryote' for a bin, from data/bin_domain.tsv (scripts/bin_domains.py: Tiara length-weighted
    fractions; eukarya >= 0.5 -> eukaryote, organelle >= 0.5 -> organelle). Unknown bins and '' count as prokaryote."""
    global _BD
    if _BD is None:
        _BD = pd.read_csv(D + 'bin_domain.tsv', sep='\t', index_col=0)['call']
    if not isinstance(b, str) or not b: return 'prokaryote'
    return _BD.get(b if b.startswith('dastool-') else 'dastool-' + b, 'prokaryote')


def is_prokaryote(b):
    return bin_domain(b) == 'prokaryote'


_EG = _EGI = _BEG = None
def _euk_groups():
    """Eukaryotic genome groups (scripts/bin_domains.py copies them from euks/mags): contig -> E## and group info."""
    global _EG, _EGI, _BEG
    if _EG is None:
        c = pd.read_csv(D + 'contig_euk_group.tsv.gz', sep='\t')
        _EG = dict(zip(c.contig, c.euk_group))
        _EGI = pd.read_csv(D + 'euk_groups.tsv', sep='\t', keep_default_na=False).set_index('group')
        _BEG = pd.read_csv(D + 'bin_domain.tsv', sep='\t', index_col=0).get('euk_group', pd.Series(dtype=str)).dropna().to_dict()
    return _EG, _EGI, _BEG


def euk_group(contig):
    """E## eukaryotic genome group of a contig, '' if none."""
    return _euk_groups()[0].get(contig, '')


def euk_group_info(gid):
    """Row of data/euk_groups.tsv for a group (status, network_ok, bp, purity, taxonomy_hint, ...) as a dict; {} if unknown."""
    I = _euk_groups()[1]
    return I.loc[gid].to_dict() if gid in I.index else {}


def is_euk_group(x):
    """True for an E## eukaryotic genome group id."""
    return isinstance(x, str) and x in _euk_groups()[1].index


def euk_group_label(gid):
    """'Eukaryote E02 (Ochromonadales)', 'Eukaryote E05', or 'Eukaryote E21 (mixed group)' for doubtful groups."""
    i = euk_group_info(gid)
    if not i: return ''
    h = i.get('taxonomy_hint', '') or ''
    tag = h if h else ('mixed group' if i.get('status') == 'doubtful' else '')
    return f'Eukaryote {gid}' + (f' ({tag})' if tag else '')


_VG = _VGI = None
def _virus_groups():
    """Virus groups (V##; euks/mags groups >= 50% geNomad-viral, copied by scripts/bin_domains.py): contig -> V## and info."""
    global _VG, _VGI
    if _VG is None:
        c = pd.read_csv(D + 'contig_virus_group.tsv.gz', sep='\t')
        _VG = dict(zip(c.contig, c.vid))
        _VGI = pd.read_csv(D + 'virus_groups.tsv', sep='\t', keep_default_na=False).set_index('vid')
    return _VG, _VGI


def virus_group(contig):
    """V## virus group of a contig, '' if none."""
    return _virus_groups()[0].get(contig, '')


def virus_group_info(vid):
    I = _virus_groups()[1]
    return I.loc[vid].to_dict() if vid in I.index else {}


def is_virus_group(x):
    return isinstance(x, str) and x in _virus_groups()[1].index


def virus_group_label(vid):
    """'Virus V01' (provenance only: V## groups are not used as labels on the site; see the viral genome-group rule in _hull_mags)."""
    return f'Virus {vid}' if is_virus_group(vid) else ''


EM = '/data/dev/rnaseq-exploration/euks/mags/out/'
def group_expression():
    """Groups (E## and V##) x RNA libraries: reads on the group's contigs per million contig-assigned reads.
    From euks/mags/out/group_expr_samples.tsv (expr.py), with ids renamed through virus_ids.tsv (e.g. E09 -> V01); groups
    registered from kernel blobs (E25-E30, V02) take the blob's abundance (da_blobs_abundance_rna.tsv, same metric) via
    blob_registry.tsv. Columns follow the libraries in those files (QC-excluded libraries already removed)."""
    import re
    E = pd.read_csv(EM + 'group_expr_samples.tsv', sep='\t', index_col=0)
    ren = pd.read_csv(EM + 'virus_ids.tsv', sep='\t').set_index('from_group').vid.to_dict()
    E = E.rename(index=ren)
    B = pd.read_csv(EM + 'da_blobs_abundance_rna.tsv', sep='\t', index_col=0)
    for x in pd.read_csv(EM + 'blob_registry.tsv', sep='\t').itertuples():
        m = re.search(r'new (?:eukaryote|virus) ([EV]\d\d)', x.decision)
        if m and m.group(1) not in E.index and x.blob in B.index: E.loc[m.group(1)] = B.loc[x.blob].reindex(E.columns)
    return E


def group_da():
    """Per-group differential abundance (RNA and DNA) from euks/mags/out/da_blobs.tsv (da.py), mapped from blobs to groups
    through blob_registry.tsv (new groups, and blobs that are 'already in' a group). Long table: group, data, contrast,
    log2_effect, q."""
    import re
    R = pd.read_csv(EM + 'blob_registry.tsv', sep='\t')
    ren = pd.read_csv(EM + 'virus_ids.tsv', sep='\t').set_index('from_group').vid.to_dict()
    mp = {}
    for x in R.itertuples():
        m = re.search(r'(?:new (?:eukaryote|virus)|already in) ([EV]\d\d)', x.decision)
        if m: mp[x.blob] = ren.get(m.group(1), m.group(1))
    A = pd.read_csv(EM + 'da_blobs.tsv', sep='\t')
    A = A[A.blob.isin(mp)].assign(group=lambda d: d.blob.map(mp))
    return A[['group', 'blob', 'data', 'contrast', 'log2_effect', 'q']]


_HM = _HMI = _HMT = None
def _hull_mags():
    """Hull-curated MAGs (EE## eukaryotic, VV## viral; euks/mags hull_mags.py, copied by scripts/bin_domains.py): contig -> MAG,
    MAG info (data/hull_mags.tsv), and the contigs Tiara calls eukaryotic.
    Viral genome-group rule (VV## and V##): >= 50% of bp in geNomad Nucleocytoviricota or another giant-virus / eukaryotic-virus
    class. Tailed-phage (Caudoviricetes) clusters are many phages, not genomes; bin_domains.py leaves them out of data/, so the
    former VV01/V01 (97% Caudoviricetes) has no genome label and its contigs keep the per-contig vOTU handling. No VV## exist now."""
    global _HM, _HMI, _HMT
    if _HM is None:
        c = pd.read_csv(D + 'contig_hull_mag.tsv.gz', sep='\t')
        _HM = dict(zip(c.contig, c.mag)); _HMT = set(c.contig[c.tiara == 'eukarya'])
        _HMI = pd.read_csv(D + 'hull_mags.tsv', sep='\t', keep_default_na=False).set_index('mag')
    return _HM, _HMI, _HMT


def hull_mag(contig):
    """EE##/VV## hull MAG of a contig, '' if none."""
    return _hull_mags()[0].get(contig, '')


def hull_mag_info(mid):
    """Row of data/hull_mags.tsv (source_group, status, network_ok, caution, bp, bp_gained, bp_lost, ...) as a dict; {} if unknown."""
    I = _hull_mags()[1]
    return I.loc[mid].to_dict() if mid in I.index else {}


def is_hull_mag(x):
    """True for an EE## or VV## hull MAG id."""
    return isinstance(x, str) and x in _hull_mags()[1].index


def is_euk_hull_mag(x):
    return is_hull_mag(x) and x.startswith('EE')


def is_virus_hull_mag(x):
    return is_hull_mag(x) and x.startswith('VV')


def hull_mag_label(mid):
    """'Virus VV##', 'Eukaryote EE05 (Minidiscus)', 'Eukaryote EE04', or 'Eukaryote EE08 (Oligohymenophorea, mixed group)'.
    The name is the identify.py taxonomy call and is shown only at high or medium confidence; doubtful MAGs say 'mixed group'."""
    i = hull_mag_info(mid)
    if not i: return ''
    if mid.startswith('VV'): return f'Virus {mid}'
    tax = i.get('taxonomy', '') if i.get('taxonomy_confidence', '') in ('high', 'medium') else ''
    tag = ', '.join(x for x in (tax, 'mixed group' if i.get('status') == 'doubtful' else '') if x)
    return f'Eukaryote {mid}' + (f' ({tag})' if tag else '')


def hull_mag_expression():
    """EE##/P## (and the retired VV01) x RNA libraries: reads on the MAG's contigs per million contig-assigned reads (euks/mags hull_expr.py; the
    same metric as group_expression() and hull_mags.tsv expr_ppm_median). QC-excluded libraries already removed."""
    return pd.read_csv(EM + 'hull_mags_expr_samples.tsv', sep='\t', index_col=0)


def hull_mag_da():
    """Per-MAG differential abundance (RNA and DNA) of the EE##/VV## hull MAGs (euks/mags hull_mags.py via dacore, the same tests as
    group_da()). Long table: group, data, contrast, log2_effect, q."""
    A = pd.read_csv(EM + 'hull_mags_da.tsv', sep='\t').rename(columns={'blob': 'group'})
    A = A[A.group.map(is_hull_mag)]                                      # wired MAGs only (drops the retired VV01)
    return A[['group', 'data', 'contrast', 'log2_effect', 'q']]


_OG = _OGI = None
def _organelle_groups():
    """Organelle groups (P## plastid, M## mitochondrial; >= 50% Tiara organelle, euks/mags register_blobs.py): contig -> id, info."""
    global _OG, _OGI
    if _OG is None:
        c = pd.read_csv(D + 'contig_organelle_group.tsv.gz', sep='\t')
        _OG = dict(zip(c.contig, c.oid))
        _OGI = pd.read_csv(D + 'organelle_groups.tsv', sep='\t', keep_default_na=False).set_index('oid')
    return _OG, _OGI


def organelle_group(contig):
    """P##/M## organelle group of a contig, '' if none."""
    return _organelle_groups()[0].get(contig, '')


def organelle_group_info(oid):
    I = _organelle_groups()[1]
    return I.loc[oid].to_dict() if oid in I.index else {}


def is_organelle_group(x):
    return isinstance(x, str) and x in _organelle_groups()[1].index


def organelle_group_label(oid):
    """'Plastid P01' or 'Mitochondrion M01'."""
    if not is_organelle_group(oid): return ''
    return ('Mitochondrion ' if oid.startswith('M') else 'Plastid ') + oid


def contig_is_plastid(contig):
    """True for a contig in a P## plastid group (the contig-level counterpart of is_plastid for unbinned contigs)."""
    return organelle_group(contig).startswith('P')


def contig_host_label(contig, fallback):
    """Host label for an UNBINNED contig, in this order: 'Virus VV##' if it is in a viral hull MAG (a virus, not a host; none at
    present, see _hull_mags);
    'Plastid P01' if it is in a P## organelle group (plastid DNA, never a bacterial host); its eukaryotic hull MAG (EE##) if it
    has one; else `fallback` (e.g. a Kaiju label). Contigs of doubtful (mixed) EE MAGs only take the label if Tiara calls the
    contig itself eukaryotic. Binned contigs are labelled through their bin. E##/V## groups are provenance only."""
    h = hull_mag(contig)
    if h.startswith('VV'): return hull_mag_label(h)
    o = organelle_group(contig)
    if o: return organelle_group_label(o)
    if not h: return fallback
    if hull_mag_info(h).get('status') == 'doubtful' and contig not in _hull_mags()[2]: return fallback
    return hull_mag_label(h)


_BHM = None
def _bin_hull_mag():
    global _BHM
    if _BHM is None:
        _BHM = pd.read_csv(D + 'bin_domain.tsv', sep='\t', index_col=0).get('hull_mag', pd.Series(dtype=str)).dropna().to_dict()
    return _BHM


def special_label(b):
    """Label that overrides GTDB/Kaiju for non-prokaryote bins: the plastid label (plastid_label) if there is one, else
    'Plastid/organelle (Tiara)' for organelle bins; for eukaryotic nuclear bins the label of the eukaryotic genome group
    holding most of the bin (hull_mag_label, e.g. comebin_3937 -> 'Eukaryote EE05'), or 'Eukaryote (Tiara)'; '' for prokaryotes."""
    if not isinstance(b, str) or not b: return ''
    p = plastid_label(b)
    if p: return p
    d = bin_domain(b)
    if d == 'organelle': return 'Plastid/organelle (Tiara)'
    if d == 'eukaryote':
        k = b if b.startswith('dastool-') else 'dastool-' + b
        g = _bin_hull_mag().get(k, '')
        return hull_mag_label(g) if g else 'Eukaryote (Tiara)'
    return ''


def is_plastid(b):
    """Plastid/organelle genome: plastid_label rule (plastid_bins.tsv) or Tiara organelle call (bin_domain)."""
    return isinstance(b, str) and bool(b) and (bool(plastid_label(b)) or bin_domain(b) == 'organelle')
