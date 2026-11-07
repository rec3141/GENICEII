"""Shared loaders for the CMO RNA-seq exploration (marine coassembly expression + sample metadata)."""
import re
import numpy as np, pandas as pd

D = '/data/dev/rnaseq-exploration/data/'
META = '/data/dev/genice-sample-metadata/outputs/cmo_rna_metadata.tsv'
# oil addition: 2024-25 diesel on DAY05 16:30 (after the DAY05 11:00 samples); 2025-26 136 L oil on DAY01
OIL_DAY = {'CMO 2024-25': 5.5, 'CMO 2025-26': 1}
CORR = '/data/dev/genice-sample-metadata/sources/sample_corrections.tsv'


def corrections(dataset):
    """Label corrections / QC flags (rna|metagenome|nanopore) from the cross-platform identity check, indexed by library id."""
    c = pd.read_csv(CORR, sep='\t', comment='#', dtype=str).fillna('')
    return c[c.dataset == dataset].set_index('id')


def meta():
    m = pd.read_csv(META, sep='\t', dtype=str).set_index('sample')
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
    return pd.read_csv(D + f'marine_coassembly.{name}.tsv', sep='\t', index_col=list(range(index_cols)))


def clr(x, pseudo=0.5):
    l = np.log(x + pseudo)
    return l - l.mean(axis=0)


def genes():
    import scipy.sparse as sp
    return sp.load_npz(D + 'gene_counts.npz').tocsr(), pd.read_pickle(D + 'genemap.pkl')


def samples():
    return pd.read_csv(D + 'marine_coassembly.samples.tsv', sep='\t')['sample'].tolist()


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
