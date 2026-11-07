"""Assemble site/data.json for the CMO diesel website from out/ tables."""
import sys, json, os; sys.path.insert(0, '.')
from common import *
O = '../out/'; SITE = '../site/'
m = meta(); s = table('sample_summary'); S = list(s.index); m = m.loc[S]
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
q = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin')
def r(x, n=3):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
def tax(b):
    t = gt.get(b.replace('dastool-', ''), '')
    d = {k: v for k, v in (x.split('__', 1) for x in t.split(';') if '__' in x)}
    return d
def label(b):
    if plastid_label(b): return plastid_label(b)
    d = tax(b)
    for k in ('s', 'g', 'f', 'o', 'c', 'p'):
        if d.get(k) and (k in 'sg' or not kaiju_bin_label(b)): return (d[k] if k in 'sg' else f"{d[k]} ({ {'f':'family','o':'order','c':'class','p':'phylum'}[k]})")
    return kaiju_bin_label(b) or 'unclassified'
pairs = s.reads_raw / 2
samples = []
for i in S:
    x = m.loc[i]
    samples.append({'id': i, 'name': x.sample_name, 'season': x.season, 'pool': x.group, 'medium': x.medium,
                    'env': x.env_medium if isinstance(x.env_medium, str) else '', 'day': r(x.day, 0), 'phase': x.phase,
                    'date': x.collection_date if isinstance(x.collection_date, str) else '',
                    'qc': {'rrna_filtered': r((s.reads_raw[i] - s.reads_norrna[i]) / s.reads_raw[i] * 100, 1),
                           'rrna_escaped': r(s.pairs_rrna_features[i] / pairs[i] * 100, 1), 'cds': r(s.pairs_cds[i] / pairs[i] * 100, 1),
                           'other_mapped': r(max(0, s.reads_mapped[i] / 2 - s.pairs_cds[i] - s.pairs_rrna_features[i]) / pairs[i] * 100, 1),
                           'unmapped': r((s.reads_nohuman[i] - s.reads_mapped[i]) / s.reads_raw[i] * 100, 1)}})
pca = pd.read_csv(O + 'pca_mags.tsv', sep='\t', index_col=0)
for d in samples: d['pc'] = [r(pca.loc[d['id'], 'PC1'], 2), r(pca.loc[d['id'], 'PC2'], 2)]
# MAG responses
resp = pd.read_csv(O + 'mag_oil_response.tsv', sep='\t')
baci = pd.read_csv(O + 'mag_baci.tsv', sep='\t')
tpm = table('bin_tpm'); bc = table('bin_counts')
dna = pd.read_csv(O + 'dna_rel_abundance.tsv', sep='\t', index_col=0); dinfo = pd.read_csv(O + 'dna_samples.tsv', sep='\t', index_col=0)
sel = set(resp[(resp.q < 0.05) & (resp.dCLR.abs() > 2.5)].bin) | set(baci[(baci.q < 0.1) & (baci.interaction.abs() > 2)].bin)
extra = json.load(open(O + 'pathways/key_bins.json')) if os.path.exists(O + 'pathways/key_bins.json') else []
sel |= set(extra)
sel = [b for b in sel if b in tpm.index]
mags = {}
for b in sel:
    d = tax(b); rr = resp[resp.bin == b]; bb = baci[baci.bin == b]
    mags[b.replace('dastool-', '')] = {
        'label': label(b), 'tax': {k: d.get(k, '') for k in 'pcofgs'}, 'size': int(q.loc[b, 'size']) if b in q.index else None,
        'completeness': r(q.loc[b, 'SCG_completeness'], 0) if b in q.index else None, 'redundancy': r(q.loc[b, 'SCG_redundancy'], 0) if b in q.index else None,
        'resp': {f"{a.season} {a.medium}": {'dclr': r(a.dCLR, 2), 'q': r(a.q, 4), 'oiled': r(a.tpm_oiled, 1), 'ctl': r(a.tpm_ctl, 1)} for a in rr.itertuples()},
        'baci': {a.season: {'i': r(a.interaction, 2), 'q': r(a.q, 3)} for a in bb.itertuples()},
        'tpm': [r(v, 1) for v in tpm.loc[b, S]],
        'dna': [r(v, 3) for v in dna.loc[b]] if b in dna.index else None}
dna_samples = [{'key': k, 'name': x.sample_name, 'season': x.season, 'pool': x.group, 'medium': x.medium, 'day': r(x.day, 0)} for k, x in dinfo.loc[dna.columns].iterrows()]
data = {'samples': samples, 'mags': mags, 'dna_samples': dna_samples,
        'unbinned_share': [r(v, 1) for v in (bc.loc['unbinned', S] / bc[S].sum() * 100)]}
for f in ['pathways/site_pathways.json']:
    if os.path.exists(O + f): data['pathways'] = json.load(open(O + f))
def fam(b):
    if b == 'unbinned': return 'unbinned'
    t = gt.get(b.replace('dastool-', ''), ''); d = dict(x.split('__', 1) for x in t.split(';') if '__' in x)
    if 'Cycloclasticus' in kaiju_bin_label(b): return 'Cycloclasticaceae'
    return d.get('f') or 'other'
KEEP = ['Cycloclasticaceae', 'Rhodobacteraceae', 'Sulfurimonadaceae', 'Porticoccaceae', 'Alteromonadaceae', 'Pseudohongiellaceae', 'Flavobacteriaceae']
fb = bc.drop(index='unbinned').groupby(lambda b: fam(b) if fam(b) in KEEP else 'other').sum()[S]
fb = fb / fb.sum() * 100
data['families'] = {'order': KEEP + ['other'], 'pct': {k: [r(v, 2) for v in fb.loc[k]] for k in KEEP + ['other']}}
kc = pd.read_csv(D + 'ko_cpm.tsv', sep='\t', index_col=0)
MK = [('porA', 'K00169', 'Pyruvate:ferredoxin oxidoreductase'), ('hydA', 'K05927', 'NiFe hydrogenase, small subunit'), ('aclA', 'K15230', 'ATP-citrate lyase (reverse TCA)'),
      ('fccB', 'K17229', 'Sulfide dehydrogenase (flavocytochrome c)'), ('soxB', 'K17224', 'Thiosulfate oxidation SoxB'), ('dsrA', 'K11180', 'Dissimilatory sulfite reductase (sulfate reduction)'),
      ('cydA', 'K00425', 'Cytochrome bd oxidase (high O2 affinity)'), ('ccoN', 'K00404', 'cbb3 oxidase (high O2 affinity)'), ('napA', 'K02567', 'Periplasmic nitrate reductase'), ('amoA/pmoA', 'K10944', 'Cu membrane monooxygenase A')]
grp = (m.season + '|' + m.group + '|' + m.medium + '|' + m.phase)
data['markers'] = [{'id': a, 'ko': k, 'name': n, 'groups': {g: r(kc.loc[k, idx].mean(), 1) for g, idx in grp.groupby(grp).groups.items()
                    if k in kc.index and g.split('|')[1] in ('CTL', 'EXP') and g.split('|')[2] in ('SW', 'SI')}} for a, k, n in MK]
eu = pd.read_csv(O + 'euk_rrna_response.tsv', sep='\t').fillna({'taxon': 'unassigned eukaryotes'})
eu = eu[eu.taxon.str.strip() != '']
data['euk'] = [{'taxon': t, **{f'{x.season}': {'oil': r(x.oiled, 2), 'ctl': r(x.ctl, 2), 'lfc': r(x.lfc, 2)} for x in d.itertuples()}} for t, d in eu.groupby('taxon')]
# cross-platform identity check (swapcheck/identity.py): original label, corrected label / QC flag, diesel-pool signature and control markers
idt = pd.read_csv('../swapcheck/identity.tsv', sep='\t', dtype={'corrected': str, 'flag': str}).fillna({'corrected': '', 'flag': ''})
data['identity'] = [{'platform': x.platform, 'library': x.library, 'label': x.label, 'corrected': x.corrected, 'flag': x.flag,
                     'diesel': r(x.diesel, 3), 'control': r(x.control, 3)} for x in idt.itertuples()]
if os.path.exists(O + 'mg_site.json'): data['mg'] = json.load(open(O + 'mg_site.json'))
from site_text import FINDINGS, CAVEATS, FOOTER
data.update(findings=FINDINGS, caveats=CAVEATS, footer=FOOTER)
json.dump(data, open(SITE + 'data.json', 'w'), separators=(',', ':'))
print('mags', len(mags), 'bytes', os.path.getsize(SITE + 'data.json'))
