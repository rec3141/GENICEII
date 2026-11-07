"""IBP mobility and sharing: MGE context of each IBP gene; IBP pairs shared across genera/domains at high identity with flank check
(blastn of +-20 kb windows: an identical block in divergent backgrounds = transfer); eukaryote IBPs vs their closest bacterial IBP.
-> ../ibp_context.tsv, ../ibp_pairs.tsv, ../ibp_euk_nearest.tsv"""
import sys, re, gzip, subprocess, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/cold/scripts')
from icecommon import *
H = '/data/dev/rnaseq-exploration/hgt/'; G = C_DIR + 'grex/'
I = pd.read_csv(C_DIR + 'ibp_genes.tsv', sep='\t', keep_default_na=False)
cg = pd.read_csv(G + 'ibp_contig_genes.tsv.gz', sep='\t', keep_default_na=False)
co = cg.set_index('gene_id')[['start', 'end', 'strand']]
I = I.join(co, on='gene_id')
# ---- MGE context (same rules as the HGT fork, hgt/scripts/s3_context.py)
p = pd.read_csv(H + 'mge/plasmid_summary.tsv', sep='\t')
v = pd.read_csv(H + 'mge/virus_summary.tsv', sep='\t'); v = v[~v.seq_name.str.contains(r'\|provirus')]
pv = pd.read_csv(H + 'mge/provirus.tsv', sep='\t'); gi = pd.read_csv(H + 'mge/genomic_islands.tsv', sep='\t'); it = pd.read_csv(H + 'mge/integrons.tsv', sep='\t')
I['plasmid_any'] = I.contig.isin(p.seq_name)
I['plasmid_strict'] = I.contig.isin(p[(p.plasmid_score >= 0.9) & ((p.n_hallmarks >= 1) | p.conjugation_genes.notna())].seq_name)
I['phage_contig'] = I.contig.isin(v[(v.n_hallmarks >= 2) & (v.virus_score >= 0.95)].seq_name)
def within(df, c, s, e, pad=0):
    grp = {k: d[[s, e]].values for k, d in df.groupby(c)}
    return [bool(((r[:, 0] - pad <= b) & (r[:, 1] + pad >= a)).any()) if (r := grp.get(ct)) is not None else False for ct, a, b in zip(I.contig, I.start, I.end)]
I['provirus'] = within(pv, 'source_seq', 'start', 'end'); I['genomic_island'] = within(gi, 'contig', 'start', 'end')
ib = it.groupby(['ID_integron', 'ID_replicon']).agg(start=('pos_beg', 'min'), end=('pos_end', 'max')).reset_index()
I['integron'] = within(ib, 'ID_replicon', 'start', 'end', pad=1000)
rx_is = re.compile(r'transposase|insertion element|\bIS\d|IS[0-9]+ family|integrase|recombinase|resolvase|invertase', re.I)
hk = cg['product'].str.contains(r'RecA|Xer[CD]|Holliday|Ruv[ABCX]|RecB|RecG|RecN|RecO|RecR|DNA repair', case=False, regex=True)
cg['is'] = cg['product'].str.contains(rx_is) & ~hk
isd = cg[cg['is']].groupby('contig')[['start', 'end']].apply(lambda d: d.values).to_dict()
dist = []
for ct, a, b in zip(I.contig, I.start, I.end):
    r = isd.get(ct); dist.append(np.nan if r is None else np.maximum(0, np.maximum(r[:, 0] - b, a - r[:, 1])).min())
I['IS_dist'] = dist; I['IS_within_10kb'] = I.IS_dist <= 10000
I['mge_borne'] = I[['plasmid_strict', 'phage_contig', 'provirus', 'integron', 'genomic_island']].any(axis=1)
I['mobile_assoc'] = I.mge_borne | I.IS_within_10kb
gcc = pd.read_csv(G + 'ibp_contig_gc.tsv', sep='\t'); gcc.columns = ['contig'] + list(gcc.columns[1:])
gcol = [c for c in gcc.columns if 'gc' in c.lower()][0]
I['contig_gc'] = I.contig.map(gcc.set_index('contig')[gcol]).astype(float)
clen = cg.groupby('contig').end.max(); I['contig_len_approx'] = I.contig.map(clen)
# neighbours (+-3 genes) for genomic context
nb = []
for _, r in I.iterrows():
    d = cg[cg.contig == r.contig].reset_index(drop=True); k = d.index[d.gene_id == r.gene_id]
    nb.append('; '.join(d.loc[max(0, k[0] - 3):k[0] + 3, 'product']) if len(k) else '')
I['neighbourhood'] = nb
I.to_csv(C_DIR + 'ibp_context.tsv', sep='\t', index=False)
cols = ['plasmid_any', 'plasmid_strict', 'phage_contig', 'provirus', 'genomic_island', 'integron', 'IS_within_10kb', 'mge_borne', 'mobile_assoc']
print('IBP genes with each MGE flag (%):'); print((I.groupby('domain')[cols].mean() * 100).round(1).to_string())
w = I.cpm_SI
print('share of ice IBP transcripts on mobile-associated genes: %.1f%%, MGE-borne %.1f%%' % (w[I.mobile_assoc].sum() / w.sum() * 100, w[I.mge_borne].sum() / w.sum() * 100))
# ---- pairs
A = pd.read_csv(G + 'ibp_ava.m8', sep='\t', header=None, names=['q', 't', 'fid', 'alnlen', 'qcov', 'tcov', 'evalue', 'bits'])
A = A[A.q < A.t]
info = I.set_index('gene_id')
for s in ['q', 't']:
    for c in ['contig', 'bin', 'genus', 'family', 'phylum', 'domain', 'cpm_SI']: A[f'{s}_{c}'] = A[s].map(info[c])
A = A[A.q_contig != A.t_contig]
gq = A.q_genus.str.replace(' (Kaiju)', '', regex=False); gt = A.t_genus.str.replace(' (Kaiju)', '', regex=False)
A['cross'] = np.select([(A.q_domain == 'eukaryote') != (A.t_domain == 'eukaryote'), A.q_phylum != A.t_phylum, gq != gt], ['cross-domain', 'cross-phylum', 'cross-genus'], 'same genus')
A.loc[(gq.isin(['unassigned', ''])) | (gt.isin(['unassigned', ''])), 'cross'] = A.loc[(gq.isin(['unassigned', ''])) | (gt.isin(['unassigned', ''])), 'cross'].replace({'cross-genus': 'unassigned', 'cross-phylum': 'unassigned'})
hi = A[(A.fid >= 0.9) & (A.qcov >= 0.8) & (A.tcov >= 0.8)].copy()
print('pairs >=90% aa / 80% cov on different contigs:', len(hi)); print(hi.cross.value_counts().to_string())
# flank check
seq, cur = {}, None
for line in open(G + 'ibp_regions.fna'):
    if line[0] == '>': cur = line[1:].split()[0]; seq[cur] = []
    else: seq[cur].append(line.strip())
seq = {k: ''.join(v) for k, v in seq.items()}
reg = pd.read_csv(G + 'ibp_regions.tsv', sep='\t', header=None, names=['region', 'gene']).set_index('gene').region
out = []
tmp = C_DIR + 'tmp_'
for _, r in hi[hi.cross != 'same genus'].iterrows():
    rq, rt = reg.get(r.q), reg.get(r.t)
    if rq not in seq or rt not in seq: continue
    open(tmp + 'q.fa', 'w').write(f'>q\n{seq[rq]}\n'); open(tmp + 't.fa', 'w').write(f'>t\n{seq[rt]}\n')
    o = subprocess.run(['blastn', '-query', tmp + 'q.fa', '-subject', tmp + 't.fa', '-outfmt', '6 pident length qstart qend sstart send', '-evalue', '1e-10'], capture_output=True, text=True).stdout.strip()
    hs = pd.DataFrame([l.split() for l in o.split('\n') if l], columns=['pid', 'len', 'qs', 'qe', 'ss', 'se']).astype(float) if o else pd.DataFrame(columns=['pid', 'len', 'qs', 'qe', 'ss', 'se'])
    qlen = len(seq[rq]); cov = np.zeros(qlen, bool)
    for _, x in hs[hs.pid >= 95].iterrows(): cov[int(x.qs) - 1:int(x.qe)] = True
    gene_nt = hs[(hs.pid >= 0)].sort_values('len', ascending=False).head(1)
    lo = hs[hs.pid < 95]
    out.append(dict(q=r.q, t=r.t, aa_ident=r.fid, cross=r.cross, q_host=f'{r.q_genus} [{r.q_domain}]', t_host=f'{r.t_genus} [{r.t_domain}]',
                    q_bin=r.q_bin, t_bin=r.t_bin, q_contig=r.q_contig, t_contig=r.t_contig, q_cpm_SI=r.q_cpm_SI, t_cpm_SI=r.t_cpm_SI,
                    top_nt_pid=float(gene_nt.pid.iloc[0]) if len(gene_nt) else np.nan, shared95_bp=int(cov.sum()), region_q=qlen, region_t=len(seq[rt]),
                    frac_region_shared=round(cov.mean(), 3), background_bp=int(lo.len.sum()), background_pid=round(float((lo.pid * lo.len).sum() / lo.len.sum()), 1) if len(lo) else np.nan))
for f in ('q.fa', 't.fa'):
    if os.path.exists(tmp + f): os.remove(tmp + f)
F = pd.DataFrame(out)
def call(x):
    short = min(x.region_q, x.region_t)
    if x.frac_region_shared > 0.9 and x.region_q > 20000 and x.region_t > 20000: return 'same lineage / redundant contig'
    if x.shared95_bp >= 0.8 * short: return 'contained (shorter window fully shared; same genome likely)'
    if x.shared95_bp < 500: return 'unresolved'
    if x.shared95_bp and x.background_bp < 2000 and x.frac_region_shared < 0.6: return 'HGT: shared block in unrelated backgrounds'
    if x.shared95_bp and x.background_pid == x.background_pid and x.background_pid < 90: return 'HGT-like: shared block, divergent background'
    return 'unresolved'
if len(F):
    F['interpretation'] = F.apply(call, axis=1)
F.to_csv(C_DIR + 'ibp_pairs.tsv', sep='\t', index=False)
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 40)
if len(F): print(F.interpretation.value_counts().to_string()); print(F.drop(columns=['q', 't']).round(2).to_string())
# eukaryote IBPs: nearest bacterial IBP
e = A[(A.cross == 'cross-domain')].copy()
e['euk'] = np.where(e.q_domain == 'eukaryote', e.q, e.t); e['bac'] = np.where(e.q_domain == 'eukaryote', e.t, e.q)
e['bac_host'] = np.where(e.q_domain == 'eukaryote', e.t_genus, e.q_genus)
best = e.sort_values('bits', ascending=False).drop_duplicates('euk')[['euk', 'bac', 'bac_host', 'fid', 'alnlen', 'qcov', 'tcov']]
euk = I[I.domain == 'eukaryote'][['gene_id', 'contig', 'plen', 'cpm_SI', 'cpm_SW']].merge(best, left_on='gene_id', right_on='euk', how='left')
euk.to_csv(C_DIR + 'ibp_euk_nearest.tsv', sep='\t', index=False)
print('eukaryote IBPs:', len(euk), ' with a bacterial IBP hit (>=30% aa, >=50% cov):', euk.bac.notna().sum())
print(euk.fid.describe().round(3).to_string())
print(euk.groupby('bac_host').size().sort_values().tail(10).to_string())
