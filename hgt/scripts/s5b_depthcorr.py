"""Add DNA depth-profile correlation (276 nanopore samples) + gene nucleotide identity to the shared pairs."""
import sys, subprocess, os; sys.path.insert(0, '.')
from hcdef import *
from scipy.stats import spearmanr
H = '/data/dev/rnaseq-exploration/hgt/'
a = pd.read_csv(H + 'shared_pairs.tsv', sep='\t', keep_default_na=False)
dep = pd.read_csv(H + 'grex/hc_contig_depths.tsv.gz', sep='\t', index_col=0)
dep = dep[[c for c in dep.columns if c.endswith('.bam')]]
L = np.log1p(dep)
ok = a.q_contig.isin(L.index) & a.t_contig.isin(L.index)
rho = np.full(len(a), np.nan); ratio = rho.copy()
for i in np.where(ok)[0]:
    x, y = L.loc[a.q_contig.iat[i]].values, L.loc[a.t_contig.iat[i]].values
    s = (x > 0) | (y > 0)
    if s.sum() >= 10: rho[i] = spearmanr(x[s], y[s]).statistic
    ratio[i] = (dep.loc[a.q_contig.iat[i]].sum() + 1) / (dep.loc[a.t_contig.iat[i]].sum() + 1)
a['depth_rho'] = rho; a['depth_ratio'] = ratio
# gene nucleotide identity via blastn on the extracted gene regions
reg = pd.read_csv(H + 'grex/regions.tsv', sep='\t', header=None, names=['region', 'gene'])
r2g = dict(zip(reg.region, reg.gene))
seqs, cur = {}, None
for line in open(H + 'grex/hc_regions.fna'):
    if line[0] == '>': cur = r2g.get(line[1:].split()[0]); seqs[cur] = []
    else: seqs[cur].append(line.strip())
need = set(a.q) | set(a.t)
with open(H + 'tmp_pair_genes.fna', 'w') as f:
    for gname in need:
        if gname in seqs: f.write(f'>{gname}\n{"".join(seqs[gname])}\n')
out = subprocess.run(['blastn', '-query', H + 'tmp_pair_genes.fna', '-subject', H + 'tmp_pair_genes.fna', '-outfmt', '6 qseqid sseqid pident length qlen slen',
                      '-max_target_seqs', '100000', '-evalue', '1e-20'], capture_output=True, text=True).stdout
b = pd.DataFrame([l.split('\t') for l in out.strip().split('\n')], columns=['q', 't', 'pident', 'len', 'qlen', 'slen'])
b[['pident', 'len', 'qlen', 'slen']] = b[['pident', 'len', 'qlen', 'slen']].astype(float)
b = b[b.len >= 0.8 * b[['qlen', 'slen']].min(axis=1)].sort_values('len', ascending=False).drop_duplicates(['q', 't'])
nt = {(q, t): p for q, t, p in zip(b.q, b.t, b.pident)}
a['nt_ident'] = [nt.get((q, t), nt.get((t, q), np.nan)) for q, t in zip(a.q, a.t)]
os.remove(H + 'tmp_pair_genes.fna')
a.to_csv(H + 'shared_pairs.tsv', sep='\t', index=False)
x = a[a.differ_at.isin(['phylum', 'class', 'order', 'family'])]
print('cross-family pairs', len(x), 'with depth rho', x.depth_rho.notna().sum())
print(x.depth_rho.describe().round(2).to_string())
print('rho<0.3 (independent hosts):', (x.depth_rho < 0.3).sum(), ' rho>0.7 (same genome likely):', (x.depth_rho > 0.7).sum())
print('nt identity cross-family:', x.nt_ident.describe().round(2).to_string())
y = a[a.differ_at == '']
print('same/unknown-taxon pairs rho:', y.depth_rho.describe().round(2).to_string())
