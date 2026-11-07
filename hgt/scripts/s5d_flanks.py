"""Flank analysis for HGT candidates: size of the shared >=98% nt block around the gene, background similarity outside it,
and mobile genes inside/at the edges of the block."""
import sys, subprocess, re, os; sys.path.insert(0, '.')
from hcdef import *
H = '/data/dev/rnaseq-exploration/hgt/'
c = pd.read_csv(H + 'hgt_candidates.tsv', sep='\t')
reg = pd.read_csv(H + 'cand_regions.tsv', sep='\t', header=None, names=['region', 'gene']).set_index('gene').region
seq, cur = {}, None
for line in open(H + 'grex/cand_regions.fna'):
    if line[0] == '>': cur = line[1:].split()[0]; seq[cur] = []
    else: seq[cur].append(line.strip())
seq = {k: ''.join(v) for k, v in seq.items()}
cg = pd.read_csv(H + 'grex/hc_contig_genes.tsv.gz', sep='\t', keep_default_na=False)
rx_is = re.compile(r'transposase|insertion element|\bIS\d|IS[0-9]+ family|integrase|resolvase|invertase|recombinase', re.I)
def genes_in(contig, a, b):
    x = cg[(cg.contig == contig) & (cg.end >= a) & (cg.start <= b)]
    return x
out = []
for _, r in c.drop_duplicates(['q_contig', 't_contig']).iterrows():
    rq, rt = reg[r.q], reg[r.t]
    with open(H + 'tq.fa', 'w') as f: f.write(f'>q\n{seq[rq]}\n')
    with open(H + 'tt.fa', 'w') as f: f.write(f'>t\n{seq[rt]}\n')
    o = subprocess.run(['blastn', '-query', H + 'tq.fa', '-subject', H + 'tt.fa', '-outfmt', '6 pident length qstart qend sstart send', '-evalue', '1e-10'],
                       capture_output=True, text=True).stdout.strip()
    hs = pd.DataFrame([l.split() for l in o.split('\n') if l], columns=['pid', 'len', 'qs', 'qe', 'ss', 'se']).astype(float) if o else pd.DataFrame(columns=['pid', 'len', 'qs', 'qe', 'ss', 'se'])
    q0 = int(rq.split(':')[1].split('-')[0]); qlen = len(seq[rq]); tlen = len(seq[rt])
    hi = hs[hs.pid >= 98]
    cov = np.zeros(qlen, bool)
    for _, x in hi.iterrows(): cov[int(x.qs) - 1:int(x.qe)] = True
    # block = contiguous covered run containing the gene (gene sits at ~20 kb offset)
    gs = cg.set_index('gene_id').loc[r.q]; g_off = int(gs.start) - q0
    i0 = j0 = min(max(g_off, 0), qlen - 1)
    if cov[i0]:
        while i0 > 0 and cov[i0 - 1:i0].all() or (i0 > 0 and cov[max(0, i0 - 200):i0].any()): i0 -= 1
        while j0 < qlen - 1 and (cov[j0 + 1] or cov[j0 + 1:j0 + 200].any()): j0 += 1
    block = j0 - i0 + 1 if cov[min(max(g_off, 0), qlen - 1)] else 0
    lo = hs[(hs.pid < 98)]
    bg_len = int(lo.len.sum()); bg_pid = float((lo.pid * lo.len).sum() / lo.len.sum()) if len(lo) else np.nan
    inb = genes_in(r.q_contig, q0 + i0, q0 + j0)
    edge = genes_in(r.q_contig, q0 + i0 - 3000, q0 + j0 + 3000)
    out.append(dict(q=r.q, t=r.t, q_contig=r.q_contig, t_contig=r.t_contig, level=r.level, category=r.q_category,
                    q_host=f'{r.q_genus or "?"} ({r.q_family or "?"})', t_host=f'{r.t_genus or "?"} ({r.t_family or "?"})', nt_ident=r.nt_ident, depth_rho=r.depth_rho,
                    region_q=qlen, region_t=tlen, shared_block_bp=block, frac_region_shared=round(cov.mean(), 3), background_aln_bp=bg_len,
                    background_pid=round(bg_pid, 1) if bg_pid == bg_pid else np.nan, block_genes=len(inb),
                    block_IS=int(inb['product'].str.contains(rx_is).sum()), edge_IS=int(edge['product'].str.contains(rx_is).sum()),
                    block_products='; '.join(inb['product'].head(25))))
for f in ('tq.fa', 'tt.fa'): os.remove(H + f)
F = pd.DataFrame(out)
# interpretation
def call(x):
    if x.frac_region_shared > 0.9 and x.region_q > 20000 and x.region_t > 20000: return 'same lineage / redundant contig (whole window identical)'
    if x.shared_block_bp and x.background_aln_bp < 2000 and x.frac_region_shared < 0.6: return 'HGT: identical block in unrelated backgrounds'
    if x.shared_block_bp and x.background_pid == x.background_pid and x.background_pid < 90: return 'HGT-like: identical block, divergent background'
    return 'unresolved (short contig or partial alignment)'
F['interpretation'] = F.apply(call, axis=1)
F.to_csv(H + 'hgt_candidates_flanks.tsv', sep='\t', index=False)
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 70)
print(F.interpretation.value_counts().to_string())
print(F[['category', 'q_host', 't_host', 'nt_ident', 'depth_rho', 'region_q', 'region_t', 'shared_block_bp', 'frac_region_shared', 'background_aln_bp', 'background_pid', 'block_genes', 'block_IS', 'edge_IS', 'interpretation']].round(2).to_string())
