"""Candidate hosts for the diesel-responsive vOTUs: Spearman correlation across all 130 RNA libraries between each vOTU's transcripts
and each MAG's TPM (data bin_tpm), top 3 MAGs; plus direct links (vOTU binned in a MAG). Correlation is only a hint at hosts.
Out: out/votu_hosts.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr
V = '/data/dev/rnaseq-exploration/viruses/'
U = pd.read_pickle(V + 'out/votu_rna_cpm.pkl'); R = pd.read_csv(V + 'out/votu_rna.tsv', sep='\t', index_col=0)
B = table('bin_tpm')[U.columns]; B = B[B.max(1) >= 500]
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).set_index('user_genome')['classification']
def lab(b):
    if special_label(b): return special_label(b)
    t = [x.split('__', 1)[1] for x in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in x and x.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b
sel = set(R[(R['oil_2025-26'] >= 20)].sort_values('baci_2025-26', ascending=False).head(15).index) | \
      set(R[(R['oil_2024-25'] >= 20)].sort_values('lfc_2024-25', ascending=False).head(15).index) | set(R.sort_values('mean_cpm', ascending=False).head(15).index)
Br = B.rank(axis=1).values; Bz = (Br - Br.mean(1, keepdims=True)) / Br.std(1, keepdims=True)
rows = []
for v in sel:
    r = U.loc[v].rank().values; z = (r - r.mean()) / (r.std() or 1)
    rho = Bz @ z / len(z)
    o = np.argsort(-rho)[:3]
    rows.append(dict(votu=v, binned=R.loc[v, 'bin'] if isinstance(R.loc[v, 'bin'], str) else '', **{f'host{i+1}': f'{lab(B.index[j])} ({B.index[j].replace("dastool-", "")})' for i, j in enumerate(o)},
                     **{f'rho{i+1}': round(float(rho[j]), 2) for i, j in enumerate(o)}))
H = pd.DataFrame(rows).set_index('votu').join(R[['cls', 'length', 'checkv_quality', 'oil_2024-25', 'ctl_2024-25', 'oil_2025-26', 'ctl_2025-26', 'baci_2025-26']])
H.to_csv(V + 'out/votu_hosts.tsv', sep='\t')
pd.set_option('display.width', 280); pd.set_option('display.max_colwidth', 45)
print(H.sort_values('oil_2024-25', ascending=False)[['length', 'oil_2024-25', 'ctl_2024-25', 'oil_2025-26', 'ctl_2025-26', 'host1', 'rho1', 'host2', 'rho2']].round(1).to_string())
