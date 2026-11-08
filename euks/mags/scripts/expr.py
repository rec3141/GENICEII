"""Expression of the final eukaryotic genome groups (out/contig_euk_group.tsv) (full-depth run).
(a) contig-level read pairs (data/marine_coassembly.contig_counts.mtx.gz; rows = contigs in ascending id order) summed per group, as a
share of all contig-assigned reads per sample (per million) -- not limited to Bakta CDS, so intron-containing genes count;
(b) CDS reads from gene_counts via genemap (Bakta calls on eukaryotic contigs are poor; rough).
Out: out/group_expr_samples.tsv, out/group_expr_summary.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
import scipy.io, scipy.sparse as sp
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
S = samples(); m = meta().loc[S]
R = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t').rename(columns={'euk_group': 'group'})
C = scipy.io.mmread(D + 'marine_coassembly.contig_counts.mtx.gz').tocsr()      # contigs x all_samples()
ALL = all_samples(); keep = [ALL.index(x) for x in S]; C = C[:, keep]          # drop QC-excluded libraries
tot = np.asarray(C.sum(0)).ravel()
# contig ids are not contiguous; matrix rows follow the assembly's contig order = ids sorted ascending (1,549,905 contigs)
allids = np.sort(pd.read_pickle(O + 'contigs_full.pkl').id.values); pos = pd.Series(np.arange(len(allids)), index=allids)
assert len(allids) == C.shape[0]
idx = pos[R.contig].values
grp = R.group.values
A = sp.csr_matrix((np.ones(len(idx)), (pd.factorize(grp)[0], idx)), shape=(len(pd.unique(grp)), C.shape[0]))
names = pd.factorize(grp)[1]
E = pd.DataFrame((A @ C).toarray() / tot * 1e6, index=names, columns=S)
# CDS-based
X, G = genes(); cds = (G['type'].values == 'CDS'); ctot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
c2g = dict(zip(R.contig, R.group)); gg = G['contig'].map(c2g).values
sel = np.where(cds & pd.notna(gg))[0]
codes, gn = pd.factorize(gg[sel])
A2 = sp.csr_matrix((np.ones(len(sel)), (codes, sel)), shape=(len(gn), X.shape[0]))
E2 = pd.DataFrame((A2 @ X).toarray() / ctot * 1e6, index=gn, columns=S)
E.to_csv(O + 'group_expr_samples.tsv', sep='\t')
core = m[m.group.isin(['CTL', 'EXP']) & m.medium.isin(['SW', 'SI'])]
def med(E, s, md, ph, g=None):
    ix = core.index[(core.season == s) & (core.medium == md) & (core.phase == ph) & ((core.group == g) if g else True)]
    return E[ix].median(1) if len(ix) else pd.Series(np.nan, index=E.index)
out = pd.DataFrame(index=E.index)
out['contig_ppm_all_med'] = E[core.index].median(1); out['cds_ppm_all_med'] = E2.reindex(E.index)[core.index].median(1)
for s in ('2024-25', '2025-26'):
    for md in ('SW', 'SI'):
        out[f'oil_{s}_{md}'] = med(E, s, md, 'oiled'); out[f'ctl_{s}_{md}'] = med(E, s, md, 'control')
    out[f'pre_exp_{s}'] = med(E, s, 'SW', 'pre', 'EXP'); out[f'pre_ctl_{s}'] = med(E, s, 'SW', 'pre', 'CTL')
    out[f'log2_{s}_SW'] = np.log2((out[f'oil_{s}_SW'] + .5) / (out[f'ctl_{s}_SW'] + .5))
out['baci_2025-26'] = np.log2((out['oil_2025-26_SW'] + .5) / (out['pre_exp_2025-26'] + .5)) - np.log2((out['ctl_2025-26_SW'] + .5) / (out['pre_ctl_2025-26'] + .5))
sw = core.index[core.medium == 'SW']; si = core.index[core.medium == 'SI']
out['ice_vs_water_log2'] = np.log2((E[si].median(1) + .5) / (E[sw].median(1) + .5))
out.to_csv(O + 'group_expr_summary.tsv', sep='\t')
pd.set_option('display.width', 250)
print(out[['contig_ppm_all_med', 'cds_ppm_all_med', 'ctl_2025-26_SW', 'oil_2025-26_SW', 'pre_ctl_2025-26', 'pre_exp_2025-26', 'baci_2025-26', 'log2_2024-25_SW', 'ice_vs_water_log2']].round(1).sort_values('contig_ppm_all_med', ascending=False).to_string())
