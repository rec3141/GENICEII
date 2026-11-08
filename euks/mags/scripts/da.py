"""Differential abundance of the kernel blobs (kernel.py) in RNA and DNA.
Blob abundance per library: RNA = summed per-library contig reads per million (128 analysed RNA libraries);
DNA = summed share of the library's depth x length (CMO nanopore libraries present in mapping/depths.txt).
Contrasts (each blob, each data type; log2 with pseudocount = 1% of the blob's median non-zero value):
  diesel_vs_control_<season>: diesel-pool vs control-pool libraries after dosing (water + ice), Mann-Whitney; effect = log2 ratio of medians
  baci_2025-26 (water): (post - pre) in the diesel pool minus (post - pre) in the control, on log2 medians; p = Mann-Whitney of
      post-dosing log values centred on each pool's own pre-dosing median
  ice_vs_water: paired libraries (same season, pool, day; RNA pairs from ice/scripts/contrasts.py), Wilcoxon signed-rank
BH across blobs within each contrast x data type. One pool per treatment: 'diesel' is a pool effect.
Out: out/da_blobs.tsv (long), out/da_blobs_wide.tsv, out/da_blobs_abundance_{rna,dna}.tsv"""
import pandas as pd, numpy as np, sys
from scipy.stats import mannwhitneyu, wilcoxon
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import meta, nanopore_samples, bh, OIL_DAY
import contrasts as C
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
K = pd.read_csv(O + 'kernel_blobs.tsv', sep='\t'); KM = pd.read_csv(O + 'kernel_blob_contigs.tsv', sep='\t')
RNA = pd.read_pickle(O + 'kernel_rna_cpm.pkl'); DNA = pd.read_pickle(O + 'kernel_dna_share.pkl')
# blob names and calls
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
dom = pd.read_csv('/data/dev/rnaseq-exploration/data/bin_domain.tsv', sep='\t').set_index('bin').call
EG = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t').set_index('contig').euk_group
VG = pd.read_csv(O + 'contig_virus_group.tsv', sep='\t').set_index('contig').vid
info = []
for i, x in K.iterrows():
    ids = KM.contig[KM.seed == x.seed]; g = T.loc[ids]; L = g.length.sum()
    if L == 0: info.append(dict(seed=x.seed, blob=f'B{i + 1:02d}', call='none')); continue
    eg = g.index.map(lambda c: EG.get(c, '')); vg = g.index.map(lambda c: VG.get(c, ''))
    egs = g.length[eg != ''].groupby(eg[eg != '']).sum(); vgs = g.length[vg != ''].groupby(vg[vg != '']).sum()
    db = g.bin.fillna(''); prok = g.length[db.map(lambda b: dom.get(b, '') == 'prokaryote').values].sum() / L
    pb = g.length[db != ''].groupby(db[db != '']).sum().sort_values(ascending=False)
    if x.virus_bp >= 0.5: call = 'virus' + (f' (= {vgs.idxmax()} {vgs.max() / L:.0%})' if len(vgs) else ' (new)')
    elif x.tiara_euk >= 0.5: call = 'eukaryote' + (f' (= {egs.idxmax()} {egs.max() / L:.0%})' if len(egs) and egs.max() / L >= 0.3 else ' (new)') + (' diffuse' if x.seed in ('9,-48', '-2,-38') else '')
    elif prok >= 0.5: call = f'prokaryote MAG {pb.index[0].replace("dastool-", "")} ({prok:.0%})'
    else: call = f'prokaryote, partly binned ({prok:.0%} in MAGs' + (f'; {pb.index[0].replace("dastool-", "")}' if len(pb) else '') + ')'
    info.append(dict(seed=x.seed, blob=f'B{i + 1:02d}', call=call, bp=int(L), gc=x.gc, tiara_euk=x.tiara_euk, virus_bp=x.virus_bp))
I = pd.DataFrame(info)
# abundance matrices
ab_r = pd.DataFrame({b: RNA.loc[KM.contig[KM.seed == s]].sum() for s, b in zip(I.seed, I.blob) if (KM.seed == s).any()}).T
ab_d = pd.DataFrame({b: DNA.reindex(KM.contig[KM.seed == s]).fillna(0).sum() for s, b in zip(I.seed, I.blob) if (KM.seed == s).any()}).T
ab_r.to_csv(O + 'da_blobs_abundance_rna.tsv', sep='\t'); ab_d.to_csv(O + 'da_blobs_abundance_dna.tsv', sep='\t')
# metadata
m = meta()
mr = m.loc[ab_r.columns, ['season', 'group', 'medium', 'phase', 'day']]
ns = nanopore_samples(); di = pd.read_csv('/data/dev/rnaseq-exploration/out/dna_samples.tsv', sep='\t')
di['sample_name'] = di.sample_name.astype(str); dmap = di.drop_duplicates('sample_name').set_index('sample_name')
md = []
for rb in ab_d.columns:
    sn = ns.sample_name[ns.run_barcode == rb]; sn = sn.iloc[0] if len(sn) else ''
    r = dmap.loc[sn] if sn in dmap.index else None
    if r is None or pd.isna(r.day) or r.group not in ('CTL', 'EXP'): md.append(dict(lib=rb, season=None, group=None, medium=None, day=np.nan, phase='other')); continue
    od = OIL_DAY['CMO ' + r.season]
    ph = 'pre' if r.day < 0 or (r.group == 'EXP' and r.day < od) and r.day < 0 else ('oiled' if (r.group == 'EXP' and r.day >= od) else ('control' if r.group == 'CTL' and r.day >= 0 else 'pre'))
    md.append(dict(lib=rb, season=r.season, group=r.group, medium=r.medium, day=r.day, phase=ph))
md = pd.DataFrame(md).set_index('lib')
print('RNA libraries', ab_r.shape[1], '; DNA libraries', ab_d.shape[1], 'of which CTL/EXP with day', int(md.season.notna().sum()))
def lg(v):
    nz = v[v > 0]; pc = 0.01 * (nz.median() if len(nz) else 1e-9); return np.log2(v + pc)
res = []
def add(blob, dt, con, eff, p, n1, n2):
    res.append(dict(blob=blob, data=dt, contrast=con, log2_effect=eff, p=p, n_a=n1, n_b=n2))
for dt, A, M in (('RNA', ab_r, mr), ('DNA', ab_d, md)):
    for b in A.index:
        v = lg(A.loc[b])
        for s in ('2024-25', '2025-26'):
            o = M.index[(M.season == s) & (M.phase == 'oiled') & M.medium.isin(['SW', 'SI'])]; c = M.index[(M.season == s) & (M.phase == 'control') & M.medium.isin(['SW', 'SI'])]
            if len(o) >= 3 and len(c) >= 3:
                add(b, dt, f'diesel_vs_control_{s}', v[o].median() - v[c].median(), mannwhitneyu(v[o], v[c]).pvalue, len(o), len(c))
        s = '2025-26'; sw = (M.season == s) & (M.medium == 'SW')
        pe, pc_, o, c = [M.index[sw & (M.phase == 'pre') & (M.group == g)] for g in ('EXP', 'CTL')] + [M.index[sw & (M.phase == ph)] for ph in ('oiled', 'control')]
        if min(len(pe), len(pc_)) >= 2 and min(len(o), len(c)) >= 3:
            eff = (v[o].median() - v[pe].median()) - (v[c].median() - v[pc_].median())
            add(b, dt, 'baci_2025-26_SW', eff, mannwhitneyu(v[o] - v[pe].median(), v[c] - v[pc_].median()).pvalue, len(o), len(c))
        # ice vs water pairs
        if dt == 'RNA':
            P = [dict(p, a=[x for x in p['a'] if x in A.columns], b=[x for x in p['b'] if x in A.columns]) for p in C.PI]
            P = [p for p in P if p['a'] and p['b']]
            dl = [v[p['a']].mean() - v[p['b']].mean() for p in P]
        else:
            dl = []
            for (s_, g_, d_), h in M[M.medium.isin(['SW', 'SI'])].groupby(['season', 'group', 'day']):
                if (h.medium == 'SI').any() and (h.medium == 'SW').any(): dl.append(v[h.index[h.medium == 'SI']].mean() - v[h.index[h.medium == 'SW']].mean())
        dl = np.array(dl)
        if len(dl) >= 5: add(b, dt, 'ice_vs_water', float(np.mean(dl)), wilcoxon(dl).pvalue if np.any(dl != 0) else 1.0, int((dl > 0).sum()), len(dl))
R = pd.DataFrame(res)
R['q'] = np.nan
for (dt, con), g in R.groupby(['data', 'contrast']): R.loc[g.index, 'q'] = bh(g.p.values)
R = R.merge(I, on='blob', how='left'); R.to_csv(O + 'da_blobs.tsv', sep='\t', index=False)
Wd = R.pivot_table(index=['blob', 'seed', 'call'], columns=['data', 'contrast'], values='log2_effect').round(2)
Wq = R.pivot_table(index=['blob', 'seed', 'call'], columns=['data', 'contrast'], values='q')
Wd.to_csv(O + 'da_blobs_wide.tsv', sep='\t')
pd.set_option('display.width', 300); pd.set_option('display.max_colwidth', 40)
print(I.to_string(index=False))
fmt = lambda e, q: '' if pd.isna(e) else f'{e:+.1f}' + ('**' if q < 0.01 else '*' if q < 0.05 else '')
show = pd.DataFrame({f'{d}:{c}': [fmt(Wd.loc[i, (d, c)], Wq.loc[i, (d, c)]) for i in Wd.index] for d, c in Wd.columns}, index=[i[0] + ' ' + i[1] for i in Wd.index])
print(show.to_string())
