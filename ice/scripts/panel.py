"""Curated ice-physiology panel: KO sets summed per sample (corrected KO CPM), paired ice-water and top-bottom deltas,
plus the main taxa carrying each category's transcripts in ice. Out: ice/out/panel.tsv, ice/out/panel_hosts.tsv"""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
import contrasts as C
TAX = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', index_col=0)['classification']
PANEL = {
 'Osmolyte uptake (opu/pro ABC, BetT, ProP)': 'K05845 K05846 K05847 K02000 K02001 K02002 K03451 K05020 K03762 K03607',
 'Glycine betaine synthesis (betA/betB)': 'K00108 K00130 K17755',
 'Ectoine synthesis (ectABC/D)': 'K06718 K00836 K06720 K10674',
 'Polyamine synthesis': 'K01585 K01581 K01480 K00797 K01611 K12251 K10536 K27683 K27502',
 'Peroxiredoxins and catalase': 'K24119 K24126 K03387 K04063 K03781 K11065 K03386',
 'Dps and osmY stress proteins': 'K04047 K04065',
 'Cold-shock proteins (cspA)': 'K03704',
 'DEAD-box RNA helicases (csdA/deaD, srmB, rhlB, rhlE)': 'K05592 K05590 K03732 K11927',
 'Ribosome hibernation (raiA, hpf, rmf)': 'K05809 K05808 K03812',
 'Branched-chain fatty acid start (bkd, fabH)': 'K00166 K00167 K09699 K00648',
 'Unsaturated fatty acids (fabA/fabB, desaturases)': 'K01716 K00647 K00507 K10255 K10256',
 'Light harvesting and PSII (LHC, psbA)': 'K08912 K08913 K08914 K08915 K08916 K08917 K02703 K02704 K02705',
 'RuBisCO (rbcL/rbcS)': 'K01601 K01602',
 'MSHA pili and adhesins': 'K12276 K12277 K12278 K12279 K12280 K12281 K12282 K12283 K12284 K12285 K12286 K12287 K21449',
 'Gas vesicles (gvpA)': 'K23262',
 'Cyanophycin and PHA storage': 'K03802 K13282 K03821 K22881',
 'Flagellar motility (bacterial)': 'K02406 K02390 K02556 K02557 K02397 K02396',
 'Methylotrophy (xoxF, exaA, mxaF)': 'K23995 K00114 K14028',
 'DMSP and methanethiol (dmdA, ddd, mddA)': 'K17486 K16953 K21310 K16954',
 'Siderophore uptake (fiu, cirA, tonB)': 'K16090 K16089 K03832',
}
K = pd.read_csv(D + 'ko_cpm.tsv', sep='\t', index_col=0)
S = pd.DataFrame({k: K.reindex(v.split()).fillna(0).sum() for k, v in PANEL.items()}).T
L = np.log2(S + 1)
rows = []
for name, P in (('icewater', C.PI), ('topbot', C.PT)):
    r, Dm = C.contrast(L, P)
    r.columns = [f'{name}_{c}' for c in r.columns]; rows.append(r)
    Dm.to_csv(C.O + f'panel_{name}_pairs.tsv', sep='\t')
R = pd.concat(rows, axis=1)
R['mean_cpm'] = S.mean(axis=1)
m = C.m
for lay in ('SW', 'TOP', 'BOT'):
    R[f'cpm_{lay}'] = S[m.index[(m.layer == lay) & (m.season == '2024-25')]].median(axis=1)
R.to_csv(C.O + 'panel.tsv', sep='\t')
pd.set_option('display.width', 250)
print(R[['icewater_delta', 'icewater_n_up', 'icewater_n', 'icewater_q', 'topbot_delta', 'topbot_n_up', 'topbot_n', 'topbot_q', 'cpm_SW', 'cpm_TOP', 'cpm_BOT']].round(2).to_string())

# who carries each category in ice: gene-level counts -> bins (corrected KO in genemap)
X, G = genes()
cds = (G.type == 'CDS').values
lib = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()  # CDS reads per library, as in hc.py
samp = samples(); ice = [samp.index(s) for s in m.index[m.medium == 'SI'] if s in samp]
cpm_ice = np.asarray(X[:, ice].multiply(1 / np.maximum(lib[ice], 1)).sum(axis=1)).ravel() * 1e6 / len(ice)  # mean CPM over ice libraries
host = []
for k, v in PANEL.items():
    ks = set(v.split())
    sel = G.KO.fillna('').str.split(',').apply(lambda x: bool(ks.intersection(x))).values & cds
    d = pd.DataFrame({'bin': G.bin[sel].fillna('unbinned').values, 'cpm': cpm_ice[sel]}).groupby('bin').cpm.sum().sort_values(ascending=False)
    tot = d.sum()
    for b, c in d.head(5).items():
        lab = special_label(b) or (TAX.get(b.replace('dastool-', ''), '') if b != 'unbinned' else 'unbinned')
        host.append(dict(category=k, bin=b, share=c / tot if tot else 0, taxon=str(lab).split(';')[-1] if ';' in str(lab) else lab))
H = pd.DataFrame(host); H.to_csv(C.O + 'panel_hosts.tsv', sep='\t', index=False)
print(H.groupby('category').head(2).round(2).to_string())
