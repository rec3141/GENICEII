"""Expression of metal and toxicity gene sets (CPM of CDS read pairs, corrected sample labels).
Out (metals/out/): cat_cpm.tsv (category x sample, confident tier), cat_cpm_all.tsv (all members), diesel.tsv (diesel vs control and
before-after-control-impact per season x medium), icewater.tsv (+ per-pair tables), hosts.tsv (category x bin CPM in diesel / control
seawater after dosing), degraders.tsv (within-host CPM of each category in the main oil degraders, by pool and phase)."""
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/metals/scripts')
from common import *
import scipy.sparse as sp
import contrasts as C
from msets import SETS
O = '/data/dev/rnaseq-exploration/metals/out/'
M = pd.read_csv(O + 'members.tsv', sep='\t', keep_default_na=False)
# confident: product evidence, or a KO on a gene whose Bakta product is not a bare hypothetical/domain-of-unknown-function
hyp = M['product'].str.contains(r'^hypothetical protein$|^DUF\d+ domain-containing protein$|^Uncharacterized protein', case=False, regex=True)
M['confident'] = M.evidence.str.contains('product') | ~hyp
M.to_csv(O + 'members.tsv', sep='\t', index=False)
X, G = genes(); S = samples(); m = meta().loc[S]
cds = G['type'].values == 'CDS'
tot = np.asarray(X[np.where(cds)[0]].sum(0)).ravel()
cats = [s[0] for s in SETS]; grp = {s[0]: s[1] for s in SETS}
def cat_matrix(sub):
    ci = sub.category.map({c: i for i, c in enumerate(cats)}).values
    A = sp.csr_matrix((np.ones(len(sub)), (ci, sub.row.values)), shape=(len(cats), X.shape[0]))
    return pd.DataFrame((A @ X).toarray() / tot * 1e6, index=cats, columns=S)
E = cat_matrix(M[M.confident]); E.to_csv(O + 'cat_cpm.tsv', sep='\t')
Ea = cat_matrix(M); Ea.to_csv(O + 'cat_cpm_all.tsv', sep='\t')

# diesel vs control, raw medians after dosing and BACI against each pool's own pre-dosing samples
rows = []
for s in ('2024-25', '2025-26'):
    for med in ('SW', 'SI'):
        q = (m.season == s) & (m.medium == med)
        oil, ctl = m.index[q & (m.phase == 'oiled')], m.index[q & (m.phase == 'control')]
        epre, cpre = m.index[q & (m.group == 'EXP') & (m.phase == 'pre')], m.index[q & (m.group == 'CTL') & (m.phase == 'pre')]
        if not len(oil) or not len(ctl): continue
        for c in cats:
            r = dict(category=c, group=grp[c], season=s, medium=med, n_oil=len(oil), n_ctl=len(ctl), n_pre_exp=len(epre), n_pre_ctl=len(cpre),
                     oil=E.loc[c, oil].median(), ctl=E.loc[c, ctl].median())
            r['log2_oil_ctl'] = np.log2((r['oil'] + 0.5) / (r['ctl'] + 0.5))
            if len(epre) and len(cpre):
                r['pre_exp'] = E.loc[c, epre].median(); r['pre_ctl'] = E.loc[c, cpre].median()
                r['log2_baci'] = r['log2_oil_ctl'] - np.log2((r['pre_exp'] + 0.5) / (r['pre_ctl'] + 0.5))
            rows.append(r)
DD = pd.DataFrame(rows); DD.to_csv(O + 'diesel.tsv', sep='\t', index=False)

# ice vs water and top vs bottom ice (paired, as in the Ice section)
L = np.log2(E + 1)
out = []
for name, P in (('icewater', C.PI), ('topbot', C.PT)):
    r, Dm = C.contrast(L, P)
    r.columns = [f'{name}_{k}' for k in r.columns]; out.append(r)
    Dm.to_csv(O + f'{name}_pairs.tsv', sep='\t')
IW = pd.concat(out, axis=1)
for lay in ('SW', 'TOP', 'BOT'):
    IW[f'cpm_{lay}'] = E[C.m.index[(C.m.layer == lay) & (C.m.season == '2024-25')]].median(axis=1)
IW.to_csv(O + 'icewater.tsv', sep='\t')

# hosts: category x bin CPM, diesel vs control seawater after dosing (per season)
bins = G['bin'].fillna('').values
Mc = M[M.confident].copy(); Mc['bin'] = bins[Mc.row.values]
hrows = []
for s in ('2024-25', '2025-26'):
    q = (m.season == s) & (m.medium == 'SW')
    for lab, idx in (('oil', m.index[q & (m.phase == 'oiled')]), ('ctl', m.index[q & (m.phase == 'control')])):
        cols = [S.index(i) for i in idx]
        cp = np.asarray(X[Mc.row.values][:, cols].multiply(1 / tot[cols]).mean(axis=1)).ravel() * 1e6
        t = pd.DataFrame({'category': Mc.category.values, 'bin': Mc.bin.values, 'cpm': cp}).groupby(['category', 'bin']).cpm.sum()
        hrows.append(t.rename(f'{lab}_{s}'))
Hs = pd.concat(hrows, axis=1).fillna(0); Hs.to_csv(O + 'hosts.tsv', sep='\t')

# within-host expression in the main oil degraders: category reads per million host CDS reads, per sample
DEG = {'dastool-metabat_9342': 'Cycloclasticus metabat_9342', 'dastool-lorbin_184': 'Cycloclasticus lorbin_184', 'dastool-lorbin_369': 'Oleispira lorbin_369',
       'dastool-lorbin_029': 'Parasphingorhabdus lorbin_029', 'dastool-lorbin_135': 'Cognaticolwellia lorbin_135', 'dastool-semibin_042': 'Colwellia semibin_042',
       'dastool-lorbin_032': 'Sulfitobacter lorbin_032', 'dastool-semibin_088': 'Sulfurimonas semibin_088',
       'dastool-vamb_4149': 'UBA9145 vamb_4149 (Pseudohongiellaceae)', 'dastool-comebin_2576': 'Porticoccaceae 50-400-T64 comebin_2576'}
ub, bi = np.unique(bins, return_inverse=True)
ind = sp.csr_matrix((cds.astype(float), (bi, np.arange(len(bins)))), shape=(len(ub), len(bins)))
BC = pd.DataFrame((ind @ X).toarray(), index=ub, columns=S)
drows = []
for b, lab in DEG.items():
    if b not in BC.index: continue
    hostreads = BC.loc[b]
    sub = Mc[Mc.bin == b]
    if not len(sub): continue
    ci = sub.category.map({c: i for i, c in enumerate(cats)}).values
    A = sp.csr_matrix((np.ones(len(sub)), (ci, sub.row.values)), shape=(len(cats), X.shape[0]))
    R = pd.DataFrame((A @ X).toarray(), index=cats, columns=S)
    ok = hostreads >= 5000
    W = (R.loc[:, ok] / hostreads[ok] * 1e6)
    for c in cats:
        for (s, g, ph), idx in m.loc[ok[ok].index].groupby(['season', 'group', 'phase']).groups.items():
            if g == 'other': continue
            v = W.loc[c, idx]
            drows.append(dict(bin=b, label=lab, category=c, season=s, pool=g, phase=ph, n=len(idx), median=v.median(), genes=int((sub.category == c).sum())))
pd.DataFrame(drows).to_csv(O + 'degraders.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
print(DD[DD.medium == 'SW'].pivot_table(index='category', columns='season', values=['log2_oil_ctl', 'log2_baci']).round(2).to_string())
print(IW[['icewater_delta', 'icewater_n_up', 'icewater_n', 'icewater_q', 'topbot_delta', 'topbot_n_up', 'topbot_n', 'topbot_q', 'cpm_SW', 'cpm_TOP', 'cpm_BOT']].round(2).to_string())
