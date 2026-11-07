"""Per-MAG transcriptome PCA and state shifts.
For each MAG with enough reads, profile = its own CDS counts / its own total (composition drops out), CLR, PCA on genes detected in >=50% of
kept samples. Then: discrete states (GMM BIC, k-means silhouette on PC1-3), depth artefact check, PERMANOVA (top PCs, permutations within medium)
for diesel-after-dosing vs control-after-dosing and vs all reference samples, BACI-style shift ratio (distance from own pool's pre-dosing centroid),
ice vs water, single changepoint per pool x season x medium on PC1, functional interpretation (within-genome category shares vs PC scores),
growth index (ribosomal-protein share), rarefaction robustness for key MAGs.
Out: states/out/{mag_states.tsv, scores.pkl, functions.tsv, loadings.pkl, changepoints.tsv, rarefy.tsv, growth.pkl}"""
import sys, os, re; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
from scipy.stats import spearmanr
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
O = '/data/dev/rnaseq-exploration/states/out/'; E0 = '/data/dev/rnaseq-exploration/'
MIN_READS, MIN_SAMPLES, NPC_TEST, NPERM = 3000, 20, 5, 999
# PCA axes are fitted on samples with >= FIT_READS reads on the genome; shallower samples (MIN_READS..FIT_READS) are projected
# onto those axes, so sampling noise in shallow libraries cannot create or rotate axes (they showed up as large PCA outliers).
FIT_READS, MIN_FIT = 10000, 15
KEY = ['lorbin_184', 'metabat_9342', 'metabat_5432', 'comebin_3787_sub', 'comebin_2576', 'semibin_088', 'vamb_549', 'metabat_5376', 'vamb_4149', 'lorbin_095', 'metabat_9833']
rng = np.random.default_rng(1)

X, G = genes(); S = samples(); m = meta().loc[S]
core = m.group.isin(['CTL', 'EXP']).values & m.medium.isin(['SW', 'SI']).values
cs = np.where(core)[0]; Sc = [S[i] for i in cs]; mc = m.iloc[cs]
cds = np.where((G['type'].values == 'CDS') & G['bin'].notna().values & (G['bin'].values != ''))[0]
Xc = X[cds][:, cs].tocsr()
bins = G['bin'].values[cds]
codes, ub = pd.factorize(bins)
import scipy.sparse as sp
A = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(ub), len(codes)))
depth = pd.DataFrame((A @ Xc).toarray(), index=ub, columns=Sc)
npass = (depth >= MIN_READS).sum(1)
sel = npass[npass >= MIN_SAMPLES].index.tolist()
print(f'{len(ub)} binned MAGs; {len(sel)} with >= {MIN_SAMPLES} core samples at >= {MIN_READS} CDS reads', flush=True)
order = {b: np.where(codes == k)[0] for k, b in enumerate(ub)}

gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
def blab(b):
    if plastid_label(b): return plastid_label(b)
    t = [y.split('__', 1)[1] for y in str(gt.get(b.replace('dastool-', ''), '')).split(';') if '__' in y and y.split('__', 1)[1]]
    return (t[-1] if t else '') or kaiju_bin_label(b) or b.replace('dastool-', '')

# function categories (global gene positions)
prod = G['product'].fillna('').values; ko = G['KO'].fillna('').values
def rx(p, flags=re.I): r = re.compile(p, flags); return np.array([bool(r.search(x)) for x in prod])
ribo = rx(r'\b(30S|50S)? ?ribosomal protein\b') & ~rx(r'methyltransferase|acetyltransferase|modifying|maturation|kinase|synthase|-L-|isomerase|binding factor|hydrolase|protease|silencing')
hcm = pd.read_pickle(E0 + 'out/hc_members.pkl'); hcset = np.zeros(len(G), bool); hcset[np.concatenate([np.asarray(v) for v in hcm.values()])] = True
ng = pd.read_csv(E0 + 'nitrogen/out/genes.tsv', sep='\t', keep_default_na=False, usecols=['gene_id']); nmem = pd.read_pickle(E0 + 'nitrogen/out/members.pkl')
NS = ['amtB (ammonium transporter)', 'glnB/glnK (PII nitrogen regulators)', 'ntrBC (Ntr two-component system)', 'glnA (glutamine synthetase)',
      'urtABCDE (urea ABC transporter)', 'nitrate/nitrite transporters (nrtABC, narK, nasA-type)']
nids = set(ng.gene_id.values[np.unique(np.concatenate([nmem[k] for k in NS]))])
nset = G['gene_id'].isin(nids).values
CATS = {
 'ribosomal proteins (growth)': ribo,
 'translation factors': rx(r'elongation factor (Tu|G|Ts|P)\b|translation initiation factor'),
 'hydrocarbon degradation': hcset,
 'N acquisition': nset,
 'phage shock / envelope stress': rx(r'phage shock protein|\bPsp[ABCDEF]\b|envelope stress|sigma-E|RpoE\b'),
 'cold shock': rx(r'cold[- ]shock'),
 'chaperones / proteases': rx(r'chaperon|\bDnaK\b|\bDnaJ\b|GroEL|GroES|\bHtpG\b|\bClpB\b|small heat shock|\bIbp[AB]\b|heat shock protein|\bHslU\b|\bHslV\b|\bLon\b protease|ATP-dependent protease'),
 'oxidative stress (ROS)': rx(r'catalase|superoxide dismutase|peroxiredoxin|alkyl hydroperoxide|glutathione peroxidase|thioredoxin reductase|\bOxyR\b|\bKatG\b|\bKatE\b|\bAhpC\b|\bSodA\b|\bSodB\b'),
 'motility / chemotaxis': rx(r'flagell|\bFli[A-Z]\b|\bFlg[A-Z]\b|\bFlh[A-Z]\b|chemotaxis|\bChe[ABRWYZ]\b|methyl-accepting|type IV pil|twitching'),
 'transporters': rx(r'ABC transporter|TRAP transporter|TonB|transporter|permease|porin'),
 'storage (PHA, polyP, glycogen)': rx(r'poly\(3-hydroxy|polyhydroxy|PHA synth|\bPhaC\b|phasin|polyphosphate kinase|exopolyphosphatase|\bPpk\b|glycogen synthase|\bGlgA\b|glycogen branching|glycogen phosphorylase'),
 'CRISPR / defence': rx(r'CRISPR|\bCas\d|restriction endonuclease|abortive infection'),
}
cats_loc = {k: v[cds] for k, v in CATS.items()}

def clr_profile(C):
    """C: genes x samples counts. log2 within-genome CPM with a depth-scaled prior (edgeR-style: prior 0.5 reads at mean depth),
    then centred per sample (CLR-like). The scaled prior keeps zeros in shallow samples from pulling them together."""
    N = C.sum(0); pr = 0.5 * N / N.mean()
    L = np.log2((C + pr[None, :]) / (N + 2 * pr)[None, :] * 1e6); return L - L.mean(0, keepdims=True)
def pca(L):
    Z = L - L.mean(1, keepdims=True)
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    ev = s ** 2 / (s ** 2).sum()
    return (Vt.T * s), U, ev           # scores samples x pcs, loadings genes x pcs
def permanova(Y, labels, strata, nperm=NPERM):
    """Pseudo-F / R2 on Euclidean distances among rows of Y for a two-group factor; permute within strata."""
    labels = np.asarray(labels); strata = np.asarray(strata)
    def stat(lb):
        gm = Y.mean(0); sst = ((Y - gm) ** 2).sum(); ssw = sum(((Y[lb == g] - Y[lb == g].mean(0)) ** 2).sum() for g in np.unique(lb))
        ssa = sst - ssw; k = len(np.unique(lb)); n = len(lb)
        return (ssa / (k - 1)) / (ssw / (n - k)), ssa / sst
    F0, R2 = stat(labels); cnt = 0
    for _ in range(nperm):
        lp = labels.copy()
        for s in np.unique(strata):
            ix = np.where(strata == s)[0]; lp[ix] = rng.permutation(lp[ix])
        if stat(lp)[0] >= F0: cnt += 1
    return F0, R2, (cnt + 1) / (nperm + 1)
def changepoint(days, v):
    """best single split of per-day means (>=2 days each side): returns (day after which it shifts, effect in pooled SD, before, after)"""
    d = pd.Series(v, index=days).groupby(level=0).mean().sort_index()
    if len(d) < 4: return None
    best = None
    for i in range(2, len(d) - 1):
        a, b = d.iloc[:i], d.iloc[i:]
        sse = ((a - a.mean()) ** 2).sum() + ((b - b.mean()) ** 2).sum()
        if best is None or sse < best[0]: best = (sse, i, a.mean(), b.mean())
    sse, i, ma, mb = best; sd = np.sqrt(sse / max(len(d) - 2, 1)) or 1e-9
    return float(d.index[i]), float((mb - ma) / sd), float(ma), float(mb)

rows, scores, funcs, loads, cps, growth, rare = [], {}, [], {}, [], {}, []
med = mc.medium.values; ses = mc.season.values; ph = mc.phase.values; pool = mc.group.values; day = mc.day.values.astype(float)
for b in sel:
    gi = order[b]; keep = (depth.loc[b].values >= MIN_READS)
    C = Xc[gi][:, keep].toarray().astype(float)
    det = ((C > 0).mean(1) >= 0.9) & (C.mean(1) >= 5)          # robustly detected genes only
    if det.sum() > 3000: det[np.where(det)[0][np.argsort(-C[det].mean(1))[3000:]]] = False
    if det.sum() < 100: continue
    Cd = C[det]; L = clr_profile(Cd)
    dep = depth.loc[b].values[keep]
    fit = dep >= FIT_READS
    if fit.sum() < MIN_FIT: fit = np.ones_like(fit)           # too few deep samples: fit on all
    mu = L[:, fit].mean(1, keepdims=True); _, U, ev = pca(L[:, fit])
    sc = (L - mu).T @ U                                        # fitted samples: identical to PCA scores; others projected
    k = np.where(keep)[0]
    r1 = spearmanr(sc[:, 0], np.log(dep))[0]; r2 = spearmanr(sc[:, 1], np.log(dep))[0]
    # discrete states on PC1-3
    Y3 = sc[:, :3]
    bic = [GaussianMixture(n, covariance_type='full', random_state=0, n_init=3).fit(Y3).bic(Y3) for n in (1, 2, 3, 4)]
    sil = {}
    for n in (2, 3, 4):
        lab = KMeans(n, n_init=10, random_state=0).fit_predict(Y3); sil[n] = silhouette_score(Y3, lab)
    kbest = int(np.argmin(bic)) + 1; nbest = max(sil, key=sil.get)
    labs = KMeans(nbest, n_init=10, random_state=0).fit_predict(Y3)
    Yt = sc[:, :NPC_TEST]
    out = dict(bin=b, taxon=blab(b), n=int(keep.sum()), n_fit=int(fit.sum()), genes=int(det.sum()), pc1=ev[0], pc2=ev[1], pc3=ev[2], depth_med=float(np.median(dep)),
               r_depth_pc1=r1, r_depth_pc2=r2, gmm_k=kbest, gmm_dbic=bic[0] - min(bic), km_k=nbest, sil=sil[nbest])
    mk, sk, pk, pl, dk = med[k], ses[k], ph[k], pool[k], day[k]
    for s in ('2024-25', '2025-26'):
        o = (sk == s) & (pk == 'oiled'); c = (sk == s) & (pk == 'control'); pr = (sk == s) & (pk == 'pre')
        if o.sum() >= 4 and c.sum() >= 4:
            ix = o | c; F, R2, p = permanova(Yt[ix], np.where(o[ix], 'o', 'c'), mk[ix])
            out.update({f'F_{s}': F, f'R2_{s}': R2, f'p_{s}': p})
            ix2 = o | c | pr; F2, R22, p2 = permanova(Yt[ix2], np.where(o[ix2], 'o', 'r'), mk[ix2])
            out.update({f'R2ref_{s}': R22, f'pref_{s}': p2})
            # shift ratio: water samples, distance from own pool's pre-dosing centroid (2025-26 has 5 pre per pool)
            pe = pr & (pl == 'EXP') & (mk == 'SW'); pc_ = pr & (pl == 'CTL') & (mk == 'SW')
            if pe.sum() >= 1 and pc_.sum() >= 1:
                de = np.linalg.norm(Yt[o & (mk == 'SW')] - Yt[pe].mean(0), axis=1); dc = np.linalg.norm(Yt[c & (mk == 'SW')] - Yt[pc_].mean(0), axis=1)
                if len(de) and len(dc): out[f'shift_{s}'] = np.median(de) / max(np.median(dc), 1e-9)
            # centroid separation in units of within-group spread (effect size)
            cen = np.linalg.norm(Yt[o].mean(0) - Yt[c].mean(0)); spr = np.sqrt((((Yt[o] - Yt[o].mean(0)) ** 2).sum(1).mean() + ((Yt[c] - Yt[c].mean(0)) ** 2).sum(1).mean()) / 2)
            out[f'sep_{s}'] = cen / spr
        # ice vs water within season (after dosing, both pools), strata = pool
        iw = (sk == s) & np.isin(pk, ['oiled', 'control'])
        if (iw & (mk == 'SI')).sum() >= 4 and (iw & (mk == 'SW')).sum() >= 4:
            F, R2, p = permanova(Yt[iw], mk[iw], pl[iw]); out.update({f'R2ice_{s}': R2, f'pice_{s}': p})
        # changepoints on PC1 per pool x medium
        for P in ('EXP', 'CTL'):
            for M in ('SW', 'SI'):
                ix = (sk == s) & (pl == P) & (mk == M)
                if ix.sum() >= 5:
                    c_ = changepoint(dk[ix], sc[ix, 0])
                    if c_: cps.append(dict(bin=b, season=s, pool=P, medium=M, day=c_[0], effect=c_[1], before=c_[2], after=c_[3], n=int(ix.sum())))
    rows.append(out)
    scores[b] = {'samples': [Sc[i] for i in k], 'pc': sc[:, :3], 'ev': ev[:3], 'cluster': labs, 'depth': dep}
    # functions: within-genome share per sample vs PC1/PC2
    tot = C.sum(0)
    for cname, msk in cats_loc.items():
        mm = msk[gi]
        if mm.sum() < 3: continue
        shr = C[mm].sum(0) / tot
        if (shr > 0).sum() < 5: continue
        funcs.append(dict(bin=b, category=cname, n_genes=int(mm.sum()), mean_share=float(shr.mean()),
                          r_pc1=spearmanr(shr, sc[:, 0])[0], r_pc2=spearmanr(shr, sc[:, 1])[0]))
        if cname.startswith('ribosomal'): growth[b] = dict(zip([Sc[i] for i in k], shr))
    gid = G.gene_id.values[cds[gi]][det]; pr_ = prod[cds[gi]][det]
    lo = {}
    for j in (0, 1):
        o_ = np.argsort(U[:, j]); lo[j] = {'neg': [(gid[i], pr_[i], float(U[i, j])) for i in o_[:12]], 'pos': [(gid[i], pr_[i], float(U[i, j])) for i in o_[::-1][:12]]}
    loads[b] = lo
    # rarefaction robustness for key MAGs: binomial thinning to the smallest kept depth
    if True:
        tgt = dep.min(); Cr = rng.binomial(Cd.astype(int), np.minimum(1, tgt / C.sum(0))[None, :]).astype(float)
        scr, _, evr = pca(clr_profile(Cr))
        r = abs(spearmanr(scr[:, 0], sc[:, 0])[0]); r2_ = abs(spearmanr(scr[:, 1], sc[:, 1])[0])
        rr = dict(bin=b, depth=float(tgt), r_pc1=r, r_pc2=r2_, pc1_rare=evr[0], r_depth_pc1_rare=spearmanr(scr[:, 0], np.log(dep))[0])
        o = (sk == '2025-26') & (pk == 'oiled'); c = (sk == '2025-26') & (pk == 'control')
        if o.sum() >= 4 and c.sum() >= 4 and b.replace('dastool-', '') in KEY:
            ix = o | c; F, R2, p = permanova(scr[ix, :NPC_TEST], np.where(o[ix], 'o', 'c'), mk[ix], nperm=499); rr.update(R2_2025_26_rare=R2, p_2025_26_rare=p)
        rare.append(rr)
    print(b, out['taxon'], out['n'], out['genes'], round(ev[0], 3), round(r1, 2), out.get('R2_2025-26'), out.get('p_2025-26'), flush=True)

R = pd.DataFrame(rows)
for s in ('2024-25', '2025-26'):
    for c in ('p', 'pref', 'pice'):
        col = f'{c}_{s}'
        if col in R: ok = R[col].notna(); R.loc[ok, col.replace('p', 'q', 1)] = bh(R.loc[ok, col].values)
R.to_csv(O + 'mag_states.tsv', sep='\t', index=False)
pd.DataFrame(funcs).to_csv(O + 'functions.tsv', sep='\t', index=False)
pd.DataFrame(cps).to_csv(O + 'changepoints.tsv', sep='\t', index=False)
pd.DataFrame(rare).to_csv(O + 'rarefy.tsv', sep='\t', index=False)
pd.to_pickle(scores, O + 'scores.pkl'); pd.to_pickle(loads, O + 'loadings.pkl'); pd.to_pickle(growth, O + 'growth.pkl')
depth.loc[sel].to_pickle(O + 'depth.pkl')
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
print(R.sort_values('R2_2025-26', ascending=False)[['bin', 'taxon', 'n', 'genes', 'pc1', 'r_depth_pc1', 'gmm_k', 'sil', 'R2_2025-26', 'q_2025-26', 'shift_2025-26', 'R2_2024-25', 'q_2024-25']].head(40).round(3).to_string())
print(pd.DataFrame(rare).round(3).to_string())
