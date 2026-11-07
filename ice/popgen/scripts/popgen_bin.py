"""popgen_bin.py BIN : from pile/BIN.npz -> pop/BIN.json (diversity, FST contrasts, pN/pS) and pop/BIN.genes.tsv (per-gene FST ice vs water)."""
import sys, json, gzip, subprocess, numpy as np, pandas as pd, os
D = '/data/dev/rnaseq-exploration/ice/popgen'
b = sys.argv[1]; os.makedirs(f'{D}/pop', exist_ok=True)
z = np.load(f'{D}/pile/{b}.npz', allow_pickle=True)
S = pd.read_csv(f'{D}/samples_map.tsv', sep='\t', index_col=0).loc[z['samples']]
C = z['counts'].astype(np.int32)            # sites x samples x 4
contig, pos = z['contig'], z['pos']
gn = list(z['gnames']); cov_g = z['cov_g']; TH = list(z['th'])
MIN = 10

def pooled(mask): return C[:, mask, :].sum(1)
def pi_site(a):
    d = a.sum(1); p = a / np.maximum(d, 1)[:, None]
    return np.where(d > 1, d / np.maximum(d - 1, 1) * (1 - (p ** 2).sum(1)), 0), d
def fst(a1, a2, mind=MIN):
    p1, d1 = pi_site(a1); p2, d2 = pi_site(a2)
    ok = (d1 >= mind) & (d2 >= mind)
    f1 = a1[ok] / d1[ok, None]; f2 = a2[ok] / d2[ok, None]
    hb = 1 - (f1 * f2).sum(1); hw = (p1[ok] + p2[ok]) / 2
    num, den = (hb - hw), hb
    # major-allele (consensus) differences with both groups >= 80% major
    fixed = int(((f1.argmax(1) != f2.argmax(1)) & (f1.max(1) >= .8) & (f2.max(1) >= .8)).sum())
    return dict(fst=float(num.sum() / den.sum()) if den.sum() > 0 else None, sites=int(ok.sum()), fixed=fixed), ok, num, den

lay = S.layer.values; sea = S.season.values; pool = S.pool.values; day = S.day.values
ice = lay != 'SW'
out = {'bin': b, 'n_sites_candidate': int(len(pos)), 'groups': {}, 'contrasts': {}, 'pairs': []}
# diversity per group (pooled)
for g in gn:
    se, l = g.split('|')
    m = (sea == se) & ((lay == 'SW') if l == 'SW' else (ice if l == 'ICE' else (lay == l)))
    if not m.any(): continue
    ps, d = pi_site(pooled(m)); den = cov_g[gn.index(g), TH.index(MIN)]
    f = pooled(m) / np.maximum(d, 1)[:, None]; minor = 1 - f.max(1)
    snv = int(((d >= MIN) & (minor >= .05)).sum())
    out['groups'][g] = dict(n=int(m.sum()), covered_bp=int(den), pi=float(ps[d >= MIN].sum() / den) if den > 1e5 else None,
                            snv_per_kb=float(snv / den * 1e3) if den > 1e5 else None, reads=int(C[:, m, :].sum()))
# group contrasts
def M(**kw):
    m = np.ones(len(S), bool)
    for k, v in kw.items(): m &= np.isin({'sea': sea, 'lay': lay, 'pool': pool}[k], v)
    return m
gene_inputs = {}
for se in ('2024-25', '2025-26'):
    W = M(sea=se, lay=['SW']); I = M(sea=se) & ice
    tests = {'ice_vs_water': (I, W), 'top_vs_bot': (M(sea=se, lay=['TOP']), M(sea=se, lay=['BOT'])),
             'top_vs_water': (M(sea=se, lay=['TOP']), W), 'bot_vs_water': (M(sea=se, lay=['BOT']), W),
             'water_ctl_vs_exp': (W & M(pool=['CTL']), W & M(pool=['EXP'])), 'ice_ctl_vs_exp': (I & M(pool=['CTL']), I & M(pool=['EXP'])),
             'water_early_vs_late': (W & (day <= 12), W & (day > 12)), 'ice_early_vs_late': (I & (day <= 12), I & (day > 12))}
    for k, (m1, m2) in tests.items():
        if m1.sum() == 0 or m2.sum() == 0: continue
        r, ok, num, den = fst(pooled(m1), pooled(m2))
        if r['sites'] >= 50: out['contrasts'][f'{se}|{k}'] = r
        if k == 'ice_vs_water': gene_inputs[se] = (ok, num, den)
    # paired same pool-day: ice sample(s) vs water sample
    for (pl, dy), sub in S[(sea == se)].groupby(['pool', 'day']):
        mi = (sea == se) & (pool == pl) & (day == dy) & ice; mw = (sea == se) & (pool == pl) & (day == dy) & (lay == 'SW')
        if mi.any() and mw.any():
            r, *_ = fst(pooled(mi), pooled(mw), mind=6)
            if r['sites'] >= 50: out['pairs'].append(dict(season=se, pool=pl, day=int(dy), kind='ice_vs_water', **r))
        for l1, l2 in (('TOP', 'BOT'),):
            m1 = (sea == se) & (pool == pl) & (day == dy) & (lay == l1); m2 = (sea == se) & (pool == pl) & (day == dy) & (lay == l2)
            if m1.any() and m2.any():
                r, *_ = fst(pooled(m1), pooled(m2), mind=6)
                if r['sites'] >= 50: out['pairs'].append(dict(season=se, pool=pl, day=int(dy), kind='top_vs_bot', **r))
    # water-vs-water null: consecutive water samples in the same pool, and CTL vs EXP water at the same day
    ws = S[(sea == se) & (lay == 'SW')].sort_values('day')
    for pl, sub in ws.groupby('pool'):
        idx = [S.index.get_loc(x) for x in sub.index]
        for i1, i2 in zip(idx[:-1], idx[1:]):
            m1 = np.zeros(len(S), bool); m1[i1] = True; m2 = np.zeros(len(S), bool); m2[i2] = True
            r, *_ = fst(pooled(m1), pooled(m2), mind=6)
            if r['sites'] >= 50: out['pairs'].append(dict(season=se, pool=pl, day=int(day[i2]), kind='water_vs_water_next', **r))
# per-gene FST ice vs water + pN/pS
cds = pd.read_csv(f'{D}/ref/binned_cds.tsv.gz', sep='\t', header=None, names=['contig', 'start', 'end', 'strand', 'gene'])
cds = cds[cds.contig.isin(set(pd.read_csv(f'{D}/ref/contig_bin_len.tsv', sep='\t').query('bin == "dastool-' + b + '"').contig))]
rows = []
site_key = pd.DataFrame({'contig': contig, 'pos': pos})
for se, (ok, num, den) in gene_inputs.items():
    sk = site_key[ok].assign(num=num, den=den)
    for c, sub in sk.groupby('contig'):
        g = cds[cds.contig == c]
        if g.empty: continue
        st, en = g.start.values, g.end.values
        j = np.searchsorted(st, sub.pos.values, side='right') - 1
        inside = (j >= 0) & (sub.pos.values <= en[np.clip(j, 0, None)])
        t = sub[inside].assign(gene=g.gene.values[j[inside]])
        for gid, x in t.groupby('gene'):
            rows.append((se, gid, len(x), x.num.sum(), x.den.sum()))
G = pd.DataFrame(rows, columns=['season', 'gene', 'sites', 'num', 'den'])
G['fst'] = G.num / G.den.where(G.den > 0)
G.to_csv(f'{D}/pop/{b}.genes.tsv', sep='\t', index=False)
# pN/pS (Nei-Gojobori style) on all samples pooled, and per season x ice/water
CODE = {a + b_ + c: aa for (a, b_, c), aa in zip([(x, y, z_) for x in 'TCAG' for y in 'TCAG' for z_ in 'TCAG'],
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
comp = str.maketrans('ACGTacgt', 'TGCAtgca')
def nsites(cod):
    n = 0.
    for i in range(3):
        for x in 'ACGT':
            if x == cod[i]: continue
            m = cod[:i] + x + cod[i + 1:]
            if CODE.get(m, '?') != CODE.get(cod, '?') and CODE.get(m) != '*': n += 1 / 3
    return n
cont = sorted(set(cds.contig))
seq = {}
if cont:
    fa = subprocess.run(['samtools', 'faidx', f'{D}/ref/allbins.fa'] + cont, capture_output=True, text=True).stdout
    for blk in fa.split('>')[1:]:
        h, *s = blk.split('\n'); seq[h.split()[0]] = ''.join(s).upper()
ACGT = 'ACGT'
NSc = {c_: nsites(c_) for c_ in CODE}
def pnps(mask):
    a = pooled(mask); d = a.sum(1); f = a / np.maximum(d, 1)[:, None]
    sort = np.argsort(-a, 1); maj, mi = sort[:, 0], sort[:, 1]
    poly = (d >= MIN) & (f[np.arange(len(f)), mi] >= .05) & (a[np.arange(len(a)), mi] >= 2)
    Nsnv = Ssnv = 0; Nsite = Ssite = 0.
    polyset = {}
    for c_, p_, m_, n_ in zip(contig[poly], pos[poly], maj[poly], mi[poly]): polyset[(c_, p_)] = (ACGT[m_], ACGT[n_])
    for r in cds.itertuples():
        s_ = seq.get(r.contig)
        if s_ is None: continue
        L = (r.end - r.start + 1) // 3 * 3
        nt = s_[r.start - 1:r.start - 1 + L] if r.strand == '+' else s_[r.end - L:r.end][::-1].translate(comp)
        for k in range(0, L, 3):
            cod = nt[k:k + 3]
            if cod not in NSc: continue
            ns = NSc[cod]; Nsite += ns; Ssite += 3 - ns
        hits = [(p_, v) for (c_, p_), v in polyset.items() if c_ == r.contig and r.start <= p_ <= r.end]
        for p_, (ma, mn) in hits:
            off = (p_ - r.start) if r.strand == '+' else (r.end - p_)
            k = off // 3 * 3
            if k + 3 > L: continue
            cod = list(nt[k:k + 3])
            if r.strand == '-': ma, mn = ma.translate(comp), mn.translate(comp)
            c1 = cod.copy(); c1[off % 3] = ma; c2 = cod.copy(); c2[off % 3] = mn
            a1, a2 = CODE.get(''.join(c1)), CODE.get(''.join(c2))
            if a1 is None or a2 is None: continue
            if a1 == a2: Ssnv += 1
            else: Nsnv += 1
    return dict(N=Nsnv, S=Ssnv, pNpS=float((Nsnv / Nsite) / (Ssnv / Ssite)) if Ssnv and Nsite else None)
out['pnps'] = {}
for se in ('2024-25', '2025-26'):
    for nm, m in (('water', M(sea=se, lay=['SW'])), ('ice', M(sea=se) & ice)):
        if m.any() and C[:, m, :].sum() > 0: out['pnps'][f'{se}|{nm}'] = pnps(m)
json.dump(out, open(f'{D}/pop/{b}.json', 'w'), indent=1)
print(b, json.dumps(out['contrasts'])[:600])
