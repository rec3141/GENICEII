"""EE## (eukaryotic) and VV## (viral) MAGs from the curated hull polygons in euks/mags/hulls.tsv (picked visually; never redrawn).
Membership: every contig >= 1 kb whose t-SNE point lies inside a polygon (matplotlib.path, boundary radius 0). Eukaryotic hulls
exclude contigs Tiara calls bacteria or archaea (counted). Hulls whose source group is no longer an E/V group, viral hulls
that are >= 50% Tiara organelle, and viral hulls with < 50% of bp in geNomad Nucleocytoviricota or a giant/eukaryotic-virus class
(tailed-phage clusters are not genomes) are skipped (no VV). NB: out/ tables written before that last rule still list VV01;
scripts/bin_domains.py applies the same rule when wiring the site. A contig inside two polygons goes to the polygon whose source group
already holds it, else to the nearer polygon centroid. IDs (shared registry out/id_registry.tsv, key 'hull:<hull id>'): a source
group with one hull -> EE<n>/VV<n> (E05 -> EE05, V01 -> VV01); with several hulls the largest by bp takes EE<n> and the others
take the next number above every E/EE in use.
Summaries: size, GC, N50, Tiara/Whokaryote purity, geNomad virus share, t-SNE RMS, expression (contig reads per million, 128
RNA libraries), BACI/ice from dacore contrasts on RNA and DNA (CMO nanopore depth x length share), bp gained (contigs not in the
source group) and lost (source-group contigs outside the hull).
Out: out/hull_mags.tsv, out/contig_hull_mag.tsv, out/hull_mags_da.tsv, out/hulls/verify/<id>.png"""
import pandas as pd, numpy as np, sys, gzip, scipy.io, re, os
from matplotlib.path import Path
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
from common import samples, all_samples, nanopore_samples, D as DD
from eukshared import load_ids, save_ids
import dacore
O = '/data/dev/rnaseq-exploration/euks/mags/out/'; HF = '/data/dev/rnaseq-exploration/euks/mags/hulls.tsv'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id'); T = T[T.tsne_x.notna() & (T.length >= 1000)]
H = pd.read_csv(HF, sep='\t', comment='#')
_valid = set(pd.read_csv(O + 'contig_euk_group.tsv', sep='\t').euk_group) | set(pd.read_csv(O + 'contig_virus_group.tsv', sep='\t').vid)
skipped_hulls = sorted(set(H.id[~H.source_group.isin(_valid)]))
H = H[H.source_group.isin(_valid)]     # e.g. V02 retired as a plastid: no VV hull
CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t').set_index('contig').euk_group
CV = pd.read_csv(O + 'contig_virus_group.tsv', sep='\t').set_index('contig').vid
src_of = pd.concat([CE, CV])
xy = T[['tsne_x', 'tsne_y']].values
GIANT_EUK_VIRUS = {'Megaviricetes', 'Maveriviricetes', 'Pokkesviricetes', 'Polintoviricetes', 'Alsuviricetes', 'Pisoniviricetes', 'Revtraviricetes'}
inside, cent, excl = {}, {}, {}
for hid, h in H.groupby('id', sort=False):
    h = h.sort_values('vertex_order'); P = Path(h[['x', 'y']].values); m = P.contains_points(xy)
    ids = T.index[m]; src = h.source_group.iloc[0]
    if src.startswith('E'):
        bad = T.tiara.reindex(ids).isin(['bacteria', 'archaea']).values; excl[hid] = (int(bad.sum()), int(T.length[ids[bad]].sum())); ids = ids[~bad]
    else: excl[hid] = (0, 0)
    if not src.startswith('E'):
        tl_ = T.length[ids][T.tiara.reindex(ids).notna().values].sum(); org_ = T.length[ids][(T.tiara.reindex(ids) == 'organelle').values].sum() / tl_ if tl_ else 0
        if org_ >= 0.5: skipped_hulls.append(f'{hid} (organelle {org_:.0%})'); continue
        # viral genome-group rule (2026-11-08): >= 50% of bp geNomad Nucleocytoviricota or a giant/eukaryotic-virus class;
        # tailed-phage (Caudoviricetes) clusters are many phages clustering on composition, not one genome (V01 hull -> retired VV01)
        gv_ = T.genomad_phylum.reindex(ids).eq('Nucleocytoviricota').values | T.genomad_class.reindex(ids).isin(GIANT_EUK_VIRUS).values
        gsh_ = T.length[ids][gv_].sum() / max(T.length[ids].sum(), 1)
        if gsh_ < 0.5: skipped_hulls.append(f'{hid} (giant/eukaryotic-virus share {gsh_:.0%}: phage cluster)'); continue
    inside[hid] = set(ids); cent[hid] = h[['x', 'y']].values.mean(0)
hsrc = H.drop_duplicates('id').set_index('id').source_group; hsrc = hsrc[hsrc.index.isin(inside)]
# overlaps
owner = {}
cnt = {}
for hid, s in inside.items():
    for c in s: cnt.setdefault(c, []).append(hid)
for c, hs in cnt.items():
    if len(hs) == 1: owner[c] = hs[0]; continue
    held = [h for h in hs if src_of.get(c, '') == hsrc[h]]
    if held: owner[c] = held[0]
    else:
        p = T.loc[c, ['tsne_x', 'tsne_y']].values.astype(float); owner[c] = min(hs, key=lambda h: np.hypot(*(p - cent[h])))
n_overlap = sum(len(v) > 1 for v in cnt.values())
M = pd.Series(owner, name='hull')
# IDs
IDS = load_ids()
def num(v): return int(re.sub(r'\D', '', v))
used = {num(v) for v in IDS.values() if re.fullmatch(r'E{1,2}\d+', str(v))} | {num(g) for g in CE.unique()}
hbp = M.groupby(M).apply(lambda s: T.length[s.index].sum())
hid2mag = {}
for src, hs in hsrc.groupby(hsrc):
    hs = sorted(hs.index, key=lambda h: -hbp.get(h, 0)); pre = 'EE' if src.startswith('E') else 'VV'
    for k, h in enumerate(hs):
        key = f'hull:{h}'
        if key not in IDS:
            if k == 0: IDS[key] = f'{pre}{num(src):02d}'
            else:
                n = max(used | {num(v) for v in IDS.values() if str(v).startswith(pre)}) + 1; IDS[key] = f'{pre}{n:02d}'
            used.add(num(IDS[key]))
        hid2mag[h] = IDS[key]
save_ids(IDS)
CM = pd.DataFrame({'contig': M.index, 'mag': M.map(hid2mag).values, 'hull': M.values, 'source_group': M.map(hsrc).values})
CM['in_source'] = [src_of.get(c, '') == s for c, s in zip(CM.contig, CM.source_group)]
CM['tiara'] = T.tiara.reindex(CM.contig).fillna('').values; CM['gc'] = T.gc.reindex(CM.contig).values; CM['length'] = T.length.reindex(CM.contig).values
assert not CM.contig.duplicated().any()
CM.sort_values(['mag', 'contig']).to_csv(O + 'contig_hull_mag.tsv', sep='\t', index=False)
# expression (RNA) and DNA abundance
Cm = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsr(); ALL = all_samples(); S = samples(); kc = [ALL.index(x) for x in S]
allids = np.sort(pd.read_pickle(O + 'contigs_full.pkl').id.values); pos = pd.Series(np.arange(len(allids)), index=allids)
tot = np.asarray(Cm[:, kc].sum(0)).ravel()
ab_r = pd.DataFrame({mg: np.asarray(Cm[pos[g.contig].values][:, kc].sum(0)).ravel() / tot * 1e6 for mg, g in CM.groupby('mag')}, index=S).T
ns = nanopore_samples(); hdr = gzip.open(O + '../dna/depths.txt.gz', 'rt').readline().rstrip('\n').split('\t')
use = [c for c in hdr[3:] if not c.endswith('-var') and c.replace('.sorted.bam', '') in set(ns.run_barcode)]; idx = [hdr.index(c) for c in use]
c2m = dict(zip(CM.contig, CM.mag)); mags = sorted(CM.mag.unique()); mi = {m: i for i, m in enumerate(mags)}
acc = np.zeros((len(mags), len(use))); totals = np.zeros(len(use))
with gzip.open(O + '../dna/depths.txt.gz', 'rt') as fh:
    next(fh)
    for line in fh:
        f = line.split('\t'); L = float(f[1]); v = np.array([float(f[i]) for i in idx]) * L; totals += v
        m = c2m.get(f[0])
        if m: acc[mi[m]] += v
ab_d = pd.DataFrame(acc / totals, index=mags, columns=[c.replace('.sorted.bam', '') for c in use])
R = dacore.run(ab_r, ab_d); R.to_csv(O + 'hull_mags_da.tsv', sep='\t', index=False)
def eff(m, d, c):
    x = R[(R.blob == m) & (R.data == d) & (R.contrast == c)]
    return (round(x.log2_effect.iloc[0], 2), x.q.iloc[0]) if len(x) else (np.nan, np.nan)
# summaries
SG = pd.concat([CE.rename('g'), CV.rename('g')]).reset_index().rename(columns={'index': 'contig'})
rows = []
for mg, g0 in CM.groupby('mag'):
    s = T.loc[g0.contig]; L = s.length.sum(); src = g0.source_group.iloc[0]; hid = g0.hull.iloc[0]
    srcC = SG.contig[SG.g == src]; srcC = srcC[srcC.isin(T.index)]
    gained = s.length[~s.index.isin(srcC)].sum(); lost = T.length[srcC[~srcC.isin(s.index)]].sum()
    if src == 'E21': lost = np.nan     # E21 is split into several hulls; per-hull loss is not meaningful
    tl = s.length[s.tiara.notna()].sum(); wl = s.length[s.whokaryote.notna()].sum()
    srt = np.sort(s.length.values)[::-1]; n50 = int(srt[np.searchsorted(np.cumsum(srt), L / 2)])
    cx, cy = np.average(s.tsne_x, weights=s.length), np.average(s.tsne_y, weights=s.length)
    expr = float(np.median(ab_r.loc[mg])); te = s.length[s.tiara == 'eukarya'].sum() / tl if tl else np.nan
    vir = s.length[s.genomad_domain == 'Viruses'].sum() / L
    b = eff(mg, 'RNA', 'baci_2025-26_SW'); iw = eff(mg, 'RNA', 'ice_vs_water'); bd = eff(mg, 'DNA', 'diesel_vs_control_2024-25')
    is_e = mg.startswith('EE')
    ok = bool(is_e and te >= 0.9 and L >= 4e6 and expr >= 10)
    rows.append(dict(mag=mg, hull=hid, source_group=src, contigs=len(s), bp=int(L), gc=round(np.average(s.gc, weights=s.length), 1), n50=n50,
                     tiara_euk=round(te, 3), tiara_bac=round(s.length[s.tiara == 'bacteria'].sum() / tl, 3) if tl else np.nan,
                     who_euk=round(s.length[s.whokaryote == 'eukaryote'].sum() / wl, 3) if wl else np.nan, virus_bp=round(vir, 3),
                     tsne_rms=round(np.sqrt(np.average((s.tsne_x - cx) ** 2 + (s.tsne_y - cy) ** 2, weights=s.length)), 2), cx=round(cx, 1), cy=round(cy, 1),
                     bp_gained=int(gained), bp_lost=lost, excluded_prok_contigs=excl[hid][0], excluded_prok_bp=excl[hid][1],
                     expr_ppm_median=round(expr, 1), baci_2025_rna=b[0], baci_2025_rna_q=b[1], ice_vs_water_rna=iw[0], ice_vs_water_rna_q=iw[1],
                     dvc_2024_dna=bd[0], dvc_2024_dna_q=bd[1], network_ok=ok,
                     note=H[(H.id == hid) & H.note.notna()].note.iloc[0] if (H[(H.id == hid)].note.notna()).any() else ''))
S2 = pd.DataFrame(rows).sort_values('mag'); S2.to_csv(O + 'hull_mags.tsv', sep='\t', index=False)
# verification panels
os.makedirs(O + 'hulls/verify', exist_ok=True)
for _, x in S2.iterrows():
    h = H[H.id == x.hull].sort_values('vertex_order'); v = h[['x', 'y']].values
    x0, x1 = v[:, 0].min() - 5, v[:, 0].max() + 5; y0, y1 = v[:, 1].min() - 5, v[:, 1].max() + 5
    b = T[T.tsne_x.between(x0, x1) & T.tsne_y.between(y0, y1)]
    fig, ax = plt.subplots(figsize=(6, 6)); ax.scatter(b.tsne_x, b.tsne_y, s=.6, c='#ced4da', lw=0)
    src = SG.contig[SG.g == x.source_group]; sc = b[b.index.isin(src)]; ax.scatter(sc.tsne_x, sc.tsne_y, s=2, c='#fd7e14', lw=0, alpha=.5, label=f'source {x.source_group}')
    mm = b[b.index.isin(CM.contig[CM.mag == x.mag])]; ax.scatter(mm.tsne_x, mm.tsne_y, s=2, c='#1c7ed6', lw=0, alpha=.6, label=f'{x.mag} members')
    ax.plot(np.r_[v[:, 0], v[0, 0]], np.r_[v[:, 1], v[0, 1]], c='k', lw=1.2)
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1); ax.set_aspect('equal'); ax.grid(alpha=.3); ax.legend(fontsize=7, loc='lower left')
    ax.set_title(f'{x.mag} ({x.hull}): {x.bp / 1e6:.1f} Mbp; +{x.bp_gained / 1e6:.1f} gained', fontsize=9)
    plt.tight_layout(); plt.savefig(O + f'hulls/verify/{x.mag}.png', dpi=80); plt.close()
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 30)
print('contigs in >1 polygon:', n_overlap, '| hulls skipped:', skipped_hulls)
print(S2[['mag', 'hull', 'source_group', 'bp', 'gc', 'tiara_euk', 'who_euk', 'virus_bp', 'tsne_rms', 'bp_gained', 'bp_lost', 'excluded_prok_bp', 'expr_ppm_median', 'baci_2025_rna', 'ice_vs_water_rna', 'network_ok']].to_string(index=False))
