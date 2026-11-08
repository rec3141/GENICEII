"""Register kernel blobs (kernel.py) as eukaryotic (E##) or viral (V##) groups, by rule. Run after finalize.py and da.py.
IDs come from the shared registry (eukshared.allocate, key 'blob:<seed>'); the eukshared split test is applied as well.
eukaryote: Tiara-eukarya >= 0.5 and geNomad virus < 0.5; skipped if >= 30% of its bp is already in an E## group; 'fragment' if < 2 Mbp;
  'diffuse' if Tiara-eukarya contigs >= 2 kb in the annulus r..2r around the blob centroid (r = 90th-pct member distance), not in
  the blob, outweigh the blob AND follow the blob's expression/coverage profile (median r >= 0.6): a disc cut from a larger
  co-varying cloud, not a genome with an expression boundary. Otherwise a new E## (next free number).
virus: geNomad virus >= 0.5; skipped if >= 30% already in a V## group; else next free V##.
prokaryote: no ID (reported only).
Appends to out/contig_euk_group.tsv, out/euk_groups.tsv, out/contig_virus_group.tsv, out/virus_groups.tsv; writes out/blob_registry.tsv.
Idempotent: previously registered blob contigs (source 'kernel_blob:...') are removed before re-adding."""
import pandas as pd, numpy as np
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
K = pd.read_csv(O + 'kernel_blobs.tsv', sep='\t'); KM = pd.read_csv(O + 'kernel_blob_contigs.tsv', sep='\t')
DA = pd.read_csv(O + 'da_blobs.tsv', sep='\t')
CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t'); CE = CE[~CE.source_bins.fillna('').str.startswith('kernel_blob:')]
EG = pd.read_csv(O + 'euk_groups.tsv', sep='\t'); EG = EG[~EG.note.fillna('').str.startswith('kernel blob')]
CV = pd.read_csv(O + 'contig_virus_group.tsv', sep='\t'); CV = CV[~CV.source_bins.fillna('').astype(str).str.startswith('kernel_blob:')]
VG = pd.read_csv(O + 'virus_groups.tsv', sep='\t'); VG = VG[~VG.note.fillna('').str.startswith('kernel blob')]
X = pd.read_csv(O + 'da_blobs_abundance_rna.tsv', sep='\t', index_col=0)
RN = pd.read_pickle(O + 'kernel_rna_cpm.pkl'); DN = pd.read_pickle(O + 'kernel_dna_share.pkl')
lR = np.log1p(RN * 10); lD = np.log1p(DN * 1e6)
def prof_corr(ids_in, ids_out):
    # median Pearson r of annulus contigs (with signal) against the blob's bp-weighted mean log profile, RNA and DNA averaged
    out = []
    for M, okf in ((lR, RN.sum(1) * 0 + 1), (lD, None)):
        a = [c for c in ids_in if c in M.index]; b = [c for c in ids_out if c in M.index]
        if len(a) < 5 or len(b) < 5: continue
        w = T.length[a]; p = M.loc[a].mul(w, axis=0).sum() / w.sum()
        Mb = M.loc[b]; Mb = Mb[Mb.sum(1) > 0]
        if len(Mb) < 5: continue
        Mc = Mb.sub(Mb.mean(1), axis=0); pc = p - p.mean()
        out.append(float(np.median((Mc @ pc) / (np.sqrt((Mc ** 2).sum(1)) * np.sqrt((pc ** 2).sum()) + 1e-12))))
    return np.mean(out) if out else np.nan
import sys; sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
from eukshared import load_ids, save_ids, allocate, split_parts
IDS = load_ids()
reg, addC, addG, addCV, addV = [], [], [], [], []
blob_of = DA.drop_duplicates('seed').set_index('seed').blob
for _, x in K.iterrows():
    ids = KM.contig[KM.seed == x.seed].tolist(); g = T.loc[ids]; L = g.length.sum(); b = blob_of.get(x.seed, '')
    if L == 0: reg.append(dict(blob=b, seed=x.seed, decision='empty')); continue
    inE = g.length[g.index.isin(CE.contig)].sum() / L; inV = g.length[g.index.isin(CV.contig)].sum() / L
    cx, cy = np.average(g.tsne_x, weights=g.length), np.average(g.tsne_y, weights=g.length)
    r = np.quantile(np.hypot(g.tsne_x - cx, g.tsne_y - cy), 0.9); dd = np.hypot(T.tsne_x - cx, T.tsne_y - cy)
    annC = T[(dd > r) & (dd <= 2 * r) & (T.tiara == 'eukarya') & (T.length >= 2000) & ~T.index.isin(ids)]
    ann = annC.length.sum(); ann_r = prof_corr(ids, list(annC.index))
    expr = float(X.loc[b].median()) if b in X.index else 0.0
    tl = g.length[g.tiara.notna()].sum(); te = g.length[g.tiara == 'eukarya'].sum() / tl if tl else 0
    common = dict(contigs=len(g), bp=int(L), gc=round(np.average(g.gc, weights=g.length), 1), n50=int(np.sort(g.length.values)[::-1][np.searchsorted(np.cumsum(np.sort(g.length.values)[::-1]), L / 2)]),
                  depth_med=round(g.depth.median(), 2), tiara_euk=round(te, 3), tsne_rms=x.rms, cx=round(cx, 1), cy=round(cy, 1), expr_ppm_median=round(expr, 1))
    if x.virus_bp >= 0.5:
        if inV >= 0.3: dec = 'virus, already in ' + ','.join(sorted(set(CV.vid[CV.contig.isin(ids)])))
        else:
            vid = allocate(IDS, f'blob:{x.seed}', 'V', taken=list(VG.vid)); dec = f'new virus {vid}'
            addV.append(dict(vid=vid, note=f'kernel blob {b} at ({x.seed})', parts=b, status='virus', network_ok=False, virus_bp=x.virus_bp, **common))
            addCV += [dict(contig=c, vid=vid, source_bins=f'kernel_blob:{b}', tiara=T.tiara.get(c, ''), gc=T.gc[c], length=T.length[c]) for c in ids]
    elif x.tiara_euk >= 0.5:
        if inE >= 0.3: dec = 'eukaryote, already in ' + ','.join(sorted(set(CE.euk_group[CE.contig.isin(ids)])))
        elif L < 2e6: dec = 'eukaryote fragment (< 2 Mbp), not registered'
        elif ann > L and ann_r >= 0.6: dec = f'diffuse eukaryotic cloud (surrounding {ann / 1e6:.0f} Mbp with the same expression profile, r {ann_r:.2f}), not registered'
        else:
            parts, sst = split_parts(g)
            if len(parts) > 1: dec = f'split by rule into {len(parts)} (not expected for kernel blobs; Ashman D {sst.get("ashman_D")})'
            eid = allocate(IDS, f'blob:{x.seed}', 'E', taken=list(EG.group)); ok = te >= 0.9 and L >= 4e6 and expr >= 10 and x.rms <= 12
            dec = f'new eukaryote {eid}' + (f'; split test D {sst.get("ashman_D")}' if sst else '')
            addG.append(dict(group=eid, note=f'kernel blob {b} at ({x.seed})', parts=b, status='network' if ok else ('little CMO RNA' if expr < 10 else ('<4 Mbp' if L < 4e6 else 'doubtful')),
                             network_ok=ok, who_euk=round(g.length[g.whokaryote == 'eukaryote'].sum() / max(g.length[g.whokaryote.notna()].sum(), 1), 3), **common))
            addC += [dict(contig=c, euk_group=eid, source_bins=f'kernel_blob:{b}', tiara=T.tiara.get(c, ''), gc=T.gc[c], length=T.length[c]) for c in ids if c not in set(CE.contig)]
    else: dec = 'prokaryote, no ID'
    reg.append(dict(blob=b, seed=x.seed, decision=dec, bp=int(L), tiara_euk=round(te, 2), virus_bp=x.virus_bp, annulus_euk_bp=int(ann), annulus_profile_r=round(ann_r, 2) if ann_r == ann_r else None, expr_ppm_median=round(expr, 1)))
pd.concat([CE, pd.DataFrame(addC)]).to_csv(O + 'contig_euk_group.tsv', sep='\t', index=False)
pd.concat([EG, pd.DataFrame(addG)]).to_csv(O + 'euk_groups.tsv', sep='\t', index=False)
pd.concat([CV, pd.DataFrame(addCV)]).to_csv(O + 'contig_virus_group.tsv', sep='\t', index=False)
pd.concat([VG, pd.DataFrame(addV)]).to_csv(O + 'virus_groups.tsv', sep='\t', index=False)
save_ids(IDS)
R = pd.DataFrame(reg); R.to_csv(O + 'blob_registry.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(R.to_string(index=False))
