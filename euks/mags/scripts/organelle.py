"""EXPLORATORY, NOT PART OF THE PIPELINE (see NOTES.md "Pending"): it rewrites the organelle tables with a flawed split
(the parent ID goes to the largest bin and the unbinned remainder is lumped into it). Fix before use.
Organelle (P/M) groups: genome count from single-copy markers, split into one P## per plastid, host links, whole-assembly check.
Markers (Bakta products on group contigs): plastid psbA, psbD, rbcL, psaA, psaB, atpB, petB; mito cox1, cob.
  genomes implied = median copy number across the markers present.
Split rule: a P/M group made of several organelle bins (data/bin_domain.tsv 'organelle', or plastid_bins.tsv) is split into one
  group per bin when the group implies > 1 genome and each bin is >= 50 kb and carries >= 1 marker; contigs in no such bin stay in
  the parent ID only if the parent keeps a bin, else they are listed as 'unassigned'. Profiles check the split: median Spearman of
  contig RNA profiles within bins vs between bins (a clean split has within > between).
  The largest bin keeps the parent ID (e.g. P01); others get IDs from the shared registry (key 'orgsplit:<parent>:<bin>').
Host link: Spearman rho of each plastid's summed RNA profile (128 libraries) against every E## group's; top match reported.
Lineage per marker copy: the plastid 16S host of its bin (plastid_bins.tsv, SILVA best hit) when present.
Whole-assembly check: every Tiara 'organelle' contig >= 1 kb, where it sits (P/M group, organelle bin, other bin, unbinned).
Out: out/organelle_markers.tsv, out/organelle_groups.tsv and out/contig_organelle_group.tsv (rewritten with split groups),
     out/organelle_host_links.tsv, out/organelle_contigs_all.tsv"""
import pandas as pd, numpy as np, sys, scipy.io, re
from scipy.stats import spearmanr
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/euks/mags/scripts')
from common import samples, all_samples, D as DD
from eukshared import load_ids, save_ids, allocate
O = '/data/dev/rnaseq-exploration/euks/mags/out/'
T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
OG = pd.read_csv(O + 'organelle_groups.tsv', sep='\t'); CO = pd.read_csv(O + 'contig_organelle_group.tsv', sep='\t')
GM = pd.read_pickle(DD + 'genemap.pkl')[['gene_id', 'contig', 'product']]
dom = pd.read_csv(DD + 'bin_domain.tsv', sep='\t').set_index('bin').call
pl = pd.read_csv(DD + 'plastid_bins.tsv', sep='\t').set_index('bin')
MK = {'psbA': r'Photosystem II protein D1|\bpsbA\b', 'psbD': r'Photosystem II D2|\bpsbD\b', 'rbcL': r'Ribulose bisphosphate carboxylase large', 'psaA': r'P700 chlorophyll a apoprotein A1|\bpsaA\b',
      'psaB': r'P700 chlorophyll a apoprotein A2|\bpsaB\b', 'atpB': r'ATP synthase subunit beta|\batpB\b', 'petB': r'[Cc]ytochrome b6\b|\bpetB\b', 'cox1': r'[Cc]ytochrome c oxidase subunit 1\b|\bcox1\b', 'cob': r'[Cc]ytochrome b\b(?!6)|\bcob\b'}
def host(b):
    if b in pl.index:
        h = pl.host.get(b); return h if isinstance(h, str) and h and not h.startswith(('uncultured', 'environmental')) else 'plastid, host unidentified'
    return ''
# RNA profiles per contig
Cm = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsr(); ALL = all_samples(); S = samples(); kc = [ALL.index(x) for x in S]
allids = np.sort(T.index.values); pos = pd.Series(np.arange(len(allids)), index=allids); tot = np.asarray(Cm[:, kc].sum(0)).ravel()
def rna(ids): return pd.DataFrame(Cm[pos[list(ids)].values][:, kc].toarray() / tot * 1e6, index=list(ids), columns=S)
IDS = load_ids(); mk_rows, new_groups, new_contigs, links = [], [], [], []
CE = pd.read_csv(O + 'contig_euk_group.tsv', sep='\t')
Eprof = {g: rna(d.contig).sum() for g, d in CE.groupby('euk_group')}
for _, og in OG.iterrows():
    ids = CO.contig[CO.oid == og.oid].tolist(); g = T.loc[ids]
    gm = GM[GM.contig.isin(ids)].copy(); gm['bin'] = T.bin.reindex(gm.contig).fillna('').values
    counts = {}
    for k, rx in MK.items():
        hit = gm[gm['product'].fillna('').str.contains(rx, regex=True)]; counts[k] = len(hit)
        for _, h in hit.iterrows(): mk_rows.append(dict(oid=og.oid, marker=k, gene=h.gene_id, contig=h.contig, bin=h.bin.replace('dastool-', ''), host_16S=host(h.bin)))
    pres = [counts[k] for k in ('psbA', 'psbD', 'rbcL', 'psaA', 'psaB', 'atpB', 'petB') if counts[k] > 0]
    n_gen = int(np.median(pres)) if pres else 0
    bins = g.bin.fillna('')
    orgbins = [b for b in bins.unique() if b and (dom.get(b, '') == 'organelle' or b in pl.index)]
    keep = []
    for b in orgbins:
        bb = g[bins == b]; nm = sum(1 for r in mk_rows if r['oid'] == og.oid and r['bin'] == b.replace('dastool-', ''))
        if bb.length.sum() >= 50000 and nm >= 1: keep.append((b, bb.length.sum(), nm))
    keep.sort(key=lambda t: -t[1]); kb = {t[0] for t in keep}
    # profile check
    R = np.log1p(rna(ids) * 10); within, between = [], []
    lab = bins.reindex(R.index).values
    rr = spearmanr(R.T.values).correlation if len(R) > 2 else np.array([[1]])
    for i in range(len(R)):
        for j in range(i + 1, len(R)):
            if lab[i] and lab[j] and lab[i] in kb and lab[j] in kb: (within if lab[i] == lab[j] else between).append(rr[i, j])
    w, bt = (np.median(within) if within else np.nan), (np.median(between) if between else np.nan)
    split = n_gen > 1 and len(keep) > 1
    parts = [(b, list(g.index[bins == b])) for b, _, _ in keep] if split else [('', ids)]
    if split:
        rest = [c for c in ids if bins[c] not in kb]
        if rest: parts[0] = (parts[0][0], parts[0][1] + rest)      # unbinned / non-organelle-bin contigs stay with the largest part (parent ID)
    for j, (b, pid) in enumerate(parts):
        oid = og.oid if j == 0 else allocate(IDS, f'orgsplit:{og.oid}:{b}', og.oid[0], taken=list(OG.oid))
        s = T.loc[pid]; L = s.length.sum(); prof = rna(pid).sum()
        best = sorted(((e, spearmanr(prof, p).correlation) for e, p in Eprof.items()), key=lambda t: -(t[1] if t[1] == t[1] else -9))[:2]
        cnt = {k: sum(1 for r in mk_rows if r['oid'] == og.oid and r['contig'] in set(pid) and r['marker'] == k) for k in MK}
        new_groups.append(dict(oid=oid, note=(og.note if j == 0 else f'split from {og.oid} (organelle bin {b.replace("dastool-", "")})') + (f'; split into {len(parts)} by organelle bins' if j == 0 and split else ''),
                               kind=og.kind, organelle_bin=b.replace('dastool-', ''), host_16S=host(b) if b else '', contigs=len(s), bp=int(L), gc=round(np.average(s.gc, weights=s.length), 1),
                               tiara_organelle=round(s.length[s.tiara == 'organelle'].sum() / max(s.length[s.tiara.notna()].sum(), 1), 3),
                               markers=' '.join(f'{k}:{v}' for k, v in cnt.items() if v), genomes_implied=int(np.median([v for k, v in cnt.items() if k not in ('cox1', 'cob') and v > 0])) if any(cnt[k] for k in cnt if k not in ('cox1', 'cob')) else 0,
                               expr_ppm_median=round(float(prof.median()), 1), host_E_top=best[0][0], host_E_rho=round(best[0][1], 2), host_E_second=f'{best[1][0]} {best[1][1]:.2f}',
                               profile_r_within_bins=round(w, 2) if w == w else None, profile_r_between_bins=round(bt, 2) if bt == bt else None,
                               genomes_implied_parent=n_gen, status='organelle (no hull yet)'))
        new_contigs += [dict(contig=c, oid=oid, source_bins=T.bin.get(c, '') if isinstance(T.bin.get(c, ''), str) else '', tiara=T.tiara.get(c, ''), gc=T.gc[c], length=T.length[c]) for c in pid]
save_ids(IDS)
pd.DataFrame(mk_rows).to_csv(O + 'organelle_markers.tsv', sep='\t', index=False)
NG = pd.DataFrame(new_groups); NG.to_csv(O + 'organelle_groups.tsv', sep='\t', index=False)
NC = pd.DataFrame(new_contigs); assert not NC.contig.duplicated().any(); NC.to_csv(O + 'contig_organelle_group.tsv', sep='\t', index=False)
# whole-assembly check
A = T[(T.tiara == 'organelle') & (T.length >= 1000)].copy()
A['where'] = np.where(A.index.isin(NC.contig), 'P/M group', np.where(A.bin.map(lambda b: dom.get(b, '') == 'organelle' if isinstance(b, str) else False), 'organelle bin (not in a P/M group)',
                     np.where(A.bin.notna() & (A.bin != ''), 'other (non-organelle) bin', 'unbinned')))
A[['length', 'gc', 'bin', 'tsne_x', 'tsne_y', 'where']].to_csv(O + 'organelle_contigs_all.tsv', sep='\t')
pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 50)
print(NG.drop(columns=['note', 'status']).to_string(index=False))
print(pd.DataFrame(mk_rows).groupby(['marker', 'bin']).size().unstack(fill_value=0).to_string())
print(A.groupby('where').length.agg(['size', 'sum']).assign(Mbp=lambda d: (d['sum'] / 1e6).round(2)).drop(columns='sum').to_string())
ob = A[A['where'] == 'organelle bin (not in a P/M group)'].groupby('bin').length.sum().sort_values(ascending=False)
print('organelle bins outside P/M groups:', len(ob), (ob / 1e3).round(0).head(25).to_dict())
