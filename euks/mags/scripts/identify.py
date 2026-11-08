"""Best-effort taxonomy of the hull-curated eukaryotic MAGs (EE##) from evidence that already exists (2026-11-08).
Run after hull_mags.py. Usage: identify.py [stage ...]  stages: proteins rrna ee08 summarize (default: all; intermediate stages
are cached in out/identify/; delete a file to redo its stage). Python /data/miniforge3/bin/python3.

Inputs
- out/identify/grex/: copies of danaSeq mag_marine outputs, fetched once with scripts/fetch_grex.sh (one ssh; mag_marine is read-only):
  eukaryotic/marferret/marferret_proteins.tsv (MetaEuk gene models vs MarFERReT v1.1.1, best hit per protein; the pipeline left the
  taxonomy column empty, so best_hit -> NCBI taxid comes from the local MarFERReT.v1.1.1.taxonomies.tab.gz and the BLAST taxdump);
  taxonomy/rrna/rrna_genes.tsv (barrnap loci, vsearch best hit against SILVA, identity and SILVA path).
- out/contig_hull_mag.tsv (members), out/contigs_full.pkl (Kaiju, Tiara), RNA contig counts (EE08 stage).

Rules
1. Protein lineage: every MetaEuk protein on a MAG's contigs with a MarFERReT best hit votes for its hit's lineage, weighted by
   alignment length. At each rank (phylum/class/order/family/genus) the call is the top taxon; support = its share of the votes
   that reach that rank (n = proteins reaching it; for kingdom, of all votes, since most protist lineages have no kingdom). The protein call is the deepest rank with support >= 0.5 and n >= MIN_N (20);
   when the median protein identity is < 50% (distant references) the call goes no deeper than class.
2. 18S/28S: rRNA genes >= 400 bp whose SILVA path starts with Eukaryota and is not Chloroplast/Mitochondria are nuclear; plastid and
   mito rRNA are counted separately. Per contig the best gene (SSU before LSU, then identity) represents the contig. The MAG's rRNA
   hit = the best SSU contig (highest identity among genes >= 600 bp, else any); rrna_agree = share of nuclear rRNA contigs whose
   SILVA path contains a name from the protein call's lineage (SILVA synonyms in SYN). rRNA counts only with >= 2 nuclear rRNA
   contigs (single rRNA contigs are often misbinned: multi-copy operons assemble and bin poorly).
3. Fallback when no protein rank passes: the top protein taxon at the deepest rank with support >= 0.3 is accepted if >= 50% of the
   nuclear rRNA contigs (here >= 1 suffices: an agreeing contig corroborates) agree with it; otherwise the deepest SILVA name shared by >= 50% of >= 2 nuclear rRNA contigs; else the
   protein top phylum/class as a guess.
4. Confidence: high = protein call with support >= 0.6 and (rrna_agree >= 0.5 or median protein identity >= 90%, i.e. a close
   reference genome); medium = protein call (support >= 0.5) not contradicted by rRNA (rrna_agree >= 0.5 or < 2 nuclear rRNA
   contigs), or support >= 0.6 from >= 300 proteins against disagreeing rRNA contigs, or the
   rule-3 protein+rRNA agreement, or an rRNA-only call from >= 3 agreeing contigs; low = everything else (best guess only).
5. Kaiju (contigs_full.pkl) is a cross-check only: bp-weighted top phylum among Kaiju-classified contigs (RefSeq has few protists).
6. Notes: known cross-links (config CROSS) and the EE08 RNA analysis (stage ee08).
Completeness: not assessed here; BUSCO is being run separately into out/busco/ and the `completeness` column is left empty.
Out: out/ee_taxonomy.tsv; out/identify/{ee_marferret_proteins.tsv.gz, rrna_ee.tsv, protein_votes.tsv, ee08_contigs.tsv}"""
import sys, os, gzip
import numpy as np, pandas as pd
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
O = '/data/dev/rnaseq-exploration/euks/mags/out/'; W = O + 'identify/'; G = W + 'grex/'; os.makedirs(W, exist_ok=True)
MARF_TAX = '/data/scratch/refdbs/marferret_db/MarFERReT.v1.1.1.taxonomies.tab.gz'; TAXDUMP = '/data/scratch/refdbs/blast/'
CALL_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']    # NCBI ranks; kingdom = Metazoa/Fungi/Viridiplantae
MIN_N, MIN_RRNA, BIG_N, LOW_PID = 20, 2, 300, 50
SYN = {'Bacillariophyta': ['Diatomea'], 'Coscinodiscophyceae': ['Mediophyceae', 'Coscinodiscophytina'], 'Phaeophyceae': ['Phaeophyceae'],
       'Choanoflagellata': ['Choanoflagellata', 'Choanoflagellida'], 'Hexanauplia': ['Copepoda', 'Hexanauplia']}
CROSS = {'EE05': 'plastid metabat_7332_sub in P01 (16S: Thalassiosira pseudonana) is linked to EE05 by RNA profile',
         'EE07': 'host candidate of giant virus contig_0000564 (rho 0.87) and virophage contig_0049668 (rho 0.74)'}


def members():
    CM = pd.read_csv(O + 'contig_hull_mag.tsv', sep='\t'); return CM[CM.mag.str.startswith('EE')][['contig', 'mag']]


def need(f):
    if not os.path.exists(f): sys.exit(f'missing {f}: run scripts/fetch_grex.sh first')


# ---------------- NCBI taxonomy ----------------
_TAX = None
def taxdump():
    global _TAX
    if _TAX is None:
        par, rank, name, merged = {}, {}, {}, {}
        with open(TAXDUMP + 'nodes.dmp') as fh:
            for l in fh:
                f = l.split('\t|\t'); par[int(f[0])] = int(f[1]); rank[int(f[0])] = f[2]
        with open(TAXDUMP + 'names.dmp') as fh:
            for l in fh:
                f = l.split('\t|\t')
                if f[3].startswith('scientific name'): name[int(f[0])] = f[1]
        with open(TAXDUMP + 'merged.dmp') as fh:
            for l in fh:
                f = l.split('\t|\t'); merged[int(f[0])] = int(f[1].split('\t')[0])
        _TAX = (par, rank, name, merged)
    return _TAX


def named(t):
    par, rank, name, merged = taxdump(); t = merged.get(t, t); d = {}
    while t in par and t != 1:
        r = rank.get(t)
        if r in CALL_RANKS + ['superkingdom', 'domain'] and r not in d: d[r] = name.get(t)
        if par[t] == t: break
        t = par[t]
    return d


# ---------------- stages ----------------
def proteins():
    out = W + 'ee_marferret_proteins.tsv.gz'
    if os.path.exists(out): return
    need(G + 'marferret_proteins.tsv')
    M = members().set_index('contig').mag
    P = pd.read_csv(G + 'marferret_proteins.tsv', sep='\t', usecols=['protein_id', 'best_hit', 'pident', 'aln_length', 'bitscore', 'pfam'])
    P['contig'] = P.protein_id.str.split('|').str[1]; P['mag'] = P.contig.map(M); P = P[P.mag.notna()].copy()
    want = set(P.best_hit); tx = {}
    with gzip.open(MARF_TAX, 'rt') as fh:
        next(fh)
        for l in fh:
            f = l.split('\t')
            if f[1] in want: tx[f[1]] = int(f[2])
    P['taxid'] = P.best_hit.map(tx)
    lin = {t: named(int(t)) for t in P.taxid.dropna().unique()}
    for r in ['superkingdom'] + CALL_RANKS: P[r] = [lin.get(t, {}).get(r) if t == t else None for t in P.taxid]
    P.to_csv(out, sep='\t', index=False, compression='gzip'); print('EE proteins with MarFERReT hits:', len(P))


def rrna():
    out = W + 'rrna_ee.tsv'
    if os.path.exists(out): return
    need(G + 'rrna_genes.tsv')
    M = members().set_index('contig').mag
    R = pd.read_csv(G + 'rrna_genes.tsv', sep='\t'); R['mag'] = R.contig_id.map(M); R = R[R.mag.notna() & R.taxonomy.notna()].copy()
    R['sub'] = R.rrna_type.map({'18S_rRNA': 'SSU', '16S_rRNA': 'SSU', '28S_rRNA': 'LSU', '23S_rRNA': 'LSU'}); R = R[R['sub'].notna()]
    R['kind'] = ['plastid' if 'Chloroplast' in x else 'mito' if 'Mitochondria' in x else 'nuclear' if x.startswith('Eukaryota') else 'prok'
                 for x in R.taxonomy]
    R[['mag', 'contig_id', 'start', 'end', 'rrna_type', 'sub', 'kind', 'gene_length', 'vsearch_identity', 'vsearch_coverage', 'taxonomy']] \
        .to_csv(out, sep='\t', index=False)
    print(R.groupby(['mag', 'kind']).size().unstack(fill_value=0).to_string())


def ee08():
    out = W + 'ee08_contigs.tsv'
    if os.path.exists(out): return
    import scipy.io
    from common import samples, all_samples, D as DD
    T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id'); M = members(); ids = M.contig[M.mag == 'EE08'].values
    C = scipy.io.mmread(DD + 'marine_coassembly.contig_counts.mtx.gz').tocsr(); S = samples(); ALL = all_samples(); kc = [ALL.index(x) for x in S]
    pos = pd.Series(np.arange(len(T)), index=np.sort(T.index.values)); R = C[pos[ids].values][:, kc].toarray()
    X = T.loc[ids, ['length', 'depth', 'gc', 'tiara', 'whokaryote', 'kaiju_domain', 'kaiju_phylum', 'kaiju_genus', 'genomad_domain', 'rrna_domain',
                    'rrna_class', 'tsne_x', 'tsne_y', 'bin']].copy()
    X['rna_reads'] = R.sum(1); X['rna_libs'] = (R > 0).sum(1); X = X.sort_values('rna_reads', ascending=False)
    X['rna_cum_share'] = (X.rna_reads.cumsum() / X.rna_reads.sum()).round(4)
    X.index.name = 'contig'; X.reset_index().to_csv(out, sep='\t', index=False)


def votes(P):
    rows = []
    for m, g in P.groupby('mag'):
        d = {'mag': m, 'proteins': len(g), 'pid_median': round(g.pident.median(), 1)}
        for r in CALL_RANKS:
            v = g.dropna(subset=[r]).groupby(r).aln_length.sum().sort_values(ascending=False)
            den = g.aln_length.sum() if r == 'kingdom' else v.sum()     # kingdom: most protists have none, so share of all votes
            d[r] = v.index[0] if len(v) else None; d[r + '_support'] = round(v.iloc[0] / den, 3) if len(v) else 0.0
            d[r + '_n'] = int(g[r].notna().sum()); d[r + '_2nd'] = f'{v.index[1]} {v.iloc[1] / v.sum():.2f}' if len(v) > 1 else ''
        rows.append(d)
    return pd.DataFrame(rows).set_index('mag')


BAD = ('uncultured', 'Incertae Sedis', 'metagenome', 'environmental', 'unidentified')
def silva_names(path):
    return [x for x in path.split(';')[1:] if x and not any(b in x for b in BAD)]


def short_path(path, k=3):
    n = silva_names(path); return ';'.join(n[-k:]) if n else path


def summarize():
    for s in ('proteins', 'rrna', 'ee08'): globals()[s]()
    P = pd.read_csv(W + 'ee_marferret_proteins.tsv.gz', sep='\t'); V = votes(P); V.to_csv(W + 'protein_votes.tsv', sep='\t')
    R = pd.read_csv(W + 'rrna_ee.tsv', sep='\t'); T = pd.read_pickle(O + 'contigs_full.pkl').set_index('id')
    H = pd.read_csv(O + 'hull_mags.tsv', sep='\t').set_index('mag'); M = members()
    rows = []
    for m in sorted(M.mag.unique()):
        h = H.loc[m]; v = V.loc[m] if m in V.index else None
        n = R[(R.mag == m) & (R.kind == 'nuclear') & (R.gene_length >= 400)].copy()
        n['ssu'] = n['sub'] == 'SSU'
        per = n.sort_values(['ssu', 'vsearch_identity'], ascending=False).drop_duplicates('contig_id')
        ssu = n[n.ssu]; ssu_l = ssu[ssu.gene_length >= 600]; best = (ssu_l if len(ssu_l) else ssu).sort_values('vsearch_identity', ascending=False)
        lsu = n[~n.ssu].sort_values('vsearch_identity', ascending=False)
        d = dict(mag=m, Mbp=round(h.bp / 1e6, 2), contigs=int(h.contigs), tiara_euk=h.tiara_euk,
                 rrna_18S=f'{short_path(best.taxonomy.iloc[0])} ({best.vsearch_identity.iloc[0]:.1f}%, {best.gene_length.iloc[0]} bp)' if len(best) else '',
                 rrna_28S=f'{short_path(lsu.taxonomy.iloc[0])} ({lsu.vsearch_identity.iloc[0]:.1f}%, {lsu.gene_length.iloc[0]} bp)' if len(lsu) else '',
                 rrna_nuclear_contigs=len(per),
                 rrna_organelle_genes=int(((R.mag == m) & R.kind.isin(['plastid', 'mito'])).sum()))
        # protein call
        call_rank = None
        if v is not None:
            d.update(proteins=int(v.proteins), protein_pid_median=v.pid_median)
            deepest = CALL_RANKS if v.pid_median >= LOW_PID else CALL_RANKS[:CALL_RANKS.index('class') + 1]
            for r in CALL_RANKS:
                d[f'protein_{r}'] = f'{v[r]} {v[r + "_support"]:.2f}' if v[r] else ''
                if r in deepest and v[r] and v[r + '_support'] >= 0.5 and v[r + '_n'] >= MIN_N: call_rank = r
        def lineage_names(upto):
            out = []
            for r in CALL_RANKS[:CALL_RANKS.index(upto) + 1]:
                x = v[r] if v is not None else None
                if x: out += [x] + SYN.get(x, [])
            return set(out)
        def agree(names, k=MIN_RRNA):
            if len(per) < k: return np.nan
            return round(np.mean([bool(names & set(silva_names(p))) for p in per.taxonomy]), 2)
        conf, rule = 'low', ''
        if call_rank:
            sup = v[call_rank + '_support']; ag = agree(lineage_names(call_rank)); call = f'{v[call_rank]} ({call_rank})'
            d['rrna_agree'] = ag
            if sup >= 0.6 and ((ag == ag and ag >= 0.5) or v.pid_median >= 90):
                conf, rule = 'high', f'protein {call_rank} support {sup:.2f} >= 0.6 and ' + (f'rRNA agrees ({ag:.2f})' if ag == ag and ag >= 0.5 else f'close reference (median id {v.pid_median}%)')
            elif ag != ag or ag >= 0.5:
                conf, rule = 'medium', f'protein {call_rank} support {sup:.2f} >= 0.5; ' + (f'rRNA uninformative ({len(per)} nuclear contigs)' if ag != ag else f'rRNA agrees ({ag:.2f})')
            elif sup >= 0.6 and v[call_rank + '_n'] >= BIG_N:
                conf, rule = 'medium', f'protein {call_rank} support {sup:.2f} from {v[call_rank + "_n"]} proteins outweighs {len(per)} disagreeing rRNA contigs'
            else:
                conf, rule = 'low', f'protein {call_rank} support {sup:.2f} but rRNA disagrees ({ag:.2f})'
        else:
            call = None
            if v is not None:
                for r in CALL_RANKS[::-1]:
                    if v[r] and v[r + '_support'] >= 0.3 and v[r + '_n'] >= MIN_N:
                        ag = agree(set([v[r]] + SYN.get(v[r], [])), k=1)
                        if ag == ag and ag >= 0.5:
                            call, conf, rule = f'{v[r]} ({r})', 'medium', f'protein {r} support {v[r + "_support"]:.2f} >= 0.3 and rRNA agrees ({ag:.2f}, {len(per)} contigs)'
                            d['rrna_agree'] = ag; break
            if call is None and len(per) >= MIN_RRNA:
                names = [silva_names(p) for p in per.taxonomy]; shared = None
                for k in range(1, max(len(x) for x in names) + 1):
                    c = pd.Series([';'.join(x[:k]) for x in names if len(x) >= k]).value_counts()
                    if len(c) and c.iloc[0] >= 0.5 * len(names): shared = (c.index[0].split(';')[-1], int(c.iloc[0]))
                    else: break
                if shared and shared[0] not in ('SAR', 'Amorphea', 'Obazoa'):
                    call = f'{shared[0]} (rRNA)'; conf = 'medium' if shared[1] >= 3 else 'low'
                    rule = f'rRNA only: {shared[1]} of {len(per)} nuclear rRNA contigs share {shared[0]}'
            if call is None and v is not None:
                r = max(['kingdom', 'phylum', 'class'], key=lambda r: (v[r + '_support'] if v[r + '_n'] >= MIN_N else 0))
                call, rule = f'{v[r]}? ({r})', f'best guess: protein {r} support {v[r + "_support"]:.2f} < 0.5'
        d.update(call=call or 'unknown', confidence=conf, rule=rule)
        # Kaiju
        x = T.loc[M.contig[M.mag == m]]; k = x[x.kaiju_phylum.notna()]
        if len(k):
            kp = k.groupby('kaiju_phylum').length.sum().sort_values(ascending=False)
            d['kaiju'] = f'{kp.index[0]} {kp.iloc[0] / kp.sum():.2f} ({k.length.sum() / x.length.sum():.0%} of bp classified)'
        d['completeness'] = ''
        d['notes'] = CROSS.get(m, '')
        rows.append(d)
    S = pd.DataFrame(rows)
    # EE08: what carries its RNA
    E = pd.read_csv(W + 'ee08_contigs.tsv', sep='\t'); top = E[E.rna_cum_share.shift(fill_value=0) < 0.95]
    pf = P[P.contig.isin(top.contig)].pfam.fillna('').str.split(';').explode()
    fam = {'PF00091': 'tubulin', 'PF03953': 'tubulin', 'PF00240': 'ubiquitin', 'PF00183': 'HSP90'}
    hits = sorted({fam[p.split('.')[0]] for p in pf if p.split('.')[0] in fam})
    msg = (f'{len(top)} contigs ({top.length.sum() / 1e3:.0f} kb) carry 95% of its RNA, in {int(top.rna_libs.min())}-{int(top.rna_libs.max())} of 128 libraries; '
           f'their proteins: {", ".join(hits) or "none recognised"} (conserved, highly expressed eukaryotic genes; likely cross-mapped reads '
           f'from many eukaryotes, so the RNA profile is not the genome\'s)')
    S.loc[S.mag == 'EE08', 'notes'] = msg
    cols = ['mag', 'Mbp', 'contigs', 'tiara_euk', 'call', 'confidence', 'rule', 'rrna_18S', 'rrna_28S', 'rrna_nuclear_contigs', 'rrna_agree',
            'rrna_organelle_genes', 'proteins', 'protein_pid_median'] + [f'protein_{r}' for r in CALL_RANKS] + ['kaiju', 'completeness', 'notes']
    S = S.reindex(columns=cols); S.to_csv(O + 'ee_taxonomy.tsv', sep='\t', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    print(S[['mag', 'Mbp', 'call', 'confidence', 'rrna_18S', 'protein_class', 'protein_genus', 'rrna_agree']].to_string(index=False))


if __name__ == '__main__':
    for s in sys.argv[1:] or ['proteins', 'rrna', 'ee08', 'summarize']: globals()[s]()
