"""Assemble website-ready cyclo.json."""
import sys, json, re
sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
C = '/data/dev/rnaseq-exploration/cyclo/'
CY = ['dastool-lorbin_184', 'dastool-metabat_9342', 'dastool-metabat_5432', 'dastool-comebin_3787_sub', 'dastool-lorbin_722_sub']
S = samples(); m = meta().loc[S]
# --- genome stats
summ = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin')
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t').set_index('user_genome')
ck = {'dastool-metabat_9342': (99.99, 0.83), 'dastool-metabat_5432': (100.0, 12.84), 'dastool-lorbin_184': (94.66, 10.21)}
sk = pd.read_csv(C + 'skani_gtdb.tsv', sep='\t')
sk['q'] = sk.Query_file.str.extract(r'(dastool-[^/]+)\.fa')[0]
def fa(b):
    seqs = {}; k = None
    for l in open(C + f'bins/{b}.fa'):
        if l[0] == '>': k = l[1:].split()[0]; seqs[k] = []
        else: seqs[k].append(l.strip())
    return {k: ''.join(v) for k, v in seqs.items()}
genes = pd.read_csv(C + 'cyclo_genes_mge.tsv', sep='\t', keep_default_na=False)
stats = []
for b in CY:
    s = fa(b); L = sorted((len(v) for v in s.values()), reverse=True); tot = sum(L)
    n50 = next(x for x, c in zip(L, np.cumsum(L)) if c >= tot / 2)
    gc = sum(v.upper().count('G') + v.upper().count('C') for v in s.values()) / tot
    g = gt.loc[b.replace('dastool-', '')] if b.replace('dastool-', '') in gt.index else gt.loc[b] if b in gt.index else None
    best = sk[sk.q == b].sort_values('ANI', ascending=False).head(3)
    stats.append(dict(bin=b, size_mb=round(tot / 1e6, 3), contigs=len(L), n50_kb=round(n50 / 1e3, 1), gc=round(gc * 100, 1),
                      longest_kb=round(L[0] / 1e3, 1), scg_completeness=float(summ.loc[b, 'SCG_completeness']), scg_redundancy=float(summ.loc[b, 'SCG_redundancy']),
                      checkm2_completeness=ck.get(b, (None, None))[0], checkm2_contamination=ck.get(b, (None, None))[1],
                      n_cds=int(((genes.cbin == b) & (genes['type'] == 'CDS')).sum()),
                      gtdb_r226=None if g is None else g['classification'].split(';')[-1], gtdb_r226_ani=None if g is None else str(g['closest_genome_ani']),
                      skani_r232=[dict(ref=r.Ref_file.replace('_genomic.fna.gz', ''), ani=r.ANI, af_ref=r.Align_fraction_ref, af_query=r.Align_fraction_query, name=r.Ref_name.split(' isolate ')[-1][:40]) for r in best.itertuples()]))
tri = pd.read_csv(C + 'skani_triangle.tsv', sep='\t', skiprows=1, header=None)
names = [x.split('/')[-1].replace('.fa', '') for x in tri[0]]
ani = pd.DataFrame(tri.iloc[:, 1:].values, index=names, columns=names)
# --- HC inventory
cats = [c for c in genes.hc.unique() if c]
inv = genes[genes.cbin.isin(CY + ['unbinned_kaiju_cyclo'])].groupby(['cbin', 'hc']).size().unstack(0).fillna(0).astype(int)
inv = inv.loc[[c for c in inv.index if c]]
ch = pd.read_csv(C + 'canthyd_best.tsv', sep='\t', header=None, names=['gene_id', 'hmm', 'score']).merge(genes[['gene_id', 'cbin', 'contig', 'product']], on='gene_id')
chi = ch.groupby(['hmm', 'cbin']).size().unstack(1).fillna(0).astype(int)
# --- MAG TPM time courses
tpm = table('bin_tpm'); bc = table('bin_counts')
share = bc.loc[CY] / bc.sum() * 100
core = m.design.isin(['CTL_SW', 'CTL_SI', 'EXP_SW', 'EXP_SI'])
tc = []
for s in S:
    if not core[s] and not str(m.loc[s, 'sample_name']).startswith('CMO_EXP_DAY6'): pass
    r = m.loc[s]
    tc.append(dict(sample=s, name=r.sample_name, season=r.season, pool=r.group, medium=r.medium, day=None if np.isnan(r.day) else int(r.day),
                   phase=r.phase, **{b.replace('dastool-', ''): round(float(share.loc[b, s]), 4) for b in CY}))
# --- operons (from operons.json) + per-operon CPM time course
ops = json.load(open(C + 'operons.json'))
cpm = np.load(C + 'cyclo_cpm.npy')
ga = pd.read_csv(C + 'cyclo_genes_annot.tsv', sep='\t', keep_default_na=False)
rowof = dict(zip(ga.gene_id, range(len(ga))))
mg = genes.set_index('gene_id')
for o in ops:
    rows = [rowof[g['id']] for g in o['genes'] if g['cat']]
    v = cpm[rows].sum(0)
    o['cpm_timecourse'] = {s: round(float(x), 2) for s, x in zip(S, v) if x > 0}
    o['plasmid_score'] = float(mg.loc[o['genes'][0]['id'], 'plasmid_score']) if o['genes'][0]['id'] in mg.index else None
    o['genomad_plasmid'] = bool(mg.loc[o['genes'][0]['id'], 'genomad_plasmid'])
    for g in o['genes']:
        if g['id'] in mg.index:
            g['island'] = bool(mg.loc[g['id'], 'island']); g['mob_near'] = bool(mg.loc[g['id'], 'mob_within5'])
ops = sorted(ops, key=lambda o: -o['pairs_oiled_total'])
ind = pd.read_csv(C + 'induction.tsv', sep='\t')
dr = pd.read_csv(C + 'dna_rna_cyclo.tsv', sep='\t')
pc = pd.read_csv(C + 'plasmid_copy.tsv', sep='\t')
out = dict(
    bins=CY, stats=stats, ani_matrix=dict(names=names, values=ani.round(2).values.tolist()),
    hc_inventory=dict(categories=inv.index.tolist(), bins=inv.columns.tolist(), counts=inv.values.tolist()),
    canthyd=dict(hmms=chi.index.tolist(), bins=chi.columns.tolist(), counts=chi.values.tolist(), note='hmmsearch --cut_nc vs CANT-HYD; best hit per protein'),
    canthyd_hits=ch[['gene_id', 'hmm', 'score', 'cbin', 'contig', 'product']].to_dict('records'),
    operons=ops, mag_timecourse=tc,
    induction=ind[ind.n_pre > 0].round(2).replace({np.nan: None}).to_dict('records'),
    dna=dr.round(4).replace({np.nan: None}).to_dict('records'),
    plasmid_copy=pc.to_dict('records'),
    oil_day=OIL_DAY)
json.dump(out, open(C + 'cyclo.json', 'w'), separators=(',', ':'), default=lambda x: None if isinstance(x, float) and np.isnan(x) else x)
import os; print(os.path.getsize(C + 'cyclo.json') / 1e6, 'MB'); print(pd.DataFrame(stats).drop(columns='skani_r232').to_string())
