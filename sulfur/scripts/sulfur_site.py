"""Second pass + site/sulfur.json: Kaiju genera for unbinned genes, dsr/apr direction by host phylum, time courses per category,
Sulfurimonas genomes (inventory, RNA, DNA, ice, pools). Inputs: sulfur/out from sulfur.py; data/ and out/ tables; ice/ outputs."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts'); sys.path.insert(0, '/data/dev/rnaseq-exploration/ice/scripts')
from common import *
O = '/data/dev/rnaseq-exploration/sulfur/out/'; E0 = '/data/dev/rnaseq-exploration/'
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
summ = pd.read_csv(O + 'cat_summary.tsv', sep='\t', index_col=0)
E = pd.read_csv(O + 'cat_cpm.tsv', sep='\t', index_col=0)
G = pd.read_csv(O + 'genes.tsv', sep='\t', keep_default_na=False)
M = pd.read_pickle(O + 'gene_cpm.pkl')
mem = pd.read_pickle(O + 'cat_host_cpm.pkl')['members']
S = list(E.columns); m = meta().loc[S]

# Kaiju genus/phylum for unbinned genes' contigs
need = set(G.contig[G.bin == ''])
kj = {}
with open(D + 'kaiju_contigs.tsv') as fh:
    next(fh)
    for line in fh:
        c, st, tid, lin = line.rstrip('\n').split('\t')
        if c in need: kj[c] = [x.strip() for x in lin.split(';')]
def kgen(c):
    l = kj.get(c)
    if not l: return 'unbinned, unclassified'
    g = l[5] if len(l) > 5 and l[5] not in ('', 'NA') else (l[4] if len(l) > 4 else l[-1])
    return f'unbinned ({g}, Kaiju)'
def kphy(c):
    l = kj.get(c); return l[1] if l and len(l) > 1 else ''
_bid = G.bin.str.replace('dastool-', '', regex=False)
G['host'] = np.where(G.bin == '', [contig_host_label(c, kgen(c)) for c in G.contig], np.where(G.taxon == _bid, G.taxon, G.taxon + ' ' + _bid))   # unbinned contigs in a eukaryotic genome group take its label
G['phylum2'] = np.where(G.bin == '', G.contig.map(kphy), G.phylum)
REDUCERS = 'Desulfobacterota|Thermodesulfobacteriota|Desulfuromonadota|Nitrospirota|Bacillota'

oil = list(m.index[m.phase == 'oiled']); ctl = list(m.index[m.phase == 'control'])
cats = list(summ.index)
J = {'categories': [], 'groups': [], 'timecourse': {}, 'samples': []}
for gname in dict.fromkeys(summ.group): J['groups'].append(gname)
for c in cats:
    if summ.loc[c, 'n_genes'] == 0: continue   # e.g. dsrD: no KO and no matching product
    x = summ.loc[c]; idx = mem[c]
    sub = G.iloc[idx]; mm = M.iloc[idx]
    hosts = {}
    for w, cols in (('diesel', oil), ('control', ctl)):
        v = mm[cols].sum(1).groupby(sub.host.values).sum(); t = v.sum()
        hosts[w] = [{'host': h, 'share': r(s / t, 3)} for h, s in v.sort_values(ascending=False).head(5).items()] if t > 0 else []
    J['categories'].append({'name': c, 'group': x.group, 'n_genes': int(x.n_genes), 'mean': r(x.mean_cpm, 1),
        'lfc24': r(x['log2_2024-25_SW']), 'lfc25': r(x['log2_2025-26_SW']), 'baci25': r(x['baci_2025-26']),
        'oil24': r(x['oil_2024-25_SW'], 1), 'ctl24': r(x['ctl_2024-25_SW'], 1), 'oil25': r(x['oil_2025-26_SW'], 1), 'ctl25': r(x['ctl_2025-26_SW'], 1),
        'pre_exp25': r(x['pre_exp_2025-26'], 1), 'pre_ctl25': r(x['pre_ctl_2025-26'], 1),
        'ice': r(x.ice_log2), 'ice_up': int(x.ice_up), 'ice_n': int(x.ice_n), 'top': r(x.top_log2), 'top_up': int(x.top_up), 'top_n': int(x.top_n),
        'hosts': hosts})
    J['timecourse'][c] = [r(v, 2) for v in E.loc[c, S]]
J['samples'] = [{'id': s, 'name': m.sample_name[s], 'season': m.season[s], 'pool': m.group[s], 'medium': m.medium[s],
                 'day': None if pd.isna(m.day[s]) else int(m.day[s])} for s in S]

# dissimilatory sulfur direction: share of dsrAB / aprAB / qmo transcripts in reducer phyla, by pool and season
dirn = []
for c in ('dsrAB (dissimilatory sulfite reductase)', 'aprAB (APS reductase)', 'qmoABC (reductive-pathway electron transfer)'):
    idx = mem[c]; sub = G.iloc[idx]; mm = M.iloc[idx]
    red = sub.phylum2.str.contains(REDUCERS, regex=True).values; unk = (sub.phylum2 == '').values
    for s in ('2024-25', '2025-26'):
        for w, ph in (('diesel', 'oiled'), ('control', 'control')):
            cols = list(m.index[(m.season == s) & (m.phase == ph)]); t = mm[cols].sum().sum() / max(len(cols), 1)
            if t <= 0: continue
            a = mm[cols][red].sum().sum() / len(cols); u = mm[cols][unk].sum().sum() / len(cols)
            dirn.append({'gene': c.split(' (')[0], 'season': s, 'pool': w, 'cpm': r(t, 2), 'reductive': r(a / t, 3), 'unknown': r(u / t, 3), 'oxidative': r(1 - (a + u) / t, 3)})
J['direction'] = dirn

# Sulfurimonas genomes
SB = ['dastool-semibin_088', 'dastool-lorbin_167', 'dastool-comebin_1140_sub']
q = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin')
inv = []
for b in SB:
    gb = G[G.bin == b]; row = {'bin': b.replace('dastool-', ''), 'size': int(q.loc[b, 'size']), 'compl': int(q.loc[b, 'SCG_completeness']), 'genes': {}}
    for c in cats:
        n = int(np.isin(np.array(mem[c]), np.where(G.bin == b)[0]).sum())
        if n: row['genes'][c] = n
    inv.append(row)
J['sulfurimonas'] = {'inventory': inv}
tpm = table('bin_tpm')
J['sulfurimonas']['tpm'] = {b.replace('dastool-', ''): [r(v, 0) for v in tpm.loc[b, S]] for b in SB if b in tpm.index}
# DNA: Illumina competitive mapping (% of reads) and nanopore (% of binned coverage)
pa = pd.read_csv(E0 + 'ice/popgen/abund_pct.tsv', sep='\t', index_col=0)
sm = pd.read_csv(E0 + 'ice/popgen/samples_map.tsv', sep='\t', index_col=0); sm = sm[~sm.index.isin(excluded_names('metagenome'))]
pa = pa[[c for c in pa.columns if c in sm.index]]
J['sulfurimonas']['dna_illumina'] = {'samples': [{'name': s, 'season': sm.season[s], 'pool': sm.pool[s], 'layer': sm.layer[s], 'day': int(sm.day[s])} for s in pa.columns],
                                     'pct': {b.replace('dastool-', ''): [r(v, 3) for v in pa.loc[b]] for b in SB if b in pa.index}}
dr = pd.read_csv(E0 + 'out/dna_rel_abundance.tsv', sep='\t', index_col=0); di = pd.read_csv(E0 + 'out/dna_samples.tsv', sep='\t', index_col=0)
dk = [k for k in dr.columns if k in di.index and di.loc[k, 'group'] in ('CTL', 'EXP') and pd.notna(di.loc[k, 'day'])]
J['sulfurimonas']['dna_nanopore'] = {'samples': [{'season': di.season[k], 'pool': di.group[k], 'medium': di.medium[k], 'day': int(di.day[k])} for k in dk],
                                     'pct': {b.replace('dastool-', ''): [r(v, 3) for v in dr.loc[b, dk]] for b in SB if b in dr.index}}
iw = pd.read_csv(E0 + 'ice/out/mag_icewater.tsv', sep='\t', index_col=0); tb = pd.read_csv(E0 + 'ice/out/mag_topbot.tsv', sep='\t', index_col=0)
ab = pd.read_csv(E0 + 'ice/popgen/abundance_icewater_dna.tsv', sep='\t', index_col=0)
ps = pd.read_csv(E0 + 'ice/popgen/popgen_summary.tsv', sep='\t').set_index('bin')
baci = pd.read_csv(E0 + 'out/mag_baci.tsv', sep='\t')
J['sulfurimonas']['stats'] = []
for b in SB:
    bb = baci[(baci.bin == b)]
    J['sulfurimonas']['stats'].append({'bin': b.replace('dastool-', ''), 'rna_ice': r(iw.delta.get(b)), 'rna_ice_up': int(iw.n_up.get(b, 0)), 'rna_ice_n': int(iw.n.get(b, 0)),
        'rna_top': r(tb.delta.get(b)), 'dna_ice': r(ab.log2_ice_water.get(b)), 'dna_ice_up': int(ab.n_pairs_up.get(b, 0)) if b in ab.index else None, 'dna_ice_n': int(ab.n_pairs.get(b, 0)) if b in ab.index else None,
        # pool and ice/water FST from 2025-26 (both pools well covered; 2024-25 control ice lost most coverage when the day-18 core was excluded)
        'fst_pool_water': r(ps['fst|2025-26|water_ctl_vs_exp'].get(b), 3) if b in ps.index and ps['sites|2025-26|water_ctl_vs_exp'].get(b, 0) >= 1000 else None,
        'fst_pool_ice': r(ps['fst|2025-26|ice_ctl_vs_exp'].get(b), 3) if b in ps.index and ps['sites|2025-26|ice_ctl_vs_exp'].get(b, 0) >= 1000 else None,
        'fst_ice_water': r(ps['fst|2025-26|ice_vs_water'].get(b), 3) if b in ps.index and ps['sites|2025-26|ice_vs_water'].get(b, 0) >= 1000 else None,
        'baci24': r(bb[bb.season == '2024-25'].interaction.mean()), 'baci25': r(bb[bb.season == '2025-26'].interaction.mean()), 'baci25_q': r(bb[bb.season == '2025-26'].q.min(), 4)})
json.dump(J, open(E0 + 'site/sulfur.json', 'w'), separators=(',', ':'), allow_nan=False)
pd.set_option('display.width', 250)
print(pd.DataFrame(dirn).to_string()); print(json.dumps(J['sulfurimonas']['stats'], indent=0)); print(inv)
for c in J['categories'][:0]: pass
print(os.path.getsize(E0 + 'site/sulfur.json') // 1000, 'kB')
