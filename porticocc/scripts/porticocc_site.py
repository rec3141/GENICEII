"""site/porticocc.json: Porticoccaceae 50-400-T64 (comebin_2576) page. Inputs: porticocc/out (porticocc.py, cummo_blast.py), data/ and out/ tables, ice/out, islands/out, hgt/mge."""
import sys, json, os; sys.path.insert(0, '/data/dev/rnaseq-exploration/scripts')
from common import *
B = '/data/dev/rnaseq-exploration/'; O = B + 'porticocc/out/'
def r(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)
BIN = 'dastool-comebin_2576'; CMP = ['dastool-lorbin_184', 'dastool-metabat_9342']
S = samples(); m = meta().loc[S]; tpm = table('bin_tpm')
J = {'samples': [{'id': s, 'name': m.sample_name[s], 'season': m.season[s], 'pool': m.group[s], 'medium': m.medium[s],
                  'day': None if pd.isna(m.day[s]) else int(m.day[s])} for s in S]}
J['tpm'] = {b.replace('dastool-', ''): [r(v, 1) for v in tpm.loc[b, S]] for b in [BIN] + CMP}
# genome stats
q = pd.read_csv(D + 'summary.tsv', sep='\t').set_index('bin').loc[BIN]
gt = pd.read_csv(D + 'gtdbtk_taxonomy.tsv', sep='\t', usecols=[0, 1]).drop_duplicates('user_genome').set_index('user_genome')['classification']
baci = pd.read_csv(B + 'out/mag_baci.tsv', sep='\t'); bb = baci[baci.bin == BIN]
iw = pd.read_csv(B + 'ice/out/mag_icewater.tsv', sep='\t', index_col=0).loc[BIN]; tb = pd.read_csv(B + 'ice/out/mag_topbot.tsv', sep='\t', index_col=0).loc[BIN]
core = m[m.group.isin(['CTL', 'EXP']) & (m.medium == 'SW')]
med = lambda s, sel: r(tpm.loc[BIN, core.index[(core.season == s) & sel]].median(), 0)
J['stats'] = {'size': int(q['size']), 'compl': r(q.SCG_completeness, 0), 'contam': r(q.get('SCG_redundancy', np.nan), 0), 'lineage': str(gt.get('comebin_2576', '')),
              'baci': {s: {'est': r(bb[bb.season == s].interaction.mean()), 'q': r(bb[bb.season == s].q.min(), 4)} for s in ('2024-25', '2025-26')},
              'ice': {'delta': r(iw.delta), 'n_up': int(iw.n_up), 'n': int(iw.n), 'q': r(iw.q, 3)}, 'top': {'delta': r(tb.delta), 'n_up': int(tb.n_up), 'n': int(tb.n)},
              'tpm': {s: {'pre_exp': med(s, (core.phase == 'pre') & (core.group == 'EXP')), 'pre_ctl': med(s, (core.phase == 'pre') & (core.group == 'CTL')),
                          'oil': med(s, core.phase == 'oiled'), 'ctl': med(s, core.phase == 'control')} for s in ('2024-25', '2025-26')},
              'tpm_max': r(tpm.loc[BIN, S].max(), 0)}
# DNA: Illumina k-mer screen and nanopore
k = pd.read_csv(B + 'out/mg_kmer.tsv', sep='\t'); k = k[k.pool.isin(['CTL', 'EXP']) & k.day.notna() & k.medium.isin(['SW', 'SI'])]
J['mg'] = [{'season': x.season, 'pool': x.pool, 'medium': x.medium, 'day': int(x.day), 'v': {b: r(getattr(x, b), 4) for b in ('comebin_2576', 'lorbin_184', 'metabat_9342')}} for x in k.itertuples()]
dr = pd.read_csv(B + 'out/dna_rel_abundance.tsv', sep='\t', index_col=0); di = pd.read_csv(B + 'out/dna_samples.tsv', sep='\t', index_col=0)
dk = [c for c in dr.columns if c in di.index and di.loc[c, 'group'] in ('CTL', 'EXP') and pd.notna(di.loc[c, 'day'])]
J['nanopore'] = [{'season': di.season[c], 'pool': di.group[c], 'medium': di.medium[c], 'day': int(di.day[c]), 'v': {b.replace('dastool-', ''): r(dr.loc[b, c], 3) for b in [BIN] + CMP}} for c in dk]
# toolkit
T = pd.read_csv(O + 'toolkit.tsv', sep='\t')
J['toolkit'] = [{'cls': x.toolkit, 'n': int(x.n), 'oil25': r(x.oil25, 1), 'ctl25': r(x.ctl25, 1), 'sh_oil': r(x.pshare_oil25), 'sh_bg': r(x.pshare_bg25), 'sh_oil24': r(x.pshare_oil24), 'lfc': r(x.log2_pshare25)}
                for x in T.itertuples() if x.n > 0]
g = pd.read_csv(O + 'genes_expr.tsv', sep='\t', keep_default_na=False)
for c in g.columns:
    if c.startswith(('cpm_', 'share_', 'pshare_', 'log2')): g[c] = pd.to_numeric(g[c], errors='coerce')
keep = g[(g.toolkit != '') & ~g.toolkit.isin(['motility / chemotaxis', 'terminal oxidases', 'sulfur / hydrogen'])].sort_values(['toolkit', 'cpm_oil_2025-26'], ascending=[True, False])
keep = keep.groupby('toolkit').head(8)
J['genes'] = [{'id': x.gene_id, 'contig': x.contig, 'p': x['product'], 'ko': x.KO, 'cls': x.toolkit, 'oil25': r(x['cpm_oil_2025-26'], 1), 'ctl25': r(x['cpm_ctl_2025-26'], 1),
               'sh_oil': r(x['pshare_oil_2025-26'], 2), 'sh_bg': r(x['pshare_bg_2025-26'], 2)} for _, x in keep.iterrows()]
J['top'] = [{'id': x.gene_id, 'p': x['product'], 'g': x.gene, 'sh_oil': r(x['pshare_oil_2025-26'], 2), 'sh_bg': r(x['pshare_bg_2025-26'], 2), 'cls': x.toolkit}
            for _, x in g.sort_values('pshare_oil_2025-26', ascending=False).head(25).iterrows()]
op = g[g.contig == 'contig_0004457']
i0 = op.index[op.gene_id == 'JPFHOD_001079468'][0]
J['operon'] = [{'id': x.gene_id, 'p': x['product'], 'ko': x.KO, 'oil25': r(x['cpm_oil_2025-26'], 1)} for _, x in op.loc[i0 - 7:i0 + 8].iterrows()]
# CuMMO placement (cummo_blast.py)
if os.path.exists(O + 'cummo_blast.json'):
    C = json.load(open(O + 'cummo_blast.json')); C['best'] = {k: v for k, v in C['best'].items() if v == v}; b = C['best']
    # whole-genome comparison with the 50-400-T64 type MAG (GCA_002162755.1), 1,020-bp fragments, blastn (ANIb-like); see NOTES.md
    C['ani'] = {'acc': 'GCA_002162755.1', 'ani': 80.8, 'af': 0.60}
    C['short'] = f"its closest characterised relatives are the propane/ethane monooxygenases of Cycloclasticus, at about {b['p_cyclo']:.0f}% identity"
    C['summary'] = (f"The CuMMO subunit A of comebin_2576 is {b['p_type']:.0f}% identical to the one in the 50-400-T64 type genome (GCA_002162755.1). "
        "That genome is a MAG from a laboratory simulation of the Deepwater Horizon deep-sea oil plume, where the fastest-growing bacteria used short-chain alkanes "
        "(Hu et al. 2017, <a href=\"https://doi.org/10.1073/pnas.1703424114\">PNAS</a>). The two genomes share about 81% average nucleotide identity over 60% of comebin_2576, "
        "so they are different species of the same genus. The type genome also carries the CuMMO and several cyclohexanone monooxygenases. "
        f"Outside the genus, the closest references are the Cycloclasticus CuMMOs that oxidise propane, ethane and butane ({b['p_cyclo']:.0f}% at most). "
        f"Next come methane monooxygenases ({b['p_pmo']:.0f}%) and bacterial ammonia monooxygenases ({b['p_amo']:.0f}%); the archaeal amoA is unrelated (24%). "
        f"So this is a separate, Porticoccaceae lineage of the CuMMO family, most likely a short-chain or sub-terminal alkane monooxygenase, and not a nitrifier or methanotroph enzyme. "
        f"The Cycloclasticus metabat_9342 CuMMO is divergent too: about {b['c_cyclo']:.0f}% identical to the known Cycloclasticus alkane CuMMOs and {b['c_pmo']:.0f}% to gammaproteobacterial pmoA. "
        "Identities this low cannot place it confidently within the family.")
    J['cummo'] = C
# islands / MGEs
I = pd.read_csv(B + 'islands/out/islands.tsv', sep='\t', keep_default_na=False, low_memory=False); I = I[I.bin == BIN]
G0 = pd.read_pickle(D + 'genemap.pkl'); cs = set(G0.contig[G0.bin == BIN])
pl = pd.read_csv(B + 'hgt/mge/plasmid_summary.tsv', sep='\t'); pl = pl[pl.seq_name.isin(cs)]
J['mge'] = {'islands': [{'id': x.island_id, 'contig': x.contig, 'len': int(x.len), 'n': int(x.n_cds), 'hc': int(x.hc_genes)} for x in I.itertuples()],
            'plasmid_like': [{'contig': x.seq_name, 'len': int(x.length), 'score': r(x.plasmid_score), 'conj': isinstance(x.conjugation_genes, str)} for x in pl.itertuples()],
            'cummo_contig_mobile': bool(('contig_0004457' in set(I.contig)) or ('contig_0004457' in set(pl.seq_name)))}
json.dump(J, open(B + 'site/porticocc.json', 'w'), separators=(',', ':'), allow_nan=False)
print(os.path.getsize(B + 'site/porticocc.json') // 1000, 'kB'); print(json.dumps(J['stats'], indent=0)[:1500])
